"""Multi-axis co-resistance — AMR determinant x plasmid Inc-type (Family C multi-axis deepening, 2026-07-11).

C-deep's cross-axis was AMR-only (AMRFinder virulence 0/845). This adds the PLASMID axis from the local
finder sweep (`scripts/plasmid_axis_sweep.py` -> `data/processed/plasmid_axis_cache.json`): the canonical
co-resistance vehicle. Two questions, WITHIN organism (Enterobacterales), de-confounded + clonality-proxied:

  1. LINKAGE: which AMR determinants / drug-classes travel on which plasmid Inc-type (lift = P(replicon|det)/
     P(replicon)) — the "resistance rides which backbone" map.
  2. IMPUTATION: can a genome's plasmid Inc-type be imputed from its resistance profile, and vice versa
     (held-out OOF-AUC 5-fold logistic, CI>0.5)? — cross-axis "predict the unobserved axis".

Reuses the C-core harness (`determinant_cooccurrence`: harvest, linkage_for_determinant, lift_table,
dedup_profiles). Verdict PASS_MULTIAXIS_LINKAGE iff >= PASS_FRACTION of testable (determinant|replicon) cells
are imputable across axes. Honest: cohorts are drug-selected; plasmid axis is Enterobacterales-only
(enterobacteriales.fsa); PlasmidFinder presence is a genotype axis (not a wet plasmid readout). Frozen surface
READ-only. Run AFTER the plasmid sweep completes.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from datetime import date as _date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import determinant_cooccurrence as dcc  # noqa: E402

PLASMID_CACHE = REPO / "data" / "processed" / "plasmid_axis_cache.json"
MIN_GENOMES = dcc.MIN_GENOMES
MIN_SUPPORT = dcc.MIN_SUPPORT
PASS_FRACTION = 0.5
SEED = 0
REP_PREFIX = "plasmid:"   # namespace Inc-type features so they never collide with determinant tokens


def load_plasmid(cache: Path):
    d = json.loads(cache.read_text(encoding="utf-8"))
    return {acc: v for acc, v in d.items()}


def build_joint(acc_org, acc_dets, plasmid):
    """Per organism -> (accs, feature matrix over determinants + plasmid:replicon, names, is_replicon mask)."""
    per_org = {}
    orgs = [o for o, n in Counter(acc_org.values()).most_common() if n >= MIN_GENOMES]
    for org in orgs:
        accs = [a for a, o in acc_org.items() if o == org and a in plasmid]
        if len(accs) < MIN_GENOMES:
            continue
        feats = {}
        for a in accs:
            s = set(acc_dets.get(a, set()))
            s |= {REP_PREFIX + r for r in plasmid[a]["replicons"]}
            feats[a] = s
        counts = Counter()
        for a in accs:
            counts.update(feats[a])
        names = sorted(n for n, c in counts.items() if c >= 1)
        nidx = {n: j for j, n in enumerate(names)}
        X = np.zeros((len(accs), len(names)), float)
        for i, a in enumerate(accs):
            for n in feats[a]:
                X[i, nidx[n]] = 1.0
        per_org[org] = (accs, X, names, nidx)
    return per_org


def run(acc_dets_repo=REPO, cache=PLASMID_CACHE, dedup=False, boot=dcc.BOOT, seed=SEED):
    acc_org_amr, acc_dets, _ = dcc.harvest(acc_dets_repo)
    plasmid = load_plasmid(cache)
    # organism from the plasmid cache (Enterobacterales), fall back to AMR harvest
    acc_org = {a: plasmid[a]["organism"] for a in plasmid}
    per_org = build_joint(acc_org, acc_dets, plasmid)

    results = {}
    tot_test = tot_link = 0
    for org, (accs, X, names, nidx) in per_org.items():
        is_rep = np.array([n.startswith(REP_PREFIX) for n in names])
        support = X.sum(axis=0)
        testable = [j for j in range(len(names))
                    if MIN_SUPPORT <= int(support[j]) <= len(accs) - MIN_SUPPORT]
        Xd = X
        note = None
        if dedup:
            keep = dcc.dedup_profiles(X, accs)
            Xd = X[keep]
            support = Xd.sum(axis=0)
            testable = [j for j in range(len(names))
                        if MIN_SUPPORT <= int(support[j]) <= Xd.shape[0] - MIN_SUPPORT]
            note = f"profile-deduped {len(accs)}->{len(keep)}"
        # CROSS-AXIS imputation: only score cells whose target uses the OTHER axis as predictors matters most;
        # we report per-feature imputability but headline the CROSS pairs (replicon<->determinant).
        per_feat = {}
        for j in testable:
            r = dcc.linkage_for_determinant(Xd, j, seed=seed, boot=boot)
            if r:
                per_feat[names[j]] = {"auc": r["auc"], "ci_lo": r["ci_lo"], "imputable": r["linked"],
                                      "is_replicon": bool(is_rep[j])}
        n_test = len(per_feat)
        n_link = sum(1 for v in per_feat.values() if v["imputable"])
        tot_test += n_test; tot_link += n_link
        # lift: determinant -> replicon (which resistance rides which backbone)
        rep_names = [n for n in names if n.startswith(REP_PREFIX)]
        det_targets = [n for n, c in Counter({names[j]: int(support[j]) for j in testable
                                              if not is_rep[j]}).most_common(8)]
        lift = dcc.lift_table(Xd, names, nidx, det_targets, min_cooc=dcc.MIN_COOC, top=5)
        # keep only replicon partners in the lift (the cross-axis map)
        cross_lift = {t: [x for x in rows if x["det"].startswith(REP_PREFIX)][:4] for t, rows in lift.items()}
        results[org] = {"organism": org, "n_genomes": Xd.shape[0], "n_replicons": len(rep_names),
                        "n_testable": n_test, "n_imputable": n_link,
                        "fraction_imputable": round(n_link / n_test, 3) if n_test else 0.0,
                        "note": note, "per_feature": per_feat, "determinant_to_plasmid_lift": cross_lift}
    frac = (tot_link / tot_test) if tot_test else 0.0
    verdict = ("PASS_MULTIAXIS_LINKAGE" if (tot_test and frac >= PASS_FRACTION)
               else ("FAIL_AXES_INDEPENDENT" if tot_test else "NO_TESTABLE"))
    return {
        "artifact": "coresistance_multiaxis", "schema": "coresistance-multiaxis-v1",
        "question": "Do AMR determinants and plasmid Inc-types co-occur / impute each other, within organism? "
                    "(which resistance rides which plasmid backbone)",
        "substrate": "AMR determinants (AMRFinder) + plasmid Inc-replicons (PlasmidFinder blastn, local sweep) "
                     "over Enterobacterales cohort genomes",
        "prereg": {"MIN_GENOMES": MIN_GENOMES, "MIN_SUPPORT": MIN_SUPPORT, "PASS_FRACTION": PASS_FRACTION,
                   "BOOT": boot, "seed": seed, "deconfound": "within-organism + optional dedup clonality proxy"},
        "deduped": dedup, "verdict": verdict,
        "total_testable": tot_test, "total_imputable": tot_link, "fraction_imputable": round(frac, 3),
        "honest_caveats": [
            "Plasmid axis is Enterobacterales-only (enterobacteriales.fsa); other organisms excluded.",
            "PlasmidFinder replicon PRESENCE is a genotype axis (blastn), not a wet plasmid-content readout.",
            "Cohorts are drug-R/S-selected; within-organism de-confound + dedup clonality proxy; associational.",
        ],
        "n_plasmid_genomes": len(plasmid),
        "per_organism": results,
    }


def render_md(res, generated):
    L = [f"# Multi-axis co-resistance — AMR determinant × plasmid Inc-type ({generated})", "",
         f"**Verdict: {res['verdict']}** — {res['total_imputable']}/{res['total_testable']} within-organism "
         f"features (determinant OR Inc-type) impute from the joint set (fraction {res['fraction_imputable']}; "
         f"bar {res['prereg']['PASS_FRACTION']}). {'(deduped)' if res['deduped'] else '(raw)'}", "",
         f"{res['question']} Substrate: {res['substrate']}.", "",
         "| organism | genomes | Inc-types | testable | **imputable** | frac | note |",
         "|---|---|---|---|---|---|---|"]
    for r in res["per_organism"].values():
        L.append(f"| {r['organism']} | {r['n_genomes']} | {r['n_replicons']} | {r['n_testable']} | "
                 f"**{r['n_imputable']}** | {r['fraction_imputable']} | {r.get('note') or ''} |")
    L += ["", "## Which resistance rides which plasmid backbone (determinant → Inc-type, by lift)"]
    for r in res["per_organism"].values():
        L.append(f"### {r['organism']}")
        for det, rows in r["determinant_to_plasmid_lift"].items():
            if rows:
                terms = ", ".join(f"{x['det'].replace('plasmid:','')}(lift {x['lift']})" for x in rows)
                L.append(f"- **{det}** → {terms}")
    L += ["", "## Honest caveats"] + [f"- {c}" for c in res["honest_caveats"]]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dedup-profiles", action="store_true")
    ap.add_argument("--boot", type=int, default=None)
    ap.add_argument("--cache", type=Path, default=PLASMID_CACHE)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    if not a.cache.exists():
        print(f"ERROR: plasmid axis cache absent at {a.cache} (run scripts/plasmid_axis_sweep.py first)",
              file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    res = run(cache=a.cache, dedup=a.dedup_profiles, boot=a.boot or dcc.BOOT)
    tag = "_dedup" if a.dedup_profiles else ""
    out = a.out or (REPO / "wiki" / f"coresistance_multiaxis_result_{today}{tag}.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(res, today), encoding="utf-8")
    print(render_md(res, today))
    print(f"\n[wrote {out} + .md]  verdict={res['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
