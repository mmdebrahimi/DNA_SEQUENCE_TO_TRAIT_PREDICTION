"""Co-resistance imputation — Family C deep-dive of the genome world-model plan (2026-07-11).

C-core showed determinant-level linkage exists. This deepens it to the PHENOTYPE-MAPPED, ACTIONABLE
capability the world model is for: **impute a genome's resistance to a drug CLASS it wasn't tested for,
from its determinants for the OTHER classes** — within organism. This is the "predict the unobserved part
of the resistance profile" payload (the joint structure made predictive), plus the co-resistance NETWORK
(which drug classes predict which).

WHY class-level (not determinant-level again): the AMRFinder `Class` field maps each determinant to the
drug(s) it confers, so a genome's class-presence vector IS its (determinant-implied) resistance profile.
Predicting a held-out class = imputing unobserved resistance — directly phenotype-relevant, unlike C-core's
determinant-level linkage.

Honest scope discovered from the cached data (2026-07-11): the virulence/stress CROSS-axis is NOT available
(0/845 genomes carry AMRFinder VIRULENCE determinants — the runs were AMR-only; STRESS is sparse), so this
deep-dive is AMR-class only. The full multi-axis Bakta/plasmid/serotype genome-map sweep (days of Docker)
remains the named Databricks-scale follow-on.

METHOD (reuses the C-core harness `determinant_cooccurrence`):
  * Per genome, drug-CLASS presence = union of the (split-on-'/') AMRFinder Classes of its AMR determinants.
  * Per organism (>= MIN_GENOMES): for each class with both-class support (>= MIN_SUPPORT present & absent),
    5-fold stratified OOF logistic imputing its presence from the OTHER classes; held-out AUC + bootstrap CI.
    IMPUTABLE iff AUC 95% CI lower > 0.5. Same within-organism de-confound + optional profile-dedup clonality
    proxy as C-core.
  * Co-resistance NETWORK: per ordered class pair (A -> B) lift = P(B|A)/P(B) (the imputation structure).

PRE-REGISTERED: PASS_CORESISTANCE_IMPUTABLE iff >= PASS_FRACTION (=0.5) of (organism, class) cells are
IMPUTABLE; else FAIL_CLASSES_INDEPENDENT. Frozen AMR surface READ-only.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from datetime import date as _date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import determinant_cooccurrence as dcc  # noqa: E402  (reuse C-core harness)

MIN_GENOMES = dcc.MIN_GENOMES
MIN_SUPPORT = dcc.MIN_SUPPORT
PASS_FRACTION = 0.5
MIN_COOC = dcc.MIN_COOC
SEED = 0


def genome_classes(acc_dets, det_class, acc):
    """The set of drug CLASSES a genome carries a determinant for (compound 'A/B' classes split)."""
    out = set()
    for d in acc_dets[acc]:
        cls = det_class.get(d)
        if cls:
            for c in cls.split("/"):
                c = c.strip()
                if c:
                    out.add(c)
    return out


def class_matrix_for_organism(acc_org, acc_dets, det_class, org):
    accs = sorted(a for a, o in acc_org.items() if o == org)
    per = {a: genome_classes(acc_dets, det_class, a) for a in accs}
    counts = Counter()
    for a in accs:
        counts.update(per[a])
    classes = sorted(counts)
    cidx = {c: j for j, c in enumerate(classes)}
    X = np.zeros((len(accs), len(classes)), float)
    for i, a in enumerate(accs):
        for c in per[a]:
            X[i, cidx[c]] = 1.0
    testable = [c for c in classes if MIN_SUPPORT <= int(X[:, cidx[c]].sum()) <= len(accs) - MIN_SUPPORT]
    return X, classes, cidx, accs, testable


def run_organism(acc_org, acc_dets, det_class, org, dedup=False, seed=SEED, boot=dcc.BOOT):
    X, classes, cidx, accs, testable = class_matrix_for_organism(acc_org, acc_dets, det_class, org)
    note = None
    if dedup:
        keep = dcc.dedup_profiles(X, accs)
        X = X[keep]
        note = f"profile-deduped {len(accs)}->{len(keep)} genomes (clonality proxy)"
        testable = [c for c in classes if MIN_SUPPORT <= int(X[:, cidx[c]].sum()) <= X.shape[0] - MIN_SUPPORT]
    per_class = {}
    for c in testable:
        r = dcc.linkage_for_determinant(X, cidx[c], seed=seed, boot=boot)
        if r:
            per_class[c] = {"auc": r["auc"], "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"],
                            "n_present": r["n_present"], "imputable": r["linked"]}
    n_test = len(per_class)
    n_imp = sum(1 for r in per_class.values() if r["imputable"])
    net = dcc.lift_table(X, classes, cidx, testable, min_cooc=MIN_COOC, top=5)
    return {"organism": org, "n_genomes": X.shape[0], "n_classes": len(classes),
            "n_testable": n_test, "n_imputable": n_imp,
            "fraction_imputable": round(n_imp / n_test, 3) if n_test else 0.0,
            "note": note, "per_class_imputation": per_class, "coresistance_network": net}


def run_all(repo=REPO, dedup=False, seed=SEED, boot=dcc.BOOT):
    acc_org, acc_dets, det_class = dcc.harvest(repo)
    orgs = [o for o, n in Counter(acc_org.values()).most_common() if n >= MIN_GENOMES]
    per_org = {o: run_organism(acc_org, acc_dets, det_class, o, dedup=dedup, seed=seed, boot=boot) for o in orgs}
    tot_test = sum(r["n_testable"] for r in per_org.values())
    tot_imp = sum(r["n_imputable"] for r in per_org.values())
    frac = (tot_imp / tot_test) if tot_test else 0.0
    verdict = ("PASS_CORESISTANCE_IMPUTABLE" if (tot_test and frac >= PASS_FRACTION)
               else ("FAIL_CLASSES_INDEPENDENT" if tot_test else "NO_TESTABLE_CLASSES"))
    return {
        "artifact": "coresistance_imputation_world_model",
        "schema": "coresistance-imputation-v1",
        "question": "Can a genome's resistance to a drug CLASS be imputed from its OTHER-class determinants, "
                    "within organism? (the 'predict the unobserved resistance profile' world-model payload)",
        "substrate": "cached AMRFinder AMR determinants -> drug-CLASS presence (self-distillation); AMR-only "
                     "(0/845 genomes carry VIRULENCE determinants -> no cross-axis; STRESS sparse)",
        "prereg": {"MIN_GENOMES": MIN_GENOMES, "MIN_SUPPORT": MIN_SUPPORT, "PASS_FRACTION": PASS_FRACTION,
                   "BOOT": boot, "seed": seed,
                   "imputable_rule": "held-out OOF-AUC 95% CI lower > 0.5 (5-fold stratified logistic on other classes)",
                   "deconfound": "within-organism + optional profile-dedup clonality proxy"},
        "deduped": dedup,
        "verdict": verdict, "total_testable": tot_test, "total_imputable": tot_imp,
        "fraction_imputable": round(frac, 3),
        "honest_caveats": [
            "Cross-axis (virulence/stress) NOT available from cached AMR-only runs; AMR drug-class only.",
            "Cohorts are drug-R/S-SELECTED -> the co-resistance structure reflects the curated cohorts.",
            "Within-organism de-confounds species; dedup is a crude clonality proxy (Mash-collapse is the follow-on).",
            "Class presence is DETERMINANT-implied (AMRFinder Class), not measured MIC per drug — a genotype axis.",
            "Self-distillation from our own caller: associational, NOT causal.",
        ],
        "per_organism": per_org,
    }


def render_md(res, generated):
    L = [f"# Co-resistance imputation — can unobserved drug-class resistance be imputed? ({generated})", "",
         f"**Verdict: {res['verdict']}** — {res['total_imputable']}/{res['total_testable']} (organism, drug-class) "
         f"cells are IMPUTABLE from the other classes (fraction {res['fraction_imputable']}; PASS bar "
         f"{res['prereg']['PASS_FRACTION']}). {'(profile-deduped)' if res['deduped'] else '(raw)'}", "",
         f"{res['question']}", "",
         "IMPUTABLE = a genome's held-out drug-class presence is predicted from its OTHER-class determinants, "
         "within organism, OOF-AUC 95% CI lower > 0.5. This is the 'impute the unobserved resistance' payload.", "",
         "| organism | genomes | classes | testable | **imputable** | frac | note |",
         "|---|---|---|---|---|---|---|"]
    for r in res["per_organism"].values():
        L.append(f"| {r['organism']} | {r['n_genomes']} | {r['n_classes']} | {r['n_testable']} | "
                 f"**{r['n_imputable']}** | {r['fraction_imputable']} | {r.get('note') or ''} |")
    L += ["", "## Per-class imputation AUC (held-out; how well each drug class is imputed from the others)"]
    for r in res["per_organism"].values():
        rows = sorted(r["per_class_imputation"].items(), key=lambda kv: -kv[1]["auc"])
        line = ", ".join(f"{c}({m['auc']}{'*' if m['imputable'] else ''})" for c, m in rows[:8])
        L.append(f"- **{r['organism']}**: {line}")
    L += ["", "## Co-resistance NETWORK (drug-class A → classes it predicts, by lift)"]
    for r in res["per_organism"].values():
        L.append(f"### {r['organism']}")
        for a, rows in r["coresistance_network"].items():
            if rows:
                terms = ", ".join(f"{x['det']}(lift {x['lift']})" for x in rows[:4])
                L.append(f"- **{a}** → {terms}")
    L += ["", "## Honest caveats"] + [f"- {c}" for c in res["honest_caveats"]]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dedup-profiles", action="store_true")
    ap.add_argument("--boot", type=int, default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    import glob
    if not glob.glob(str(REPO / "data" / "raw" / "*" / "amrfinder_runs" / "*" / "main.tsv")):
        print("ERROR: no cached AMRFinder runs (gitignored)", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    res = run_all(REPO, dedup=a.dedup_profiles, boot=a.boot or dcc.BOOT)
    tag = "_dedup" if a.dedup_profiles else ""
    out = a.out or (REPO / "wiki" / f"coresistance_imputation_result_{today}{tag}.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(res, today), encoding="utf-8")
    print(render_md(res, today))
    print(f"\n[wrote {out} + .md]  verdict={res['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
