"""Determinant co-occurrence / linkage world model + phenotype->genotype inverter — Family C of the
genome world-model plan (2026-07-11).

QUESTION: is there JOINT structure (co-resistance LINKAGE) in a genome's AMR determinants BEYOND base
rates and BEYOND species — i.e. do a genome's OTHER determinants predict a held-out determinant, within
organism, better than chance? The deployed decoder scores each determinant independently; this models the
"what-goes-with-what" joint distribution the decoder throws away, and inverts it to the phenotype->genotype
("vice-versa") direction the user asked for.

SUBSTRATE (data in hand — NO new Docker/network): the CACHED AMRFinder `main.tsv` determinant calls across
every `data/raw/*/amrfinder_runs/` cohort (845 distinct genomes; self-distillation from our own decoder's
determinant caller). Element symbol = feature token; Class = drug it confers; Type=AMR only for the linkage
core. Organism = the cohort's organism (the de-confound axis).

DE-CONFOUND (the project's standing discipline — beat WITHIN organism, never conflate species/lineage):
  * Linkage is tested WITHIN each organism with >= MIN_GENOMES genomes (species held constant), so any
    predictive power is genuine within-species co-resistance linkage, not "Klebsiella != Acinetobacter".
  * Clonality proxy: a `--dedup-profiles` pass keeps ONE genome per unique determinant-profile (a crude
    clonal-expansion control — identical determinant sets are likely the same clone re-sampled). Reported
    alongside the raw number so clonal inflation is DISCLOSED, not hidden. Full Mash-clonality correction
    is the follow-on (needs the assemblies + Docker).

PRE-REGISTERED DESIGN (derived; plan `plans/Genome_World_Model_Creative_Data_Reuse_Plan_2026-07-11.md`):
  * Testable determinant (within organism) = both-class support: present in >= MIN_SUPPORT AND absent in
    >= MIN_SUPPORT genomes (=10; safe for 5-fold stratification).
  * C2 linkage: for each testable determinant, 5-fold stratified OOF logistic predicting its presence from
    ALL OTHER determinants (within-organism). Held-out AUC; paired bootstrap (B=BOOT) CI. LINKED = AUC 95%
    CI lower bound > 0.5.
  * Verdict PASS_LINKAGE_STRUCTURE iff >= PASS_FRACTION (=0.5) of testable determinants are LINKED (there IS
    joint co-resistance structure beyond base rates, within organism); else FAIL_CONDITIONALLY_INDEPENDENT.
  * C3a "vice-versa" LIFT table: per target determinant, top co-occurring determinants by lift =
    P(d2|d1)/P(d2), within-organism, co-occurring in >= MIN_COOC genomes.
  * C3b phenotype->genotype inversion: per AMRFinder drug Class, the determinants ranked by within-organism
    frequency (given a target resistance, which determinant sets produce it).
  * verify-in-batch: top-lift pairs inspected for known cassettes (sul1+aadA/aac, blaCTX-M co-members, QRDR
    gyrA+parC).

Honest caveats: the cohorts are DRUG-R/S-SELECTED (resistance-enriched), so co-occurrence reflects the
curated cohorts + selection, NOT a random population sample; within-organism de-confounds species but the
dedup is only a crude clonality proxy. A PASS here (linkage is real, expected) demonstrates the world model
captures real joint structure; a FAIL would mean within-organism conditional independence. NOT a sequence
embedding. Frozen AMR surface is READ-only / untouched.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict

REPO = Path(__file__).resolve().parent.parent

MIN_GENOMES = 60
MIN_SUPPORT = 10          # determinant both-class floor (>=10 present & >=10 absent; safe for 5-fold)
CV_FOLDS = 5
BOOT = 1000
LINKED_AUC = 0.5
PASS_FRACTION = 0.5
MIN_COOC = 8              # lift pair floor
SEED = 0
_DRUGS = ("ciprofloxacin", "ceftriaxone", "gentamicin", "tetracycline", "meropenem", "oxacillin")


def _norm(p: str) -> str:
    return p.replace(os.sep, "/")


def organism_of(main_tsv_norm: str) -> str:
    b = main_tsv_norm.split("/amrfinder_runs/")[0].split("/data/raw/")[-1]
    for drug in _DRUGS:
        b = re.sub(rf"_{drug}(_indep)?$", "", b)
    return re.sub(r"_(provdisjoint|indep|cipro|xsource)$", "", b)


def harvest(repo: Path = REPO):
    """Return (acc_org: acc->organism, acc_dets: acc->set(AMR Element symbol), det_class: symbol->drug Class)."""
    acc_org: dict[str, str] = {}
    acc_dets: dict[str, set[str]] = defaultdict(set)
    det_class: dict[str, Counter] = defaultdict(Counter)
    for mt in glob.glob(str(repo / "data" / "raw" / "*" / "amrfinder_runs" / "*" / "main.tsv")):
        m = _norm(mt)
        org = organism_of(m)
        acc = m.split("/amrfinder_runs/")[1].split("/main.tsv")[0]
        acc_org.setdefault(acc, org)
        try:
            with open(mt, encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if (row.get("Type") or "").upper() != "AMR":
                        continue
                    sym = (row.get("Element symbol") or "").strip()
                    if not sym:
                        continue
                    acc_dets[acc].add(sym)
                    cls = (row.get("Class") or "").strip().upper()
                    if cls:
                        det_class[sym][cls] += 1
        except (OSError, csv.Error):
            continue
    det_class_1 = {s: c.most_common(1)[0][0] for s, c in det_class.items() if c}
    return acc_org, acc_dets, det_class_1


def matrix_for_organism(acc_org, acc_dets, org, min_support=MIN_SUPPORT):
    accs = sorted(a for a, o in acc_org.items() if o == org)
    counts = Counter()
    for a in accs:
        counts.update(acc_dets[a])
    dets = sorted(d for d, n in counts.items() if n >= 1)
    didx = {d: j for j, d in enumerate(dets)}
    X = np.zeros((len(accs), len(dets)), float)
    for i, a in enumerate(accs):
        for d in acc_dets[a]:
            X[i, didx[d]] = 1.0
    testable = [d for d in dets if min_support <= int(X[:, didx[d]].sum()) <= len(accs) - min_support]
    return X, dets, didx, accs, testable


def dedup_profiles(X, accs):
    """Keep one row per unique determinant-profile (crude clonality proxy). Returns kept row indices."""
    seen: dict[tuple, int] = {}
    keep = []
    for i in range(X.shape[0]):
        key = tuple(np.flatnonzero(X[i]).tolist())
        if key not in seen:
            seen[key] = i
            keep.append(i)
    return keep


def _auc(scores, labels):
    scores = np.asarray(scores, float); labels = np.asarray(labels, bool)
    pos, neg = scores[labels], scores[~labels]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    order = np.argsort(allv, kind="mergesort"); sv = allv[order]
    ranks = np.empty(len(allv), float); i = 0
    while i < len(allv):
        j = i
        while j < len(allv) and sv[j] == sv[i]:
            j += 1
        ranks[order[i:j]] = (i + j - 1) / 2.0 + 1.0
        i = j
    u = ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2.0
    return u / (len(pos) * len(neg))


def linkage_for_determinant(X, target_j, seed=SEED, boot=BOOT):
    """OOF logistic predicting determinant target_j from ALL OTHERS; held-out AUC + bootstrap CI."""
    y = (X[:, target_j] > 0).astype(int)
    others = np.delete(X, target_j, axis=1)
    n_min = int(min(y.sum(), len(y) - y.sum()))
    folds = min(CV_FOLDS, n_min)
    if folds < 2:
        return None
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    clf = LogisticRegression(max_iter=2000, C=1.0, solver="liblinear")
    oof = cross_val_predict(clf, others, y, cv=skf, method="predict_proba")[:, 1]
    auc = _auc(oof, y.astype(bool))
    rng = np.random.default_rng(seed)
    yb = y.astype(bool); n = len(y)
    aucs = np.empty(boot, float)
    for b in range(boot):
        idx = rng.integers(0, n, n)
        aucs[b] = _auc(oof[idx], yb[idx])
    aucs = aucs[~np.isnan(aucs)]
    lo, hi = (np.percentile(aucs, [2.5, 97.5]) if len(aucs) else (float("nan"), float("nan")))
    return {"auc": round(float(auc), 4), "ci_lo": round(float(lo), 4), "ci_hi": round(float(hi), 4),
            "n_present": int(y.sum()), "linked": bool(lo > LINKED_AUC)}


def lift_table(X, dets, didx, targets, min_cooc=MIN_COOC, top=6):
    n = X.shape[0]
    out = {}
    for t in targets:
        if t not in didx:
            continue
        ti = didx[t]
        pt = X[:, ti] > 0
        rows = []
        for d in dets:
            if d == t:
                continue
            dj = didx[d]
            pd = X[:, dj] > 0
            cooc = int(np.sum(pt & pd))
            if cooc < min_cooc:
                continue
            p_d = pd.mean()
            p_d_given_t = cooc / max(int(pt.sum()), 1)
            lift = p_d_given_t / p_d if p_d > 0 else float("inf")
            rows.append({"det": d, "cooc": cooc, "lift": round(float(lift), 2),
                         "p_given_target": round(float(p_d_given_t), 3)})
        rows.sort(key=lambda r: -r["lift"])
        out[t] = rows[:top]
    return out


def drug_inversion(acc_org, acc_dets, det_class, org, top=8):
    """Phenotype->genotype: per AMRFinder drug Class, determinants ranked by within-organism frequency."""
    accs = [a for a, o in acc_org.items() if o == org]
    by_class: dict[str, Counter] = defaultdict(Counter)
    for a in accs:
        for d in acc_dets[a]:
            cls = det_class.get(d)
            if cls:
                by_class[cls][d] += 1
    return {cls: [{"det": d, "n": n} for d, n in c.most_common(top)] for cls, c in sorted(by_class.items())}


def run_organism(acc_org, acc_dets, det_class, org, dedup=False, seed=SEED, boot=BOOT):
    X, dets, didx, accs, testable = matrix_for_organism(acc_org, acc_dets, org)
    note = None
    if dedup:
        keep = dedup_profiles(X, accs)
        X = X[keep]
        note = f"profile-deduped {len(accs)}->{len(keep)} genomes (clonality proxy)"
        # recompute testable on deduped matrix
        testable = [d for d in dets if MIN_SUPPORT <= int(X[:, didx[d]].sum()) <= X.shape[0] - MIN_SUPPORT]
    per_det = {}
    for d in testable:
        r = linkage_for_determinant(X, didx[d], seed=seed, boot=boot)
        if r:
            per_det[d] = r
    n_test = len(per_det)
    n_linked = sum(1 for r in per_det.values() if r["linked"])
    targets = [d for d, _ in Counter({d: int(X[:, didx[d]].sum()) for d in testable}).most_common(6)]
    return {"organism": org, "n_genomes": X.shape[0], "n_determinants": len(dets),
            "n_testable": n_test, "n_linked": n_linked,
            "fraction_linked": round(n_linked / n_test, 3) if n_test else 0.0,
            "note": note, "per_determinant": per_det,
            "lift_table": lift_table(X, dets, didx, targets),
            "phenotype_to_genotype_inversion": drug_inversion(acc_org, acc_dets, det_class, org)}


def run_all(repo: Path = REPO, dedup=False, seed=SEED, boot=BOOT):
    acc_org, acc_dets, det_class = harvest(repo)
    oc = Counter(acc_org.values())
    orgs = [o for o, n in oc.most_common() if n >= MIN_GENOMES]
    per_org = {o: run_organism(acc_org, acc_dets, det_class, o, dedup=dedup, seed=seed, boot=boot) for o in orgs}
    tot_test = sum(r["n_testable"] for r in per_org.values())
    tot_linked = sum(r["n_linked"] for r in per_org.values())
    frac = (tot_linked / tot_test) if tot_test else 0.0
    verdict = ("PASS_LINKAGE_STRUCTURE" if (tot_test > 0 and frac >= PASS_FRACTION)
               else ("FAIL_CONDITIONALLY_INDEPENDENT" if tot_test > 0 else "NO_TESTABLE_DETERMINANTS"))
    return {
        "artifact": "determinant_cooccurrence_world_model",
        "schema": "determinant-cooccurrence-v1",
        "question": "Is there within-organism co-resistance LINKAGE (do a genome's other determinants predict "
                    "a held-out determinant beyond base rates)? + the phenotype->genotype inversion.",
        "substrate": "cached AMRFinder AMR-determinant calls across data/raw/*/amrfinder_runs (self-distillation)",
        "prereg": {"MIN_GENOMES": MIN_GENOMES, "MIN_SUPPORT": MIN_SUPPORT, "CV_FOLDS": CV_FOLDS, "BOOT": boot,
                   "LINKED_AUC": LINKED_AUC, "PASS_FRACTION": PASS_FRACTION, "MIN_COOC": MIN_COOC, "seed": seed,
                   "deconfound": "within-organism (species held constant) + optional profile-dedup clonality proxy",
                   "linked_rule": "held-out OOF-AUC 95% CI lower bound > 0.5 (5-fold stratified logistic on other determinants)"},
        "deduped": dedup,
        "n_genomes_total": len(acc_org),
        "verdict": verdict,
        "total_testable": tot_test, "total_linked": tot_linked, "fraction_linked": round(frac, 3),
        "honest_caveats": [
            "Cohorts are DRUG-R/S-SELECTED (resistance-enriched) -> co-occurrence reflects the curated cohorts + "
            "selection, NOT a random population sample. Linkage is 'within these cohorts'.",
            "Within-organism de-confounds SPECIES; --dedup-profiles is only a CRUDE clonality proxy (identical "
            "determinant sets). Full Mash-clonality correction is the follow-on (needs the assemblies + Docker).",
            "AMRFinder point-mutation determinants (gyrA_S83L) are chromosomal/organism-specific; acquired genes "
            "(sul1, tet(A)) are mobile/plasmid. Linkage mixes both mechanisms.",
            "A PASS (linkage exists) is the expected + valid finding — it demonstrates the world model captures "
            "real joint structure; a FAIL would mean within-organism conditional independence.",
            "Self-distillation from our own AMRFinder caller: associational, NOT a causal claim.",
        ],
        "per_organism": per_org,
    }


def render_md(res, generated):
    L = [f"# Determinant co-occurrence / linkage world model + phenotype->genotype inverter ({generated})", "",
         f"**Verdict: {res['verdict']}** — {res['total_linked']}/{res['total_testable']} within-organism "
         f"testable determinants are LINKED (fraction {res['fraction_linked']}; PASS bar "
         f"{res['prereg']['PASS_FRACTION']}). {'(profile-deduped clonality proxy)' if res['deduped'] else '(raw)'}", "",
         f"{res['question']}", "",
         "Substrate: cached AMRFinder determinant calls (self-distillation). LINKED = a determinant's held-out "
         "presence is predicted from the OTHER determinants, WITHIN organism, with OOF-AUC 95% CI lower > 0.5.", "",
         "| organism | genomes | testable | **linked** | frac | note |",
         "|---|---|---|---|---|---|"]
    for r in res["per_organism"].values():
        L.append(f"| {r['organism']} | {r['n_genomes']} | {r['n_testable']} | **{r['n_linked']}** | "
                 f"{r['fraction_linked']} | {r.get('note') or ''} |")
    L += ["", "## C3a — 'vice-versa' co-occurrence LIFT (what determinants travel together; verify-in-batch)"]
    for r in res["per_organism"].values():
        L.append(f"### {r['organism']}")
        for t, rows in r["lift_table"].items():
            if rows:
                terms = ", ".join(f"{x['det']}(lift {x['lift']}, n{x['cooc']})" for x in rows[:4])
                L.append(f"- **{t}** → {terms}")
    L += ["", "## C3b — phenotype→genotype inversion (given a drug Class, the determinants that confer it)"]
    for r in res["per_organism"].values():
        inv = r["phenotype_to_genotype_inversion"]
        L.append(f"### {r['organism']}")
        for cls, rows in list(inv.items())[:8]:
            terms = ", ".join(f"{x['det']}({x['n']})" for x in rows[:5])
            L.append(f"- **{cls}** ← {terms}")
    L += ["", "## Honest caveats"] + [f"- {c}" for c in res["honest_caveats"]]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dedup-profiles", action="store_true", help="clonality proxy: 1 genome per determinant-profile")
    ap.add_argument("--boot", type=int, default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    boot = a.boot or BOOT
    if not glob.glob(str(REPO / "data" / "raw" / "*" / "amrfinder_runs" / "*" / "main.tsv")):
        print("ERROR: no cached AMRFinder runs under data/raw/*/amrfinder_runs (gitignored)", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    res = run_all(REPO, dedup=a.dedup_profiles, boot=boot)
    tag = "_dedup" if a.dedup_profiles else ""
    out = a.out or (REPO / "wiki" / f"determinant_cooccurrence_result_{today}{tag}.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(res, today), encoding="utf-8")
    print(render_md(res, today))
    print(f"\n[wrote {out} + .md]  verdict={res['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
