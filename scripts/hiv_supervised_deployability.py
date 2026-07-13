"""Deployability stress-test of the supervised HIV blind-spot rescue (2026-07-12).

`hiv_supervised_vs_catalog.py` found a supervised full-sequence model scores 0.889 AUROC on the NNRTI
catalog's blind spot (vs ESM zero-shot 0.449) — BUT that was plain 5-fold CV on the deployed dataset,
which (a) does not group by patient (the .Full metadata shows 1.35x patient duplication incl. ONE patient
with 508 isolates -> likely leakage inflation) and (b) is all in-distribution. This runs the honest
out-of-distribution / de-leaked tests using the richer `NNRTI_DataSet.Full.txt` (Subtype / PtID / RefID):

  (A) PATIENT-GROUPED CV (GroupKFold by PtID) — removes patient leakage; still same subtypes/studies.
  (B) LEAVE-ONE-STUDY-OUT (GroupKFold by RefID) — train on some studies, predict held-out studies =
      "predict a NEW lab's isolates". The deployment-relevant test.
  (C) CROSS-SUBTYPE (train on subtype B, predict non-B) — genuine OOD (underpowered on the blind spot;
      reported honestly).

For each split we report the OVERALL AUROC and the BLIND-SPOT AUROC (catalog-negative test isolates),
vs the mutation-burden baseline + shuffled null. PASS bar (same as the prior): blind-spot AUROC >= 0.65
AND > burden AND null < 0.55.

VERDICT:
  DEPLOYABLE_HOLDS_OOD        — leave-study-out blind-spot PASSES (survives predicting new studies).
  DE_LEAKED_HOLDS_NOT_OOD     — patient-grouped PASSES but leave-study-out fails (real but in-distribution).
  IN_DISTRIBUTION_INFLATED    — collapses once patient leakage is removed (the 0.889 was leakage).

Frozen decoder surface untouched (research-only). Reuses the prior script's one-hot builder + catalog.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import random
import statistics as st
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_val_predict

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

spec = importlib.util.spec_from_file_location(
    "hiv_supervised_vs_catalog", REPO / "scripts" / "hiv_supervised_vs_catalog.py")
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
ev = S.ev
import dna_decode.data.hiv_amr as H  # noqa: E402

FULL = REPO / "data" / "raw" / "hiv" / "NNRTI_DataSet.Full.txt"
DRUG, CUTOFF, SEED, C = "EFV", 3.0, 0, 1.0
majors = H.NNRTI_RT_MAJOR_DRMS


def load_full():
    rows = list(csv.DictReader(open(FULL, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    have = [r for r in rows if r.get(DRUG) not in ("NA", "", "-", None)]
    return have, pcols


def _fit_predict_grouped(X, y, groups):
    n_splits = min(5, len(set(groups)))
    return cross_val_predict(LogisticRegression(max_iter=2000, C=C, solver="liblinear"), X, y,
                             groups=groups, cv=GroupKFold(n_splits=n_splits), method="predict_proba")[:, 1]


def _blindspot_metrics(y, scores, neg_mask, burden, seed=SEED):
    ys = [y[i] for i in range(len(y)) if neg_mask[i]]
    ss = [scores[i] for i in range(len(y)) if neg_mask[i]]
    bs = [burden[i] for i in range(len(y)) if neg_mask[i]]
    if len(set(ys)) < 2:
        return {"n": len(ys), "R": int(sum(ys)), "auroc": None, "burden": None, "null": None,
                "pass": False, "note": "degenerate (single-class blind-spot test set)"}
    a = ev.auroc(ys, ss)
    ab = ev.auroc(ys, bs)
    rng = random.Random(seed)
    nulls = []
    for _ in range(200):
        yy = ys[:]
        rng.shuffle(yy)
        nulls.append(ev.auroc(yy, ss))
    an = st.median(nulls)
    return {"n": len(ys), "R": int(sum(ys)), "auroc": round(a, 4), "burden": round(ab, 4),
            "null": round(an, 4), "pass": bool(a >= 0.65 and a > ab and an < 0.55)}


def main():
    prot = ev.rt_protein()
    have, pcols = load_full()
    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in have if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(have):
            drifted.add(p)

    def catalog_call(r):
        return any(f"{H._RT_WT.get(p,'?')}{p}{aa}" in majors for p, aa in ev.isolate_muts(r, pcols))

    X, feats = S.build_onehot(have, pcols, drifted, prot)
    y = np.array([1 if float(r[DRUG]) >= CUTOFF else 0 for r in have])
    cat = np.array([1 if catalog_call(r) else 0 for r in have])
    neg = (cat == 0)
    burden = [len(ev.isolate_muts(r, pcols)) for r in have]
    ptid = np.array([r["PtID"] for r in have])
    refid = np.array([r["RefID"] for r in have])
    subtype = np.array([r["Subtype"] for r in have])
    print(f"cohort n={len(have)}  R={int(y.sum())}  features={len(feats)}  "
          f"patients={len(set(ptid))}  studies={len(set(refid))}  non-B={int((subtype!='B').sum())}")
    print(f"catalog-negative overall: n={int(neg.sum())} R={int(y[neg].sum())}")

    results = {}

    # (A) patient-grouped CV
    oof_pt = _fit_predict_grouped(X, y, ptid)
    results["A_patient_grouped"] = {
        "overall_auroc": round(ev.auroc(y.tolist(), oof_pt.tolist()), 4),
        "blind_spot": _blindspot_metrics(y.tolist(), oof_pt.tolist(), neg.tolist(), burden)}

    # (B) leave-one-study-out CV
    oof_ref = _fit_predict_grouped(X, y, refid)
    results["B_leave_study_out"] = {
        "overall_auroc": round(ev.auroc(y.tolist(), oof_ref.tolist()), 4),
        "blind_spot": _blindspot_metrics(y.tolist(), oof_ref.tolist(), neg.tolist(), burden)}

    # (C) cross-subtype: train B, predict non-B
    B = subtype == "B"
    clf = LogisticRegression(max_iter=2000, C=C, solver="liblinear").fit(X[B], y[B])
    pnb = clf.predict_proba(X[~B])[:, 1]
    ynb, negnb = y[~B].tolist(), neg[~B].tolist()
    burden_nb = [burden[i] for i in range(len(have)) if not B[i]]
    results["C_cross_subtype_trainB_testNonB"] = {
        "n_test_nonB": int((~B).sum()),
        "overall_auroc": round(ev.auroc(ynb, pnb.tolist()), 4) if len(set(ynb)) > 1 else None,
        "blind_spot": _blindspot_metrics(ynb, pnb.tolist(), negnb, burden_nb)}

    a_leak = results["A_patient_grouped"]["blind_spot"]["auroc"]
    b_ood = results["B_leave_study_out"]["blind_spot"]["auroc"]
    b_pass = results["B_leave_study_out"]["blind_spot"]["pass"]
    a_pass = results["A_patient_grouped"]["blind_spot"]["pass"]
    verdict = ("DEPLOYABLE_HOLDS_OOD" if b_pass
               else ("DE_LEAKED_HOLDS_NOT_OOD" if a_pass else "IN_DISTRIBUTION_INFLATED"))

    res = {
        "artifact": "hiv_supervised_deployability", "schema": "hiv-supervised-deployability-v1",
        "date": str(_date.today()), "drug": DRUG, "cutoff_fold": CUTOFF, "seed": SEED,
        "in_distribution_prior": {"blind_spot_auroc": 0.889, "note": "plain 5-fold CV, no patient grouping"},
        "cohort": {"n": len(have), "R": int(y.sum()), "patients": len(set(ptid)), "studies": len(set(refid)),
                   "non_B": int((subtype != "B").sum()), "catalog_negative_n": int(neg.sum())},
        "results": results,
        "verdict": verdict,
        "honest_note": (
            "Compares the 0.889 in-distribution blind-spot AUROC against (A) patient-grouped de-leaking, "
            "(B) leave-one-study-out OOD (the deployment test), (C) cross-subtype. The blind-spot subsets in "
            "(C) are small (non-B=~291 total) -> often single-class/underpowered, reported as such. The "
            "load-bearing numbers are (A) [is 0.889 patient-leakage?] and (B) [does it survive predicting a "
            "new study?]."),
    }
    for k, v in results.items():
        bs = v["blind_spot"]
        print(f"\n{k}: overall AUROC {v.get('overall_auroc')}")
        print(f"   blind-spot: n={bs['n']} R={bs['R']}  AUROC={bs['auroc']}  "
              f"burden={bs.get('burden')}  null={bs.get('null')}  PASS={bs['pass']}"
              + (f"  [{bs['note']}]" if bs.get('note') else ""))
    print(f"\nin-distribution prior blind-spot AUROC = 0.889 (no patient grouping)")
    print(f"VERDICT: {verdict}")

    out = REPO / "wiki" / f"hiv_supervised_deployability_{_date.today()}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[wrote {out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
