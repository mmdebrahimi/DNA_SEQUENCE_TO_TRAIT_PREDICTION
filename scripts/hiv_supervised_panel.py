"""Is the supervised blind-spot rescue GENERAL across the HIV RT drug panel, or an EFV quirk? (2026-07-12)

Recommendation 1: run the supervised-full-sequence-vs-catalog blind-spot test + the deployability splits
(patient-grouped de-leaking + leave-one-study-out OOD) across ALL 11 HIV RT drugs — 5 NNRTIs (mutant-level
catalog) + 6 NRTIs (position-based catalog) — on the free independent Stanford PhenoSense fold-change label.

Per drug: overall AUROC + BLIND-SPOT AUROC (catalog-negative test isolates) under BOTH the patient-grouped
and the leave-study-out splits, vs mutation-burden + shuffled-null. PASS = blind-spot AUROC >= 0.65 AND
> burden AND null < 0.55, under the LEAVE-STUDY-OUT split (the deployment-relevant one).

VERDICT:
  GENERAL_RESCUE      — the deployment-split blind-spot rescue PASSES on a majority of powered drugs.
  EFV_QUIRK           — it passes on <= 1 drug (EFV only-ish).
  MIXED               — passes on some but not a clear majority.

Reuses the one-hot builder, catalog calls, and deployability splits from the prior two scripts. Per-drug
cutoffs: NRTI from hiv_nrti_validate.NRTI_LOWER_CUTOFF; NNRTI = illustrative fold>=3 (the class uses a
cutoff-free AUC primarily — 3.0 is a sensitivity cutoff, fine for a supervised-vs-catalog DELTA). Frozen
decoder surface untouched (research-only).
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "scripts"))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


D = _load("hiv_supervised_deployability", REPO / "scripts" / "hiv_supervised_deployability.py")
S = D.S            # hiv_supervised_vs_catalog
ev = D.ev          # hiv_esm_vs_catalog
import dna_decode.data.hiv_amr as H  # noqa: E402
import hiv_nrti_validate as NV  # noqa: E402

RAW = REPO / "data" / "raw" / "hiv"
NNRTI = {"EFV": 3.0, "NVP": 3.0, "ETR": 3.0, "RPV": 3.0, "DOR": 3.0}
# NRTI reuse-column -> (drug name, cutoff) from the validated per-drug lower cutoffs
NRTI = {NV._NRTI_COL[d]: NV.NRTI_LOWER_CUTOFF[d] for d in NV.NRTI_DRUGS}   # {'3TC':5.0,'ABC':2.0,...}
SEED = 0


def _prep(dataset):
    prot = ev.rt_protein()
    rows = list(csv.DictReader(open(RAW / dataset, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in rows if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(rows):
            drifted.add(p)
    return prot, rows, pcols, drifted


def _nnrti_call(r, pcols):
    return any(f"{H._RT_WT.get(p,'?')}{p}{aa}" in H.NNRTI_RT_MAJOR_DRMS for p, aa in ev.isolate_muts(r, pcols))


def _nrti_call(r, pcols):
    return any(p in H.NRTI_MAJOR_POSITIONS for p, _ in ev.isolate_muts(r, pcols))


def run_drug(col, cutoff, prot, rows, pcols, drifted, catalog_call):
    have = [r for r in rows if r.get(col) not in ("NA", "", "-", None)]
    if len(have) < 200:
        return {"drug": col, "powered": False, "note": f"n={len(have)}<200"}
    X, feats = S.build_onehot(have, pcols, drifted, prot)
    y = np.array([1 if float(r[col]) >= cutoff else 0 for r in have])
    if len(set(y.tolist())) < 2:
        return {"drug": col, "powered": False, "note": "single-class at cutoff"}
    cat = np.array([1 if catalog_call(r, pcols) else 0 for r in have])
    neg = (cat == 0)
    burden = [len(ev.isolate_muts(r, pcols)) for r in have]
    ptid = np.array([r.get("PtID", str(i)) for i, r in enumerate(have)])
    refid = np.array([r.get("RefID", "NA") for r in have])
    out = {"drug": col, "powered": True, "n": len(have), "R": int(y.sum()), "cutoff": cutoff,
           "catalog_negative_n": int(neg.sum()), "catalog_negative_R": int(y[neg].sum())}
    for tag, groups in (("patient_grouped", ptid), ("leave_study_out", refid)):
        oof = D._fit_predict_grouped(X, y, groups)
        out[tag] = {"overall_auroc": round(ev.auroc(y.tolist(), oof.tolist()), 4),
                    "blind_spot": D._blindspot_metrics(y.tolist(), oof.tolist(), neg.tolist(), burden)}
    return out


def main():
    results = []
    # NNRTI panel
    prot, rows, pcols, drifted = _prep("NNRTI_DataSet.Full.txt")
    for col, cut in NNRTI.items():
        print(f"[panel] NNRTI {col} ...", flush=True)
        results.append({**run_drug(col, cut, prot, rows, pcols, drifted, _nnrti_call), "class": "NNRTI"})
    # NRTI panel
    prot, rows, pcols, drifted = _prep("NRTI_DataSet.Full.txt")
    for col, cut in NRTI.items():
        print(f"[panel] NRTI {col} ...", flush=True)
        results.append({**run_drug(col, cut, prot, rows, pcols, drifted, _nrti_call), "class": "NRTI"})

    powered = [r for r in results if r.get("powered")]
    def _pass(r):  # deployment-split blind-spot pass
        bs = r["leave_study_out"]["blind_spot"]
        return bool(bs.get("pass"))
    n_pass = sum(1 for r in powered if _pass(r))
    testable = [r for r in powered if r["leave_study_out"]["blind_spot"].get("auroc") is not None]
    verdict = ("GENERAL_RESCUE" if testable and n_pass >= max(2, (len(testable) + 1) // 2)
               else ("EFV_QUIRK" if n_pass <= 1 else "MIXED"))
    res = {
        "artifact": "hiv_supervised_panel", "schema": "hiv-supervised-panel-v1", "date": str(_date.today()),
        "question": "Is the supervised blind-spot rescue GENERAL across the HIV RT drug panel (deployment "
                    "split = leave-one-study-out), or an EFV quirk?",
        "pass_bar": "leave-study-out blind-spot AUROC >= 0.65 AND > burden AND null < 0.55",
        "n_drugs": len(results), "n_powered": len(powered), "n_testable_blindspot": len(testable),
        "n_pass_deployment_split": n_pass, "verdict": verdict,
        "honest_caveats": [
            "NNRTI cutoff = illustrative fold>=3 (class primary metric is cutoff-free AUC); NRTI = validated "
            "per-drug lower cutoffs. The headline is the supervised-vs-catalog+null DELTA at a fixed cutoff.",
            "NRTI catalog is POSITION-based (any major-position mutation) -> its 'blind spot' is smaller/"
            "different from the NNRTI mutant-level catalog; a NRTI blind-spot subset can be small/underpowered.",
            "Same scope as the EFV result: supervised needs training labels; in-distribution to the Stanford "
            "knowledge base; the leave-study-out split is the deployment-relevant generalization test.",
        ],
        "per_drug": results,
    }
    out = REPO / "wiki" / f"hiv_supervised_panel_{_date.today()}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("\n=== HIV supervised panel — leave-study-out (deployment) blind-spot ===")
    print(f"{'drug':>5} {'class':>6} {'n':>5} {'cat-neg':>8} {'blind_auroc':>12} {'burden':>7} {'PASS':>5}")
    for r in results:
        if not r.get("powered"):
            print(f"{r['drug']:>5} {r.get('class',''):>6}  {r.get('note','')}")
            continue
        bs = r["leave_study_out"]["blind_spot"]
        print(f"{r['drug']:>5} {r['class']:>6} {r['n']:>5} {r['catalog_negative_n']:>8} "
              f"{str(bs.get('auroc')):>12} {str(bs.get('burden')):>7} {str(bs.get('pass')):>5}")
    print(f"\nVERDICT: {verdict}  ({n_pass}/{len(testable)} powered+testable drugs pass the deployment split)")
    print(f"[wrote {out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
