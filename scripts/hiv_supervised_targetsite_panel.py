"""Extend the supervised blind-spot complement to the PI (protease) + INSTI (integrase) genes (2026-07-12).

Same method as the RT panel (hiv_supervised_panel.py), DIFFERENT proteins: a supervised full-sequence
one-hot model over the protease / integrase per-position residues vs the position-based Stanford HIVDB
catalog, on the free PhenoSense fold-change label, de-confounded by patient-grouped + leave-one-study-out.

Per drug: overall AUROC + BLIND-SPOT AUROC (isolates the catalog calls susceptible — no major-position
mutation) under BOTH splits. PASS = blind-spot AUROC >= 0.65 AND > mutation-burden AND null < 0.55, under
LEAVE-STUDY-OUT (the deployment split).

VERDICT:
  GENERALIZES_TO_TARGETSITE_GENES  — deployment-split rescue PASSES on a majority of powered PI+INSTI drugs.
  RT_SPECIFIC                      — passes on <= 1 (the RT win didn't transfer to other HIV genes).
  MIXED                            — some but not a clear majority.

Frozen decoder surface + hiv_amr catalog untouched (READ-only). Reuses the one-hot builder + deployability
splits from the supervised-vs-catalog scripts + the committed PR/IN HXB2 references.
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
S, ev = D.S, D.ev
import dna_decode.data.hiv_amr as H  # noqa: E402

RAW = REPO / "data" / "raw" / "hiv"
REF = REPO / "data" / "hiv_ref"
CUTOFF, SEED = 3.0, 0
# (class, dataset, gene-CDS ref, drug reuse-columns, major positions)
GENES = {
    "PI": {"dataset": "PI_DataSet.Full.txt", "ref": "HIV1_PR_HXB2_cds.fna",
           "cols": ("FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"),
           "positions": set(H.PI_CLASS.positions)},
    "INSTI": {"dataset": "INI_DataSet.Full.txt", "ref": "HIV1_IN_HXB2_cds.fna",
              "cols": ("RAL", "EVG", "DTG", "BIC", "CAB"),
              "positions": set(H.INSTI_CLASS.positions)},
}


def translate(ref_path):
    seq = "".join(l.strip() for l in open(ref_path) if not l.startswith(">"))
    return "".join(ev.CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def prep(dataset, prot):
    rows = list(csv.DictReader(open(RAW / dataset, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in rows if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(rows):
            drifted.add(p)
    return rows, pcols, drifted


def run_drug(col, prot, rows, pcols, drifted, positions):
    have = [r for r in rows if r.get(col) not in ("NA", "", "-", None)]
    if len(have) < 150:
        return {"drug": col, "powered": False, "note": f"n={len(have)}<150"}
    X, feats = S.build_onehot(have, pcols, drifted, prot)
    y = np.array([1 if float(r[col]) >= CUTOFF else 0 for r in have])
    if min(int(y.sum()), len(y) - int(y.sum())) < 15:
        return {"drug": col, "powered": False, "note": f"minority class {min(int(y.sum()), len(y)-int(y.sum()))}<15"}
    cat = np.array([1 if any(p in positions for p, _ in ev.isolate_muts(r, pcols)) else 0 for r in have])
    neg = (cat == 0)
    burden = [len(ev.isolate_muts(r, pcols)) for r in have]
    ptid = np.array([r.get("PtID", str(i)) for i, r in enumerate(have)])
    refid = np.array([r.get("RefID", "NA") for r in have])
    out = {"drug": col, "powered": True, "n": len(have), "R": int(y.sum()),
           "catalog_negative_n": int(neg.sum()), "catalog_negative_R": int(y[neg].sum())}
    for tag, groups in (("patient_grouped", ptid), ("leave_study_out", refid)):
        oof = D._fit_predict_grouped(X, y, groups)
        out[tag] = {"overall_auroc": round(ev.auroc(y.tolist(), oof.tolist()), 4),
                    "blind_spot": D._blindspot_metrics(y.tolist(), oof.tolist(), neg.tolist(), burden)}
    return out


def main():
    results = []
    for cls, spec in GENES.items():
        prot = translate(REF / spec["ref"])
        rows, pcols, drifted = prep(spec["dataset"], prot)
        print(f"[targetsite] {cls}: gene {len(prot)} aa; {len(rows)} isolates", flush=True)
        for col in spec["cols"]:
            print(f"[targetsite] {cls} {col} ...", flush=True)
            results.append({**run_drug(col, prot, rows, pcols, drifted, spec["positions"]), "class": cls})

    powered = [r for r in results if r.get("powered")]
    testable = [r for r in powered if r["leave_study_out"]["blind_spot"].get("auroc") is not None]
    n_pass = sum(1 for r in testable if r["leave_study_out"]["blind_spot"].get("pass"))
    verdict = ("GENERALIZES_TO_TARGETSITE_GENES" if testable and n_pass >= max(2, (len(testable) + 1) // 2)
               else ("RT_SPECIFIC" if n_pass <= 1 else "MIXED"))
    res = {"artifact": "hiv_supervised_targetsite_panel", "schema": "hiv-supervised-targetsite-panel-v1",
           "date": str(_date.today()),
           "question": "Does the supervised blind-spot rescue extend from RT to the PI (protease) + INSTI "
                       "(integrase) genes (deployment split = leave-one-study-out)?",
           "cutoff": CUTOFF, "pass_bar": "leave-study-out blind-spot >= 0.65 AND > burden AND null < 0.55",
           "n_drugs": len(results), "n_powered": len(powered), "n_testable": len(testable),
           "n_pass_deployment_split": n_pass, "verdict": verdict,
           "honest_caveats": [
               "PI + INSTI catalogs are POSITION-based -> blind spot = no major-position mutation (like NRTI). "
               "Cutoff = illustrative fold>=3 (no per-drug clinical PI/INSTI cutoff in-repo); headline is the "
               "supervised-vs-catalog+null DELTA.",
               "2nd-gen INSTIs (DTG/BIC/CAB) have few resistant isolates -> often under-powered, reported honestly.",
               "Supervised needs the free label to train; in-distribution to the Stanford knowledge base; "
               "leave-study-out is the deployment-relevant generalization test.",
           ],
           "per_drug": results}
    out = REPO / "wiki" / f"hiv_supervised_targetsite_panel_{_date.today()}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("\n=== PI + INSTI supervised panel — leave-study-out (deployment) blind-spot ===")
    print(f"{'drug':>5} {'class':>6} {'n':>5} {'cat-neg':>8} {'blind_auroc':>12} {'burden':>7} {'PASS':>5}")
    for r in results:
        if not r.get("powered"):
            print(f"{r['drug']:>5} {r.get('class',''):>6}  {r.get('note','')}"); continue
        bs = r["leave_study_out"]["blind_spot"]
        print(f"{r['drug']:>5} {r['class']:>6} {r['n']:>5} {r['catalog_negative_n']:>8} "
              f"{str(bs.get('auroc')):>12} {str(bs.get('burden')):>7} {str(bs.get('pass')):>5}")
    print(f"\nVERDICT: {verdict}  ({n_pass}/{len(testable)} powered+testable pass the deployment split)")
    print(f"[wrote {out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
