"""Validate the S. pneumoniae β-lactam PBP→MIC engine vs WET-LAB measured AST (GPS Poland cohort).

Pipeline per isolate: GPS PBP types (pbp1a/pbp2b/pbp2x, MOESM4) -> APT -> CDC Ref_PBPtype_MIC lookup -> MIC
-> R/S via pneumo_breakpoints -> vs measured MIC-derived R/S (MOESM3). Reports per (drug, breakpoint-context).

HONESTY: the engine (PBP-type->MIC lookup) is CDC's deterministic method (faithful-to-tool); the PBP TYPES
are GPS's (our-PBP-typer swap deferred). Label = wet-lab measured AST (independent). No-call when the PBP
type is novel ('NEW') or absent from the CDC table.

    uv run python scripts/pneumo_betalactam_validate.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from dna_decode.data.pneumo_breakpoints import classify                            # noqa: E402
from dna_decode.organism_rules.pneumo_betalactam import load_pbp_mic_table, predict_rs  # noqa: E402

GPS = Path("D:/dna_decode_cache/pneumo_gps")
DB = REPO / "data" / "pneumo_betalactam_db" / "Ref_PBPtype_MIC.csv"
# measured MIC column (MOESM3) per drug + breakpoint contexts to score.
MEAS_COL = {"penicillin": "Penicillin", "meropenem": "Meropenem", "ceftriaxone": "Ceftriaxone",
            "cefotaxime": "Cefotaxime"}
CONTEXTS = {"penicillin": ("meningitis", "non_meningitis"), "meropenem": ("meningitis",),
            "ceftriaxone": ("meningitis", "non_meningitis"), "cefotaxime": ("meningitis", "non_meningitis")}


def _mic(v):
    v = (v or "").strip().replace("<=", "").replace(">=", "").replace("<", "").replace(">", "")
    try:
        return float(v)
    except ValueError:
        return None


def main(argv=None) -> int:
    if not DB.exists():
        print(f"ERROR: CDC PBP-MIC table not found at {DB} (copy from a GlobalPneumoSeq/spn-pbp-amr clone)")
        return 2
    table = load_pbp_mic_table(DB)
    m3 = {r["ERR"].strip(): r for r in csv.DictReader(open(GPS / "sd_3.csv", encoding="utf-8")) if r.get("ERR")}
    m4 = {r["Sample_ID"].strip(): r for r in csv.DictReader(open(GPS / "sd_4.csv", encoding="utf-8")) if r.get("Sample_ID")}
    both = set(m3) & set(m4)

    out = {"schema": "pneumo-betalactam-pbp-mic-validation-v1",
           "engine": "CDC Ref_PBPtype_MIC lookup (GlobalPneumoSeq/spn-pbp-amr) -> pneumo_breakpoints",
           "genotype_source": "GPS pipeline PBP types (our-PBP-typer-from-genome swap deferred)",
           "label": "wet-lab measured AST (GPS Poland, MIC -> R/S via CLSI pneumo breakpoints)",
           "n_isolates": len(both), "cells": {}}
    for drug, contexts in CONTEXTS.items():
        for ctx in contexts:
            tp = tn = fp = fn = nocall = 0
            for k in both:
                meas = classify(drug, ctx, _mic(m3[k].get(MEAS_COL[drug])))
                if meas not in ("R", "S"):
                    continue
                r = predict_rs(table, m4[k].get("pbp1a"), m4[k].get("pbp2b"), m4[k].get("pbp2x"), drug, ctx)
                pred = r["prediction"]
                pred = "R" if pred in ("R", "I") else "S" if pred == "S" else None
                if pred is None:
                    nocall += 1
                    continue
                if meas == "R" and pred == "R": tp += 1
                elif meas == "S" and pred == "S": tn += 1
                elif meas == "S" and pred == "R": fp += 1
                else: fn += 1
            n = tp + tn + fp + fn
            out["cells"][f"{drug}@{ctx}"] = {
                "n": n, "nocall": nocall, "accuracy": round((tp + tn) / n, 3) if n else None,
                "sensitivity": round(tp / (tp + fn), 3) if (tp + fn) else None,
                "specificity": round(tn / (tn + fp), 3) if (tn + fp) else None,
                "tp": tp, "fp": fp, "tn": tn, "fn": fn}
            print(f"{drug}@{ctx}: n={n} nocall={nocall} acc={out['cells'][f'{drug}@{ctx}']['accuracy']} "
                  f"sens={out['cells'][f'{drug}@{ctx}']['sensitivity']} spec={out['cells'][f'{drug}@{ctx}']['specificity']} "
                  f"TP{tp} FP{fp} TN{tn} FN{fn}")
    (REPO / "wiki" / "pneumo_betalactam_pbp_validation.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("-> wiki/pneumo_betalactam_pbp_validation.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
