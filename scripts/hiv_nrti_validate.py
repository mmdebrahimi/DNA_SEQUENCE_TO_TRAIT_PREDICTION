"""Validate the v0 HIV NRTI cell vs HIVDB PhenoSense + the underlying-tool OLS baseline (Wave B, class 2).

Second HIV drug class. Same circularity-safe contract as the NNRTI validation (label = PhenoSense
fold-change, NOT HIVDB Sierra), but the v0 NRTI catalog is POSITION-BASED (Stanford published NRTI majors
as positions only: 41,65,70,74,75,151,184,210,215). So it deliberately OVER-CALLS (T215 revertants / V75
polymorphisms) — the per-drug specificity hit is the honest v0 finding (motivates a mutant-specific v0.1).

Reuses the NNRTI machinery wholesale: load_rows/_parse_fold/_auc_rank (validate) + _build_design_matrix/
_confusion/_auc (baseline OLS reimplementation of DRMcv.R). Per-drug clinical lower cutoffs sourced from
DRMcv.R lines 173-178: 3TC=5, ABC=2, AZT=3, D4T=1.5, DDI=1.5, TDF=1.5.

DATA: data/raw/hiv/NRTI_DataSet.txt (gitignored; cite Rhee 2003).
"""
from __future__ import annotations

import json
import statistics
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_predict

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.hiv_amr import NRTI_MAJOR_POSITIONS, NRTI_RT_WT, call_nrti_from_observed
from scripts.hiv_nnrti_baseline import CV_FOLDS, _auc, _build_design_matrix, _confusion, _position_columns
from scripts.hiv_nnrti_validate import _auc_rank, _parse_fold, load_rows

DEFAULT_DATA = REPO / "data" / "raw" / "hiv" / "NRTI_DataSet.txt"
NRTI_DRUGS = ["lamivudine", "abacavir", "zidovudine", "stavudine", "didanosine", "tenofovir"]
_NRTI_COL = {"lamivudine": "3TC", "abacavir": "ABC", "zidovudine": "AZT",
             "stavudine": "D4T", "didanosine": "DDI", "tenofovir": "TDF"}
# Clinical lower cutoffs (susceptible/non-susceptible) — DRMcv.R confusion-matrix lines 173-178.
NRTI_LOWER_CUTOFF = {"lamivudine": 5.0, "abacavir": 2.0, "zidovudine": 3.0,
                     "stavudine": 1.5, "didanosine": 1.5, "tenofovir": 1.5}


def _observed_nrti(row: dict[str, str]) -> set[str]:
    """<wt><pos><mut> substitutions at the NRTI major positions from the P-columns ('-'=consensus)."""
    out: set[str] = set()
    for pos in NRTI_MAJOR_POSITIONS:
        cell = (row.get(f"P{pos}") or "").strip()
        if cell in ("", "-", ".", "NA"):
            continue
        wt = NRTI_RT_WT[pos]
        for aa in cell:
            if aa.isalpha() and aa != wt:
                out.add(f"{wt}{pos}{aa}")
    return out


def run(path: Path = DEFAULT_DATA) -> dict:
    rows = load_rows(path)
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
    pos_cols = _position_columns(header)
    X_all, feat_names = _build_design_matrix(rows, pos_cols)

    per_drug = {}
    for drug in NRTI_DRUGS:
        col = _NRTI_COL[drug]
        cutoff = NRTI_LOWER_CUTOFF[drug]
        keep, folds, catalog_R, fold_R, fold_S = [], [], [], [], []
        for i, row in enumerate(rows):
            fold = _parse_fold(row.get(col, ""))
            if fold is None or fold <= 0:
                continue
            keep.append(i); folds.append(fold)
            call_R = call_nrti_from_observed(drug, {"RT": _observed_nrti(row)}).prediction == "R"
            catalog_R.append(call_R)
            (fold_R if call_R else fold_S).append(fold)
        if len(keep) < 30:
            per_drug[drug] = {"n": len(keep), "note": "too few isolates"}
            continue
        fold = np.array(folds)
        actual_R = fold >= cutoff
        catalog_R = np.array(catalog_R, dtype=bool)
        # OLS underlying-tool baseline (5-fold CV OOF predicted log10 fold)
        Xd = X_all[np.array(keep)]
        oof_fold = np.power(10.0, cross_val_predict(LinearRegression(), Xd, np.log10(fold), cv=CV_FOLDS))
        ols_R = oof_fold >= cutoff
        cat = _confusion(catalog_R, actual_R)
        ols = {**_confusion(ols_R, actual_R), "auc": _auc(oof_fold, actual_R)}
        per_drug[drug] = {
            "n_isolates": len(keep),
            "prevalence_R_at_cutoff": round(float(actual_R.mean()), 3),
            "clinical_lower_cutoff_fold": cutoff,
            "catalog": {**cat, "auc_call_separates_fold": (
                            round(_auc_rank(fold_R, fold_S), 4) if (fold_R and fold_S) else None),
                        "median_fold_R": round(statistics.median(fold_R), 2) if fold_R else None,
                        "median_fold_S": round(statistics.median(fold_S), 2) if fold_S else None},
            "ols_baseline": ols,
            "delta_ols_minus_catalog_balacc": (
                round(ols["balanced_accuracy"] - cat["balanced_accuracy"], 3)
                if (cat["balanced_accuracy"] is not None and ols["balanced_accuracy"] is not None) else None),
        }
    return {
        "artifact": "hiv_nrti_v0_validation", "schema": "hiv-nrti-validation-v0",
        "label_source": "Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra)",
        "caller": "dna_decode.data.hiv_amr v0 (POSITION-BASED NRTI major-position catalog)",
        "underlying_tool": "Stanford DRMcv.R OLS, reimplemented (sklearn); 5-fold CV",
        "dataset": str(path), "n_isolates": len(rows), "n_ols_features_min10": len(feat_names),
        "honest_caveats": [
            "v0 NRTI is POSITION-BASED -> over-calls T215 revertants / V75 polymorphisms (spec hit, esp AZT/D4T)",
            "no Subtype column -> per-subtype transfer check is v0.1 (unfiltered set)",
            "fold>=lower-cutoff binarization uses DRMcv.R's per-drug clinical cutoffs",
            "v0.1: mutant-specific NRTI catalog (data-derived OLS coefficients / sourced SDRM list)",
        ],
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from Hedlin/Stanford DRMcv.R",
        "per_drug": per_drug,
    }


def render_md(result: dict, generated: str) -> str:
    lines = [f"# HIV NRTI v0 cell - validation vs HIVDB PhenoSense + OLS baseline ({generated})", ""]
    lines.append(f"Caller = {result['caller']} (POSITION-BASED - deliberately over-calls).")
    lines.append(f"Label = {result['label_source']}. {result['n_isolates']} isolates; "
                 f"OLS features = {result['n_ols_features_min10']}.")
    lines.append("")
    lines.append("| Drug | n | prev R | cutoff | catalog sens/spec/**balacc** (AUCsep) | OLS **balacc** (AUC) | d(OLS-cat) |")
    lines.append("|---|---|---|---|---|---|---|")
    for drug, m in result["per_drug"].items():
        if "catalog" not in m:
            lines.append(f"| {drug} | {m.get('n')} | — | — | {m.get('note','')} | — | — |")
            continue
        c, o = m["catalog"], m["ols_baseline"]
        lines.append(f"| {drug} | {m['n_isolates']} | {m['prevalence_R_at_cutoff']} | "
                     f"{m['clinical_lower_cutoff_fold']} | "
                     f"{c['sens']}/{c['spec']}/**{c['balanced_accuracy']}** ({c['auc_call_separates_fold']}) | "
                     f"**{o['balanced_accuracy']}** ({o['auc']}) | {m['delta_ols_minus_catalog_balacc']} |")
    lines.append("")
    lines.append("## Honest caveats")
    for c in result["honest_caveats"]:
        lines.append(f"- {c}")
    lines.append("")
    lines.append(f"Citation: {result['citation']}.")
    return "\n".join(lines)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args(argv)
    if not args.data.exists():
        print(f"ERROR: dataset not found at {args.data}", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    result = run(args.data)
    out_md = args.out_md or (REPO / "wiki" / f"hiv_nrti_v0_validation_{today}.md")
    out_md.write_text(render_md(result, today), encoding="utf-8")
    out_md.with_suffix(".json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(render_md(result, today))
    print(f"\n[wrote {out_md} + .json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
