"""Validate the v0 HIV PI/INSTI/CAI catalogs vs the UNDERLYING TOOL — the OLS regression baseline.

The project's load-bearing discipline (validate-the-wrapper-vs-the-tool): a policy layer over a curated
tool must be validated against NAIVE use of the underlying tool on the same data — in-cohort accuracy alone
only proves the tool works, not that the layer ADDS VALUE. NNRTI + NRTI already have this baseline; this
extends it to the three target-site classes added 2026-06-22.

  - WRAPPER         = dna_decode.data.hiv_amr (the deterministic position-based PI/INSTI or mutant-level CAI catalog).
  - UNDERLYING TOOL = Stanford DRMcv.R-style OLS: log10(fold-change) ~ binary mutation-presence features,
    5-fold cross-validated (reimplemented in Python/sklearn — R not installed). REUSES the exact
    `_build_design_matrix` / `_confusion` / `_auc` machinery from scripts/hiv_nnrti_baseline (gene-generic).

CUTOFF (honest): a UNIFORM illustrative fold>=3 boundary (the same one scripts/hiv_targetsite_validate
already uses + flags), NOT a per-drug clinical breakpoint. Both the catalog and the OLS are scored at the
SAME cutoff + the SAME isolates, so the DELTA (the headline) is a fair wrapper-vs-tool comparison; only the
absolute calibration would need per-drug clinical cutoffs (a v0.1 item — those weren't sourced here, no
fabrication). For INSTIs especially the clinical lower cutoff is drug-specific (<3 for 1st-gen) so the
absolute balaccs are illustrative; the delta is the signal.

DATA: data/raw/hiv/{PI,INI,CAI}_DataSet.txt (gitignored; cite Rhee 2003).
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_predict

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.hiv_amr import call_hiv_observed  # noqa: E402
from scripts.hiv_nnrti_baseline import (  # noqa: E402 (reuse the gene-generic OLS machinery)
    CV_FOLDS, _auc, _build_design_matrix, _confusion, _position_columns,
)
from scripts.hiv_nnrti_validate import _parse_fold, load_rows  # noqa: E402
from scripts.hiv_targetsite_validate import (  # noqa: E402
    ILLUSTRATIVE_FOLD_CUTOFF, _CLASS_SPECS, _observed_mutations,
)


def run(label: str, path: Path) -> dict:
    cls, _, drug_cols = _CLASS_SPECS[label]
    cutoff = ILLUSTRATIVE_FOLD_CUTOFF
    rows = load_rows(path)
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
    pos_cols = _position_columns(header)
    X_all, feat_names = _build_design_matrix(rows, pos_cols)

    per_drug = {}
    for drug, col in drug_cols.items():
        keep, folds, catalog_R = [], [], []
        for i, row in enumerate(rows):
            fold = _parse_fold(row.get(col, ""))
            if fold is None or fold <= 0:
                continue
            keep.append(i)
            folds.append(fold)
            obs = _observed_mutations(row, cls)
            catalog_R.append(call_hiv_observed(drug, {cls.gene: obs}).prediction == "R")
        if len(keep) < 30:
            per_drug[drug] = {"n": len(keep), "note": "too few isolates"}
            continue
        idx = np.array(keep)
        fold = np.array(folds)
        actual_R = fold >= cutoff
        catalog_R = np.array(catalog_R, dtype=bool)
        if actual_R.all() or not actual_R.any():
            per_drug[drug] = {"n_isolates": len(keep),
                              "note": "single-class at illustrative cutoff (no S/R contrast)",
                              "catalog": _confusion(catalog_R, actual_R)}
            continue
        Xd = X_all[idx]
        y = np.log10(fold)
        oof_logfold = cross_val_predict(LinearRegression(), Xd, y, cv=CV_FOLDS)
        oof_fold = np.power(10.0, oof_logfold)
        ols_R = oof_fold >= cutoff
        per_drug[drug] = {
            "n_isolates": len(keep),
            "prevalence_R_at_cutoff": round(float(actual_R.mean()), 3),
            "illustrative_lower_cutoff_fold": cutoff,
            "catalog": _confusion(catalog_R, actual_R),
            "ols_baseline": {**_confusion(ols_R, actual_R), "auc": _auc(oof_fold, actual_R)},
        }
        cb = per_drug[drug]["catalog"]["balanced_accuracy"]
        ob = per_drug[drug]["ols_baseline"]["balanced_accuracy"]
        per_drug[drug]["delta_ols_minus_catalog_balacc"] = (
            round(ob - cb, 3) if (cb is not None and ob is not None) else None)
    return {
        "artifact": f"hiv_{label.lower()}_validate_vs_underlying_tool",
        "schema": "hiv-targetsite-baseline-v0",
        "drug_class": label, "gene": cls.gene,
        "call_mode": "mutant-level" if cls.major_drms is not None else "position-based",
        "wrapper": f"dna_decode.data.hiv_amr ({label} catalog)",
        "underlying_tool": "Stanford DRMcv.R-style OLS (log10 fold ~ binary mutation presence), 5-fold CV, "
                           "Python/sklearn reimpl (shared machinery with the NNRTI baseline)",
        "cutoff_note": "UNIFORM illustrative fold>=3 (NOT per-drug clinical); both models scored at the same "
                       "cutoff so the DELTA is the fair wrapper-vs-tool signal; absolute calibration needs "
                       "per-drug clinical cutoffs (v0.1, not sourced here)",
        "n_features_min10": len(feat_names),
        "interpretation": "small delta => the deterministic catalog MATCHES the full regression (adds "
                          "interpretability, not error); large positive delta => per-drug/mutant signal the "
                          "catalog misses (bounds a v0.1 refinement)",
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; method Stanford DRMcv.R",
        "per_drug": per_drug,
    }


def render_md(result: dict, generated: str) -> str:
    L = [f"# HIV {result['drug_class']} v0 — catalog vs underlying tool (OLS) ({generated})", ""]
    L.append(f"Wrapper = {result['wrapper']} ({result['call_mode']}). Gene = {result['gene']}.")
    L.append(f"Underlying tool = {result['underlying_tool']}.")
    L.append(f"Cutoff: {result['cutoff_note']}. OLS features (>=10 isolates) = {result['n_features_min10']}.")
    L.append("")
    L.append("| Drug | n | prev R | catalog sens/spec/**balacc** | OLS sens/spec/**balacc** (AUC) | d(OLS-cat) |")
    L.append("|---|---|---|---|---|---|")
    for drug, m in result["per_drug"].items():
        if "catalog" not in m:
            L.append(f"| {drug} | {m.get('n')} | — | {m.get('note','')} | — | — |")
            continue
        c = m["catalog"]
        if "ols_baseline" not in m:
            L.append(f"| {drug} | {m['n_isolates']} | — | {c['sens']}/{c['spec']}/**{c['balanced_accuracy']}** "
                     f"| {m.get('note','')} | — |")
            continue
        o = m["ols_baseline"]
        L.append(f"| {drug} | {m['n_isolates']} | {m['prevalence_R_at_cutoff']} | "
                 f"{c['sens']}/{c['spec']}/**{c['balanced_accuracy']}** | "
                 f"{o['sens']}/{o['spec']}/**{o['balanced_accuracy']}** ({o['auc']}) | "
                 f"{m['delta_ols_minus_catalog_balacc']} |")
    L.append("")
    L.append(f"**Interpretation:** {result['interpretation']}.")
    L.append("")
    L.append(f"Citation: {result['citation']}.")
    return "\n".join(L)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--class", dest="cls", required=True, choices=sorted(_CLASS_SPECS))
    ap.add_argument("--data", type=Path, default=None)
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args(argv)
    _, fname, _ = _CLASS_SPECS[args.cls]
    data = args.data or (REPO / "data" / "raw" / "hiv" / fname)
    if not data.exists():
        print(f"ERROR: dataset not found at {data}", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    result = run(args.cls, data)
    out_md = args.out_md or (REPO / "wiki" / f"hiv_{args.cls.lower()}_baseline_vs_ols_{today}.md")
    out_md.write_text(render_md(result, today), encoding="utf-8")
    out_md.with_suffix(".json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(render_md(result, today))
    print(f"\n[wrote {out_md} + .json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
