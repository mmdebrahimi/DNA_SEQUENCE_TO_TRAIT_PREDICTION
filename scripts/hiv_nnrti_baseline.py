"""Validate the v0 HIV NNRTI catalog vs the UNDERLYING TOOL — a faithful Python reimplementation of
Stanford's DRMcv.R cross-validated OLS regression baseline (Wave B, "validate-the-wrapper-vs-the-tool").

The project's standing discipline: a policy layer over a curated tool must be validated against NAIVE use
of the underlying tool on the same data — in-cohort accuracy alone only proves the tool works, not that the
layer adds value. Here:
  - WRAPPER  = dna_decode.data.hiv_amr (the deterministic "any major NNRTI DRM present -> R" catalog).
  - UNDERLYING TOOL = Stanford's DRMcv.R: ordinary least squares of log10(fold-change) ~ binary
    mutation-presence features, 5-fold cross-validated. Reimplemented here in Python (sklearn) — R is not
    installed; the reimplementation mirrors DRMcv.R: design matrix = 1 if the isolate carries amino-acid `aa`
    at RT position `p` (features kept iff present in >= min_muts isolates), Y = log10(fold), lm + k-fold CV.

CLINICAL CUTOFF (sourced from DRMcv.R lines 179-182, not invented): for EFV/NVP/ETR/RPV the lower
fold-resistance cutoff is 3 (susceptible < 3; intermediate 3-10; high-level > 10). v0 binary label =
fold >= 3 (non-susceptible). DOR is newer than the 2014 script -> 3 reused + flagged.

HEADLINE: per drug, the catalog's balanced accuracy vs the OLS baseline's, at the SAME clinical cutoff +
the SAME isolates. A small delta = the simple deterministic catalog matches the full regression (the
honest "the wrapper adds interpretability, not error"); a large delta on ETR/RPV/DOR = the per-drug signal
the class-level v0 misses (bounds + motivates the v0.1 per-drug catalog). The OLS per-mutation coefficients
are ALSO the data-derived per-drug DRM signal for that v0.1 catalog (provenance-clean, no web search).

DATA: data/raw/hiv/NNRTI_DataSet.txt (gitignored; cite Rhee 2003).
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

from dna_decode.data.hiv_amr import call_from_observed_substitutions
from scripts.hiv_nnrti_validate import (
    DEFAULT_DATA, DRUGS, _DRUG_COL, _observed_rt_mutations, _parse_fold, load_rows,
)

# Sourced from DRMcv.R confusion-matrix cutoffs (lower bound = susceptible/non-susceptible boundary).
NNRTI_CLINICAL_LOWER_CUTOFF = {"efavirenz": 3.0, "nevirapine": 3.0, "etravirine": 3.0,
                               "rilpivirine": 3.0, "doravirine": 3.0}  # DOR reused (not in 2014 script)
MIN_MUTS = 10        # DRMcv.R default: drop a mutation feature present in < 10 isolates
CV_FOLDS = 5         # DRMcv.R default nfold


def _position_columns(header: list[str]) -> list[str]:
    return [h for h in header if len(h) > 1 and h[0] == "P" and h[1:].isdigit()]


def _build_design_matrix(rows, pos_cols):
    """Binary mutation-presence design matrix (mirrors DRMcv.R buildX + the min_muts filter).

    Feature (p, aa) = 1 iff the isolate carries amino acid `aa` at RT position p (`-`=consensus). Keeps
    only features present in >= MIN_MUTS isolates. Returns (X[n_rows, n_feats], feature_names)."""
    counts: dict[str, int] = {}
    per_row_feats: list[set[str]] = []
    for row in rows:
        feats = set()
        for c in pos_cols:
            cell = (row.get(c) or "").strip()
            if cell in ("", "-", ".", "NA"):
                continue
            pos = c[1:]
            for aa in cell:
                if aa.isalpha():
                    feats.add(f"{pos}{aa}")
        per_row_feats.append(feats)
        for f in feats:
            counts[f] = counts.get(f, 0) + 1
    feature_names = sorted(f for f, n in counts.items() if n >= MIN_MUTS)
    fidx = {f: j for j, f in enumerate(feature_names)}
    X = np.zeros((len(rows), len(feature_names)), dtype=float)
    for i, feats in enumerate(per_row_feats):
        for f in feats:
            j = fidx.get(f)
            if j is not None:
                X[i, j] = 1.0
    return X, feature_names


def _confusion(pred_R, actual_R):
    tp = int(np.sum(pred_R & actual_R)); fp = int(np.sum(pred_R & ~actual_R))
    tn = int(np.sum(~pred_R & ~actual_R)); fn = int(np.sum(~pred_R & actual_R))
    sens = tp / (tp + fn) if (tp + fn) else None
    spec = tn / (tn + fp) if (tn + fp) else None
    acc = (tp + tn) / (tp + fp + tn + fn) if (tp + fp + tn + fn) else None
    bal = (sens + spec) / 2 if (sens is not None and spec is not None) else None
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "sens": round(sens, 3) if sens is not None else None,
            "spec": round(spec, 3) if spec is not None else None,
            "accuracy": round(acc, 3) if acc is not None else None,
            "balanced_accuracy": round(bal, 3) if bal is not None else None}


def _auc(score, label_R):
    pos = score[label_R]; neg = score[~label_R]
    if len(pos) == 0 or len(neg) == 0:
        return None
    allv = np.concatenate([pos, neg])
    order = np.argsort(allv, kind="mergesort")
    ranks = np.empty(len(allv), dtype=float)
    i = 0
    sorted_v = allv[order]
    while i < len(allv):
        j = i
        while j < len(allv) and sorted_v[j] == sorted_v[i]:
            j += 1
        ranks[order[i:j]] = (i + j - 1) / 2.0 + 1.0
        i = j
    rank_sum_pos = ranks[:len(pos)].sum()
    u = rank_sum_pos - len(pos) * (len(pos) + 1) / 2.0
    return round(u / (len(pos) * len(neg)), 4)


def run(path: Path = DEFAULT_DATA, seed: int = 0) -> dict:
    rows = load_rows(path)
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
    pos_cols = _position_columns(header)
    X_all, feat_names = _build_design_matrix(rows, pos_cols)

    per_drug = {}
    for drug in DRUGS:
        cutoff = NNRTI_CLINICAL_LOWER_CUTOFF[drug]
        col = _DRUG_COL[drug]
        keep, folds, catalog_R = [], [], []
        for i, row in enumerate(rows):
            fold = _parse_fold(row.get(col, ""))
            if fold is None or fold <= 0:
                continue
            keep.append(i)
            folds.append(fold)
            obs = _observed_rt_mutations(row)
            catalog_R.append(call_from_observed_substitutions(drug, {"RT": obs}).prediction == "R")
        if len(keep) < 30:
            per_drug[drug] = {"n": len(keep), "note": "too few isolates"}
            continue
        idx = np.array(keep)
        fold = np.array(folds)
        actual_R = fold >= cutoff
        catalog_R = np.array(catalog_R, dtype=bool)
        Xd = X_all[idx]
        y = np.log10(fold)
        # 5-fold cross-validated OOF prediction of log10(fold) (the underlying-tool baseline).
        oof_logfold = cross_val_predict(LinearRegression(), Xd, y, cv=CV_FOLDS)
        oof_fold = np.power(10.0, oof_logfold)
        ols_R = oof_fold >= cutoff
        per_drug[drug] = {
            "n_isolates": len(keep),
            "prevalence_R_at_cutoff": round(float(actual_R.mean()), 3),
            "clinical_lower_cutoff_fold": cutoff,
            "catalog": _confusion(catalog_R, actual_R),
            "ols_baseline": {**_confusion(ols_R, actual_R), "auc": _auc(oof_fold, actual_R)},
        }
        cb = per_drug[drug]["catalog"]["balanced_accuracy"]
        ob = per_drug[drug]["ols_baseline"]["balanced_accuracy"]
        per_drug[drug]["delta_ols_minus_catalog_balacc"] = (
            round(ob - cb, 3) if (cb is not None and ob is not None) else None)
    return {
        "artifact": "hiv_nnrti_validate_vs_underlying_tool",
        "schema": "hiv-nnrti-baseline-v0",
        "wrapper": "dna_decode.data.hiv_amr (deterministic class-level major-DRM catalog)",
        "underlying_tool": "Stanford DRMcv.R OLS (log10 fold ~ binary mutation presence), 5-fold CV, "
                           "reimplemented in Python/sklearn (R not installed)",
        "clinical_cutoff_source": "DRMcv.R confusion-matrix lower cutoffs (EFV/NVP/ETR/RPV=3; DOR reused)",
        "n_features_min10": len(feat_names),
        "interpretation": "small delta => the simple deterministic catalog MATCHES the full regression "
                          "(adds interpretability, not error); large delta on 2nd-gen NNRTIs => per-drug "
                          "signal the class-level v0 misses (bounds the v0.1 per-drug catalog)",
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; method Hedlin/Stanford DRMcv.R",
        "per_drug": per_drug,
    }


def render_md(result: dict, generated: str) -> str:
    lines = [f"# HIV NNRTI v0 — validate the catalog vs the underlying tool (OLS) ({generated})", ""]
    lines.append(f"Wrapper = {result['wrapper']}.")
    lines.append(f"Underlying tool = {result['underlying_tool']}.")
    lines.append(f"Clinical cutoff source = {result['clinical_cutoff_source']}. "
                 f"OLS features (>=10 isolates) = {result['n_features_min10']}.")
    lines.append("")
    lines.append("Both scored on the SAME isolates + the SAME clinical cutoff (fold>=lower). "
                 "balacc = balanced accuracy ((sens+spec)/2).")
    lines.append("")
    lines.append("| Drug | n | prev R | catalog sens/spec/**balacc** | OLS sens/spec/**balacc** (AUC) | d(OLS-cat) |")
    lines.append("|---|---|---|---|---|---|")
    for drug, m in result["per_drug"].items():
        if "catalog" not in m:
            lines.append(f"| {drug} | {m.get('n')} | — | {m.get('note','')} | — | — |")
            continue
        c, o = m["catalog"], m["ols_baseline"]
        lines.append(f"| {drug} | {m['n_isolates']} | {m['prevalence_R_at_cutoff']} | "
                     f"{c['sens']}/{c['spec']}/**{c['balanced_accuracy']}** | "
                     f"{o['sens']}/{o['spec']}/**{o['balanced_accuracy']}** ({o['auc']}) | "
                     f"{m['delta_ols_minus_catalog_balacc']} |")
    lines.append("")
    lines.append(f"**Interpretation:** {result['interpretation']}.")
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
    out_md = args.out_md or (REPO / "wiki" / f"hiv_nnrti_baseline_vs_ols_{today}.md")
    out_md.write_text(render_md(result, today), encoding="utf-8")
    out_md.with_suffix(".json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(render_md(result, today))
    print(f"\n[wrote {out_md} + .json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
