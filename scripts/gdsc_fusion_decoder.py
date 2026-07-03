"""GDSC fusion decoder — powers the DepMap fusion cases with a broad-coverage drug-response screen.

The PRISM 19Q4 subset (569 lines) had ZERO BCR-ABL1 cell lines (CML/heme lines absent) and only ~9 ALK-fusion
lines, so the fusion feature-match result (`depmap_fusion_methylation.py`) was directional-only. GDSC (Genomics
of Drug Sensitivity in Cancer) screens ~950 cell lines and DOES include the CML lines, so the BCR-ABL case
becomes real and the ALK case gains power.

Join chain: GDSC fitted dose-response (keyed by COSMIC_ID) -> DepMap `sample_info.csv` bridge (COSMIC_ID ->
DepMap_ID + lineage) -> the same CCLE_fusions.csv the DepMap decoder uses (keyed by DepMap_ID). Response metric
is LN_IC50 (lower = more sensitive). De-confounding by lineage reuses the promoted
`dna_decode.deconfound.group_centered_biomarker_t` (binary fusion presence -> within-lineage t): does a fusion
separate sensitivity WITHIN a lineage, not just because BCR-ABL concentrates in the (drug-sensitive) leukemia
lineage. Data on D: (gitignored). GDSC1 carries imatinib/ponatinib; GDSC2 carries dasatinib/nilotinib/crizotinib
(LN_IC50 is not comparable across the two assays, so each case names its dataset).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.deconfound import group_centered_biomarker_t  # noqa: E402  (promoted primitive)

# (drug, GDSC dataset, fusion-name substring, biomarker gene, direction). "sensitize": fusion+ => lower LN_IC50.
CASES = [
    ("imatinib", "GDSC1", "BCR--ABL1", "ABL1", "sensitize"),
    ("dasatinib", "GDSC2", "BCR--ABL1", "ABL1", "sensitize"),
    ("nilotinib", "GDSC2", "BCR--ABL1", "ABL1", "sensitize"),
    ("ponatinib", "GDSC1", "BCR--ABL1", "ABL1", "sensitize"),
    ("crizotinib", "GDSC2", "--ALK", "ALK", "sensitize"),
]


def load_gdsc(gdsc_dir: Path) -> pd.DataFrame:
    """Combined GDSC1+GDSC2 fitted dose-response, bridged to DepMap_ID + lineage. One row per (dataset, drug,
    cell). Reads the cached parquets (built once from the 21/29 MB xlsx)."""
    si = pd.read_csv(gdsc_dir / "sample_info.csv").dropna(subset=["COSMIC_ID"]).copy()
    si["COSMIC_ID"] = si["COSMIC_ID"].astype(int)
    cos2dm = si.set_index("COSMIC_ID")["DepMap_ID"].to_dict()
    cos2lin = si.set_index("COSMIC_ID")["lineage"].to_dict()
    frames = []
    for ds in ("GDSC1", "GDSC2"):
        pq = gdsc_dir / f"{ds.lower()}_fitted.parquet"
        if not pq.exists():
            continue
        g = pd.read_parquet(pq, columns=["COSMIC_ID", "DRUG_NAME", "LN_IC50", "Z_SCORE"])
        g = g.dropna(subset=["COSMIC_ID"]).copy()
        g["COSMIC_ID"] = g["COSMIC_ID"].astype(int)
        g["DATASET"] = ds
        g["drug"] = g["DRUG_NAME"].astype(str).str.lower()
        g["DepMap_ID"] = g["COSMIC_ID"].map(cos2dm)
        g["lineage"] = g["COSMIC_ID"].map(cos2lin).fillna("unknown")
        frames.append(g.dropna(subset=["DepMap_ID"]))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def fusion_ids(fusions_csv: Path, pattern: str) -> set:
    """DepMap_IDs carrying a fusion whose #FusionName contains `pattern`."""
    fu = pd.read_csv(fusions_csv, usecols=["DepMap_ID", "#FusionName"])
    return set(fu[fu["#FusionName"].astype(str).str.contains(pattern, regex=False)]["DepMap_ID"])


def _direction_ok(value: float, direction: str) -> bool:
    return (value < 0) if direction == "sensitize" else (value > 0)


def score_case(gdsc: pd.DataFrame, fus_ids: set, drug: str, dataset: str, direction: str) -> dict:
    sub = gdsc[(gdsc["drug"] == drug) & (gdsc["DATASET"] == dataset)]
    # one LN_IC50 per cell line (mean over replicate curves)
    s = sub.groupby("DepMap_ID")["LN_IC50"].mean()
    lin = sub.drop_duplicates("DepMap_ID").set_index("DepMap_ID")["lineage"].reindex(s.index)
    ids = list(s.index)
    y = s.values.astype(float)
    g = np.array([1 if c in fus_ids else 0 for c in ids], dtype=float)
    groups = lin.fillna("unknown").values
    n_pos = int(g.sum())
    pos, neg = y[g == 1], y[g == 0]
    # global effect + Mann-Whitney
    if n_pos >= 1 and (g == 0).sum() >= 1:
        u_p = float(mannwhitneyu(pos, neg, alternative="less" if direction == "sensitize" else "greater")[1])
        global_delta = float(pos.mean() - neg.mean())
    else:
        u_p, global_delta = float("nan"), float("nan")
    bt = group_centered_biomarker_t(y, g, groups)
    within_t = bt["within_lineage_t"]
    powered = n_pos >= 8
    matched = bool(_direction_ok(within_t, direction) and abs(within_t) > 1.0 and n_pos >= 4)
    return {
        "drug": drug, "dataset": dataset, "biomarker_fusion": None, "direction": direction,
        "n_cell_lines": len(ids), "n_fusion_positive": n_pos, "powered": powered,
        "global_ln_ic50_delta": round(global_delta, 2) if not np.isnan(global_delta) else None,
        "global_mannwhitney_p": round(u_p, 5) if not np.isnan(u_p) else None,
        "within_lineage_t": within_t,
        "per_lineage_delta_ln_ic50": bt["per_lineage_delta_lfc"],
        "direction_ok": bool(_direction_ok(within_t, direction)),
        "feature_match_powered": bool(matched and powered),
    }


def run(gdsc_dir: Path, fusions_csv: Path) -> dict:
    gdsc = load_gdsc(gdsc_dir)
    out = {"n_gdsc_rows": len(gdsc), "n_cell_lines_total": int(gdsc["DepMap_ID"].nunique()) if len(gdsc) else 0,
           "cases": []}
    for drug, dataset, pat, gene, direction in CASES:
        fus = fusion_ids(fusions_csv, pat)
        rec = score_case(gdsc, fus, drug, dataset, direction)
        rec["biomarker_fusion"] = pat
        rec["biomarker_gene"] = gene
        out["cases"].append(rec)
    return out


def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gdsc-dir", type=Path, default=Path("D:/dna_decode_cache/gdsc"))
    ap.add_argument("--fusions-csv", type=Path, default=Path("D:/dna_decode_cache/depmap_pilot/ccle_fusions.csv"))
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "gdsc_fusion_scores.json")
    a = ap.parse_args(argv)
    res = run(a.gdsc_dir, a.fusions_csv)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"GDSC rows={res['n_gdsc_rows']} cell_lines={res['n_cell_lines_total']}")
    for r in res["cases"]:
        mark = "POWERED-MATCH" if r["feature_match_powered"] else ("dir-ok" if r["direction_ok"] else "  -  ")
        pw = "" if r["powered"] else " [underpowered]"
        print(f"[{mark:13}] {r['drug']:11} ({r['dataset']}) x {r['biomarker_fusion']} "
              f"n+={r['n_fusion_positive']}{pw}: global dLN_IC50={r['global_ln_ic50_delta']} "
              f"(MWU p={r['global_mannwhitney_p']}) | within-lineage t={r['within_lineage_t']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
