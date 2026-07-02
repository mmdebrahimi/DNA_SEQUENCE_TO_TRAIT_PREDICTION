"""DepMap FUSION + METHYLATION modalities — the feature-match law extended to two more feature types.

The multimodal DepMap decoder proved mutation / copy-number / expression EACH carry the de-confounded signal
for the matching mechanism (BRAF-mutation->vemurafenib; EGFR-amplification/expression->erlotinib). This adds
two feature types that are INVISIBLE to those three:

  - GENE FUSION (binary presence, cell x fusion): a fusion (EML4-ALK, BCR-ABL) is not an SNV, not a copy
    gain, and does not raise single-gene expression of the partner -> ONLY a fusion-presence feature captures
    it. Test: crizotinib / alectinib / lorlatinib -> ALK fusion.
  - PROMOTER METHYLATION (continuous RRBS beta, cell x gene): methylation SILENCES a gene with no coding
    change -> ONLY methylation captures it. Test: temozolomide -> MGMT silencing (sensitize); olaparib /
    topotecan -> SLFN11 silencing (resist).

For each case the matched modality is scored de-confounded (residualize by primary tissue / lineage), and the
"wrong" modalities (mutation / copy-number / expression of the SAME biomarker gene) are scored alongside to
show the matched feature type WINS — the direct multi-modal proof that feature type must match mechanism type.

De-confounding reuses the promoted `dna_decode.deconfound` primitives: binary fusion -> within-lineage
biomarker t; continuous methylation -> clade-centered Spearman. Data on D: (gitignored-class). NOTE: PRISM
19Q4 is a 569-line SUBSET -> fusion cases are UNDERPOWERED (ALK n~9; BCR-ABL n=0, CML lines absent) and
reported as DIRECTIONAL; methylation is well-powered (~487 bridged lines).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import scripts.depmap_decoder as DD  # noqa: E402  (reuse PRISM/mutation load + gene_matrix)
from dna_decode.deconfound import (  # noqa: E402  (promoted de-confounding primitives)
    group_centered_association,
    group_centered_biomarker_t,
)

# --- fusion cases: (drug, fusion-name substring, biomarker gene, direction) ---
# direction "sensitize": fusion+ => MORE sensitive => LOWER LFC => negative within-lineage t.
FUSION_CASES = [
    ("crizotinib", "--ALK", "ALK", "sensitize"),
    ("alectinib", "--ALK", "ALK", "sensitize"),
    ("lorlatinib", "--ALK", "ALK", "sensitize"),
]
# --- methylation cases: (drug, gene, direction) ---
# "sensitize": high methylation => silenced => sensitive => negative rho (MGMT->TMZ).
# "resist":    high methylation => silenced => resistant => positive rho (SLFN11->PARPi/topo).
METHYL_CASES = [
    ("temozolomide", "MGMT", "sensitize"),
    ("olaparib", "SLFN11", "resist"),
    ("topotecan", "SLFN11", "resist"),
]
_FUSION_GENES = sorted({c[2] for c in FUSION_CASES})
_METHYL_GENES = sorted({c[1] for c in METHYL_CASES})


def load_base(data: Path):
    """PRISM LFC (depmap_id index) + treat + lineage map + ccle_name->depmap_id bridge + mutation matrix."""
    lfc, treat, cells, mut = DD.load(data)
    r2d = cells.set_index("row_name")["depmap_id"].to_dict()
    tis = cells.set_index("depmap_id")["primary_tissue"].to_dict()
    ccle2dm = cells.dropna(subset=["ccle_name"]).set_index("ccle_name")["depmap_id"].to_dict()
    lfc = lfc.rename(index=r2d)
    lfc = lfc[~lfc.index.duplicated()]
    return lfc, treat, tis, ccle2dm, mut


def fusion_presence(data: Path, keep_ids) -> pd.DataFrame:
    """cell(depmap_id) x fusion-pattern binary: 1 if the cell carries a fusion whose #FusionName contains
    the pattern. Absence of a row for a cell means no reported fusion (=> 0)."""
    fu = pd.read_csv(data / "ccle_fusions.csv", usecols=["DepMap_ID", "#FusionName"])
    fu = fu[fu["DepMap_ID"].isin(keep_ids)]
    patterns = sorted({c[1] for c in FUSION_CASES})              # unique patterns (multiple drugs share one)
    out = pd.DataFrame(0, index=sorted(keep_ids), columns=patterns, dtype="int8")
    for pat in patterns:
        hit = set(fu[fu["#FusionName"].astype(str).str.contains(pat, regex=False)]["DepMap_ID"])
        for cid in hit:
            if cid in out.index:
                out.at[cid, pat] = 1
    return out


def methylation_genes(data: Path, ccle2dm: dict, keep_ids) -> pd.DataFrame:
    """gene x cell(depmap_id) mean promoter methylation (RRBS TSS_1kb beta, 0..1). A gene may have several
    TSS_id rows -> averaged. RRBS columns are CCLE names -> bridged to depmap_id via ccle2dm."""
    path = data / "ccle_rrbs_tss_1kb.txt"
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    cell_cols = header[7:]
    # keep only cell columns that bridge to a PRISM depmap_id
    keep_cells = [c for c in cell_cols if c in ccle2dm and ccle2dm[c] in keep_ids]
    usecols = ["gene"] + keep_cells
    rows = []
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=50000, dtype=str):
        rows.append(chunk[chunk["gene"].isin(_METHYL_GENES)])
    sub = pd.concat(rows) if rows else pd.DataFrame(columns=usecols)
    for c in keep_cells:                                             # RRBS marks missing as '     NA' (spaces)
        sub[c] = pd.to_numeric(sub[c].str.strip(), errors="coerce")
    gm = sub.groupby("gene")[keep_cells].mean()                      # gene x ccle_name
    gm = gm.rename(columns={c: ccle2dm[c] for c in keep_cells})       # -> gene x depmap_id
    return gm


def _read_gene_cols(path: Path, genes) -> pd.DataFrame:
    """Read col 0 (cell-line id -> index) + the columns whose gene symbol (before ' (') is in `genes`, by
    INTEGER position (robust to the unnamed-first-column quirk). Used for the CN/expr cross-check."""
    if not path.exists():
        return pd.DataFrame()
    with open(path, encoding="utf-8") as fh:
        cols = fh.readline().rstrip("\n").split(",")
    keep = set(genes)
    pos = [0] + [i for i, c in enumerate(cols) if c.split(" (")[0] in keep]
    return pd.read_csv(path, usecols=pos, index_col=0)


def _direction_ok(value: float, direction: str) -> bool:
    return (value < 0) if direction == "sensitize" else (value > 0)


def _wrong_modality_stats(gene, lfc, cols, tis, mut, cn, ex):
    """Cross-check: does mutation / copy-number / expression of the SAME biomarker gene carry the signal?
    Returns clade-centered stats so we can show the matched modality wins."""
    y_all = lfc[cols].mean(axis=1)
    out = {}
    G = DD.gene_matrix(mut, set(lfc.index)) if mut is not None else pd.DataFrame()
    # mutation (binary) -> within-lineage biomarker t
    if gene in getattr(G, "columns", []):
        common = y_all.index.intersection(G.index)
        y = y_all.loc[common].astype(float).values
        g = G.loc[common, gene].astype(float).values
        lin = np.array([tis.get(c, "unknown") for c in common])
        out["mutation"] = {"within_lineage_t": group_centered_biomarker_t(y, g, lin)["within_lineage_t"],
                           "n_carriers": int(g.sum())}
    # copy-number + expression (continuous) -> clade-centered rho
    for mod, src in (("copy_number", cn), ("expression", ex)):
        col = next((c for c in getattr(src, "columns", []) if c.split(" (")[0] == gene), None)
        if col is None:
            continue
        common = y_all.index.intersection(src.index)
        if len(common) < 20:
            continue
        y = y_all.loc[common].astype(float).values
        x = pd.to_numeric(src.loc[common, col], errors="coerce").astype(float).values
        lin = np.array([tis.get(c, "unknown") for c in common])
        _, within, n = group_centered_association(y, x, lin)
        out[mod] = {"within_lineage_rho": None if np.isnan(within) else round(float(within), 3), "n": n}
    return out


def run(data: Path) -> dict:
    lfc, treat, tis, ccle2dm, mut = load_base(data)
    lfc = lfc[[isinstance(i, str) for i in lfc.index]]           # drop unmapped (NaN depmap_id) rows
    keep_ids = set(lfc.index)
    fus = fusion_presence(data, keep_ids)
    meth = methylation_genes(data, ccle2dm, keep_ids)
    biomarkers = set(_FUSION_GENES) | set(_METHYL_GENES)
    cn = _read_gene_cols(data / "ccle_gene_cn.csv", biomarkers)
    ex = _read_gene_cols(data / "ccle_expression.csv", biomarkers)

    out = {"n_cells": len(keep_ids), "fusion_cases": [], "methylation_cases": []}

    # --- FUSION (binary presence -> within-lineage biomarker t) ---
    for drug, pat, gene, direction in FUSION_CASES:
        cols = DD._drug_cols(treat, drug)
        if not cols or pat not in fus.columns:
            continue
        y_all = lfc[cols].mean(axis=1)
        common = y_all.index.intersection(fus.index)
        y = y_all.loc[common].astype(float).values
        g = fus.loc[common, pat].astype(float).values
        lin = np.array([tis.get(c, "unknown") for c in common])
        ok = ~np.isnan(y)
        y, g, lin = y[ok], g[ok], lin[ok]
        bt = group_centered_biomarker_t(y, g, lin)
        n_pos = int(g.sum())
        matched_ok = bool(_direction_ok(bt["within_lineage_t"], direction) and abs(bt["within_lineage_t"]) > 1.0)
        wrong = _wrong_modality_stats(gene, lfc, cols, tis, mut, cn, ex)
        out["fusion_cases"].append({
            "drug": drug, "fusion": pat, "biomarker_gene": gene, "direction": direction,
            "n_fusion_positive": n_pos, "underpowered": n_pos < 12,
            "fusion_within_lineage_t": bt["within_lineage_t"],
            "fusion_per_lineage_delta_lfc": bt["per_lineage_delta_lfc"],
            "fusion_direction_ok": matched_ok,
            "wrong_modalities": wrong,
            "feature_match": matched_ok,   # fusion carries a correct-direction de-confounded signal
        })

    # --- METHYLATION (continuous -> clade-centered Spearman) ---
    for drug, gene, direction in METHYL_CASES:
        cols = DD._drug_cols(treat, drug)
        if not cols or gene not in meth.index:
            continue
        y_all = lfc[cols].mean(axis=1)
        mser = meth.loc[gene]
        common = y_all.index.intersection(mser.index)
        y = y_all.loc[common].astype(float).values
        x = pd.to_numeric(mser.loc[common], errors="coerce").astype(float).values
        lin = np.array([tis.get(c, "unknown") for c in common])
        g_rho, within, n = group_centered_association(y, x, lin)
        matched_ok = bool(not np.isnan(within) and abs(within) > 0.1 and _direction_ok(within, direction))
        wrong = _wrong_modality_stats(gene, lfc, cols, tis, mut, cn, ex)
        out["methylation_cases"].append({
            "drug": drug, "biomarker_gene": gene, "direction": direction, "n": n,
            "methylation_global_rho": round(float(g_rho), 3) if not np.isnan(g_rho) else None,
            "methylation_within_lineage_rho": None if np.isnan(within) else round(float(within), 3),
            "methylation_direction_ok": matched_ok,
            "wrong_modalities": wrong,
            "feature_match": matched_ok,
        })
    return out


def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "depmap_fusion_methylation_scores.json")
    a = ap.parse_args(argv)
    res = run(a.data)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"cells={res['n_cells']}")
    print("\n== FUSION (binary presence; underpowered in PRISM subset) ==")
    for r in res["fusion_cases"]:
        mark = "MATCH" if r["feature_match"] else "  -  "
        up = " [UNDERPOWERED]" if r["underpowered"] else ""
        print(f"[{mark}] {r['drug']:11} x {r['fusion']} (n+={r['n_fusion_positive']}){up}: "
              f"fusion within-lineage_t={r['fusion_within_lineage_t']} | wrong={ {k: v.get('within_lineage_t', v.get('within_lineage_rho')) for k,v in r['wrong_modalities'].items()} }")
    print("\n== METHYLATION (continuous; well-powered) ==")
    for r in res["methylation_cases"]:
        mark = "MATCH" if r["feature_match"] else "  -  "
        print(f"[{mark}] {r['drug']:13} x {r['biomarker_gene']} methylation ({r['direction']}, n={r['n']}): "
              f"within-lineage_rho={r['methylation_within_lineage_rho']} (global={r['methylation_global_rho']}) | "
              f"wrong={ {k: v.get('within_lineage_rho', v.get('within_lineage_t')) for k,v in r['wrong_modalities'].items()} }")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
