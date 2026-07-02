"""DepMap cell-line drug-response decoder — the cancer analog of the yeast substrate (Track A).

Genotype (CCLE damaging/hotspot mutations, gene x cell binary) x phenotype (PRISM drug LFC, lower = more
sensitive) joined by DepMap ID, de-confounded by cell-line LINEAGE (primary tissue) -- the DepMap analog of
yeast clades. Reuses the yeast within-clade machinery. Unlike yeast (copy-number mechanisms, attribution
failed), DepMap has GENE-LEVEL pharmacogenomic mechanisms (BRAF->vemurafenib, EGFR->erlotinib) so attribution
is directly testable against biological ground truth.

Data (free figshare, DepMap 19Q4 + PRISM 19Q4) staged under --data.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts.yeast_growth_decoder import _cv_r2, _r2, _within_clade_r2  # noqa: E402  (reuse machinery)

# drug -> known biomarker gene (for the attribution ground-truth test)
TARGET_DRUGS = {"vemurafenib": "BRAF", "dabrafenib": "BRAF", "erlotinib": "EGFR",
                "gefitinib": "EGFR", "nutlin-3": "TP53", "selumetinib": "BRAF"}


def load(data: Path):
    lfc = pd.read_csv(data / "prism_lfc.csv", index_col=0)                       # cell(row_name) x compound
    treat = pd.read_csv(data / "prism_treat.csv")                                # column_name -> name/moa/target
    cells = pd.read_csv(data / "prism_cells.csv")                                # row_name -> depmap_id, primary_tissue
    mut = pd.read_csv(data / "ccle_mutations.csv",
                      usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification",
                               "isDeleterious", "isCOSMIChotspot"],
                      dtype=str, low_memory=False)                               # robust: 276MB mixed-type file
    return lfc, treat, cells, mut


def gene_matrix(mut: pd.DataFrame, keep_cells) -> pd.DataFrame:
    """cell(DepMap_ID) x gene binary: 1 if the cell has a DAMAGING or COSMIC-HOTSPOT mutation in the gene."""
    dmg = mut[(mut["isDeleterious"] == "True") | (mut["isCOSMIChotspot"] == "True")]  # string booleans
    dmg = dmg[dmg["Tumor_Sample_Barcode"].isin(keep_cells)]
    g = (dmg.groupby(["Tumor_Sample_Barcode", "Hugo_Symbol"]).size().unstack(fill_value=0) > 0).astype("int8")
    return g


def _drug_cols(treat: pd.DataFrame, drug: str) -> list[str]:
    m = treat["name"].astype(str).str.lower() == drug.lower()
    return list(treat[m]["column_name"])


def univariate_top(y: np.ndarray, G: np.ndarray, genes, n=12):
    """t-stat of phenotype (mutant vs WT) per gene; sensitivity => mutant LFC LOWER => negative t."""
    ts = np.zeros(G.shape[1])
    for j in range(G.shape[1]):
        g = G[:, j]
        mut, wt = y[g == 1], y[g == 0]
        if len(mut) < 5 or len(wt) < 5:
            continue
        se = np.sqrt(mut.std() ** 2 / len(mut) + wt.std() ** 2 / len(wt))
        ts[j] = (mut.mean() - wt.mean()) / se if se > 0 else 0
    order = np.argsort(np.abs(ts))[::-1][:n]
    return [(genes[j], round(float(ts[j]), 2)) for j in order]


def within_lineage_biomarker(y, g, lineage, min_both=4) -> dict:
    """DE-CONFOUNDED single-gene attribution: does the biomarker separate response WITHIN lineages (not just
    because the mutation concentrates in a sensitive lineage)? Center y by lineage mean, compare mut vs WT."""
    yr = y.copy()
    for L in np.unique(lineage):
        m = lineage == L
        if m.sum() > 1:
            yr[m] = y[m] - np.nanmean(y[m])
    mr, wr = yr[g == 1], yr[g == 0]
    se = np.sqrt(mr.std() ** 2 / max(len(mr), 1) + wr.std() ** 2 / max(len(wr), 1))
    t_within = float((mr.mean() - wr.mean()) / se) if se > 0 else 0.0
    per = {}
    for L in np.unique(lineage):
        m = lineage == L; gm = g[m]; ym = y[m]
        if (gm == 1).sum() >= min_both and (gm == 0).sum() >= min_both:
            per[str(L)] = round(float(ym[gm == 1].mean() - ym[gm == 0].mean()), 2)
    return {"within_lineage_t": round(t_within, 2), "per_lineage_delta_lfc": per}


def run(data: Path) -> dict:
    lfc, treat, cells, mut = load(data)
    r2d = cells.set_index("row_name")["depmap_id"].to_dict()
    tis = cells.set_index("depmap_id")["primary_tissue"].to_dict()
    lfc = lfc.rename(index=r2d)                                                  # -> DepMap_ID rows
    lfc = lfc[~lfc.index.duplicated()]
    G = gene_matrix(mut, set(lfc.index))
    common = lfc.index.intersection(G.index)
    lfc, G = lfc.loc[common], G.loc[common]
    lineage = np.array([tis.get(c, "unknown") for c in common])
    keep = G.columns[G.sum(axis=0) >= 10]                                        # genes mutated in >=10 lines
    Gm = G[keep].values.astype("float32"); genes = np.array(keep)
    out = {"n_cells": len(common), "n_genes": len(genes), "n_lineages": len(set(lineage)), "drugs": []}
    for drug, target in TARGET_DRUGS.items():
        cols = _drug_cols(treat, drug)
        if not cols:
            continue
        y = lfc[cols].mean(axis=1).values.astype(float)                          # avg LFC across doses
        ok = ~np.isnan(y)
        if ok.sum() < 100:
            continue
        yo, Go, lin = y[ok], Gm[ok], lineage[ok]
        # signal-vs-null (permute)
        top = univariate_top(yo, Go, genes)
        real_max = max(abs(t) for _, t in top)
        perm = []
        for s in range(10):
            yp = yo.copy(); np.random.default_rng(s).shuffle(yp)
            perm.append(max(abs(t) for _, t in univariate_top(yp, Go, genes, n=5)))
        # attribution: rank of the known target gene
        gi = {g: i for i, g in enumerate(genes)}
        tgt_rank = next((r for r, (g, _) in enumerate(univariate_top(yo, Go, genes, n=len(genes)), 1)
                         if g == target), None)
        # within-lineage MULTIVARIATE decoder (de-confounded; ~0 for single-gene mechanisms -- see the memo)
        within, nlin = _within_clade_r2(Go, yo, lin, min_n=25)
        naive = _cv_r2(Go, yo)
        # DE-CONFOUNDED SINGLE-GENE attribution (the right test for a single-gene biomarker)
        wl = (within_lineage_biomarker(yo, Go[:, gi[target]], lin) if target in gi else {})
        out["drugs"].append({
            "drug": drug, "target": target, "n": int(ok.sum()),
            "top_genes": top[:6], "signal_real_maxt": round(real_max, 2),
            "signal_perm_maxt": round(float(np.max(perm)), 2),
            "target_gene_rank": tgt_rank, "target_is_top1": top[0][0] == target,
            "target_in_top5": target in [g for g, _ in top[:5]],
            "target_within_lineage_t": wl.get("within_lineage_t"),
            "target_per_lineage_delta_lfc": wl.get("per_lineage_delta_lfc"),
            "within_lineage_multivariate_r2": round(within, 3), "naive_r2": round(naive, 3),
            "n_lineages_used": nlin})
    return out


def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "depmap_decoder_scores.json")
    a = ap.parse_args(argv)
    res = run(a.data)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"cells={res['n_cells']} genes={res['n_genes']} lineages={res['n_lineages']}")
    for d in res["drugs"]:
        att = "TOP1" if d["target_is_top1"] else ("top5" if d["target_in_top5"] else f"rank{d['target_gene_rank']}")
        print(f"\n{d['drug']} (target {d['target']}, n={d['n']}): ATTRIB {att} | "
              f"target within-lineage_t={d['target_within_lineage_t']} | "
              f"multivariate within-lineage r2={d['within_lineage_multivariate_r2']}")
        print("   top genes: " + ", ".join(f"{g}({t})" for g, t in d["top_genes"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
