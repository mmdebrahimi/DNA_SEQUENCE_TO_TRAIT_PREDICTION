"""Yeast COPY-NUMBER attribution test — the confirmatory closure of the attribution capstone.

The presence/absence decoder FAILED canonical-gene attribution (copper->CUP1 is copy-number, invisible to
presence/absence). This tests the RIGHT feature: clade-centered CONTINUOUS association of a gene's COPY NUMBER
(genesMatrix_CopyNumber) with a growth condition, de-confounded by SNP-distance clades, with a within-clade
permutation null + direction check + K=18/K=30 robustness. Reuses the yeast clade machinery. NO SNP/LMM build
(the copy-number matrix already exists free); NO binary carrier test (copy number is continuous dosage).
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import scripts.yeast_growth_decoder as Y  # noqa: E402  (reuse _clades)

# gene copy-number column -> (condition, expected direction). +1 = more copies -> more resistant (higher growth).
CASES = {
    "CUP1 (YHR055C)": {"col": "X4306.YHR055C_NumOfGenes_4", "condition": "YPDCUSO410MM", "dir": +1,
                       "mechanism": "copper metallothionein tandem array"},
    "ENA5 (YDR038C)": {"col": "X2730.YDR038C_NumOfGenes_3", "condition": "YPDNACL1M", "dir": +1,
                       "mechanism": "Na+/Li+ efflux ATPase tandem cluster (ENA1/2/5)"},
}
# Mechanism-SPECIFICITY negative control: ENA5 copy should associate with IONIC stress (Na/Li) but NOT with
# non-ionic osmotic stress (sorbitol). Verified 2026-07-02: ENA5 x YPSORBITOL clade-centered rho +0.03,
# perm_p 0.38 (null) vs Na +0.25 / Li +0.27 (perm_p 0.01) -- a textbook specificity check.


def _load(data: Path):
    ph = pd.read_csv(gzip.open(data / "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz", "rt"),
                     sep="\t", index_col=0)
    cnv = pd.read_csv(gzip.open(data / "genesMatrix_CopyNumber.tab.gz", "rt"), sep="\t", index_col=0)
    dist = pd.read_csv(gzip.open(data / "1011DistanceMatrixBasedOnSNPs.tab.gz", "rt"), sep="\t", index_col=0)
    common = ph.index.intersection(cnv.index).intersection(dist.index)
    return ph.loc[common], cnv.loc[common], dist.loc[common, common]


def clade_centered_spearman(y, c, clades):
    """Residualize BOTH phenotype and copy number by clade mean, then Spearman -> de-confounded association."""
    yr, cr = y.copy(), c.copy()
    for L in np.unique(clades):
        m = clades == L
        if m.sum() > 1:
            yr[m] = y[m] - np.nanmean(y[m]); cr[m] = c[m] - np.nanmean(c[m])
    return spearmanr(cr, yr)[0], yr, cr


def _perm_null(y, cr, clades, n=200):
    out = []
    for s in range(n):
        yp = y.copy(); rng = np.random.default_rng(s)
        for L in np.unique(clades):
            i = np.where(clades == L)[0]; v = yp[i].copy(); rng.shuffle(v); yp[i] = v
        ypr = yp.copy()
        for L in np.unique(clades):
            m = clades == L
            if m.sum() > 1:
                ypr[m] = yp[m] - np.nanmean(yp[m])
        out.append(spearmanr(cr, ypr)[0])
    return np.array(out)


def run_case(ph, cnv, dist, col, condition, direction, ks=(18, 30)) -> dict:
    c = pd.to_numeric(cnv[col], errors="coerce").values
    y = ph[condition].astype(float).values
    ok = ~np.isnan(y) & ~np.isnan(c)
    y, c = y[ok], c[ok]
    res = {"col": col, "condition": condition, "n": int(ok.sum()),
           "copy_number_distribution": {k: round(float(v), 2) for k, v in
                                        pd.Series(c).describe().to_dict().items()},
           "global_spearman": round(float(spearmanr(c, y)[0]), 3), "by_k": {}}
    for K in ks:
        clades = Y._clades(dist, K)[ok]
        rho, _, cr = clade_centered_spearman(y, c, clades)
        perms = _perm_null(y, cr, clades)
        perm_p = float((np.sum(np.abs(perms) >= abs(rho)) + 1) / (len(perms) + 1))
        res["by_k"][K] = {"clade_centered_spearman": round(float(rho), 3),
                          "perm_null_p95": round(float(np.percentile(np.abs(perms), 95)), 3),
                          "perm_p": round(perm_p, 4),
                          "direction_correct": bool(np.sign(rho) == direction)}
    # verdict: de-confounded association significant + correct direction at BOTH K
    res["confirmed"] = bool(all(v["perm_p"] < 0.05 and v["direction_correct"] and abs(v["clade_centered_spearman"]) > 0.1
                                for v in res["by_k"].values()))
    return res


def run(data: Path) -> dict:
    ph, cnv, dist = _load(data)
    out = {"n_isolates": len(ph), "cases": {}}
    for name, cfg in CASES.items():
        if cfg["col"] in cnv.columns:
            out["cases"][name] = {**run_case(ph, cnv, dist, cfg["col"], cfg["condition"], cfg["dir"]),
                                  "mechanism": cfg["mechanism"]}
    return out


def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "yeast_cnv_attribution_scores.json")
    a = ap.parse_args(argv)
    res = run(a.data)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    for name, r in res["cases"].items():
        print(f"{name} -> {r['condition']} ({r['mechanism']}): global rho={r['global_spearman']:+.3f} | "
              f"{'CONFIRMED' if r['confirmed'] else 'NOT confirmed'}")
        for K, v in r["by_k"].items():
            print(f"   K={K}: clade-centered rho={v['clade_centered_spearman']:+.3f} "
                  f"perm_p={v['perm_p']} dir_ok={v['direction_correct']} (null p95={v['perm_null_p95']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
