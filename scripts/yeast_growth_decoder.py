"""Lineage-DE-CONFOUNDED yeast growth decoder on the 1002 Yeast Genomes substrate.

Answers the ONLY question that matters given the project's history (embeddings 0-for-4, all structure-
learners): is the genotype->growth signal MECHANISM or POPULATION STRUCTURE? Uses the project's own 3-part
test (the embedding-niche bar), NOT a naive global CV r2.

Inputs (fetched free from 1002genomes.u-strasbg.fr; paths via --data):
  phenoMatrix_35ConditionsNormalizedByYPD.tab.gz   isolate x 35 lab growth conditions (quantitative)
  genesMatrix_PresenceAbsence.tab.gz               isolate x 7796 genes (0/1/NA)  [genotype features]
  1011DistanceMatrixBasedOnSNPs.tab.gz             isolate x isolate SNP distance [structure -> clades]

Per condition, 4 numbers:
  naive_r2        global 5-fold Ridge on gene-PA           (structure-POLLUTED; the misleading number)
  clade_only_r2   global 5-fold Ridge on clade one-hot     (pure structure baseline)
  loco_r2         leave-one-clade-out Ridge on gene-PA     (generalize ACROSS structure)
  within_clade_r2 pooled within-clade 5-fold Ridge on gene-PA  (the DE-CONFOUNDED mechanistic signal)

Verdict per condition (frozen thresholds):
  MECHANISM   within_clade_r2 >= 0.05  (real signal inside clades, can't be structure)
  STRUCTURE   within_clade_r2 < 0.05 AND naive_r2 >= 0.10 (predicts, but only the structure part)
  WEAK        neither (little signal either way)
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

REPO = Path(__file__).resolve().parent.parent
WITHIN_CLADE_MECH = 0.05
NAIVE_STRUCTURE = 0.10


def _load(data_dir: Path):
    ph = pd.read_csv(gzip.open(data_dir / "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz", "rt"),
                     sep="\t", index_col=0)
    pa = pd.read_csv(gzip.open(data_dir / "genesMatrix_PresenceAbsence.tab.gz", "rt"),
                     sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce")
    dist = pd.read_csv(gzip.open(data_dir / "1011DistanceMatrixBasedOnSNPs.tab.gz", "rt"),
                       sep="\t", index_col=0)
    common = ph.index.intersection(pa.index).intersection(dist.index)
    return ph.loc[common], pa.loc[common], dist.loc[common, common]


# De-confounding primitives PROMOTED to the installable package (2026-07-02); this script is a thin runner.
from dna_decode.deconfound import (  # noqa: E402
    cluster_from_distance as _clades,
    cv_r2 as _cv_r2,
    r2 as _r2,
    within_group_r2 as _within_clade_r2,
)


def run(data_dir: Path, conditions: list[str] | None = None, k_clades: int = 18) -> dict:
    ph, pa, dist = _load(data_dir)
    clades = _clades(dist, k_clades)
    # genotype matrix: impute NaN->0 (gene absent), drop near-invariant genes
    X = pa.values.astype("float32")
    X = np.nan_to_num(X, nan=0.0)
    var = X.var(axis=0)
    X = X[:, var > 1e-4]
    onehot = pd.get_dummies(pd.Series(clades)).values.astype("float32")
    conds = conditions or list(ph.columns)
    rows = []
    for c in conds:
        y = ph[c].astype(float).values
        ok = ~np.isnan(y)
        if ok.sum() < 100:
            continue
        Xo, yo, ch, oh = X[ok], y[ok], clades[ok], onehot[ok]
        naive = _cv_r2(Xo, yo)
        clade_only = _cv_r2(oh, yo)
        loco = _cv_r2(Xo, yo, groups=ch)
        within, n_clades_used = _within_clade_r2(Xo, yo, ch)
        verdict = ("MECHANISM" if within >= WITHIN_CLADE_MECH
                   else "STRUCTURE" if naive >= NAIVE_STRUCTURE
                   else "WEAK")
        rows.append({"condition": c, "n": int(ok.sum()), "naive_r2": round(naive, 3),
                     "clade_only_r2": round(clade_only, 3), "loco_r2": round(loco, 3),
                     "within_clade_r2": round(within, 3), "n_clades_used": n_clades_used,
                     "verdict": verdict})
    return {"n_isolates": len(ph), "k_clades": k_clades, "conditions": rows}


def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, required=True, help="dir with the 3 gz matrices")
    ap.add_argument("--conditions", type=str, default=None, help="comma-separated; default all 35")
    ap.add_argument("--k-clades", type=int, default=18)
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "yeast_growth_decoder_scores.json")
    a = ap.parse_args(argv)
    conds = a.conditions.split(",") if a.conditions else None
    res = run(a.data, conds, a.k_clades)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"{'condition':22} {'n':>4} {'naive':>7} {'clade':>7} {'loco':>7} {'within':>7}  verdict")
    for r in res["conditions"]:
        print(f"{r['condition']:22} {r['n']:>4} {r['naive_r2']:>7} {r['clade_only_r2']:>7} "
              f"{r['loco_r2']:>7} {r['within_clade_r2']:>7}  {r['verdict']}")
    mech = sum(1 for r in res["conditions"] if r["verdict"] == "MECHANISM")
    print(f"\nMECHANISM: {mech}/{len(res['conditions'])} conditions; artifact -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
