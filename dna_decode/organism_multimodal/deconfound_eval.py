"""The de-confounding evaluator: POOLED vs WITHIN-population cross-individual accuracy.

The load-bearing discipline (population-structure = the clonality analogue): a
grouping variable the genotype tracks (here: population) inflates pooled r. We report
each group's WITHIN-group Spearman vs that group's OWN null, never a single pooled
number. Pure numpy/scipy; no I/O.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import spearmanr


@dataclass(frozen=True)
class DeconfoundResult:
    pooled_rho: float
    within_rho_mean: float          # sample-size-weighted mean of per-group rho
    per_group_rho: dict[str, float]
    per_group_null_rho: dict[str, float]   # permutation null (mean |rho| over shuffles)
    per_group_n: dict[str, int]
    inflation: float                # pooled - within_mean (the confound size)


def evaluate(y_true: np.ndarray, y_pred: np.ndarray, groups: np.ndarray,
             n_perm: int = 200, seed: int = 0) -> DeconfoundResult:
    rng = np.random.default_rng(seed)
    pooled = float(spearmanr(y_true, y_pred).statistic)
    per_g, per_null, per_n = {}, {}, {}
    for g in sorted(set(groups.tolist())):
        m = groups == g
        n = int(m.sum())
        per_n[g] = n
        if n < 5 or np.std(y_pred[m]) == 0 or np.std(y_true[m]) == 0:
            per_g[g] = float("nan"); per_null[g] = float("nan"); continue
        per_g[g] = float(spearmanr(y_true[m], y_pred[m]).statistic)
        nulls = [abs(spearmanr(y_true[m], rng.permutation(y_pred[m])).statistic)
                 for _ in range(n_perm)]
        per_null[g] = float(np.mean(nulls))
    valid = [(per_g[g], per_n[g]) for g in per_g if not np.isnan(per_g[g])]
    within = (sum(r * n for r, n in valid) / sum(n for _, n in valid)) if valid else float("nan")
    return DeconfoundResult(pooled, within, per_g, per_null, per_n, pooled - within)
