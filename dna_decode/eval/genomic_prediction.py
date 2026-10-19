"""Cross-validated genomic prediction — the confound-free decoding-validation arm (Bloom-2013 yeast).

Genotype matrix X (n_individuals x n_markers) -> quantitative phenotype y. Held-out predictive skill via
K-fold CV ridge regression (the standard genomic-prediction / GBLUP-adjacent baseline). Reports predictive
Pearson r + R2 (coefficient of determination) on held-out folds, plus a LABEL-PERMUTATION null so a
positive can't be a fitting artifact. Pure sklearn/numpy; no foundation model, no network.

Substrate discipline: this arm is meaningful ONLY on a confound-free design (e.g. a single-cross segregant
panel where recombination randomizes ancestry -> no population-structure confound). On natural-population
genotypes the same code needs within-group de-confounding (see the Arabidopsis closed negative) -- pass
`groups=` to both entry points for that. The grouped path exists but its CALIBRATION status (whether it
responds differentially to real structure) is a separate question: see `wiki/f3_deconfound_calibration_*`
and `eval/regime.py`. A grouped number from an UNCALIBRATED de-confounder is not evidence of de-confounding.
"""
from __future__ import annotations

from collections.abc import Hashable, Sequence
from dataclasses import dataclass

import numpy as np


class GroupedFoldError(ValueError):
    """A `groups` argument that cannot support honest grouped CV.

    Raised rather than degraded: silently falling back to random folds would report a
    de-confounded number that was never de-confounded.
    """


@dataclass
class GPResult:
    trait: str
    n: int
    n_markers: int
    predictive_r: float          # mean held-out Pearson r (pred vs true), the genomic-prediction headline
    predictive_r2: float         # predictive_r ** 2
    r2_score: float              # held-out coefficient of determination (can be < 0)
    best_alpha: float
    null_r_mean: float           # label-permutation null: mean predictive r under shuffled y
    null_r_p95: float            # 95th percentile of the null r (positive-control threshold)
    beats_null: bool             # predictive_r > null_r_p95
    # Which fold structure produced every number above. An artifact that does not say this is
    # ambiguous between a random-CV and a grouped-CV result, which are different quantities.
    fold_scheme: str = "random"  # "random" | "grouped"
    n_groups: int | None = None  # distinct groups when fold_scheme == "grouped", else None

    def as_dict(self) -> dict:
        return {
            "trait": self.trait, "n": self.n, "n_markers": self.n_markers,
            "predictive_r": round(self.predictive_r, 4), "predictive_r2": round(self.predictive_r2, 4),
            "r2_score": round(self.r2_score, 4), "best_alpha": self.best_alpha,
            "null_r_mean": round(self.null_r_mean, 4), "null_r_p95": round(self.null_r_p95, 4),
            "beats_null": self.beats_null,
            "fold_scheme": self.fold_scheme, "n_groups": self.n_groups,
        }


def _kfold_indices(n: int, k: int, seed: int):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    return [idx[i::k] for i in range(k)]      # k interleaved test folds


def _grouped_fold_indices(groups: Sequence[Hashable], k: int, seed: int) -> list[np.ndarray]:
    """Assign WHOLE groups to folds, so no group spans train and test.

    Returns the same list-of-test-index-arrays shape `_kfold_indices` returns, which is what
    makes every downstream consumer — including the permutation null, which reuses the caller's
    `folds` object — grouped by construction rather than by a second code path.

    Deterministic given `seed`: groups are permuted once (so equal-size groups tie-break by
    seed) then stable-sorted largest-first and dropped into the lightest fold so far.
    """
    arr = np.asarray(groups, dtype=object)
    members: dict[Hashable, list[int]] = {}
    for i, g in enumerate(arr):
        members.setdefault(g, []).append(i)

    keys = list(members)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(keys))
    keys = [keys[i] for i in order]
    keys.sort(key=lambda g: -len(members[g]))          # stable: seed breaks size ties

    folds: list[list[int]] = [[] for _ in range(k)]
    for g in keys:
        lightest = min(range(k), key=lambda f: (len(folds[f]), f))
        folds[lightest].extend(members[g])
    return [np.asarray(sorted(f), dtype=int) for f in folds]


def _resolve_folds(y, groups, n_splits: int, seed: int):
    """Validate `groups` and return (folds, fold_scheme, n_groups).

    REFUSES rather than degrades on every shape that would make a grouped number dishonest.
    `groups` must already be filtered to the same rows as `y` (see `_filter_rows`).
    """
    n = len(y)
    if groups is None:
        return _kfold_indices(n, n_splits, seed), "random", None

    if len(groups) != n:
        raise GroupedFoldError(
            f"groups has length {len(groups)} but {n} rows survive NaN filtering; "
            "a misaligned groups array silently mislabels rows"
        )
    distinct = {g for g in groups}
    if len(distinct) < n_splits:
        raise GroupedFoldError(
            f"only {len(distinct)} distinct groups for {n_splits}-fold grouped CV; "
            "at least one fold would be empty"
        )
    counts = {g: sum(1 for x in groups if x == g) for g in distinct}
    biggest, biggest_n = max(counts.items(), key=lambda kv: (kv[1], str(kv[0])))
    limit = n / n_splits
    if biggest_n > limit:
        raise GroupedFoldError(
            f"group {biggest!r} holds {biggest_n}/{n} rows ({biggest_n / n:.1%}), above the "
            f"{1 / n_splits:.1%} one-fold share; grouped CV would become leave-{biggest_n / n:.0%}-out. "
            "Use a finer grouping or fewer splits."
        )
    return _grouped_fold_indices(groups, n_splits, seed), "grouped", len(distinct)


def _filter_rows(X, y, groups):
    """Drop rows with NaN y, applying the SAME mask to `groups`.

    Filtering y without filtering groups misaligns the labels against the rows — a silent
    corruption that still produces a plausible number.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = ~np.isnan(y)
    if groups is None:
        return X[keep], y[keep], None
    g = np.asarray(groups, dtype=object)
    if g.shape[0] != keep.shape[0]:
        raise GroupedFoldError(
            f"groups has length {g.shape[0]} but X/y have {keep.shape[0]} rows before NaN filtering"
        )
    return X[keep], y[keep], list(g[keep])


def _inner_split(train_idx: np.ndarray, groups, rng_seed: int):
    """Inner validation split for the alpha pick.

    With `groups`, the split is by WHOLE GROUP: a positional slice of `train_idx` can put two
    members of one group on opposite sides, so the alpha choice leaks within the fold even
    though the OUTER fold is clean.
    """
    n_in = max(1, len(train_idx) // 5)
    if groups is None:
        return train_idx[:n_in], train_idx[n_in:]

    by_group: dict[Hashable, list[int]] = {}
    for i in train_idx:
        by_group.setdefault(groups[i], []).append(int(i))
    keys = list(by_group)
    rng = np.random.default_rng(rng_seed)
    keys = [keys[i] for i in rng.permutation(len(keys))]

    va: list[int] = []
    for g in keys:
        if len(va) >= n_in and len(va) < len(train_idx):
            break
        va.extend(by_group[g])
    tr2 = [i for i in map(int, train_idx) if i not in set(va)]
    if not tr2 or not va:          # one group held everything: fall back to the positional split
        return train_idx[:n_in], train_idx[n_in:]
    return np.asarray(va, dtype=int), np.asarray(tr2, dtype=int)


def _cv_predictive_r(X, y, alphas, folds, groups=None):
    from sklearn.linear_model import Ridge

    preds = np.full(len(y), np.nan)
    chosen_alpha = []
    for f, test_idx in enumerate(folds):
        train_idx = np.concatenate([folds[j] for j in range(len(folds)) if j != f])
        Xtr, Xte = X[train_idx], X[test_idx]
        ytr = y[train_idx]
        mu, sd = Xtr.mean(0), Xtr.std(0)
        sd[sd == 0] = 1.0
        Xtr_s, Xte_s = (Xtr - mu) / sd, (Xte - mu) / sd
        ymu = ytr.mean()
        # inner 1-fold alpha pick (small grid) by train-subset validation
        best_a, best_score = alphas[0], -np.inf
        if len(alphas) > 1:
            va, tr2 = _inner_split(train_idx, groups, rng_seed=f)
            Xi, yi = (X[tr2] - mu) / sd, y[tr2]
            Xv, yv = (X[va] - mu) / sd, y[va]
            for a in alphas:
                m = Ridge(alpha=a).fit(Xi, yi - yi.mean())
                pv = m.predict(Xv) + yi.mean()
                s = -np.mean((pv - yv) ** 2)
                if s > best_score:
                    best_score, best_a = s, a
        chosen_alpha.append(best_a)
        m = Ridge(alpha=best_a).fit(Xtr_s, ytr - ymu)
        preds[test_idx] = m.predict(Xte_s) + ymu
    r = float(np.corrcoef(preds, y)[0, 1]) if np.std(preds) > 0 else 0.0
    ss_res = float(np.sum((y - preds) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2_score = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return r, r2_score, float(np.median(chosen_alpha))


def _cv_r_with(X, y, folds, fit_predict):
    """Generic held-out predictive r + R2 for any fit_predict(Xtr, ytr, Xte) -> preds callback."""
    preds = np.full(len(y), np.nan)
    for f, test_idx in enumerate(folds):
        train_idx = np.concatenate([folds[j] for j in range(len(folds)) if j != f])
        preds[test_idx] = fit_predict(X[train_idx], y[train_idx], X[test_idx])
    r = float(np.corrcoef(preds, y)[0, 1]) if np.std(preds) > 0 else 0.0
    ss_res = float(np.sum((y - preds) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return r, (1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0)


def _ridge_fit_predict(Xtr, ytr, Xte, alpha=100.0):
    from sklearn.linear_model import Ridge

    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1.0
    ymu = ytr.mean()
    m = Ridge(alpha=alpha).fit((Xtr - mu) / sd, ytr - ymu)
    return m.predict((Xte - mu) / sd) + ymu


def _gbm_fit_predict(Xtr, ytr, Xte):
    # sklearn gradient-boosted trees (no extra dep); captures nonlinearity / epistasis a linear
    # model cannot represent.
    from sklearn.ensemble import HistGradientBoostingRegressor

    m = HistGradientBoostingRegressor(max_iter=300, max_depth=4, learning_rate=0.05,
                                      l2_regularization=1.0, random_state=0)
    m.fit(Xtr, ytr)
    return m.predict(Xte)


def cv_model_gp(X, y, *, model: str = "gbm", trait: str = "trait", n_splits: int = 5,
                seed: int = 0, n_perm: int = 0, ridge_alpha: float = 100.0,
                groups: Sequence[Hashable] | None = None) -> GPResult:
    """CV genomic prediction with a pluggable model (ridge / gbm=gradient-boosted trees). `n_perm=0`
    skips the permutation null (for a head-to-head benchmark where the ridge arm already established
    the null baseline).

    `groups`: when given, folds are built by WHOLE GROUP (see `_grouped_fold_indices`). The null
    below reuses this same `folds` object, so it is grouped by construction. `groups=None` is a
    bit-for-bit no-op.
    """
    X, y, groups = _filter_rows(X, y, groups)
    n = len(y)
    if n < n_splits + 1:
        raise ValueError(f"too few samples ({n}) for {n_splits}-fold CV")
    folds, fold_scheme, n_groups = _resolve_folds(y, groups, n_splits, seed)
    fp = {"ridge": lambda a, b, c: _ridge_fit_predict(a, b, c, ridge_alpha),
          "gbm": _gbm_fit_predict}[model]
    r, r2_score = _cv_r_with(X, y, folds, fp)
    if n_perm > 0:
        rng = np.random.default_rng(seed + 1)
        nulls = np.array([_cv_r_with(X, rng.permutation(y), folds, fp)[0] for _ in range(n_perm)])
        nmean, np95 = float(nulls.mean()), float(np.percentile(nulls, 95))
    else:
        nmean, np95 = 0.0, 0.0
    return GPResult(trait=trait, n=n, n_markers=X.shape[1], predictive_r=r,
                    predictive_r2=r * r if r > 0 else 0.0, r2_score=r2_score,
                    best_alpha=(ridge_alpha if model == "ridge" else 0.0),
                    null_r_mean=nmean, null_r_p95=np95, beats_null=(r > np95 if n_perm > 0 else True),
                    fold_scheme=fold_scheme, n_groups=n_groups)


def cv_ridge_gp(X, y, *, trait: str = "trait", n_splits: int = 5,
                alphas=(1.0, 10.0, 100.0, 1000.0), seed: int = 0, n_perm: int = 100,
                groups: Sequence[Hashable] | None = None) -> GPResult:
    """K-fold CV ridge genomic prediction with a label-permutation null.

    X: (n, n_markers) genotype; y: (n,) quantitative phenotype. Rows with NaN y are dropped.

    `groups`: when given, folds are built by WHOLE GROUP, so no group spans train and test —
    the within-group de-confounding this module's docstring names as the natural-population
    requirement. The null loop below reuses the SAME `folds` object, so the null is computed
    under the grouped scheme too; a null sharing the metric's fold structure is the only kind
    that bounds the metric it is compared against. `groups=None` is a bit-for-bit no-op.
    """
    X, y, groups = _filter_rows(X, y, groups)
    n = len(y)
    if n < n_splits + 1:
        raise ValueError(f"too few samples ({n}) for {n_splits}-fold CV")

    folds, fold_scheme, n_groups = _resolve_folds(y, groups, n_splits, seed)
    r, r2_score, alpha = _cv_predictive_r(X, y, list(alphas), folds, groups)

    rng = np.random.default_rng(seed + 1)
    null_rs = []
    for _ in range(n_perm):
        yp = rng.permutation(y)
        # fixed alpha for the null (cheaper, fair) -- and the SAME folds, so a grouped run has
        # a grouped null without a second code path.
        nr, _, _ = _cv_predictive_r(X, yp, [alpha], folds, groups)
        null_rs.append(nr)
    null_rs = np.array(null_rs)
    if n_perm > 0:
        nmean, np95 = float(null_rs.mean()), float(np.percentile(null_rs, 95))
    else:  # n_perm=0: skip the null (empty percentile would crash)
        nmean, np95 = 0.0, 0.0

    return GPResult(
        trait=trait, n=n, n_markers=X.shape[1], predictive_r=r, predictive_r2=r * r if r > 0 else 0.0,
        r2_score=r2_score, best_alpha=alpha, null_r_mean=nmean, null_r_p95=np95,
        beats_null=(r > np95 if n_perm > 0 else True),
        fold_scheme=fold_scheme, n_groups=n_groups,
    )
