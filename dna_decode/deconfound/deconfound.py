"""De-confounding primitives — measure within-group signal after population/lineage structure is removed.

Canonical, data-independent implementations promoted from the 2026-07-02 research scripts. The recurring
pattern: center BOTH the phenotype and the feature by group (clade / cell-line lineage) mean so the
between-group offset (population structure) cancels, then measure the residual association. A within-group
permutation null bounds the false-positive rate.
"""
from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold


def r2(y_true, y_pred) -> float:
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def cv_r2(X, y, groups=None, alpha=10.0, seed=0) -> float:
    """5-fold CV r2. If groups given, leave-one-group-out; pooled out-of-fold r2."""
    pred = np.full(len(y), np.nan)
    if groups is not None:
        for g in np.unique(groups):
            te = groups == g; tr = ~te
            if tr.sum() < 5 or te.sum() < 1:
                continue
            m = Ridge(alpha=alpha).fit(X[tr], y[tr])
            pred[te] = m.predict(X[te])
    else:
        kf = KFold(n_splits=5, shuffle=True, random_state=seed)
        for tr, te in kf.split(X):
            m = Ridge(alpha=alpha).fit(X[tr], y[tr])
            pred[te] = m.predict(X[te])
    ok = ~np.isnan(pred)
    return r2(y[ok], pred[ok]) if ok.sum() > 10 else float("nan")


def within_group_r2(X, y, groups, min_n=30, alpha=10.0) -> tuple[float, int]:
    """Pooled out-of-fold r2 computed INSIDE each group (>=min_n), scored on GROUP-CENTERED RESIDUALS so the
    between-group mean (population structure) cancels from BOTH truth and prediction. Measures ONLY
    within-group explained variance — structure can't inflate it (the honest de-confounded metric)."""
    yr = np.full(len(y), np.nan)                       # group-centered truth
    pr = np.full(len(y), np.nan)                       # residual prediction
    used = 0
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        if len(idx) < min_n:
            continue
        used += 1
        Xg, yg = X[idx], y[idx]
        kf = KFold(n_splits=5, shuffle=True, random_state=0)
        for tr, te in kf.split(Xg):
            mu = yg[tr].mean()
            m = Ridge(alpha=alpha).fit(Xg[tr], yg[tr] - mu)   # predict within-group residual
            pr[idx[te]] = m.predict(Xg[te])
            yr[idx[te]] = yg[te] - mu                          # center test truth by TRAIN mean
    ok = ~np.isnan(pr)
    if ok.sum() <= 10:
        return float("nan"), used
    ss_res = float(np.sum((yr[ok] - pr[ok]) ** 2))
    ss_tot = float(np.sum(yr[ok] ** 2))                # residual variance (already ~zero-mean)
    return (1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0), used


def cluster_from_distance(dist, k: int) -> np.ndarray:
    """Average-linkage hierarchical clusters (K groups) from a square distance matrix (e.g. SNP distances)."""
    D = np.asarray(dist.values if hasattr(dist, "values") else dist, dtype=float)
    D = (D + D.T) / 2.0
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    return fcluster(Z, t=k, criterion="maxclust")


def group_centered_spearman(y, x, groups):
    """De-confounded CONTINUOUS association: residualize both y and x by group mean, then Spearman.
    Returns (rho, y_resid, x_resid)."""
    yr, xr = y.copy(), x.copy()
    for g in np.unique(groups):
        m = groups == g
        if m.sum() > 1:
            yr[m] = y[m] - np.nanmean(y[m]); xr[m] = x[m] - np.nanmean(x[m])
    return spearmanr(xr, yr)[0], yr, xr


def group_centered_association(y, x, groups):
    """NaN-aware de-confounded association. Returns (global_rho, within_group_rho, n_used)."""
    ok = ~np.isnan(y) & ~np.isnan(x)
    y, x, grp = y[ok], x[ok], groups[ok]
    yr, xr = y.copy(), x.copy()
    for g in np.unique(grp):
        m = grp == g
        if m.sum() > 1:
            yr[m] = y[m] - np.nanmean(y[m]); xr[m] = x[m] - np.nanmean(x[m])
    gl = spearmanr(x, y)[0]
    return (float(gl) if not np.isnan(gl) else 0.0), \
           (float(spearmanr(xr, yr)[0]) if ok.sum() > 20 else float("nan")), int(ok.sum())


def group_centered_biomarker_t(y, g, groups, min_both=4) -> dict:
    """De-confounded SINGLE-GENE (binary carrier) attribution: does the biomarker separate the phenotype
    WITHIN groups (not just because the carrier concentrates in a phenotype-extreme group)? Center y by group,
    compare carrier vs non-carrier. Returns within_group t + per-group carrier-minus-noncarrier delta."""
    yr = y.copy()
    for g_ in np.unique(groups):
        m = groups == g_
        if m.sum() > 1:
            yr[m] = y[m] - np.nanmean(y[m])
    mr, wr = yr[g == 1], yr[g == 0]
    se = np.sqrt(mr.std() ** 2 / max(len(mr), 1) + wr.std() ** 2 / max(len(wr), 1))
    t_within = float((mr.mean() - wr.mean()) / se) if se > 0 else 0.0
    per = {}
    for g_ in np.unique(groups):
        m = groups == g_; gm = g[m]; ym = y[m]
        if (gm == 1).sum() >= min_both and (gm == 0).sum() >= min_both:
            per[str(g_)] = round(float(ym[gm == 1].mean() - ym[gm == 0].mean()), 2)
    return {"within_lineage_t": round(t_within, 2), "per_lineage_delta_lfc": per}


def univariate_top(y, G, genes, n=12):
    """Top |t| binary features (mutant vs WT) associated with a continuous phenotype. G = units x features."""
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


def permutation_null(y, x_resid, groups, n=200):
    """Within-group label-permutation null for a group-centered Spearman: shuffle y INSIDE each group,
    re-center, correlate against the fixed residualized feature x_resid. Returns the array of null rhos."""
    out = []
    for s in range(n):
        yp = y.copy(); rng = np.random.default_rng(s)
        for g in np.unique(groups):
            i = np.where(groups == g)[0]; v = yp[i].copy(); rng.shuffle(v); yp[i] = v
        ypr = yp.copy()
        for g in np.unique(groups):
            m = groups == g
            if m.sum() > 1:
                ypr[m] = yp[m] - np.nanmean(yp[m])
        out.append(spearmanr(x_resid, ypr)[0])
    return np.array(out)
