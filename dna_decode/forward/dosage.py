"""Forward-cell DOSAGE head — turn a rank-score into a CALIBRATED MAGNITUDE prediction with honest intervals.

The forward variant-effect cell (`variant_effect.predict_effect`) produces a score whose RANK correlates
with the measured effect (Spearman). A decoder should say more than "ranks low" — it should predict the
actual effect MAGNITUDE with a calibrated uncertainty. This module maps any method's arbitrary-scale score
to the measured-effect scale (monotone isotonic calibrator) and wraps it in a split-conformal prediction
interval with a PRE-REGISTERED held-out coverage target.

Reuses the split-conformal definition from `scripts/hiv_quantitative_calibration._conformal_q` (J2's
Family-B MIC-calibration helper) — `conformal_q` here is the SAME finite-sample formula, verified byte-equal
in the tests (non-duplication: same math, in-package for clean deps).

THE LOAD-BEARING HONESTY RAIL (J2's lesson): split-conformal coverage holds even for a USELESS model — the
interval just widens to the MARGINAL distribution. So coverage alone does NOT prove the score is
informative. This module ALSO reports `interval_narrowing = 1 - q/marginal_q`: how much the score's
conditioning shrinks the interval vs a no-features (predict-the-mean) baseline. A calibrated + INFORMATIVE
dosage head needs BOTH nominal coverage AND meaningful narrowing.
"""
from __future__ import annotations

from dataclasses import dataclass


def conformal_q(abs_res, alpha: float) -> float:
    """Finite-sample split-conformal quantile: the ceil((m+1)(1-alpha))/m empirical quantile of |residuals|.
    Byte-equal to hiv_quantitative_calibration._conformal_q (the shared canonical definition; asserted in
    tests). abs_res is a sequence of absolute calibration residuals; alpha = miscoverage (1 - coverage)."""
    import numpy as np
    r = np.asarray(abs_res, dtype=float)
    m = len(r)
    if m == 0:
        return float("nan")
    import math
    k = math.ceil((m + 1) * (1 - alpha))
    if k > m:
        return float(np.max(r))
    return float(np.sort(r)[k - 1])


def _spearman_sign(x, y) -> int:
    import numpy as np
    x = np.asarray(x, float); y = np.asarray(y, float)

    def rank(v):
        order = np.argsort(v, kind="mergesort")
        r = np.empty(len(v)); r[order] = np.arange(len(v))
        return r
    rx, ry = rank(x), rank(y)
    c = np.corrcoef(rx, ry)[0, 1] if len(x) > 1 else 0.0
    return 1 if (c >= 0 or np.isnan(c)) else -1


@dataclass
class DosageResult:
    coverage: float           # held-out fraction of test y inside [lo, hi] (the honest coverage number)
    target: float             # pre-registered coverage target
    halfwidth: float          # conformal q — interval is point +/- q on the measured-effect scale
    marginal_halfwidth: float # conformal q of the predict-the-mean baseline (no features)
    interval_narrowing: float # 1 - halfwidth/marginal_halfwidth (informativeness; >0 = score narrows it)
    point_spearman: float     # rank corr of the calibrated point estimate vs test y (sanity)
    point_rmse: float
    n_fit: int
    n_calib: int
    n_test: int


def _isotonic_fit(fit_x, fit_y):
    from sklearn.isotonic import IsotonicRegression
    sign = _spearman_sign(fit_x, fit_y)
    iso = IsotonicRegression(increasing=(sign >= 0), out_of_bounds="clip")
    iso.fit(fit_x, fit_y)
    return iso


def dosage_intervals(fit_x, fit_y, calib_x, calib_y, test_x, coverage: float = 0.8):
    """Split-conformal dosage intervals for `test_x`.
      - fit isotonic score->effect on the FIT split,
      - conformal q from |calib_y - iso(calib_x)| at alpha = 1 - coverage on the CALIB split,
      - per-test: point = iso(test_x), interval = point +/- q.
    Returns (point[], lo[], hi[], q, marginal_q). marginal_q uses the predict-the-mean baseline (no
    features) so the caller can measure how much the score NARROWS the interval."""
    import numpy as np
    iso = _isotonic_fit(fit_x, fit_y)
    alpha = 1 - coverage
    q = conformal_q(np.abs(np.asarray(calib_y, float) - iso.predict(calib_x)), alpha)
    # marginal (no-feature) baseline: predict the FIT mean; conformal q on the calib set
    mean_y = float(np.mean(fit_y))
    marginal_q = conformal_q(np.abs(np.asarray(calib_y, float) - mean_y), alpha)
    point = iso.predict(test_x)
    return point, point - q, point + q, q, marginal_q


def evaluate_dosage(fit_x, fit_y, calib_x, calib_y, test_x, test_y, coverage: float = 0.8) -> DosageResult:
    """End-to-end: fit+calibrate on fit/calib, evaluate held-out coverage + informativeness on test."""
    import numpy as np
    test_y = np.asarray(test_y, float)
    point, lo, hi, q, marginal_q = dosage_intervals(fit_x, fit_y, calib_x, calib_y, test_x, coverage)
    covered = float(np.mean((test_y >= lo) & (test_y <= hi)))
    narrowing = float(1 - q / marginal_q) if (marginal_q and not np.isnan(marginal_q)) else float("nan")
    sp = _spearman_sign  # reuse ranker
    # point Spearman vs test y
    def spearman(a, b):
        ra, rb = _rank(a), _rank(b)
        return float(np.corrcoef(ra, rb)[0, 1]) if len(a) > 1 else float("nan")
    return DosageResult(
        coverage=round(covered, 4), target=coverage, halfwidth=round(float(q), 4),
        marginal_halfwidth=round(float(marginal_q), 4), interval_narrowing=round(narrowing, 4),
        point_spearman=round(spearman(point, test_y), 4),
        point_rmse=round(float(np.sqrt(np.mean((point - test_y) ** 2))), 4),
        n_fit=len(fit_x), n_calib=len(calib_x), n_test=len(test_x))


def _rank(v):
    import numpy as np
    v = np.asarray(v, float)
    order = np.argsort(v, kind="mergesort")
    r = np.empty(len(v)); r[order] = np.arange(len(v))
    return r
