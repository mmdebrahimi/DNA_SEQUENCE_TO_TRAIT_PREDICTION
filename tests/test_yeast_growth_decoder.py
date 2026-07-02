"""Pin the yeast decoder's de-confounding logic on synthetic data (no network/data dependency).

The REAL-data verdict lives in wiki/yeast_growth_decoder_result_2026-07-02.md; these tests pin that the
within-clade CV actually isolates within-clade signal (positive when a feature drives y inside clades;
~0 when the driver is purely between-clade structure)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.yeast_growth_decoder import _cv_r2, _r2, _within_clade_r2  # noqa: E402


def test_r2_perfect_and_null():
    y = np.array([1.0, 2, 3, 4])
    assert abs(_r2(y, y) - 1.0) < 1e-9
    assert _r2(y, np.full_like(y, y.mean())) == 0.0


def test_within_clade_detects_within_signal():
    rng = np.random.default_rng(0)
    n_per, K = 60, 4
    clades = np.repeat(np.arange(K), n_per)
    # a within-clade driver: feature g correlates with y INSIDE every clade (not just between)
    g = rng.integers(0, 2, size=n_per * K).astype(float)
    clade_mean = np.repeat(rng.normal(0, 5, K), n_per)          # big between-clade offset (structure)
    y = clade_mean + 2.0 * g + rng.normal(0, 0.3, n_per * K)     # within-clade effect of g
    X = np.column_stack([g, rng.normal(size=n_per * K)])         # g + a noise feature
    within, used = _within_clade_r2(X, y, clades, min_n=30)
    assert used == K
    assert within > 0.3                                          # within-clade signal recovered


def test_within_clade_null_on_pure_structure():
    rng = np.random.default_rng(1)
    n_per, K = 60, 4
    clades = np.repeat(np.arange(K), n_per)
    clade_mean = np.repeat(rng.normal(0, 5, K), n_per)
    y = clade_mean + rng.normal(0, 0.3, n_per * K)               # ONLY between-clade structure, no within
    X = rng.normal(size=(n_per * K, 3))                          # features unrelated to within-clade y
    within, _ = _within_clade_r2(X, y, clades, min_n=30)
    assert within < 0.1                                          # no within-clade signal to find


def test_cv_r2_runs():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(200, 5))
    y = X[:, 0] * 3 + rng.normal(0, 0.1, 200)
    assert _cv_r2(X, y) > 0.8
