"""Offline tests for the CV genomic-prediction arm (dna_decode/eval/genomic_prediction.py).

Synthetic genotype/phenotype with KNOWN structure: a heritable trait (y = X.beta + noise) must give a
high held-out predictive r that BEATS the label-permutation null; a pure-noise trait must NOT. Pins the
null-control logic so a positive can't be a fitting artifact. No network, no foundation model.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from dna_decode.eval.genomic_prediction import cv_model_gp, cv_ridge_gp  # noqa: E402


def _synthetic(n=200, m=50, n_causal=8, noise=0.5, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.integers(0, 2, size=(n, m)).astype(float) * 2 - 1   # {-1,+1} marker coding (segregant-like)
    beta = np.zeros(m)
    beta[:n_causal] = rng.normal(0, 1, n_causal)
    y = X @ beta + rng.normal(0, noise, n)
    return X, y


def test_heritable_trait_beats_null():
    X, y = _synthetic(noise=0.5)
    res = cv_ridge_gp(X, y, trait="heritable", n_perm=50, seed=1)
    assert res.predictive_r > 0.5           # real signal reconstructs well out-of-fold
    assert res.beats_null                    # and beats the permutation null
    assert res.null_r_mean < 0.25            # null is centred near 0


def test_pure_noise_does_not_beat_null():
    rng = np.random.default_rng(3)
    X = rng.integers(0, 2, size=(200, 50)).astype(float) * 2 - 1
    y = rng.normal(0, 1, 200)                # phenotype independent of genotype
    res = cv_ridge_gp(X, y, trait="noise", n_perm=50, seed=2)
    assert not res.beats_null                 # no genuine signal -> must not clear the null
    assert res.predictive_r < 0.3


def test_nan_phenotypes_dropped_and_dict():
    X, y = _synthetic(n=120)
    y[:10] = np.nan
    res = cv_ridge_gp(X, y, trait="t", n_perm=20)
    assert res.n == 110                        # NaN rows dropped
    d = res.as_dict()
    assert d["trait"] == "t" and "predictive_r" in d and "beats_null" in d


def test_xgboost_catches_epistasis_linear_ridge_misses():
    # y has a pure INTERACTION term (X0*X1) that a linear model cannot represent + one main effect.
    # XGBoost (nonlinear) should out-predict ridge on this epistatic trait.
    rng = np.random.default_rng(7)
    n, m = 400, 30
    X = rng.integers(0, 2, size=(n, m)).astype(float) * 2 - 1
    y = 2.0 * (X[:, 0] * X[:, 1]) + 0.7 * X[:, 2] + rng.normal(0, 0.4, n)
    xgb = cv_model_gp(X, y, model="gbm", n_perm=0, seed=1)
    rid = cv_model_gp(X, y, model="ridge", n_perm=0, seed=1)
    assert xgb.predictive_r > rid.predictive_r + 0.1   # nonlinear captures the epistasis ridge can't


def test_n_perm_zero_skips_null_no_crash():
    # n_perm=0 must skip the null cleanly (empty percentile used to crash) -- both models.
    X, y = _synthetic(n=120)
    for m in ("ridge", "gbm"):
        res = cv_model_gp(X, y, model=m, n_perm=0)
        assert res.null_r_p95 == 0.0 and res.beats_null is True
    rid = cv_ridge_gp(X, y, n_perm=0)          # the fixed path
    assert rid.null_r_p95 == 0.0


def test_too_few_samples_raises():
    X, y = _synthetic(n=4)
    with pytest.raises(ValueError):
        cv_ridge_gp(X, y, n_splits=5)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as e:
            print(f"FAIL {fn.__name__}: {e}")
    print("done")
