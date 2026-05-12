"""Tests for Step 9 — Baseline XGBoost classifiers."""
from __future__ import annotations

import numpy as np
import pytest

# xgboost / sklearn aren't in the global Python; tests skip if absent
xgboost = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")

from dna_decode.models.classifiers import (  # noqa: E402
    ClassifierTrainingError,
    MIN_TRAINING_SAMPLES,
    TrainedClassifier,
    XGBParams,
    aggregate_strain_features,
    feature_importance,
    predict_proba,
    train_xgboost_classifier,
)


# ---- aggregate_strain_features ----


def test_aggregate_mean_returns_per_dim_mean():
    gene_embeds = np.array([[1.0, 2.0], [3.0, 4.0]])
    out = aggregate_strain_features(gene_embeds, "mean")
    assert out.shape == (2,)
    np.testing.assert_array_almost_equal(out, [2.0, 3.0])


def test_aggregate_max_returns_per_dim_max():
    gene_embeds = np.array([[1.0, 5.0], [3.0, 4.0]])
    out = aggregate_strain_features(gene_embeds, "max")
    np.testing.assert_array_almost_equal(out, [3.0, 5.0])


def test_aggregate_mean_plus_max_concatenates():
    gene_embeds = np.array([[1.0, 2.0], [3.0, 4.0]])
    out = aggregate_strain_features(gene_embeds, "mean+max")
    assert out.shape == (4,)
    np.testing.assert_array_almost_equal(out, [2.0, 3.0, 3.0, 4.0])


def test_aggregate_nan_aware():
    gene_embeds = np.array([[1.0, 2.0], [np.nan, 4.0]])
    out = aggregate_strain_features(gene_embeds, "mean")
    np.testing.assert_array_almost_equal(out, [1.0, 3.0])  # NaN row excluded per-col


def test_aggregate_invalid_shape_raises():
    with pytest.raises(ValueError, match="2-D"):
        aggregate_strain_features(np.array([1.0, 2.0]), "mean")


def test_aggregate_unknown_aggregation_raises():
    with pytest.raises(ValueError, match="Unknown"):
        aggregate_strain_features(np.zeros((2, 2)), "median")


# ---- train_xgboost_classifier ----


def _synthetic_separable(n: int = 50, dim: int = 8, seed: int = 0):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, dim))
    y = (X[:, 0] + 0.3 * rng.standard_normal(n) > 0).astype(int)
    return X, y


def test_train_basic_round_trip():
    X, y = _synthetic_separable(50)
    clf = train_xgboost_classifier(X, y, drug_name="cipro")
    assert isinstance(clf, TrainedClassifier)
    assert clf.drug_name == "cipro"
    assert clf.feature_dim == 8
    assert clf.calibrated is True


def test_train_uncalibrated():
    X, y = _synthetic_separable(50)
    clf = train_xgboost_classifier(X, y, drug_name="cipro", calibrate=False)
    assert clf.calibrated is False


def test_train_below_min_samples_raises():
    X = np.zeros((5, 4))
    y = np.zeros(5, dtype=int)
    with pytest.raises(ClassifierTrainingError, match="samples"):
        train_xgboost_classifier(X, y, drug_name="cipro")


def test_train_shape_mismatch_raises():
    X = np.zeros((20, 4))
    y = np.zeros(15, dtype=int)
    with pytest.raises(ClassifierTrainingError, match="!="):
        train_xgboost_classifier(X, y, drug_name="cipro")


def test_train_single_class_raises():
    X = np.zeros((20, 4))
    y = np.ones(20, dtype=int)
    with pytest.raises(ClassifierTrainingError, match="single-class"):
        train_xgboost_classifier(X, y, drug_name="cipro")


def test_train_custom_params():
    X, y = _synthetic_separable(50)
    custom = XGBParams(n_estimators=50, max_depth=3)
    clf = train_xgboost_classifier(X, y, drug_name="cipro", params=custom)
    assert clf.params.n_estimators == 50
    assert clf.params.max_depth == 3


# ---- predict_proba ----


def test_predict_proba_returns_probabilities():
    X, y = _synthetic_separable(60)
    clf = train_xgboost_classifier(X, y, drug_name="cipro")
    proba = predict_proba(clf, X)
    assert proba.shape == (60,)
    assert ((proba >= 0) & (proba <= 1)).all()


def test_predict_proba_separable_data_above_random():
    """Easy synthetic signal should give predict_proba clearly above 0.5 for y=1."""
    X, y = _synthetic_separable(80, seed=1)
    clf = train_xgboost_classifier(X, y, drug_name="cipro")
    proba = predict_proba(clf, X)
    # On training data the model should be well above random
    correct = ((proba >= 0.5) == y.astype(bool)).mean()
    assert correct > 0.85


def test_predict_proba_shape_mismatch_raises():
    X, y = _synthetic_separable(30)
    clf = train_xgboost_classifier(X, y, drug_name="cipro")
    with pytest.raises(ValueError, match="feature_dim"):
        predict_proba(clf, np.zeros((5, clf.feature_dim + 1)))


def test_predict_proba_handles_nan_inputs():
    """Cache.bulk_get returns NaN for missing pairs; predict_proba should not crash."""
    X, y = _synthetic_separable(40)
    clf = train_xgboost_classifier(X, y, drug_name="cipro")
    X_with_nan = X[:5].copy()
    X_with_nan[0, 0] = np.nan
    proba = predict_proba(clf, X_with_nan)
    assert not np.any(np.isnan(proba))


# ---- feature_importance ----


def test_feature_importance_shape_matches_feature_dim():
    X, y = _synthetic_separable(50, dim=12)
    clf = train_xgboost_classifier(X, y, drug_name="cipro")
    imp = feature_importance(clf)
    assert imp.shape == (12,)


def test_feature_importance_uncalibrated():
    X, y = _synthetic_separable(50, dim=8)
    clf = train_xgboost_classifier(X, y, drug_name="cipro", calibrate=False)
    imp = feature_importance(clf)
    assert imp.shape == (8,)
    assert (imp >= 0).all()


def test_feature_importance_first_feature_dominant_in_separable_data():
    """Synthetic data has signal only in feature 0; importance should reflect that."""
    X, y = _synthetic_separable(80, dim=8)
    clf = train_xgboost_classifier(X, y, drug_name="cipro", calibrate=False)
    imp = feature_importance(clf)
    # Feature 0 should be among the top-3 by importance
    top3 = np.argsort(imp)[::-1][:3]
    assert 0 in top3.tolist()
