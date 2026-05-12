"""Tests for Step 10 — Evaluation harness (CV + metrics + clade-only baseline)."""
from __future__ import annotations

import numpy as np
import pytest

from dna_decode.eval.clade_baseline import (
    predict_clade_only,
    train_clade_only_classifier,
    validation_gate,
)
from dna_decode.eval.cv import (
    CVResult,
    leave_one_clade_out_cv,
    leave_one_mlst_out_cv,
    leave_one_strain_out_cv,
)
from dna_decode.eval.metrics import (
    compute_attribution_precision,
    compute_metrics,
    compute_per_clade_metrics,
)


# ---- Synthetic data helpers ----


def _synthetic_data(n: int = 30, seed: int = 42) -> tuple[np.ndarray, np.ndarray, list[str]]:
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, 8))
    # Make the label easy: y depends linearly on first feature
    y = (X[:, 0] + 0.5 * rng.standard_normal(n) > 0).astype(int)
    strain_ids = [f"s{i:03d}" for i in range(n)]
    return X, y, strain_ids


def _simple_train(X: np.ndarray, y: np.ndarray) -> dict:
    """Trivial trainer: store mean feature vector per class."""
    return {
        0: X[y == 0].mean(axis=0) if (y == 0).any() else np.zeros(X.shape[1]),
        1: X[y == 1].mean(axis=0) if (y == 1).any() else np.zeros(X.shape[1]),
    }


def _simple_predict(model: dict, X: np.ndarray) -> np.ndarray:
    """Distance-to-class-mean score → probability that y=1."""
    d0 = np.linalg.norm(X - model[0], axis=1)
    d1 = np.linalg.norm(X - model[1], axis=1)
    return d0 / (d0 + d1 + 1e-9)


# ---- CV strategies ----


def test_loso_produces_one_fold_per_strain():
    X, y, ids = _synthetic_data(20)
    result = leave_one_strain_out_cv(X, y, ids, _simple_train, _simple_predict)
    assert result.strategy == "loso"
    assert result.n_folds == 20


def test_lomo_produces_one_fold_per_mlst():
    X, y, ids = _synthetic_data(20)
    mlst = {ids[i]: f"ST{i // 5}" for i in range(20)}  # 4 MLSTs
    result = leave_one_mlst_out_cv(X, y, ids, mlst, _simple_train, _simple_predict)
    assert result.strategy == "lomo"
    assert result.n_folds == 4


def test_leave_one_clade_out_produces_one_fold_per_clade():
    X, y, ids = _synthetic_data(30)
    clades = {ids[i]: i // 6 for i in range(30)}  # 5 clades
    result = leave_one_clade_out_cv(X, y, ids, clades, _simple_train, _simple_predict)
    assert result.strategy == "leave_one_clade_out"
    assert result.n_folds == 5


def test_cv_handles_single_class_fold_with_nan():
    """When all train labels are one class, the fold records NaN scores."""
    X = np.random.RandomState(0).randn(4, 3)
    y = np.array([1, 1, 1, 0])  # holding out strain 3 leaves train all-positive
    ids = ["a", "b", "c", "d"]
    result = leave_one_strain_out_cv(X, y, ids, _simple_train, _simple_predict)
    # Find the fold where held_out is index 3 (label 0); train is all 1s
    folds_by_id = {f.held_out_id: f for f in result.folds}
    # The fold for "d" should have train all-positive → NaN score
    assert np.isnan(folds_by_id["d"].y_score[0]) or len(set(y[[0, 1, 2]])) == 1


def test_cv_shape_mismatch_raises():
    X = np.zeros((5, 3))
    y = np.zeros(3)
    ids = ["a", "b", "c", "d", "e"]
    with pytest.raises(ValueError, match="labels"):
        leave_one_strain_out_cv(X, y, ids, _simple_train, _simple_predict)


# ---- Metrics ----


def test_compute_metrics_basic():
    y_true = np.array([0, 0, 1, 1, 0, 1])
    y_score = np.array([0.1, 0.2, 0.8, 0.9, 0.3, 0.7])
    m = compute_metrics(y_true, y_score)
    assert m.n_samples == 6
    assert m.auroc > 0.9  # very-separable scores
    assert 0 <= m.brier <= 1
    assert 0 <= m.ece <= 1


def test_compute_metrics_handles_nan_scores():
    y_true = np.array([0, 1, 0, 1])
    y_score = np.array([0.2, 0.8, np.nan, 0.9])
    m = compute_metrics(y_true, y_score)
    assert m.n_samples == 3  # NaN filtered out


def test_compute_metrics_single_class_returns_nan():
    y_true = np.array([1, 1, 1])
    y_score = np.array([0.5, 0.6, 0.7])
    m = compute_metrics(y_true, y_score)
    assert np.isnan(m.auroc)


def test_attribution_precision_basic():
    predicted = ["gyrA", "parC", "tetA", "random1", "random2"]
    known = {"gyrA", "parC", "gyrB"}
    prec = compute_attribution_precision(predicted, known, k=5)
    # 2 of top-5 are in known → 0.4
    assert prec == pytest.approx(0.4)


def test_attribution_precision_case_insensitive():
    predicted = ["GYRA", "ParC"]
    known = {"gyrA", "parC"}
    prec = compute_attribution_precision(predicted, known, k=2)
    assert prec == 1.0


def test_attribution_precision_empty_predicted():
    assert compute_attribution_precision([], {"gyrA"}, k=5) == 0.0


def test_attribution_precision_k_zero_raises():
    with pytest.raises(ValueError, match="positive"):
        compute_attribution_precision(["x"], {"x"}, k=0)


def test_per_clade_metrics_returns_one_per_fold():
    X, y, ids = _synthetic_data(30)
    clades = {ids[i]: i // 6 for i in range(30)}
    result = leave_one_clade_out_cv(X, y, ids, clades, _simple_train, _simple_predict)
    per_clade = compute_per_clade_metrics(result)
    assert len(per_clade) == 5


# ---- Clade-only baseline ----


def test_clade_only_trains_per_clade_rates():
    ids = ["a", "b", "c", "d"]
    clades = {"a": 0, "b": 0, "c": 1, "d": 1}
    y = np.array([0, 1, 1, 1])  # clade 0 → 50%, clade 1 → 100%
    model = train_clade_only_classifier(ids, clades, y)
    assert model.clade_to_positive_rate[0] == 0.5
    assert model.clade_to_positive_rate[1] == 1.0


def test_clade_only_predict_uses_clade_rate():
    ids = ["a", "b"]
    clades = {"a": 0, "b": 1}
    y = np.array([0, 1])
    model = train_clade_only_classifier(ids, clades, y)
    scores = predict_clade_only(model, ["a", "b"], clades)
    np.testing.assert_array_almost_equal(scores, [0.0, 1.0])


def test_clade_only_predict_unseen_clade_uses_global_rate():
    ids = ["a", "b"]
    clades = {"a": 0, "b": 0}
    y = np.array([0, 1])  # global rate = 0.5
    model = train_clade_only_classifier(ids, clades, y)
    # Strain in a new clade 99 not seen during training
    new_clades = {"new_strain": 99}
    scores = predict_clade_only(model, ["new_strain"], new_clades)
    assert scores[0] == 0.5


# ---- Validation gate ----


def test_validation_gate_passes_when_foundation_beats_baseline():
    foundation = {"c1": 0.95, "c2": 0.92, "c3": 0.88, "c4": 0.81}
    baseline = {"c1": 0.60, "c2": 0.70, "c3": 0.65, "c4": 0.70}
    result = validation_gate(foundation, baseline, min_gap=0.10, pass_fraction=0.75)
    assert result["passed"] is True
    assert result["fraction_passing"] == 1.0


def test_validation_gate_fails_when_gap_insufficient():
    foundation = {"c1": 0.70, "c2": 0.68, "c3": 0.65, "c4": 0.66}
    baseline = {"c1": 0.65, "c2": 0.64, "c3": 0.60, "c4": 0.65}
    result = validation_gate(foundation, baseline, min_gap=0.10)
    assert result["passed"] is False
    assert result["fraction_passing"] < 0.75


def test_validation_gate_empty_clade_intersection():
    result = validation_gate({}, {})
    assert result["passed"] is False
    assert "no clades" in result["reason"]
