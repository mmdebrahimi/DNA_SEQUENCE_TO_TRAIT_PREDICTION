"""Step 9 — Baseline XGBoost classifiers on frozen foundation-model embeddings.

Per-drug binary R/S classification. Mean-pools per-gene embeddings within a
strain into a single feature vector, then trains an XGBoost classifier with
sigmoid probability calibration (CalibratedClassifierCV) for AUROC + Brier
honesty.

Used by Wave 3 + Wave 6 leaderboard. Foundation-model variants (Evo,
DNABERT-2, NT, GENA-LM) all feed into the same classifier code — Step 7's
wrapper interface is uniform.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(frozen=True)
class XGBParams:
    """Hyperparameters for the XGBoost classifier head."""

    n_estimators: int = 200
    max_depth: int = 6
    learning_rate: float = 0.1
    subsample: float = 0.8
    colsample_bytree: float = 0.8
    eval_metric: str = "auc"
    random_state: int = 42


@dataclass
class TrainedClassifier:
    """Trained XGBoost + optional sigmoid-calibration wrapper.

    `model` is the sklearn-API-compatible object (XGBClassifier or
    CalibratedClassifierCV depending on `calibrated`). `feature_dim` records the
    expected input dimension for shape sanity at predict time.
    """

    model: Any
    drug_name: str
    feature_dim: int
    calibrated: bool
    params: XGBParams = field(default_factory=XGBParams)


class ClassifierTrainingError(Exception):
    """Training failed (insufficient samples, single-class labels, etc.)."""


MIN_TRAINING_SAMPLES = 10


def aggregate_strain_features(
    gene_embeddings: np.ndarray,
    aggregation: str = "mean",
) -> np.ndarray:
    """Pool per-gene embeddings within a strain into a single feature vector.

    Args:
        gene_embeddings: shape (n_genes, embedding_dim).
        aggregation: 'mean' (default) | 'max' | 'mean+max' (concatenates both).

    Returns:
        1-D array of shape (embedding_dim,) for 'mean' / 'max'; (2*embedding_dim,)
        for 'mean+max'.
    """
    if gene_embeddings.ndim != 2:
        raise ValueError(
            f"gene_embeddings must be 2-D (n_genes, dim); got shape {gene_embeddings.shape}"
        )
    # NaN-aware aggregation handles bulk_get NaN-fills from cache
    if aggregation == "mean":
        return np.nanmean(gene_embeddings, axis=0).astype(np.float32)
    if aggregation == "max":
        return np.nanmax(gene_embeddings, axis=0).astype(np.float32)
    if aggregation == "mean+max":
        return np.concatenate(
            [
                np.nanmean(gene_embeddings, axis=0),
                np.nanmax(gene_embeddings, axis=0),
            ]
        ).astype(np.float32)
    raise ValueError(f"Unknown aggregation: {aggregation!r}")


def train_xgboost_classifier(
    X: np.ndarray,
    y: np.ndarray,
    drug_name: str,
    params: XGBParams | None = None,
    calibrate: bool = True,
) -> TrainedClassifier:
    """Train XGBoost on mean-pooled strain embeddings → binary R/S labels.

    Args:
        X: (n_strains, embedding_dim) feature matrix.
        y: (n_strains,) binary labels (0 = susceptible, 1 = resistant).
        drug_name: for downstream attribution + reporting.
        params: XGBParams overrides; defaults if None.
        calibrate: wrap in sklearn CalibratedClassifierCV with sigmoid (Platt
            scaling) for honest probability outputs.

    Returns:
        TrainedClassifier holding the fitted model.
    """
    params = params or XGBParams()

    if X.shape[0] < MIN_TRAINING_SAMPLES:
        raise ClassifierTrainingError(
            f"Need >= {MIN_TRAINING_SAMPLES} training samples; got {X.shape[0]}"
        )
    if X.shape[0] != len(y):
        raise ClassifierTrainingError(
            f"X.shape[0]={X.shape[0]} != len(y)={len(y)}"
        )
    if len(set(y.tolist())) < 2:
        raise ClassifierTrainingError(
            f"Drug {drug_name!r}: training labels are single-class "
            f"({set(y.tolist())}); need both 0 and 1 present"
        )

    # Import xgboost + sklearn lazily so the module imports clean without them
    try:
        import xgboost as xgb
    except ImportError as e:
        raise ClassifierTrainingError(
            "xgboost not installed; run `uv sync` to install Phase 1 deps"
        ) from e

    base = xgb.XGBClassifier(
        n_estimators=params.n_estimators,
        max_depth=params.max_depth,
        learning_rate=params.learning_rate,
        subsample=params.subsample,
        colsample_bytree=params.colsample_bytree,
        eval_metric=params.eval_metric,
        random_state=params.random_state,
        verbosity=0,
        use_label_encoder=False,
    )

    if calibrate:
        try:
            from sklearn.calibration import CalibratedClassifierCV
        except ImportError as e:
            raise ClassifierTrainingError(
                "scikit-learn not installed; run `uv sync`"
            ) from e
        # 3-fold internal CV for calibration; falls back gracefully if n is small
        cv_folds = min(3, max(2, (y == 1).sum(), (y == 0).sum()))
        model = CalibratedClassifierCV(base, cv=cv_folds, method="sigmoid")
        model.fit(X, y)
    else:
        base.fit(X, y)
        model = base

    return TrainedClassifier(
        model=model,
        drug_name=drug_name,
        feature_dim=X.shape[1],
        calibrated=calibrate,
        params=params,
    )


def predict_proba(classifier: TrainedClassifier, X: np.ndarray) -> np.ndarray:
    """Return per-sample P(resistant=1) under the trained classifier.

    Tolerates NaN rows from cache.bulk_get(fill_missing=True) — replaces them
    with the global mean so xgboost doesn't crash on NaN inputs.
    """
    if X.shape[1] != classifier.feature_dim:
        raise ValueError(
            f"X.shape[1]={X.shape[1]} != trained feature_dim={classifier.feature_dim}"
        )
    # NaN-fill with column means (per-batch); avoids xgboost NaN errors
    X = X.copy()
    nan_mask = np.isnan(X)
    if nan_mask.any():
        col_means = np.nanmean(X, axis=0)
        # Fall back to zero if a whole column is NaN
        col_means = np.where(np.isnan(col_means), 0.0, col_means)
        X[nan_mask] = np.take(col_means, np.where(nan_mask)[1])

    proba = classifier.model.predict_proba(X)
    return proba[:, 1].astype(np.float32)


def feature_importance(classifier: TrainedClassifier) -> np.ndarray:
    """Return per-feature importance scores (length = feature_dim).

    For calibrated models, averages importance across the wrapped per-fold
    XGBClassifier instances. For non-calibrated, returns the XGBClassifier's
    feature_importances_ directly.
    """
    model = classifier.model
    if classifier.calibrated:
        # CalibratedClassifierCV.calibrated_classifiers_ holds the per-fold fits
        importances = []
        for fitted in getattr(model, "calibrated_classifiers_", []):
            base_estimator = getattr(fitted, "estimator", None) or getattr(
                fitted, "base_estimator", None
            )
            if base_estimator is not None and hasattr(base_estimator, "feature_importances_"):
                importances.append(base_estimator.feature_importances_)
        if not importances:
            return np.zeros(classifier.feature_dim, dtype=np.float32)
        return np.mean(np.stack(importances), axis=0).astype(np.float32)
    return model.feature_importances_.astype(np.float32)
