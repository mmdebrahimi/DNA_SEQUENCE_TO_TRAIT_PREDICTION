"""Step 10 — Evaluation metrics.

Classification: AUROC, AUPRC, F1, accuracy at 0.5, Brier score, ECE.
Attribution: top-K overlap with known resistance loci.
Per-clade reporting: aggregate CVResult by held-out clade for the
clade-specific failure-mode check.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Iterable

import numpy as np

from dna_decode.eval.cv import CVResult, FoldResult


@dataclass(frozen=True)
class Metrics:
    """Classification metrics for a single drug × CV strategy."""

    auroc: float
    auprc: float
    f1: float
    accuracy_at_05: float
    brier: float
    ece: float  # expected calibration error (10-bin)
    n_samples: int

    def to_dict(self) -> dict[str, float | int]:
        return {
            "auroc": self.auroc,
            "auprc": self.auprc,
            "f1": self.f1,
            "accuracy_at_05": self.accuracy_at_05,
            "brier": self.brier,
            "ece": self.ece,
            "n_samples": self.n_samples,
        }


def _nan_safe_filter(
    y_true: np.ndarray, y_score: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Drop NaN-scored entries (single-class CV folds emit NaN)."""
    mask = ~np.isnan(y_score)
    return y_true[mask], y_score[mask]


def compute_metrics(y_true: np.ndarray, y_score: np.ndarray) -> Metrics:
    """Compute classification metrics. Tolerates NaN-scored single-class folds."""
    from sklearn.metrics import (
        average_precision_score,
        brier_score_loss,
        f1_score,
        roc_auc_score,
    )

    y_true_f, y_score_f = _nan_safe_filter(y_true, y_score)
    n = len(y_true_f)

    if n == 0 or len(set(y_true_f.tolist())) < 2:
        # Cannot compute AUROC / AUPRC without both classes present
        return Metrics(
            auroc=np.nan,
            auprc=np.nan,
            f1=np.nan,
            accuracy_at_05=np.nan,
            brier=np.nan,
            ece=np.nan,
            n_samples=n,
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        auroc = float(roc_auc_score(y_true_f, y_score_f))
        auprc = float(average_precision_score(y_true_f, y_score_f))

    y_pred = (y_score_f >= 0.5).astype(int)
    f1 = float(f1_score(y_true_f, y_pred, zero_division=0))
    acc = float((y_pred == y_true_f).mean())
    brier = float(brier_score_loss(y_true_f, y_score_f))
    ece = _expected_calibration_error(y_true_f, y_score_f, n_bins=10)

    return Metrics(
        auroc=auroc,
        auprc=auprc,
        f1=f1,
        accuracy_at_05=acc,
        brier=brier,
        ece=ece,
        n_samples=n,
    )


def _expected_calibration_error(
    y_true: np.ndarray, y_score: np.ndarray, n_bins: int = 10
) -> float:
    """10-bin Expected Calibration Error.

    Bins predictions by confidence; weighted-average abs(predicted_prob - true_prob).
    """
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    n = len(y_true)
    for i in range(n_bins):
        in_bin = (y_score >= bins[i]) & (y_score < bins[i + 1])
        if i == n_bins - 1:  # include 1.0 in last bin
            in_bin = (y_score >= bins[i]) & (y_score <= bins[i + 1])
        count = in_bin.sum()
        if count == 0:
            continue
        bin_acc = float(y_true[in_bin].mean())
        bin_conf = float(y_score[in_bin].mean())
        ece += (count / n) * abs(bin_conf - bin_acc)
    return float(ece)


def compute_attribution_precision(
    predicted_loci: list[str], known_loci: Iterable[str], k: int
) -> float:
    """Top-K precision of predicted loci against literature known-loci set.

    Args:
        predicted_loci: gene symbols ranked by attribution score (descending).
        known_loci: literature-curated resistance-locus set for the drug.
        k: top-K cutoff.

    Returns:
        |top-K ∩ known| / k
    """
    if k <= 0:
        raise ValueError(f"k must be positive, got {k}")
    if not predicted_loci:
        return 0.0
    top_k = {locus.lower() for locus in predicted_loci[:k]}
    known_lower = {locus.lower() for locus in known_loci}
    return len(top_k & known_lower) / k


def compute_per_clade_metrics(
    cv_result: CVResult,
) -> dict[str, Metrics]:
    """Aggregate CV metrics per held-out clade/group.

    Catches the "0.9 on 8 clades, 0.45 on 2 clades" failure mode (per
    post-tech-plan brainstorm M2 — Codex flagged that aggregate AUROC
    hides clade-specific failure).
    """
    return {
        fold.held_out_id: compute_metrics(fold.y_true, fold.y_score)
        for fold in cv_result.folds
    }


def compute_within_mlst_permutation_control(
    cv_result: CVResult,
    labels: np.ndarray,
    mlst_assignments: dict[str, str],
    seed: int = 0,
    n_permutations: int = 100,
) -> dict[str, float]:
    """Within-MLST label-shuffle control.

    Shuffles resistance labels within each MLST clade; recomputes CV metrics.
    If the original model's AUROC remains high under shuffle, it's exploiting
    clade signature rather than mechanistic resistance signal.

    Returns dict with 'observed_auroc', 'permuted_auroc_mean',
    'permuted_auroc_std', 'p_value'.
    """
    observed = compute_metrics(cv_result.all_y_true, cv_result.all_y_score).auroc
    rng = np.random.default_rng(seed)

    permuted_aurocs: list[float] = []
    # Group strain indices by MLST (need a strain_id list — derive from folds)
    # For simplicity, we permute the GLOBAL labels within each MLST group using
    # the fold's held_out_indices as a proxy for index space.
    all_indices = np.arange(len(labels))
    mlst_groups: dict[str, list[int]] = {}
    # NOTE: this helper expects the caller to pass labels indexed-aligned with
    # the strain_ids used to construct mlst_assignments. The proxy here groups
    # by MLST id derived from the assignments dict.
    for idx, strain in enumerate(mlst_assignments.keys()):
        mlst_id = mlst_assignments[strain]
        mlst_groups.setdefault(mlst_id, []).append(idx)

    for _ in range(n_permutations):
        permuted = labels.copy()
        for indices in mlst_groups.values():
            if len(indices) < 2:
                continue
            sub = np.array(indices)
            perm = rng.permutation(sub)
            permuted[sub] = labels[perm]
        # Recompute predicted scores under permutation is non-trivial; for now,
        # we report only the label-shuffle distribution shape as a sanity check.
        # A full implementation re-trains under permuted labels (expensive).
        # Phase 1 records the distribution metadata; full re-train is Phase 2.
        permuted_aurocs.append(
            float(compute_metrics(permuted, cv_result.all_y_score).auroc)
        )

    perm_arr = np.array([a for a in permuted_aurocs if not np.isnan(a)])
    if len(perm_arr) == 0:
        return {
            "observed_auroc": observed,
            "permuted_auroc_mean": float("nan"),
            "permuted_auroc_std": float("nan"),
            "p_value": float("nan"),
        }
    p_value = float((perm_arr >= observed).mean())
    return {
        "observed_auroc": float(observed),
        "permuted_auroc_mean": float(perm_arr.mean()),
        "permuted_auroc_std": float(perm_arr.std()),
        "p_value": p_value,
    }
