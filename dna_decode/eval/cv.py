"""Step 10 — Cross-validation strategies.

Three CV strategies, in order of strictness:
1. Leave-one-strain-out (LOSO): per-strain folds. Baseline.
2. Leave-one-MLST-out (LOMO): holds out entire MLST sequence-types.
3. Leave-one-Mash-clade-out: holds out Mash/ANI-distance-based phylo clusters.
   Phase 1 primary CV strategy per post-tech-plan brainstorm M2.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np


@dataclass
class FoldResult:
    """Per-fold CV output."""

    held_out_id: str  # strain_id / mlst_id / clade_id depending on strategy
    held_out_indices: list[int]
    train_indices: list[int]
    y_true: np.ndarray  # held-out true labels
    y_score: np.ndarray  # held-out predicted scores
    n_train: int
    n_test: int


@dataclass
class CVResult:
    """CV-strategy output: per-fold + aggregate."""

    strategy: str  # 'loso' | 'lomo' | 'leave_one_clade_out'
    drug: str
    folds: list[FoldResult] = field(default_factory=list)

    @property
    def n_folds(self) -> int:
        return len(self.folds)

    @property
    def all_y_true(self) -> np.ndarray:
        if not self.folds:
            return np.array([])
        return np.concatenate([f.y_true for f in self.folds])

    @property
    def all_y_score(self) -> np.ndarray:
        if not self.folds:
            return np.array([])
        return np.concatenate([f.y_score for f in self.folds])

    @property
    def strain_ids(self) -> list[str]:
        """Held-out IDs in fold order — alignment contract for downstream paired comparisons."""
        return [f.held_out_id for f in self.folds]


def _validate_inputs(
    features: np.ndarray, labels: np.ndarray, strain_ids: list[str]
) -> None:
    if features.shape[0] != len(labels):
        raise ValueError(
            f"features.shape[0]={features.shape[0]} != len(labels)={len(labels)}"
        )
    if features.shape[0] != len(strain_ids):
        raise ValueError(
            f"features.shape[0]={features.shape[0]} != len(strain_ids)={len(strain_ids)}"
        )


def leave_one_strain_out_cv(
    features: np.ndarray,
    labels: np.ndarray,
    strain_ids: list[str],
    train_fn: Callable[[np.ndarray, np.ndarray], object],
    predict_fn: Callable[[object, np.ndarray], np.ndarray],
    drug: str = "",
) -> CVResult:
    """LOSO: one fold per strain."""
    _validate_inputs(features, labels, strain_ids)
    result = CVResult(strategy="loso", drug=drug)

    for i, strain in enumerate(strain_ids):
        train_idx = [j for j in range(len(strain_ids)) if j != i]
        held_out = [i]
        if len(set(labels[train_idx].tolist())) < 2:
            # single-class fold; record empty y_score for warning downstream
            result.folds.append(
                FoldResult(
                    held_out_id=strain,
                    held_out_indices=held_out,
                    train_indices=train_idx,
                    y_true=labels[held_out],
                    y_score=np.array([np.nan]),
                    n_train=len(train_idx),
                    n_test=1,
                )
            )
            continue
        model = train_fn(features[train_idx], labels[train_idx])
        y_score = predict_fn(model, features[held_out])
        result.folds.append(
            FoldResult(
                held_out_id=strain,
                held_out_indices=held_out,
                train_indices=train_idx,
                y_true=labels[held_out],
                y_score=y_score,
                n_train=len(train_idx),
                n_test=1,
            )
        )
    return result


def _leave_one_group_out_cv(
    features: np.ndarray,
    labels: np.ndarray,
    strain_ids: list[str],
    group_assignments: dict[str, str],
    train_fn: Callable[[np.ndarray, np.ndarray], object],
    predict_fn: Callable[[object, np.ndarray], np.ndarray],
    strategy_name: str,
    drug: str = "",
    *,
    allow_unassigned: bool = False,
) -> CVResult:
    """Generic leave-one-group-out CV. Used by LOMO + leave-one-clade-out.

    Per `plans/Stage2_N150_Prep_Plan.md` Phase A.5 defensive fix: missing-
    strain assignments raise `ValueError` by default. The `__unassigned__`
    bucket was a silent failure mode -- Stage 2's Mash-clade-out CV depends
    on every strain having a clade label, and silent bucketing would warp
    fold structure (one giant fold for all unassigned strains).

    Pass `allow_unassigned=True` to recover the original silent-bucket
    behavior; this is for the rare cases where the caller has explicitly
    audited the unassigned set.

    Raises:
        ValueError: if `allow_unassigned=False` and any strain_id is missing
            from `group_assignments`, with the count + sample of missing IDs.
    """
    _validate_inputs(features, labels, strain_ids)

    # Map each strain to its group; one fold per unique group
    if not allow_unassigned:
        missing = [s for s in strain_ids if s not in group_assignments]
        if missing:
            sample = missing[:5]
            raise ValueError(
                f"_leave_one_group_out_cv ({strategy_name}): {len(missing)} strain_ids "
                f"missing from group_assignments. Sample: {sample}. "
                f"Pass `allow_unassigned=True` to bucket them into '__unassigned__' instead "
                f"(silent failure mode -- not recommended for Stage 2)."
            )

    groups = [group_assignments.get(s, "__unassigned__") for s in strain_ids]
    unique_groups = sorted(set(groups))

    result = CVResult(strategy=strategy_name, drug=drug)
    for group in unique_groups:
        held_out = [i for i, g in enumerate(groups) if g == group]
        train_idx = [i for i, g in enumerate(groups) if g != group]
        if not train_idx:
            continue
        if len(set(labels[train_idx].tolist())) < 2:
            result.folds.append(
                FoldResult(
                    held_out_id=group,
                    held_out_indices=held_out,
                    train_indices=train_idx,
                    y_true=labels[held_out],
                    y_score=np.full(len(held_out), np.nan),
                    n_train=len(train_idx),
                    n_test=len(held_out),
                )
            )
            continue
        model = train_fn(features[train_idx], labels[train_idx])
        y_score = predict_fn(model, features[held_out])
        result.folds.append(
            FoldResult(
                held_out_id=group,
                held_out_indices=held_out,
                train_indices=train_idx,
                y_true=labels[held_out],
                y_score=y_score,
                n_train=len(train_idx),
                n_test=len(held_out),
            )
        )
    return result


def leave_one_accession_out_cv(
    features: np.ndarray,
    labels: np.ndarray,
    strain_ids: list[str],
    accession_assignments: dict[str, str],
    train_fn: Callable[[np.ndarray, np.ndarray], object],
    predict_fn: Callable[[object, np.ndarray], np.ndarray],
    drug: str = "",
    *,
    allow_unassigned: bool = False,
) -> CVResult:
    """Leave-one-accession-out CV — leakage-safe when duplicate assembly_accession
    values exist (e.g., the cipro N=147 cohort had GCA_025200635.1 registered as
    two strain_ids; leave_one_strain_out_cv would leak the same genome across
    train + held-out folds).

    Per `LESSONS_LEARNED.md` 2026-05-22 "duplicate accession in LOSO cohort =
    same-genome train/test leakage by construction": every cohort builder now
    asserts uniqueness, but legacy cohorts + cross-machine pulls may still
    surface duplicates. Use this CV strategy whenever `find_duplicate_accessions`
    on the cohort returns non-empty.

    Returns CVResult.strategy = "leave_one_accession_out".
    """
    return _leave_one_group_out_cv(
        features=features,
        labels=labels,
        strain_ids=strain_ids,
        group_assignments=accession_assignments,
        train_fn=train_fn,
        predict_fn=predict_fn,
        strategy_name="leave_one_accession_out",
        drug=drug,
        allow_unassigned=allow_unassigned,
    )


def leave_one_mlst_out_cv(
    features: np.ndarray,
    labels: np.ndarray,
    strain_ids: list[str],
    mlst_assignments: dict[str, str],
    train_fn: Callable[[np.ndarray, np.ndarray], object],
    predict_fn: Callable[[object, np.ndarray], np.ndarray],
    drug: str = "",
    *,
    allow_unassigned: bool = False,
) -> CVResult:
    """LOMO: holds out entire MLST sequence-types.

    Raises ValueError if any strain_id lacks an MLST assignment (set
    `allow_unassigned=True` for the legacy silent-bucket behavior).
    """
    return _leave_one_group_out_cv(
        features, labels, strain_ids, mlst_assignments, train_fn, predict_fn, "lomo", drug,
        allow_unassigned=allow_unassigned,
    )


def leave_one_clade_out_cv(
    features: np.ndarray,
    labels: np.ndarray,
    strain_ids: list[str],
    clade_assignments: dict[str, int],
    train_fn: Callable[[np.ndarray, np.ndarray], object],
    predict_fn: Callable[[object, np.ndarray], np.ndarray],
    drug: str = "",
    *,
    allow_unassigned: bool = False,
) -> CVResult:
    """Phase 1 primary CV: holds out entire Mash/ANI-distance clades.

    Stricter than LOMO. Catches lineage memorization that low-resolution MLST
    misses (per post-tech-plan brainstorm M2).

    Raises ValueError if any strain_id lacks a clade assignment (set
    `allow_unassigned=True` for the legacy silent-bucket behavior).
    """
    str_clades = {s: str(c) for s, c in clade_assignments.items()}
    return _leave_one_group_out_cv(
        features,
        labels,
        strain_ids,
        str_clades,
        train_fn,
        predict_fn,
        "leave_one_clade_out",
        drug,
        allow_unassigned=allow_unassigned,
    )
