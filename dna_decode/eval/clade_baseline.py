"""Step 10 — Clade-only baseline classifier.

Trained on JUST one-hot clade-membership features (no sequence embeddings).
Functions as the null baseline: if the foundation-model embedding's per-clade-
held-out AUROC ≤ clade-only AUROC + 0.05, the embedding model is learning
clade signature, NOT mechanistic resistance signal.

Per post-tech-plan brainstorm M2 (Codex's M2 critique).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class CladeOnlyModel:
    """Class-conditional mean predictor over clade memberships.

    For each clade in the training set, records the empirical proportion of
    positive labels. Prediction for a strain returns its clade's positive-rate.
    Falls back to global-mean if the clade is unseen.
    """

    clade_to_positive_rate: dict[int, float]
    global_positive_rate: float


def train_clade_only_classifier(
    strain_ids: list[str],
    clade_assignments: dict[str, int],
    labels: np.ndarray,
) -> CladeOnlyModel:
    """Fit clade-only baseline on training strains.

    Args:
        strain_ids: training strain identifiers.
        clade_assignments: strain -> clade_id mapping (full pan-genome).
        labels: 0/1 array aligned with strain_ids.
    """
    if len(strain_ids) != len(labels):
        raise ValueError(f"strain_ids len {len(strain_ids)} != labels len {len(labels)}")

    clade_pos: dict[int, list[float]] = {}
    for sid, label in zip(strain_ids, labels):
        clade = clade_assignments.get(sid, -1)
        clade_pos.setdefault(clade, []).append(float(label))

    clade_to_rate = {
        clade: float(np.mean(vals)) for clade, vals in clade_pos.items()
    }
    global_rate = float(labels.mean()) if len(labels) else 0.5
    return CladeOnlyModel(
        clade_to_positive_rate=clade_to_rate,
        global_positive_rate=global_rate,
    )


def predict_clade_only(
    model: CladeOnlyModel,
    strain_ids: list[str],
    clade_assignments: dict[str, int],
) -> np.ndarray:
    """Return per-strain positive-class score under the clade-only model."""
    scores = np.empty(len(strain_ids), dtype=np.float32)
    for i, sid in enumerate(strain_ids):
        clade = clade_assignments.get(sid, -1)
        if clade in model.clade_to_positive_rate:
            scores[i] = model.clade_to_positive_rate[clade]
        else:
            scores[i] = model.global_positive_rate
    return scores


def validation_gate(
    foundation_per_clade_auroc: dict[str, float],
    clade_only_per_clade_auroc: dict[str, float],
    min_gap: float = 0.10,
    pass_fraction: float = 0.75,
) -> dict[str, object]:
    """Validation gate: foundation embedding must beat clade-only baseline.

    Phase 1 ships iff foundation model's per-clade AUROC ≥ clade-only AUROC + min_gap
    on at least `pass_fraction` of held-out clades.

    Returns dict with 'passed', 'fraction_passing', 'per_clade_gap' details.
    """
    clades = set(foundation_per_clade_auroc.keys()) & set(clade_only_per_clade_auroc.keys())
    if not clades:
        return {
            "passed": False,
            "fraction_passing": 0.0,
            "per_clade_gap": {},
            "reason": "no clades present in both per-clade dicts",
        }

    per_clade_gap = {
        c: foundation_per_clade_auroc[c] - clade_only_per_clade_auroc[c] for c in clades
    }
    passing = [c for c, gap in per_clade_gap.items() if not np.isnan(gap) and gap >= min_gap]
    fraction = len(passing) / len(clades)
    return {
        "passed": fraction >= pass_fraction,
        "fraction_passing": fraction,
        "per_clade_gap": per_clade_gap,
        "min_gap": min_gap,
        "pass_fraction_required": pass_fraction,
    }
