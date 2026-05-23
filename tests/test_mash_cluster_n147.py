"""Tests for scripts/mash_cluster_n147.py pure-logic helpers.

Mash + Docker not available on this laptop -- functional smoke-test happens
on the Precision 7780 when Codex runs the script. These tests pin the
threshold-sweep scoring + selection logic against synthetic distance matrices.
"""
from __future__ import annotations

import math

import numpy as np
import pytest


# ---- score_threshold ----


def test_score_threshold_3_balanced_clades():
    """3 clades of 4 strains each -- satisfies all criteria."""
    from scripts.mash_cluster_n147 import score_threshold

    strain_ids = [f"s_{i}" for i in range(12)]
    # Block-diagonal distance matrix: intra-block = 0.01, inter-block = 0.10
    matrix = np.full((12, 12), 0.10, dtype=np.float32)
    for block in range(3):
        for i in range(block * 4, (block + 1) * 4):
            for j in range(block * 4, (block + 1) * 4):
                matrix[i, j] = 0.0 if i == j else 0.01
    assignments = {sid: i // 4 for i, sid in enumerate(strain_ids)}
    score = score_threshold(matrix, assignments, strain_ids, threshold=0.05)
    assert score.n_clades == 3
    assert score.max_clade_fraction == pytest.approx(4 / 12)
    assert score.satisfies_min_clades
    assert score.satisfies_max_fraction
    assert score.fully_satisfied
    # variance_ratio = intra-mean (0.01) / inter-mean (0.10) = 0.10
    assert score.variance_ratio == pytest.approx(0.10, abs=1e-6)


def test_score_threshold_dominant_single_clade_fails_max_fraction():
    """One clade of 10 + 2 singletons -- max-clade fraction 0.83 > 0.60."""
    from scripts.mash_cluster_n147 import score_threshold

    strain_ids = [f"s_{i}" for i in range(12)]
    matrix = np.full((12, 12), 0.10, dtype=np.float32)
    for i in range(12):
        matrix[i, i] = 0.0
    assignments = {sid: 0 if i < 10 else (1 if i == 10 else 2) for i, sid in enumerate(strain_ids)}
    score = score_threshold(matrix, assignments, strain_ids, threshold=0.05)
    assert score.n_clades == 3
    assert score.max_clade_fraction == pytest.approx(10 / 12)
    assert score.satisfies_min_clades
    assert not score.satisfies_max_fraction
    assert not score.fully_satisfied


def test_score_threshold_single_clade_fails_min_clades():
    """All in one clade -- n_clades = 1, fails min_clades=3."""
    from scripts.mash_cluster_n147 import score_threshold

    strain_ids = [f"s_{i}" for i in range(12)]
    matrix = np.full((12, 12), 0.001, dtype=np.float32)
    for i in range(12):
        matrix[i, i] = 0.0
    assignments = {sid: 0 for sid in strain_ids}
    score = score_threshold(matrix, assignments, strain_ids, threshold=0.05)
    assert score.n_clades == 1
    assert not score.satisfies_min_clades
    assert not score.fully_satisfied
    assert math.isnan(score.variance_ratio)  # no inter-clade pairs


def test_score_threshold_empty_cohort():
    from scripts.mash_cluster_n147 import score_threshold
    score = score_threshold(np.zeros((0, 0)), {}, [], threshold=0.05)
    assert score.n_clades == 0
    assert not score.fully_satisfied


# ---- pick_best_threshold ----


def test_pick_best_threshold_picks_lowest_variance_ratio_among_qualifying():
    """Among qualifying thresholds, pick the one with lowest variance_ratio."""
    from scripts.mash_cluster_n147 import pick_best_threshold

    strain_ids = [f"s_{i}" for i in range(12)]
    # Use a block-diagonal matrix; clustering at different thresholds yields
    # different clade counts.
    matrix = np.full((12, 12), 0.20, dtype=np.float32)
    # 4 tight blocks of 3 strains each
    for block in range(4):
        for i in range(block * 3, (block + 1) * 3):
            for j in range(block * 3, (block + 1) * 3):
                matrix[i, j] = 0.0 if i == j else 0.01

    def cluster_fn(t: float) -> dict[str, int]:
        # Simulate union-find on pairs < threshold
        parent = list(range(12))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for i in range(12):
            for j in range(i + 1, 12):
                if matrix[i, j] < t:
                    ra, rb = find(i), find(j)
                    if ra != rb:
                        parent[ra] = rb
        roots: dict[int, int] = {}
        out: dict[str, int] = {}
        next_id = 0
        for i, sid in enumerate(strain_ids):
            r = find(i)
            if r not in roots:
                roots[r] = next_id
                next_id += 1
            out[sid] = roots[r]
        return out

    chosen, scores, fellback = pick_best_threshold(
        matrix, strain_ids, cluster_fn,
        candidates=(0.005, 0.02, 0.05, 0.30),
    )
    # 0.005 -> all singletons (12 clades, max_frac 1/12) -> qualifies; variance_ratio = 0 (no intra) -> NaN
    # 0.02 -> 4 clades of 3 -> max_frac 3/12 = 0.25 < 0.60 -> qualifies; var ratio low
    # 0.05 -> same as 0.02 (no inter-block pairs in [0.02, 0.05]) -> qualifies, same var ratio
    # 0.30 -> 1 clade of 12 -> fails min_clades
    assert chosen in (0.02, 0.05)  # both qualify with same variance_ratio (tied lowest)
    assert not fellback


def test_pick_best_threshold_falls_back_when_none_qualify():
    """If no threshold satisfies criteria -> return fallback."""
    from scripts.mash_cluster_n147 import pick_best_threshold

    strain_ids = [f"s_{i}" for i in range(12)]
    matrix = np.full((12, 12), 0.5, dtype=np.float32)  # all far apart
    for i in range(12):
        matrix[i, i] = 0.0

    def cluster_fn(t: float) -> dict[str, int]:
        # No pairs ever cluster -> 12 singletons
        return {sid: i for i, sid in enumerate(strain_ids)}

    chosen, scores, fellback = pick_best_threshold(
        matrix, strain_ids, cluster_fn,
        candidates=(0.02, 0.05),
        fallback=0.05,
    )
    # 12 singletons: max_clade_frac = 1/12 (passes max_fraction), n_clades=12 (passes min_clades)
    # BUT variance_ratio is NaN (no intra-clade pairs) -> doesn't qualify
    assert chosen == 0.05
    assert fellback


# ---- per_clade_label_balance ----


def test_per_clade_label_balance_counts_R_S_unknown():
    from scripts.mash_cluster_n147 import per_clade_label_balance
    assignments = {"a": 0, "b": 0, "c": 0, "d": 1, "e": 1, "f": 2}
    labels = {"a": 1, "b": 0, "c": None, "d": 1, "e": 1, "f": 0}
    balance = per_clade_label_balance(assignments, labels)
    assert balance[0] == {"R": 1, "S": 1, "unknown": 1, "n": 3}
    assert balance[1] == {"R": 2, "S": 0, "unknown": 0, "n": 2}
    assert balance[2] == {"R": 0, "S": 1, "unknown": 0, "n": 1}


def test_per_clade_label_balance_handles_missing_strain():
    """Strain in assignments but not in labels -> counted as unknown."""
    from scripts.mash_cluster_n147 import per_clade_label_balance
    assignments = {"a": 0, "b": 0}
    labels = {"a": 1}  # "b" missing
    balance = per_clade_label_balance(assignments, labels)
    assert balance[0] == {"R": 1, "S": 0, "unknown": 1, "n": 2}
