"""Unit tests for the clonality clustering + cluster-weighted metric module.

Pins the brainstorm-mandated behaviors: greedy-representative dedup resists a chain
that single-linkage would merge (C2), mixed-label clones are DISCORDANT not
majority-voted (C1), the cluster-weighted metric collapses a clone to one vote, and
the Wilson CI is well-behaved at the small/extreme N the layer produces (C3).
"""
from __future__ import annotations

import numpy as np
import pytest

from dna_decode.eval.clonality import (
    cluster_class,
    cluster_members,
    cluster_weighted_confusion,
    effective_lineage_n,
    greedy_representative_clusters_from_matrix,
    wilson_ci,
)
from dna_decode.eval.phylogeny import DistanceMatrix


def _dm(sids, m):
    return DistanceMatrix(strain_ids=list(sids), matrix=np.array(m, dtype=np.float32))


# --------------------------------------------------------------------------- #
# greedy-representative clustering (C2 — chaining-resistant)
# --------------------------------------------------------------------------- #
def test_greedy_resists_chain_single_linkage_would_merge():
    # A~B (.01), B~C (.01), A!~C (.09). Single-linkage union-find merges all three;
    # greedy-representative keeps A's cluster (A,B) separate from C.
    dm = _dm(["A", "B", "C"], [[0, .01, .09], [.01, 0, .01], [.09, .01, 0]])
    clusters = greedy_representative_clusters_from_matrix(dm, threshold=0.02)
    assert clusters == {"A": 0, "B": 0, "C": 1}


def test_greedy_is_deterministic_via_sorted_ids():
    # Representative selection walks SORTED ids, so input order can't change the result.
    m = [[0, .005, .005], [.005, 0, .005], [.005, .005, 0]]
    a = greedy_representative_clusters_from_matrix(_dm(["z", "a", "m"], m), 0.01)
    b = greedy_representative_clusters_from_matrix(_dm(["m", "z", "a"], m), 0.01)
    assert a == b
    # all within threshold of the first sorted rep ("a") -> one cluster
    assert set(a.values()) == {0}


def test_greedy_empty_matrix():
    assert greedy_representative_clusters_from_matrix(_dm([], np.zeros((0, 0))), 0.02) == {}


def test_greedy_all_singletons_when_far():
    dm = _dm(["A", "B"], [[0, 0.5], [0.5, 0]])
    assert greedy_representative_clusters_from_matrix(dm, 0.02) == {"A": 0, "B": 1}


# --------------------------------------------------------------------------- #
# cluster_class (C1 — DISCORDANT, never majority-voted)
# --------------------------------------------------------------------------- #
def test_cluster_class_pure():
    assert cluster_class([1, 1, 1]) == "R"
    assert cluster_class(["S", "S"]) == "S"


def test_cluster_class_mixed_is_discordant_not_majority():
    # 2 R + 1 S would majority-vote to R; the honest answer is DISCORDANT.
    assert cluster_class([1, 1, 0]) == "DISCORDANT"


def test_cluster_class_accepts_both_label_conventions():
    assert cluster_class([1, "R"]) == "R"
    assert cluster_class([0, "s"]) == "S"


def test_cluster_class_empty_raises():
    with pytest.raises(ValueError):
        cluster_class([])


def test_cluster_class_bad_label_raises():
    with pytest.raises(ValueError):
        cluster_class([2])


# --------------------------------------------------------------------------- #
# cluster_weighted_confusion
# --------------------------------------------------------------------------- #
def test_weighted_collapses_clone_to_one_vote():
    # 5-isolate R clone (all pred R) + 1 distinct S lineage (pred S) -> tp=1, tn=1.
    preds = {f"r{i}": "R" for i in range(5)} | {"s0": "S"}
    labels = {f"r{i}": 1 for i in range(5)} | {"s0": 0}
    clusters = {f"r{i}": 0 for i in range(5)} | {"s0": 1}
    c = cluster_weighted_confusion(preds, labels, clusters)
    assert c["tp"] == 1 and c["tn"] == 1 and c["n_scored"] == 2
    assert c["n_clusters_R"] == 1 and c["n_clusters_S"] == 1 and c["sens"] == 1.0


def test_weighted_excludes_discordant():
    preds = {"a": "R", "b": "R"}
    labels = {"a": 1, "b": 0}  # same cluster, mixed labels
    clusters = {"a": 0, "b": 0}
    c = cluster_weighted_confusion(preds, labels, clusters)
    assert c["n_discordant"] == 1 and c["n_scored"] == 0
    assert c["sens"] is None and c["spec"] is None


def test_weighted_clone_tie_prediction_abstains():
    # A 2-member R clone whose members disagree (R/S) -> tie -> cluster abstains, not scored.
    preds = {"a": "R", "b": "S"}
    labels = {"a": 1, "b": 1}
    clusters = {"a": 0, "b": 0}
    c = cluster_weighted_confusion(preds, labels, clusters)
    assert c["n_cluster_abstain"] == 1 and c["n_scored"] == 0 and c["n_clusters_R"] == 1


def test_weighted_raw_vs_weighted_divergence():
    # 10 R isolates in ONE clone (pred R) vs 1 distinct R lineage (pred S, a miss).
    # Raw sens = 10/11 = 0.909; weighted sens = 1/2 = 0.5 (clone collapses to 1 vote).
    preds = {f"r{i}": "R" for i in range(10)} | {"rx": "S"}
    labels = {f"r{i}": 1 for i in range(10)} | {"rx": 1}
    clusters = {f"r{i}": 0 for i in range(10)} | {"rx": 1}
    c = cluster_weighted_confusion(preds, labels, clusters)
    assert c["sens"] == 0.5 and c["tp"] == 1 and c["fn"] == 1


# --------------------------------------------------------------------------- #
# wilson_ci (C3)
# --------------------------------------------------------------------------- #
def test_wilson_perfect_small_n():
    lo, hi = wilson_ci(8, 8)
    assert 0.6 < lo < 0.72 and hi == 1.0  # NOT a degenerate [1.0, 1.0]


def test_wilson_zero_trials():
    assert wilson_ci(0, 0) == (0.0, 1.0)


def test_wilson_half():
    lo, hi = wilson_ci(5, 10)
    assert lo < 0.5 < hi and 0.0 <= lo and hi <= 1.0


def test_wilson_clamped_to_unit_interval():
    lo, hi = wilson_ci(0, 3)
    assert lo == 0.0 and 0.0 < hi < 1.0


# --------------------------------------------------------------------------- #
# effective_lineage_n + cluster_members
# --------------------------------------------------------------------------- #
def test_effective_lineage_n_counts_same_label_clusters():
    labels = {"r0": 1, "r1": 1, "rx": 1, "s0": 0}
    clusters = {"r0": 0, "r1": 0, "rx": 1, "s0": 2}  # 2 R lineages, 1 S lineage
    assert effective_lineage_n(clusters, labels, "R") == 2
    assert effective_lineage_n(clusters, labels, "S") == 1


def test_effective_lineage_n_excludes_discordant():
    labels = {"a": 1, "b": 0}
    clusters = {"a": 0, "b": 0}  # one discordant cluster
    assert effective_lineage_n(clusters, labels, "R") == 0
    assert effective_lineage_n(clusters, labels, "S") == 0


def test_cluster_members_sorted():
    assert cluster_members({"b": 0, "a": 0, "c": 1}) == {0: ["a", "b"], 1: ["c"]}
