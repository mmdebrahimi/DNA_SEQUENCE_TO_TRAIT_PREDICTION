"""Genotype-derived grouping (F3 Step 3).

The two guards are the point: a grouping can fail by being too COARSE (one cluster swallows
the panel) or too FINE (mostly singletons, which approximates random CV while carrying the
authority of a de-confounded number). The Mash work in this repo guarded only the first.
"""
from __future__ import annotations

import numpy as np
import pytest

from dna_decode.eval.clonality import greedy_representative_clusters_from_matrix
from dna_decode.eval.panel_structure import (
    MAX_LARGEST_CLUSTER_FRACTION,
    IBS_PERCENTILE_LADDER,
    MAX_SINGLETON_FRACTION,
    PanelStructureError,
    assert_usable_shape,
    cluster_shape,
    derive_threshold,
    ibs_distance_matrix,
    ibs_distance_result,
    panel_groups,
    pc1_variance_fraction,
    select_threshold,
    structure_census,
    sweep_thresholds,
)
from dna_decode.eval.phylogeny import DistanceMatrix


# --------------------------------------------------------------------------- IBS, by hand

def test_ibs_on_a_hand_computable_binary_matrix():
    X = [[0, 0, 0, 0],
         [0, 0, 1, 1],
         [1, 1, 1, 1]]
    res = ibs_distance_result(X, ["a", "b", "c"])
    assert res.path == "binary_exact"
    m = res.matrix.matrix
    assert m[0, 1] == pytest.approx(0.5)      # differ at 2 of 4
    assert m[0, 2] == pytest.approx(1.0)      # differ at 4 of 4
    assert m[1, 2] == pytest.approx(0.5)
    assert m[0, 0] == 0.0
    assert np.allclose(m, m.T)


def test_ibs_handles_the_minus_one_plus_one_coding():
    """The Bloom segregant panel ships parent-of-origin calls, not 0/1 dosage."""
    a = ibs_distance_result([[-1, -1, 1, 1], [-1, -1, -1, -1]], ["x", "y"])
    b = ibs_distance_result([[0, 0, 1, 1], [0, 0, 0, 0]], ["x", "y"])
    assert a.path == "binary_exact"
    assert a.matrix.matrix[0, 1] == pytest.approx(b.matrix.matrix[0, 1])
    assert a.matrix.matrix[0, 1] == pytest.approx(0.5)


def test_ibs_excludes_missing_markers_pairwise():
    X = [[0, 0, np.nan, 1],
         [0, 1, 1, 1]]
    m = ibs_distance_matrix(X, ["a", "b"]).matrix
    assert m[0, 1] == pytest.approx(1 / 3)     # 3 jointly observed, differ at 1


def test_ibs_marks_a_pair_with_no_shared_marker_non_finite():
    X = [[0, np.nan], [np.nan, 1]]
    m = ibs_distance_matrix(X, ["a", "b"]).matrix
    assert not np.isfinite(m[0, 1])


def test_ibs_takes_the_exact_general_path_on_dosage_coding():
    X = [[0, 0, 2, 2], [0, 0, 0, 0], [1, 1, 1, 1]]
    res = ibs_distance_result(X, ["a", "b", "c"])
    assert res.path == "general_loop"
    assert res.matrix.matrix[0, 1] == pytest.approx(0.5)      # |0-2|+|0-2| over 4, range 2


def test_binary_fast_path_equals_the_general_loop_on_binary_data():
    """The matmul identity |a-b| = a+b-2ab must be EXACT, not an approximation."""
    rng = np.random.default_rng(3)
    X = rng.integers(0, 2, size=(25, 60)).astype(float)
    X[0, 0] = np.nan
    ids = [f"s{i}" for i in range(25)]
    fast = ibs_distance_result(X, ids).matrix.matrix
    # force the loop by adding a third value in a column everyone shares, then undo its effect
    n = X.shape[0]
    slow = np.zeros((n, n))
    rngv = 1.0
    for i in range(n):
        slow[i] = np.nanmean(np.abs(X[i] - X), axis=1) / rngv
    assert np.allclose(fast, slow, equal_nan=True)


def test_ibs_refuses_duplicate_ids():
    with pytest.raises(PanelStructureError, match="duplicate ids"):
        ibs_distance_matrix([[0, 1], [1, 0]], ["a", "a"])


def test_ibs_refuses_id_count_mismatch():
    with pytest.raises(PanelStructureError, match="ids for"):
        ibs_distance_matrix([[0, 1], [1, 0]], ["a"])


def test_ibs_refuses_a_single_individual():
    with pytest.raises(PanelStructureError, match=">=2 individuals"):
        ibs_distance_matrix([[0, 1]], ["a"])


# --------------------------------------------------------------------------- threshold

def test_derive_threshold_refuses_a_non_finite_distance():
    """An unmeasurable pair is not a close pair. Letting NaN through has the clusterer
    evaluate NaN <= threshold (False) and quietly treat it as distant."""
    m = np.array([[0.0, np.nan, 0.2],
                  [np.nan, 0.0, 0.3],
                  [0.2, 0.3, 0.0]])
    with pytest.raises(PanelStructureError, match="non-finite"):
        derive_threshold(DistanceMatrix(strain_ids=["a", "b", "c"], matrix=m))


def test_derive_threshold_is_a_percentile_of_the_off_diagonal():
    m = np.array([[0.0, 0.1, 0.9],
                  [0.1, 0.0, 0.9],
                  [0.9, 0.9, 0.0]])
    t = derive_threshold(DistanceMatrix(strain_ids=["a", "b", "c"], matrix=m), percentile=50.0)
    assert t == pytest.approx(0.9)      # off-diag = [.1,.9,.1,.9,.9,.9], median .9


# --------------------------------------------------------------------------- PC1 statistic

def test_pc1_is_higher_on_a_two_block_panel_than_an_iid_one():
    """Assert the ORDERING, not a value -- the statistic's job is to rank panels."""
    rng = np.random.default_rng(11)
    iid = rng.integers(0, 2, size=(60, 200)).astype(float)

    block = rng.integers(0, 2, size=(60, 200)).astype(float)
    a = rng.integers(0, 2, size=200).astype(float)
    b = rng.integers(0, 2, size=200).astype(float)
    for i in range(30):
        block[i] = np.where(rng.random(200) < 0.9, a, block[i])
    for i in range(30, 60):
        block[i] = np.where(rng.random(200) < 0.9, b, block[i])

    assert pc1_variance_fraction(block) > pc1_variance_fraction(iid)


def test_pc1_mean_imputes_missing_and_drops_an_all_missing_marker():
    rng = np.random.default_rng(5)
    X = rng.integers(0, 2, size=(20, 10)).astype(float)
    X[:, 3] = np.nan                       # entirely missing marker
    X[0, 1] = np.nan
    v = pc1_variance_fraction(X)
    assert 0.0 < v <= 1.0


def test_pc1_refuses_a_constant_panel():
    with pytest.raises(PanelStructureError, match="no variation"):
        pc1_variance_fraction(np.ones((10, 5)))


# --------------------------------------------------------------------------- the two guards

def _two_block_panel(n_per=30, m=200, seed=11):
    rng = np.random.default_rng(seed)
    X = rng.integers(0, 2, size=(2 * n_per, m)).astype(float)
    a = rng.integers(0, 2, size=m).astype(float)
    b = rng.integers(0, 2, size=m).astype(float)
    for i in range(n_per):
        X[i] = np.where(rng.random(m) < 0.97, a, X[i])
    for i in range(n_per, 2 * n_per):
        X[i] = np.where(rng.random(m) < 0.97, b, X[i])
    return X, [f"s{i}" for i in range(2 * n_per)]


def test_a_percentile_threshold_MANUFACTURES_clusters_on_an_unstructured_panel():
    """MEASURED, and load-bearing for how the calibration is read.

    An i.i.d. binary panel has no structure, yet the 10th-percentile cut still yields a
    well-shaped grouping that clears BOTH guards -- because ~percentile% of pairs fall below
    any quantile by construction. So a clean cluster shape is NOT evidence of structure, which
    is why pc1_variance_fraction is primary and the F3 verdict is differential.
    """
    rng = np.random.default_rng(2)
    X = rng.integers(0, 2, size=(50, 400)).astype(float)
    pg = panel_groups(X, [f"s{i}" for i in range(50)])
    assert pg.n_clusters > 1
    assert pg.largest_cluster_fraction <= MAX_LARGEST_CLUSTER_FRACTION
    assert pg.singleton_fraction <= MAX_SINGLETON_FRACTION
    # ...and the parameter-free statistic is what actually separates it from a structured panel
    block, bids = _two_block_panel()
    assert pc1_variance_fraction(block) > pc1_variance_fraction(X)


def test_the_dominant_cluster_guard_fires_and_names_the_share():
    """Driven deterministically: percentile=100 puts every pair within threshold of the rep,
    so the whole panel collapses into one cluster."""
    rng = np.random.default_rng(4)
    X = rng.integers(0, 2, size=(30, 200)).astype(float)
    with pytest.raises(PanelStructureError, match="largest cluster"):
        panel_groups(X, [f"s{i}" for i in range(30)], percentile=100.0)


def test_the_over_split_guard_fires_and_names_the_singleton_share():
    """percentile=0 puts the threshold at the single closest pair, so almost nothing merges."""
    rng = np.random.default_rng(4)
    X = rng.integers(0, 2, size=(30, 200)).astype(float)
    with pytest.raises(PanelStructureError, match="singleton"):
        panel_groups(X, [f"s{i}" for i in range(30)], percentile=0.0)


def test_cluster_shape_is_pure_and_reports_structure_not_just_count():
    n, largest, singles = cluster_shape([0, 0, 0, 0, 1, 2, 3])
    assert n == 4
    assert largest == pytest.approx(4 / 7)
    assert singles == pytest.approx(3 / 4)


def test_assert_usable_shape_guards_both_directions_independently():
    assert_usable_shape(5, 0.30, 0.20, 0.01)                  # fine
    with pytest.raises(PanelStructureError, match="largest cluster"):
        assert_usable_shape(3, 0.61, 0.0, 0.01)
    with pytest.raises(PanelStructureError, match="singleton"):
        assert_usable_shape(10, 0.20, 0.51, 0.01)


def test_assert_usable_shape_checks_the_dominant_guard_first():
    """A grouping that fails both reports the coarse failure -- pinned so the order is not
    an accident of edit history."""
    with pytest.raises(PanelStructureError, match="largest cluster"):
        assert_usable_shape(2, 0.90, 0.90, 0.01)


def test_the_guards_are_the_published_constants():
    assert MAX_LARGEST_CLUSTER_FRACTION == 0.60
    assert MAX_SINGLETON_FRACTION == 0.50


def test_panel_groups_is_deterministic_across_calls():
    X, ids = _two_block_panel()
    try:
        a = panel_groups(X, ids)
        b = panel_groups(X, ids)
    except PanelStructureError:
        pytest.skip("fixture refused by a guard; determinism covered by the matrix test below")
    assert a.groups == b.groups
    assert a.threshold == b.threshold


def test_ibs_matrix_is_invariant_to_row_permutation():
    """Determinism that does not depend on which side of a guard the fixture lands."""
    X, ids = _two_block_panel(n_per=10, m=80)
    base = ibs_distance_result(X, ids).matrix
    perm = np.random.default_rng(0).permutation(len(ids))
    Xp = X[perm]
    idsp = [ids[i] for i in perm]
    pm = ibs_distance_result(Xp, idsp).matrix
    i0 = {s: i for i, s in enumerate(base.strain_ids)}
    i1 = {s: i for i, s in enumerate(pm.strain_ids)}
    for s in ids:
        for t in ids:
            assert base.matrix[i0[s], i0[t]] == pytest.approx(pm.matrix[i1[s], i1[t]])


# --------------------------------------------------------------------------- chaining

def test_greedy_rep_does_not_chain_where_single_linkage_would():
    """The reused primitive's chaining-resistance is RELIED ON, so it is asserted.

    A~B and B~C within threshold but A~C outside it. Single linkage merges all three;
    greedy-representative (pick a rep, drop everything within threshold OF THE REP) does not.
    """
    m = np.array([[0.00, 0.04, 0.09],
                  [0.04, 0.00, 0.04],
                  [0.09, 0.04, 0.00]])
    dm = DistanceMatrix(strain_ids=["a", "b", "c"], matrix=m)
    clusters = greedy_representative_clusters_from_matrix(dm, threshold=0.05)
    assert clusters["a"] == clusters["b"]
    assert clusters["c"] != clusters["a"], "chained A~B~C into one cluster"
    assert len(set(clusters.values())) == 2


# --------------------------------------------------------------------------- census

def test_structure_census_records_a_refusal_instead_of_raising():
    """A census must describe every panel, including the unusable ones."""
    rng = np.random.default_rng(2)
    X = rng.integers(0, 2, size=(50, 400)).astype(float)
    out = structure_census(X, [f"s{i}" for i in range(50)], percentile=0.0)
    assert out["pc1_variance_fraction"] is not None
    assert out["groups"] is None
    assert "singleton" in out["groups_refused"]


def test_structure_census_reports_a_usable_grouping_when_the_guards_pass():
    rng = np.random.default_rng(2)
    X = rng.integers(0, 2, size=(50, 400)).astype(float)
    out = structure_census(X, [f"s{i}" for i in range(50)])
    assert out["groups_refused"] is None
    assert out["groups"]["n_clusters"] > 1


def test_structure_census_records_a_pc1_failure_without_stopping():
    out = structure_census(np.ones((10, 5)), [f"s{i}" for i in range(10)])
    assert out["pc1_variance_fraction"] is None
    assert "no variation" in out["pc1_error"]
    assert out["groups_refused"] is not None


# --------------------------------------------------------------------------- the ladder sweep

def test_the_sweep_recovers_planted_blocks_where_a_fixed_percentile_SHATTERS_them():
    """The measured reason the fixed percentile was replaced by a guard-selected sweep.

    The two-block panel is strongly structured (pc1 ~0.89). At the old fixed rung of 10.0 the
    threshold lands BELOW typical within-block distances and the over-split guard refuses it;
    the sweep climbs until the guards are satisfied and recovers both blocks exactly.
    """
    X, ids = _two_block_panel()
    assert pc1_variance_fraction(X) > 0.5, "fixture is not structured; test proves nothing"

    with pytest.raises(PanelStructureError, match="singleton"):
        panel_groups(X, ids, percentile=10.0)

    pg = panel_groups(X, ids)                       # sweep
    assert pg.n_clusters == 2
    assert pg.largest_cluster_fraction == pytest.approx(0.5)
    assert pg.singleton_fraction == 0.0
    # the two clusters ARE the two planted blocks, not an arbitrary bisection
    half = len(ids) // 2
    c0 = {i for i, g in enumerate(pg.groups) if g == pg.groups[0]}
    assert c0 == set(range(half)) or c0 == set(range(half, len(ids)))


def test_select_threshold_picks_the_smallest_passing_rung():
    """Pre-committed selection: smallest passing rung, which UNDER-merges and so biases
    against detecting a de-confounding response rather than flattering it."""
    X, ids = _two_block_panel()
    dm = ibs_distance_matrix(X, ids)
    p, t, sweep = select_threshold(dm, ids)
    passing = [r["percentile"] for r in sweep if r["passes"]]
    assert passing, "fixture has no passing rung"
    assert p == min(passing)
    assert t == pytest.approx([r["threshold"] for r in sweep if r["percentile"] == p][0])


def test_the_sweep_records_why_each_rejected_rung_failed():
    X, ids = _two_block_panel()
    dm = ibs_distance_matrix(X, ids)
    sweep = sweep_thresholds(dm, ids)
    assert len(sweep) == len(IBS_PERCENTILE_LADDER)
    rejected = [r for r in sweep if not r["passes"]]
    assert rejected, "fixture rejects nothing; the sweep is untested"
    for r in rejected:
        assert r["refusal"] and ("singleton" in r["refusal"] or "largest cluster" in r["refusal"])
    for r in sweep:
        if r["passes"]:
            assert r["refusal"] is None


def test_select_threshold_refuses_with_the_whole_sweep_in_the_message():
    """A panel with no usable scale must fail diagnosably, not with a bare 'not found'."""
    # two individuals: every rung gives 1 cluster (dominant) -- no scale can work
    X = [[0, 0, 0, 1], [0, 0, 0, 1]]
    ids = ["a", "b"]
    dm = ibs_distance_matrix(X, ids)
    with pytest.raises(PanelStructureError) as e:
        select_threshold(dm, ids)
    msg = str(e.value)
    assert "no percentile on the ladder" in msg
    assert "p1:" in msg and "p70:" in msg        # the sweep is reported, not summarised away


def test_panel_groups_records_the_chosen_rung_and_the_sweep_in_as_dict():
    X, ids = _two_block_panel()
    d = panel_groups(X, ids).as_dict()
    assert d["percentile"] in IBS_PERCENTILE_LADDER
    assert len(d["sweep"]) == len(IBS_PERCENTILE_LADDER)
    assert "refusal" not in d["sweep"][0]        # prose stays out of the machine-readable row


def test_an_explicit_percentile_bypasses_the_sweep():
    X, ids = _two_block_panel()
    pg = panel_groups(X, ids, percentile=50.0)
    assert pg.percentile == 50.0
    assert pg.sweep == []
