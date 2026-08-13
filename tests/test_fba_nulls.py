"""Margin-preserving null — pure, no solver, no model.

The property that makes this null worth building is EXACT margin preservation. If a shuffle can move a
row or column sum, the null is measuring base rate again and the whole point is lost, so the margin
invariant is tested directly and on random matrices rather than assumed.
"""
from __future__ import annotations

import random

from dna_decode.fba.nulls import (
    calls_from_matrix,
    curveball_shuffle,
    margin_preserving_null,
    margins,
    matrix_from_calls,
)

CONDS = ("a", "b", "c", "d")
GENES = ["g1", "g2", "g3"]


def test_matrix_round_trips_through_the_calls_shape():
    calls = {"a": {"g1": True, "g2": False, "g3": True},
             "b": {"g1": False, "g2": True, "g3": False},
             "c": {"g1": True, "g2": True, "g3": False},
             "d": {"g1": False, "g2": False, "g3": False}}
    mat = matrix_from_calls(GENES, CONDS, calls)
    assert mat == [[1, 0, 1, 0], [0, 1, 1, 0], [1, 0, 0, 0]]
    assert calls_from_matrix(GENES, CONDS, mat) == calls


def test_a_missing_cell_reads_as_not_essential_not_as_an_error():
    """Deletion scripts omit cells when wild-type growth is zero; that must not crash the null."""
    mat = matrix_from_calls(GENES, CONDS, {"a": {"g1": True}})
    assert mat == [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]


def test_the_shuffle_preserves_BOTH_margins_exactly_on_random_matrices():
    """THE property. If a trade could move a row or column sum, the null degenerates back into the
    rate-matched one it was built to replace."""
    rng = random.Random(11)
    for trial in range(25):
        n_rows, n_cols = rng.randint(2, 12), rng.randint(2, 10)
        mat = [[rng.randint(0, 1) for _ in range(n_cols)] for _ in range(n_rows)]
        before = margins(mat)
        after = margins(curveball_shuffle(mat, seed=trial))
        assert after == before, f"margins moved on trial {trial}"


def test_the_shuffle_does_not_mutate_its_input():
    mat = [[1, 0, 0, 1], [0, 1, 1, 0]]
    snapshot = [row[:] for row in mat]
    curveball_shuffle(mat, seed=3)
    assert mat == snapshot


def test_the_shuffle_actually_moves_cells_when_it_can():
    """A null that returns the input unchanged would trivially 'fail to beat' the observation."""
    mat = [[1, 1, 0, 0], [0, 0, 1, 1]]
    got = curveball_shuffle(mat, n_trades=200, seed=5)
    assert got != mat
    assert margins(got) == margins(mat)


def test_a_matrix_with_no_swappable_pair_is_returned_intact():
    """Identical rows have no discordant columns -- there is nothing to trade, and that is not a bug."""
    mat = [[1, 0], [1, 0], [1, 0]]
    assert curveball_shuffle(mat, n_trades=50, seed=1) == mat


def test_a_single_row_is_returned_intact():
    assert curveball_shuffle([[1, 0, 1]], seed=0) == [[1, 0, 1]]


def test_the_null_is_deterministic_for_a_given_seed():
    calls = {c: {g: (hash((g, c)) % 3 == 0) for g in GENES} for c in CONDS}

    def score(cs):
        return sum(1 for c in CONDS for g in GENES if cs[c][g]) / (len(CONDS) * len(GENES))

    a = margin_preserving_null(GENES, CONDS, calls, score, n_draws=10, seed0=4)
    b = margin_preserving_null(GENES, CONDS, calls, score, n_draws=10, seed0=4)
    assert a["mean"] == b["mean"] and a["max"] == b["max"]


def test_a_score_that_depends_only_on_the_TOTAL_cannot_beat_this_null():
    """The sharpest statement of what changed. A predictor scored purely on how MANY cells it calls
    essential gets an identical score on every margin-preserving draw -- so this null gives it zero
    credit, where the rate-matched null would have given it the benefit of matching the base rate."""
    calls = {"a": {"g1": True, "g2": False, "g3": True},
             "b": {"g1": False, "g2": True, "g3": False},
             "c": {"g1": True, "g2": True, "g3": False},
             "d": {"g1": False, "g2": False, "g3": False}}

    def total_only(cs):
        return float(sum(1 for c in CONDS for g in GENES if cs[c][g]))

    got = margin_preserving_null(GENES, CONDS, calls, total_only, n_draws=25, seed0=0)
    assert got["mean"] == got["max"] == 5.0        # every draw identical -- no base-rate credit
