"""Tests for the CE baseline panel's pure parts: the two new nulls, the position replay, the bar.

The panel's trustworthiness rests on three things that can each fail silently, so each is pinned:

1. The two new nulls are PROBABILITY models. A table whose rows do not sum to 1 still produces
   finite, plausible cross-entropies.
2. `rebuild_positions` must reproduce the scorer's SHARED-RNG draw. A per-window re-seed would give
   every window identical positions and would still reconcile if the gate checked only window 0 --
   so the distinctness of consecutive windows is asserted directly.
3. `compute_verdict` must be reachable at all four outcomes and must be ORDERED. A bar that cannot
   return its own bad outcomes is not a bar.
"""
from __future__ import annotations

import math
import sys

import pytest

from scripts.glm_baseline_panel_ce import (
    PanelError,
    compute_verdict,
    comp_nll,
    fit_composition,
    fit_site_independent,
    rebuild_positions,
    site_nll,
)


def _cell(delta: float, resolved: bool) -> dict:
    return {"paired_delta_nats": delta, "resolved": resolved}


def _cells(comp: tuple, site: tuple, ctx: tuple, arm: str = "left_only") -> dict:
    return {
        f"{arm}|shuffled_null_ecoli": _cell(*comp),
        f"{arm}|site_independent_ecoli": _cell(*site),
        f"{arm}|markov_ecoli": _cell(*ctx),
    }


# --------------------------------------------------------------------------- the module stays pure


def test_importing_the_panel_does_not_load_torch():
    """The whole point of the design: it runs beside the test suite without holding a model."""
    assert "torch" not in sys.modules


# --------------------------------------------------------------------------- the two new nulls


def test_site_independent_rows_are_distributions():
    table = fit_site_independent(["ACGT" * 2, "TTTT" * 2], length=8)
    assert len(table) == 8
    for row in table:
        assert math.isclose(sum(row), 1.0, abs_tol=1e-12)
        assert all(p > 0 for p in row)  # add-one smoothing: no zero-probability base


def test_composition_is_a_distribution():
    probs = fit_composition(["ACGT", "AAAA"])
    assert math.isclose(sum(probs), 1.0, abs_tol=1e-12)
    assert probs[0] > probs[1]  # A is the most common base in that pair


def test_site_independent_uses_the_position_and_composition_does_not():
    """The two nulls differ exactly in whether position is known -- the panel's whole premise.

    The test must compare the SAME BASE at two positions where it has different frequencies. An
    earlier version compared two positions that were each fully determined (always A, always C), so
    both scored -log(1/2) and the assertion tied for a reason that had nothing to do with position.

    Here A is common at position 0 (3 of 3 sequences) and absent at position 1, so a position-aware
    model must score A differently at the two, while a composition-only model cannot.
    """
    seqs = ["AC", "AC", "AG"]
    table = fit_site_independent(seqs, length=2)
    probs = fit_composition(seqs)
    assert site_nll(table, "AA", 0) < site_nll(table, "AA", 1)   # position IS used
    assert comp_nll(probs, "AA", 0) == comp_nll(probs, "AA", 1)   # position is NOT used
    # and the position-aware score of A at 0 is strictly better than the composition-only one,
    # because position 0 carries information that pooling the window destroys
    assert site_nll(table, "AA", 0) < comp_nll(probs, "AA", 0)


def test_site_independent_nll_is_the_hand_computed_value():
    """Two sequences, add-one smoothing: counts at position 0 are A=2, C=1, G=1, T=2 over 6."""
    table = fit_site_independent(["A" + "C" * 3, "T" + "C" * 3], length=4)
    assert math.isclose(site_nll(table, "A" + "C" * 3, 0), -math.log(2 / 6), abs_tol=1e-12)


def test_a_uniform_fit_gives_exactly_ln4():
    """Equal counts at a position must put that position at chance -- the sanity anchor."""
    table = fit_site_independent(["A", "C", "G", "T"], length=1)
    assert math.isclose(site_nll(table, "A", 0), math.log(4), abs_tol=1e-12)


@pytest.mark.parametrize("bad", ["ACGN", "ACG"])
def test_site_independent_refuses_bad_input(bad):
    with pytest.raises(PanelError):
        fit_site_independent([bad], length=4)


def test_empty_fits_refuse():
    with pytest.raises(PanelError):
        fit_site_independent([], length=4)
    with pytest.raises(PanelError):
        fit_composition([])


# --------------------------------------------------------------------------- the position replay


def test_rebuild_positions_shape_and_bounds():
    posns = rebuild_positions(5, 7, seed=0, length=150)
    assert len(posns) == 5
    for pos in posns:
        assert len(pos) == len(set(pos)) == 7
        assert pos == sorted(pos)
        assert min(pos) >= 3          # ORDER: a full order-3 left context must exist
        assert max(pos) < 150


def test_rebuild_positions_differ_across_windows():
    """THE TRAP. A per-window re-seed yields identical draws everywhere and still reconciles on
    window 0 alone, which is why the gate checks every window and why this is pinned here."""
    posns = rebuild_positions(8, 10, seed=0, length=150)
    assert len({tuple(p) for p in posns}) == 8


def test_rebuild_positions_is_deterministic_and_seed_sensitive():
    assert rebuild_positions(4, 6, seed=0) == rebuild_positions(4, 6, seed=0)
    assert rebuild_positions(4, 6, seed=0) != rebuild_positions(4, 6, seed=1)


def test_rebuild_positions_is_a_prefix_stream():
    """The RNG is shared, so asking for more windows must EXTEND the stream, not redraw it.

    This is what makes the replay valid for any n_windows the artifacts happen to record.
    """
    assert rebuild_positions(3, 5, seed=0) == rebuild_positions(10, 5, seed=0)[:3]


# --------------------------------------------------------------------------- the frozen bar


def test_verdict_clears_full_panel_only_when_all_three_resolve_for_the_generator():
    v = compute_verdict(_cells((0.02, True), (0.01, True), (0.01, True)), "left_only")
    assert v["VERDICT"] == "GENERATOR_CLEARS_FULL_PANEL"


def test_verdict_beats_only_weak_nulls_when_context_null_wins():
    v = compute_verdict(_cells((0.02, True), (0.005, False), (-0.01, True)), "left_only")
    assert v["VERDICT"] == "GENERATOR_BEATS_ONLY_THE_WEAK_NULLS"
    assert v["vs_markov_ecoli"] == "null"
    assert v["vs_site_independent"] == "tie"


def test_verdict_new_null_dominates_outranks_everything_else():
    """Ordering matters: a NEW null beating the generator is the finding, even alongside a win."""
    v = compute_verdict(_cells((0.02, True), (-0.01, True), (0.01, True)), "left_only")
    assert v["VERDICT"] == "NEW_NULL_DOMINATES_GENERATOR"


def test_verdict_unresolved_when_nothing_separates():
    v = compute_verdict(_cells((0.001, False), (0.001, False), (0.001, False)), "left_only")
    assert v["VERDICT"] == "UNRESOLVED"


def test_an_unresolved_cell_never_counts_as_a_win_however_large_its_delta():
    """A big point estimate with an interval spanning zero is not evidence -- pinned, because this
    is the exact shape that lets a negative get reported as a win."""
    v = compute_verdict(_cells((0.5, False), (0.5, False), (0.5, False)), "left_only")
    assert v["VERDICT"] == "UNRESOLVED"


def test_verdict_reads_the_arm_it_is_asked_for():
    cells = {}
    cells.update(_cells((0.02, True), (0.01, True), (0.01, True), arm="bidirectional"))
    cells.update(_cells((0.02, True), (0.01, True), (-0.01, True), arm="left_only"))
    assert compute_verdict(cells, "bidirectional")["VERDICT"] == "GENERATOR_CLEARS_FULL_PANEL"
    assert compute_verdict(cells, "left_only")["VERDICT"] == "GENERATOR_BEATS_ONLY_THE_WEAK_NULLS"


def test_every_named_outcome_is_reachable():
    """A bar with an unreachable outcome is a bar that cannot return that answer."""
    from scripts.glm_baseline_panel_ce import PREREGISTERED

    reached = {
        compute_verdict(_cells((0.02, True), (0.01, True), (0.01, True)), "left_only")["VERDICT"],
        compute_verdict(_cells((0.02, True), (0.0, False), (-0.01, True)), "left_only")["VERDICT"],
        compute_verdict(_cells((-0.02, True), (0.01, True), (0.01, True)), "left_only")["VERDICT"],
        compute_verdict(_cells((0.0, False), (0.0, False), (0.0, False)), "left_only")["VERDICT"],
    }
    assert reached == set(PREREGISTERED["outcomes"])
