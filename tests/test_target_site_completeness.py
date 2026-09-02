"""The target-site completeness signal: complementary to position-novelty, and augment-only."""
from __future__ import annotations

import pytest

from dna_decode.data.target_site_completeness import (
    TARGET_SITE_COMPLETENESS, UNMEASURED_CELLS, completeness_units_for, is_measured, matching_units,
)
from dna_decode.eval.doubt import NONE, STRONG, doubt_one_line, target_site_doubt
from dna_decode.eval.position_novelty import flag_for_cell


# --- the two signals cover DIFFERENT blind spots; neither subsumes the other --------------------

def test_position_novelty_is_silent_on_a_non_catalogued_position():
    """The entire reason this signal exists. If position-novelty caught V179F, the index would be
    redundant and should be deleted."""
    assert flag_for_cell(["V179F"], "hiv-nnrti-rt").position_novel is False
    assert flag_for_cell(["K103R"], "hiv-nnrti-rt").position_novel is True


def test_the_completeness_signal_fires_where_position_novelty_cannot():
    sigs = {s.kind: s for s in target_site_doubt("efavirenz", {"RT": {"V179F"}}).signals}
    assert sigs["position_novelty"].tier == NONE
    assert sigs["target_site_completeness"].tier == STRONG


def test_position_novelty_still_fires_on_its_own_shape():
    """Regression: adding a signal must not disturb the incumbent."""
    sigs = {s.kind: s for s in target_site_doubt("efavirenz", {"RT": {"K103R"}}).signals}
    assert sigs["position_novelty"].tier != NONE
    assert sigs["target_site_completeness"].tier == NONE


# --- the human-facing line must name the signal that actually fired ------------------------------

def test_the_printed_line_reports_the_reason_of_the_signal_that_fired():
    """It hardcoded signals[0], so a STRONG tier printed alongside the OTHER signal's 'found nothing'
    prose -- a self-contradicting disclosure, the worst failure mode for a human-facing line."""
    block = target_site_doubt("efavirenz", {"RT": {"V179F"}}).as_dict()
    line = doubt_one_line(block)
    assert block["max_tier"] == STRONG
    assert "V179F" in line and "does not carry" in line
    assert "no uncatalogued substitution" not in line


# --- three states, never collapsed ---------------------------------------------------------------

def test_a_never_measured_cell_says_so_rather_than_reporting_clean():
    sigs = {s.kind: s for s in target_site_doubt("nirmatrelvir", {"Mpro": {"E166V"}}).signals}
    s = sigs["target_site_completeness"]
    assert s.tier == NONE and s.evidence["measured"] is False
    assert "NOT an absence of doubt" in s.reason


def test_measured_and_quiet_is_distinguishable_from_never_measured():
    sigs = {s.kind: s for s in target_site_doubt("efavirenz", {"RT": {"K101P"}}).signals}
    s = sigs["target_site_completeness"]
    assert s.tier == NONE and s.evidence["measured"] is True


def test_an_unsurfaced_path_mentions_both_screens():
    """The not-assessable message named only position-novelty while two screens now exist."""
    sig = target_site_doubt("efavirenz", None).signals[0]
    assert "completeness" in sig.reason and "position-novelty" in sig.reason


# --- the index's entry bar ------------------------------------------------------------------------

def test_every_listed_unit_is_pure_and_survives_familywise_correction():
    for cell, units in TARGET_SITE_COMPLETENESS.items():
        for sub, st in units.items():
            assert st["carriers_labelled_s"] == 0, f"{sub} has a susceptible carrier -- signal ENDS"
            assert st["carriers_labelled_r"] >= 5, f"{sub} underpowered"
            assert st["purity_surprise_p"] <= 0.05 / st["n_units_tested"], f"{sub} fails correction"
            assert st["artifact"].endswith(".json") and st["label"]


def test_a_listed_unit_is_not_already_representable_by_the_catalog():
    """A gap the catalog already covers is not a gap."""
    for cell, units in TARGET_SITE_COMPLETENESS.items():
        for sub in units:
            assert flag_for_cell([sub], cell).position_novel is False


def test_unmeasured_cells_are_declared_and_disjoint_from_measured_ones():
    assert UNMEASURED_CELLS
    assert not (UNMEASURED_CELLS & set(TARGET_SITE_COMPLETENESS))


def test_matching_is_case_insensitive_and_stable():
    assert [s for s, _ in matching_units(["v179f"], "hiv-nnrti-rt")] == ["V179F"]
    assert matching_units([], "hiv-nnrti-rt") == []
    assert is_measured("hiv-nnrti-rt") and not is_measured("nope")
    assert completeness_units_for("nope") == {}


def test_the_block_still_refuses_to_emit_a_call():
    """The layer's load-bearing constraint, re-checked with a second signal present."""
    for obs in ({"RT": {"V179F"}}, {"RT": {"K103R"}}, None):
        d = target_site_doubt("efavirenz", obs).as_dict()   # as_dict() runs assert_no_call
        assert "prediction" not in d and "call" not in d
