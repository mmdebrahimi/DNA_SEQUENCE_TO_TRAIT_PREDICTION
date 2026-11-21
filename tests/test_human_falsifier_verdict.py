"""The memo's verdicts must be the ones the FROZEN rule yields from the committed grid.

And the memo must not make the comparison the bar prohibits. Both are enforced here so the prose
cannot drift from the computation, and so a later edit cannot quietly introduce the easy, invalid
cross-substrate comparison.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "human_falsifier_margin.py"
GRID = ROOT / "wiki" / "human_falsifier_margin_2026-11-21.json"
VERDICT = ROOT / "wiki" / "human_falsifier_verdict_2026-11-21.json"
MEMO = ROOT / "wiki" / "human_falsifier_margin_2026-11-21.md"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_hfm", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def grid():
    return json.loads(GRID.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def memo():
    return MEMO.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def flat(memo):
    """Whitespace-collapsed memo.

    Prose assertions must not break on a line wrap: the memo genuinely contains "ONE substrate
    change, not a substrate sweep" but with a newline inside it, and an exact-substring test failed
    on content that was present. Match against this, not the raw text, for any multi-word phrase.
    """
    import re
    return re.sub(r"\s+", " ", memo)


def test_the_grid_is_complete_and_every_cell_scored(grid):
    assert grid["complete"] is True
    assert len(grid["rows"]) == 8
    assert all(r["status"] == "SCORED" for r in grid["rows"])


def test_the_uniform_control_validates_the_instrument_in_every_cell(grid):
    """This gates the whole reading per the frozen bar: a uniform control near chance would mean
    the discriminator is broken and the null's value says nothing."""
    for r in grid["rows"]:
        assert r["uniform_control"] >= 0.95, (r["cell"], r["mode"], r["uniform_control"])


def test_the_verdict_artifact_matches_a_FRESH_computation(mod, grid):
    """The committed verdict must be reproducible from the committed grid."""
    fresh = mod.compute_verdict(grid["rows"])
    stored = json.loads(VERDICT.read_text(encoding="utf-8"))
    assert fresh["by_geometry_and_mode"] == stored["by_geometry_and_mode"]
    assert fresh["adaptor_diagnostic"] == stored["adaptor_diagnostic"]


def test_the_primary_geometry_reads_WEAKER_and_the_memo_says_so(mod, grid, memo):
    v = mod.compute_verdict(grid["rows"])["by_geometry_and_mode"]
    for mode in ("kmer", "both"):
        g = v["full_230nt_as_shipped|%s" % mode]
        assert g["verdict"] == "CHAIN_IS_WEAKER_ON_HUMAN"
        assert g["cells_contradict"] is False
    assert "CHAIN_IS_WEAKER_ON_HUMAN" in memo


def test_the_secondary_geometry_reads_COMPARABLE_because_the_cells_CONTRADICT(mod, grid, memo):
    v = mod.compute_verdict(grid["rows"])["by_geometry_and_mode"]
    for mode in ("kmer", "both"):
        g = v["element_200nt_adaptors_stripped|%s" % mode]
        assert g["verdict"] == "CHAIN_IS_COMPARABLE"
        assert g["cells_contradict"] is True, "the conservative rule must fire on a real split"
        assert g["per_cell"]["HepG2"] == "CHAIN_IS_STRONGER_ON_HUMAN"
        assert g["per_cell"]["K562"] == "CHAIN_IS_WEAKER_ON_HUMAN"
    assert "cells contradict" in memo.lower()


def test_the_memo_records_the_adaptor_artifact_with_its_numbers(mod, grid, memo):
    """The load-bearing correction: the primary geometry's signal is mostly constant flanks."""
    a = mod.compute_verdict(grid["rows"])["adaptor_diagnostic"]["HepG2|both"]
    assert a["null_drop"] > 0.40, "HepG2/both must lose >0.40 of its null signal on stripping"
    assert abs(a["ceiling_move"]) < 0.05, "while its noise floor barely moves"
    assert "-0.4296" in memo or "0.4296" in memo
    assert "uninformative" in memo


def test_the_memo_does_NOT_swap_the_frozen_primary(memo, flat):
    """Swapping primary and secondary after seeing which is favourable is what a pre-registration
    exists to prevent. The memo must say it did not."""
    assert "not** swapped" in flat or "not swapped" in flat or "NOT swapped" in flat
    assert "frozen PRIMARY" in memo


def test_the_memo_makes_NO_raw_cross_substrate_comparison(memo):
    """The bar's prohibition. The bacterial 0.6679 null level must not appear as a comparator."""
    assert "0.6679" not in memo, (
        "the bacterial raw null level appears in the memo; only MARGINS may cross the substrate "
        "boundary, and 0.6679 is a level")


def test_the_memo_refuses_to_claim_a_generator_would_succeed(memo, flat):
    assert "says nothing about whether a generator would succeed" in flat.lower()
    assert "GENERator-v2-prokaryote-1.2b-base" in memo


def test_the_memo_keeps_the_cell_disagreement_unexplained(memo, flat):
    assert "UNEXPLAINED" in memo
    assert "not explained away" in flat


def test_the_memo_carries_the_structural_limits(memo, flat):
    low = flat.lower()
    assert "one substrate change, not a substrate sweep" in low
    assert "0.941" in memo and "not borrowed" in low
    assert "necessary-not-sufficient" in low


def test_the_forecast_outcome_is_reported_honestly(flat):
    """The bar forecast 0.5 for WEAKER; the secondary returned COMPARABLE, so neither."""
    assert "neither confirmed nor refuted" in flat
