"""The colour-cell family freeze: the roster guard, and proof the guard is not vacuous.

WHY THIS EXISTS. The family reached 19 CLI cells before anyone asked whether they could be validated;
`dna_decode/pigment/substrate_screen.py` found 40 of 65 loci record no causal variant at all (7 cells
record none for ANY locus -> unvalidatable as written) and 14 of the 25 that do are off-panel. The freeze
makes adding cell #20 fail loudly instead of silently.

HONESTY: this is an ATTENTION/SCOPE freeze, not enforcement -- ratifiable by deliberately editing
`FROZEN_COLOUR_ROUTES`. The point is that the edit must be conscious.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.data.colour_cell_freeze import (  # noqa: E402
    EXPECTED_TOTALS, FREEZE_DATE, FREEZE_RATIONALE, FROZEN_COLOUR_ROUTES, freeze_status, gates_for,
    screen_summaries,
)
from dna_decode.pigment.substrate_screen import collect, summarise, trait_for_species  # noqa: E402


# ------------------------------------------------------------------------ gates_for (pure)

@pytest.mark.parametrize("unrec,blocked,expect", [
    (2, 0, ("G9",)),            # alpaca: nothing recorded, nothing off-panel to speak of
    (0, 1, ("G10",)),           # buffalo: recorded but structural
    (4, 2, ("G9", "G10")),      # cat: both walls at once
    (0, 0, ()),                 # donkey / roe deer: fully SNV-tractable, no gate applies
])
def test_gates_for_counts(unrec, blocked, expect):
    assert gates_for(unrec, blocked) == expect


def test_a_fully_tractable_cell_gets_an_empty_gate_tuple_not_a_placeholder():
    """An empty tuple is the honest answer for donkey/roe deer. Inventing a gate for a cell that passes
    both would make the freeze look better-justified than it is."""
    assert gates_for(0, 0) == ()


# ------------------------------------------------------------------------ freeze_status

def test_freeze_status_accepts_both_the_bare_trait_and_the_cli_route():
    a = freeze_status("rabbitcolor")
    b = freeze_status("dna-rabbitcolor")
    assert a == b
    assert a["frozen"] is True
    assert a["screen_verdict"] == "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED"
    assert a["gates"] == ("G9",)


def test_freeze_status_on_an_unknown_trait_says_nothing_rather_than_guessing():
    """A trait outside the family is not 'unfrozen because it passed' -- the freeze has no opinion."""
    s = freeze_status("amr")
    assert s["frozen"] is False
    assert s["screen_verdict"] is None
    assert s["gates"] == ()


def test_the_two_fully_tractable_cells_are_frozen_but_carry_no_gate():
    """donkey + roe deer pass G9 and G10 cleanly; they are frozen on scope grounds, not evidence
    grounds, and the status must show that distinction rather than blurring it."""
    for trait in ("donkeycolor", "roedeercolor"):
        s = freeze_status(trait)
        assert s["frozen"] is True
        assert s["screen_verdict"] == "FULLY_SNV_TRACTABLE"
        assert s["gates"] == ()


def test_every_frozen_route_has_a_live_screen_verdict():
    """A roster entry with no derivable verdict would be a name with nothing behind it."""
    for trait in FROZEN_COLOUR_ROUTES:
        assert freeze_status(trait)["screen_verdict"] is not None, trait


# ------------------------------------------------------------------------ the roster guard

def _live_colour_traits() -> set[str]:
    """The colour family DERIVED from the live catalogs -- never a second hand-listed copy."""
    return {trait_for_species(sp) for sp in collect()}


def test_the_frozen_roster_matches_the_live_catalogs():
    """THE FREEZE GUARD. A 20th colour cell fails here.

    To ratify a new cell deliberately: screen it against G9/G10 first
    (`uv run python scripts/colour_cell_substrate_screen.py`), then add its trait to
    FROZEN_COLOUR_ROUTES. This is scope friction by design, not a hard invariant.
    """
    live = _live_colour_traits()
    assert live == set(FROZEN_COLOUR_ROUTES), (
        f"colour family drifted from the freeze roster.\nadded: {sorted(live - set(FROZEN_COLOUR_ROUTES))}\n"
        f"removed: {sorted(set(FROZEN_COLOUR_ROUTES) - live)}\n{FREEZE_RATIONALE}")


def test_the_frozen_roster_matches_the_cli_registry():
    """The other direction: a colour trait routable from the CLI but absent from the catalogs (or vice
    versa) means the family and its screen have come apart."""
    from dna_decode.cli import TRAITS
    cli_colour = {k for k in TRAITS if k.endswith("color") or k == "plumage"}
    assert cli_colour == set(FROZEN_COLOUR_ROUTES)


def test_the_freeze_guard_is_not_vacuous(monkeypatch):
    """MANDATORY. A guard that passes while checking nothing is worse than no guard -- that exact defect
    was caught twice in the session that produced this plan (a control fixture pinned at +1.000 by
    construction; a retraction guard satisfied by the quotation inside its own retraction).

    This does NOT re-implement the comparison and check the arithmetic; that would prove only that set
    inequality works. It injects a synthetic 20th cell into the LIVE catalog collection and asserts the
    REAL guard function raises -- exercising the same code path a genuine new cell would take.
    """
    import dna_decode.pigment.substrate_screen as ss

    real_collect = ss.collect

    def collect_with_a_20th_cell():
        data = dict(real_collect())
        data["unicorn"] = [{"locus": "U", "gene": "HORN", "variant_class": "UNRECORDED",
                            "snv_panel_scorable": None, "source": ""}]
        return data

    monkeypatch.setattr(ss, "collect", collect_with_a_20th_cell)
    # this module imported `collect` by name, so patch the local binding the guard actually calls
    monkeypatch.setitem(globals(), "collect", collect_with_a_20th_cell)

    assert "unicorncolor" in _live_colour_traits(), "the injection did not reach the guard's input"
    with pytest.raises(AssertionError, match="drifted from the freeze roster"):
        test_the_frozen_roster_matches_the_live_catalogs()


def test_the_freeze_declares_its_date_and_rationale():
    assert FREEZE_DATE == "2026-08-26"
    assert "unvalidatable as written" in FREEZE_RATIONALE
    assert "G9/G10" in FREEZE_RATIONALE


# ------------------------------------------------------------------------ headline counts

def test_the_headline_counts_are_pinned_in_one_place():
    """The memo and the code must not diverge silently. If a catalog gains a causal variant these MOVE --
    update EXPECTED_TOTALS and the memo in the same commit rather than loosening the test."""
    data = collect()
    tot = {"n_cells": len(data), "n_loci": 0, "n_unrecorded": 0, "n_snv_panel_blocked": 0}
    for rows in data.values():
        s = summarise(rows)
        tot["n_loci"] += s["n_loci"]
        tot["n_unrecorded"] += s["n_unrecorded"]
        tot["n_snv_panel_blocked"] += s["n_snv_panel_blocked"]
    assert tot == EXPECTED_TOTALS


def test_seven_cells_are_unscreenable():
    """The finding that makes this a curation wall rather than only a substrate wall."""
    unscreenable = {t for t, s in screen_summaries().items()
                    if s["verdict"] == "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED"}
    assert unscreenable == {"alpacacolor", "cattlecolor", "mousecolor", "pigcolor",
                            "pigeoncolor", "rabbitcolor", "sheepcolor"}


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
