"""The 19 colour CellContracts must declare their real gate and a TRUTHFUL promotion path.

TWO DEFECTS THIS PINS (both live in the repo before 2026-08-26):

  1. All 19 carried `incoming_data_gate="n/a"` while gates G9/G10 apply directly to them.
  2. SEVEN cells promised a promotion path that CANNOT WORK. `dna-rabbitcolor` read
     "MEASURED needs a free rabbit genotype+observed-colour cohort" -- but none of its 5 loci record a
     causal variant, so no cohort could ever score the rule. A contract that names a promotion path must
     name a SUFFICIENT one.

Assertions here read the LIVE `substrate_screen` derivation, never a second copy of the prose -- a
wording-only change is otherwise trivial to fake green.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.data import cell_registry as cr  # noqa: E402
from dna_decode.data.colour_cell_freeze import (  # noqa: E402
    FROZEN_COLOUR_ROUTES, freeze_status, screen_summaries,
)

_UNSCREENABLE = "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED"


def _colour_contracts() -> dict:
    return {c.route.removeprefix("dna-"): c for c in cr.cells()
            if c.route.removeprefix("dna-") in FROZEN_COLOUR_ROUTES}


def test_all_nineteen_colour_cells_have_a_contract():
    got = _colour_contracts()
    assert set(got) == set(FROZEN_COLOUR_ROUTES), sorted(set(FROZEN_COLOUR_ROUTES) - set(got))


def test_no_gated_colour_contract_still_declares_n_a():
    """The field exists for exactly this and was "n/a" on all 19 before 2026-08-26.

    SCOPED to cells a gate actually applies to. donkey + roe deer pass G9 and G10, so "n/a" is the
    CORRECT value for them -- asserting "no colour cell says n/a" would force a false gate onto the two
    cells that cleared the screen. That they were screened-and-passed is carried by the CLI validation
    string ("screened CLEAR"), not by inventing a gate id here.
    """
    stale = [t for t, c in _colour_contracts().items()
             if c.incoming_data_gate.strip() == "n/a" and freeze_status(t)["gates"]]
    assert not stale, f"colour cells with a live gate still declaring 'n/a': {sorted(stale)}"
    # and the guard must be non-vacuous: 17 of 19 DO carry a gate
    gated = [t for t in _colour_contracts() if freeze_status(t)["gates"]]
    assert len(gated) == 17, f"expected 17 gated colour cells, found {len(gated)}"


@pytest.mark.parametrize("trait", sorted(FROZEN_COLOUR_ROUTES))
def test_each_declared_gate_matches_the_live_screen_verdict(trait):
    """The declared gate must agree EXACTLY with what the screen DERIVES right now, so a future catalog
    change that flips a verdict fails HERE instead of drifting silently.

    The field's contract is a comma-separated list of BARE gate ids, or "n/a" when none applies -- the
    same shape the pre-existing cells use ('G1,G7,G8'). Reasoning and counts live in `demotion_rule`, the
    CLI `validation` string, and the memo; this field is the machine-readable projection. Writing prose
    here broke `test_incoming_gate_subset_of_known_gates`, which parses it by splitting on commas.
    """
    declared = _colour_contracts()[trait].incoming_data_gate
    gates = freeze_status(trait)["gates"]
    expected = ",".join(gates) if gates else "n/a"
    assert declared == expected, f"{trait}: screen derives {expected!r}, contract declares {declared!r}"


def test_the_declared_gates_parse_under_the_registry_wide_contract():
    """Pins the interop that the first version of this edit broke: every colour gate must survive the
    comma-split the registry-wide guard performs, against the known G1-G10 vocabulary."""
    known = {f"G{i}" for i in range(1, 11)}
    for trait, c in _colour_contracts().items():
        if c.incoming_data_gate == "n/a":
            continue
        assert {t.strip() for t in c.incoming_data_gate.split(",")} <= known, trait


def test_the_two_clear_cells_declare_n_a_rather_than_inventing_a_gate():
    """donkey + roe deer pass G9 and G10; 'n/a' is the honest value. That they were SCREENED and passed
    (rather than never screened) is carried by the CLI validation string, which says 'screened CLEAR'."""
    from dna_decode.cli import TRAITS
    for trait in ("donkeycolor", "roedeercolor"):
        assert _colour_contracts()[trait].incoming_data_gate == "n/a"
        assert "screened CLEAR" in TRAITS[trait]["validation"]


@pytest.mark.parametrize("trait", sorted(
    t for t, s in screen_summaries().items() if s["verdict"] == _UNSCREENABLE))
def test_unscreenable_cells_state_the_real_precondition(trait):
    """REGRESSION. These 7 promised a cohort would promote them. It cannot."""
    rule = _colour_contracts()[trait].demotion_rule
    assert "UNVALIDATABLE AS WRITTEN" in rule, f"{trait}: {rule[:120]!r}"
    assert "NOT SUFFICIENT" in rule, f"{trait} does not say a cohort is insufficient"
    assert "Curating the causal variants" in rule, f"{trait} does not name the real precondition"


def test_the_exact_cohort_only_wording_that_was_wrong_is_gone():
    """Pins the specific rabbit string, so the false promise cannot silently return."""
    rule = _colour_contracts()["rabbitcolor"].demotion_rule
    assert "MEASURED needs a free rabbit genotype+observed-colour cohort" not in rule


def test_partially_tractable_cells_keep_their_cohort_path():
    """The correction is SCOPED. A cell with recorded SNV loci CAN be scored by a cohort, so its
    promotion path stays -- over-applying the fix would swap one false claim for another."""
    for trait in ("catcolor", "plumage", "horsecolor"):
        assert screen_summaries()[trait]["verdict"] != _UNSCREENABLE
        assert "UNVALIDATABLE AS WRITTEN" not in _colour_contracts()[trait].demotion_rule


def test_the_dog_cell_keeps_its_measured_result_in_the_demotion_rule():
    """coatcolor already recorded the measured substrate limitation; Step 3 only added its gate."""
    rule = _colour_contracts()["coatcolor"].demotion_rule
    assert "0.994" in rule and "SUBSTRATE-LIMITED" in rule


def test_all_nineteen_remain_knowledge_baseline():
    """No tier was moved -- demoting the 7 is an open question, deliberately out of scope."""
    for trait, c in _colour_contracts().items():
        assert c.evidence_tier.name == "KNOWLEDGE_BASELINE", trait


def test_non_colour_contracts_were_not_touched_by_the_scoped_edit():
    """`incoming_data_gate="n/a"` appears 48x repo-wide; a blanket replace would have rewritten 29
    unrelated contracts. This asserts the edit stayed inside the colour family."""
    others = [c for c in cr.cells() if c.route.removeprefix("dna-") not in FROZEN_COLOUR_ROUTES]
    still_na = [c for c in others if c.incoming_data_gate.strip() == "n/a"]
    assert len(still_na) >= 25, (
        f"only {len(still_na)} non-colour contracts still declare 'n/a' — the scoped edit likely "
        f"leaked outside the colour family")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
