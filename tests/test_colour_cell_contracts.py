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


def test_no_colour_contract_still_declares_gate_n_a():
    """The field exists for exactly this and was empty on all 19."""
    stale = [t for t, c in _colour_contracts().items() if c.incoming_data_gate.strip() == "n/a"]
    assert not stale, f"colour cells still declaring no incoming data gate: {sorted(stale)}"


@pytest.mark.parametrize("trait", sorted(FROZEN_COLOUR_ROUTES))
def test_each_declared_gate_matches_the_live_screen_verdict(trait):
    """The declared gate must agree with what the screen DERIVES for that cell right now, so a future
    catalog change that flips a verdict fails HERE instead of drifting silently."""
    c = _colour_contracts()[trait]
    gates = freeze_status(trait)["gates"]
    declared = c.incoming_data_gate
    for g in gates:
        assert g in declared, f"{trait}: screen says {gates}, contract omits {g}: {declared!r}"
    if not gates:
        # donkey / roe deer pass both gates; the contract must say CLEAR, not invent a gate
        assert "CLEAR" in declared, f"{trait} passes G9+G10 but its contract does not say so: {declared!r}"
        assert "G9 causal-variant-unrecorded" not in declared
        assert "G10 variant-class-off-panel" not in declared


def test_every_declared_gate_cites_the_screen_artifact():
    for trait, c in _colour_contracts().items():
        assert "colour_cell_substrate_screen_2026-08-26" in c.incoming_data_gate, trait


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
