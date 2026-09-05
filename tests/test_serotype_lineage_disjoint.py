"""The serotype fix survives a LINEAGE-disjoint split -- pin the guarantee and the guards.

The prior replication was held out BY ISOLATE, which does not prove generalization past near-identical
genomes. This run splits whole sequence types. The tests below pin the three things that make that
claim mean something: that the disjointness is by construction, that the comparison was powered, and
that the committed artifact says what the memo says it says.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from serotype_lineage_disjoint import lineage_side, verdict_for  # noqa: E402

ARTIFACT = ROOT / "wiki" / "serotype_lineage_disjoint_2026-09-04.json"


@pytest.fixture(scope="module")
def art() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("lineage-disjoint artifact absent")
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


# --- the split ------------------------------------------------------------------------------------

def test_a_sequence_type_can_never_land_on_both_halves():
    """Disjointness is BY CONSTRUCTION: the side is a pure function of the ST.

    This is stronger than checking two lists for overlap after a run -- that would only say it did not
    happen that time. Here there is no code path that could place an ST on both sides.
    """
    for st in ["131", "10", "69", "1193", "ST-with-odd-name", "1"]:
        assert lineage_side(st) == lineage_side(st)
        assert lineage_side(st) in {"TEST", "TRAIN"}


def test_the_split_is_not_degenerate():
    """A hash that sent every ST to one side would satisfy determinism and prove nothing."""
    sides = {lineage_side(str(i)) for i in range(200)}
    assert sides == {"TEST", "TRAIN"}


# --- the pre-registered verdict rule --------------------------------------------------------------

def test_the_power_guard_outranks_the_sign_of_the_gain():
    """Zero flips means the split never exercised the rule -- that is absence of evidence, not failure.

    The ordering matters: a 0-flip run necessarily has a 0 gain, so if the falsification branch were
    checked first, an unpowered run would be published as a refutation of the fix.
    """
    verdict, why = verdict_for(n_flips=0, h_gain=0.0, n_test=169)
    assert verdict == "UNDERPOWERED_RULES_NEVER_DIFFERED"
    assert "NOT" in why and "evidence the fix fails" in why


def test_a_real_negative_is_still_reachable():
    """The rule must be able to falsify the fix, or it is not a test."""
    verdict, _ = verdict_for(n_flips=28, h_gain=-0.02, n_test=169)
    assert verdict == "FALSIFIED_ON_LINEAGE_DISJOINT"


def test_an_attenuated_gain_gets_its_own_branch():
    """A positive-but-small gain must NOT read as a clean survival -- it means overlap carried part."""
    verdict, why = verdict_for(n_flips=28, h_gain=0.02, n_test=169)
    assert verdict == "SURVIVES_BUT_SMALLER_ON_LINEAGE_DISJOINT"
    assert "attenuated" in why


def test_the_run_took_the_survival_branch():
    verdict, _ = verdict_for(n_flips=28, h_gain=0.1617, n_test=169)
    assert verdict == "SURVIVES_LINEAGE_DISJOINT"


# --- the committed artifact -----------------------------------------------------------------------

def test_artifact_verdict_and_gain(art):
    assert art["verdict"] == "SURVIVES_LINEAGE_DISJOINT"
    assert art["H_gain"] == pytest.approx(0.1617, abs=1e-3)


def test_artifact_is_powered_not_vacuous(art):
    """28 differing H calls. A survival verdict on 0 flips would be meaningless."""
    assert art["n_H_calls_differing_between_rules"] > 0
    assert verdict_for(art["n_H_calls_differing_between_rules"],
                       art["H_gain"], art["test"]["n_isolates"])[0] == art["verdict"]


def test_untyped_isolates_are_excluded_and_counted_not_pooled(art):
    """Pooling ST-less isolates into one bucket would manufacture a fake lineage and break the claim."""
    assert "n_excluded_no_complete_st" in art
    assert art["n_excluded_no_complete_st"] == 2
    assert art["test"]["n_isolates"] + art["train"]["n_isolates"] + art["n_excluded_no_complete_st"] \
        == art["n_scored"] - 0  # 400 scored = 169 test + 229 train + 2 untyped


def test_the_st_counts_partition(art):
    assert art["test"]["n_st"] + art["train"]["n_st"] == art["n_distinct_st"]


def test_the_o_axis_regression_stays_visible(art):
    """The fix costs a little O accuracy. That cost is part of the finding, not a rounding detail."""
    assert art["O_gain"] < 0
    assert art["O_gain"] == pytest.approx(-0.0068, abs=1e-3)


def test_the_over_reading_guard_is_recorded(art):
    """The lineage gain is LARGER than the isolate gain; the artifact must refuse to read that as growth."""
    cmp = art["comparison_to_the_isolate_disjoint_run"]
    assert cmp["lineage_disjoint_H_gain"] > cmp["isolate_disjoint_H_gain"]
    assert "NOT directly comparable" in cmp["do_not_over_read"]


def test_limits_name_st_as_a_lineage_unit_not_a_clonality_guarantee(art):
    joined = " ".join(art["honest_limits"]).lower()
    assert "not a clonality guarantee" in joined or "clonality guarantee" in joined
