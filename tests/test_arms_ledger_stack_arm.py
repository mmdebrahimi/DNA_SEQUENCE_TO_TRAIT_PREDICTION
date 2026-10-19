"""Guards for the declared orthogonal-modality stack arm.

The declaration must ADD an arm and change nothing else. The seven measured arms are
pinned by metric so a future edit cannot quietly restate a result, and the refusal
`kill_arm` makes on a redundancy reason is pinned too -- because withdrawing the
kill option is only safe while that refusal actually fires.

Every mutating test runs against a COPY of the ledger in `tmp_path`. A test that
killed an arm in the committed ledger would be a destructive test.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, r"C:/my_skills/skills/soraya/scripts")

import arms_ledger as al  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
LEDGER = REPO / "soraya_runs" / "arms-ledger.json"
CRITERION_PREFIX = "ceiling: best achievable prediction of bacterial promoter"
NEW_ARM = "orthogonal-modality-stack-promcalc-plus-N-channels"

# The seven measured arms, pinned by metric. Restating a measured number elsewhere
# would now have to change this line too.
MEASURED = {
    "gc-single-feature": 0.3112,
    "conv-3width-13k-params-150bp": 0.3234,
    "pyramid-4scale-1.70decades-150bp": 0.3167,
    "promoter-calculator-biophysical-346params-fwd_max_tx": 0.3968,
    "promcalc-full-decomposition-ridge-150bp": 0.3925,
    "fragment-rung-321k-rows-size-lever-LB-M9": 0.185,
    "promcalc-plus-conv-rank-average-hybrid-150bp": 0.4254,
}

EXPECTED_PLANNED = {
    "multiresolution-pyramid-4-decades-UNRUN",
    "promoteratlas-base-embeddings-frozen-150bp-padded-to-200",
    "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN",
    NEW_ARM,
}


def _live_ledger() -> al.ArmsLedger:
    return al.ArmsLedger(cwd=str(REPO))


def _criterion(ledger: al.ArmsLedger) -> str:
    return next(k for k in ledger._data["criteria"] if k.startswith(CRITERION_PREFIX))


@pytest.fixture
def sandbox(tmp_path: Path) -> Path:
    """A throwaway copy of the real ledger -- mutating tests must not touch the real one."""
    dest = tmp_path / "soraya_runs"
    dest.mkdir(parents=True)
    shutil.copy2(LEDGER, dest / "arms-ledger.json")
    return tmp_path


def test_the_stack_arm_is_declared_and_planned():
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    assert NEW_ARM in arms
    arm = arms[NEW_ARM]
    assert arm["status"] == al.STATUS_PLANNED
    assert arm["metric"] is None
    assert arm["cost_tier"] == "cheap"
    assert arm["evidence"] == "", "a planned arm has no evidence yet"


def test_the_stack_arm_records_why_it_is_worth_running():
    ledger = _live_ledger()
    arm = ledger._data["criteria"][_criterion(ledger)][NEW_ARM]
    rationale = arm["rationale"]
    # The declaration rests on the shuffled-conv control, so it must be named.
    assert "0.4254" in rationale
    assert "SHUFFLED-CONV" in rationale.upper()
    assert "-0.1232" in rationale


def test_the_withdrawn_kill_option_is_recorded_not_merely_dropped():
    ledger = _live_ledger()
    arm = ledger._data["criteria"][_criterion(ledger)][NEW_ARM]
    note = arm["kill_option_withdrawn"]
    assert "redundancy" in note.lower()
    assert "authority" in note.lower()
    assert "scope reduction" in arm["verdict_consequence"].lower()


def test_the_seven_measured_arms_are_unchanged():
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    measured = {n: a["metric"] for n, a in arms.items()
                if a["status"] == al.STATUS_MEASURED}
    assert measured == MEASURED, "a measured arm's metric moved"
    for name in MEASURED:
        assert arms[name]["evidence"].strip(), f"{name} lost its evidence string"


def test_the_arm_set_is_exactly_seven_measured_and_four_planned():
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    planned = {n for n, a in arms.items() if a["status"] == al.STATUS_PLANNED}
    assert planned == EXPECTED_PLANNED
    assert len(arms) == 11


def test_ceiling_verdict_is_blocked_untried_arm_naming_four():
    ledger = _live_ledger()
    verdict = ledger.ceiling_verdict(_criterion(ledger))
    assert verdict.verdict == "blocked:untried-arm"
    assert set(verdict.planned) == EXPECTED_PLANNED
    assert len(verdict.planned) == 4
    assert not verdict.met


def test_the_eukaryotic_arm_is_still_present():
    """Step 2's north-star step-4 cross-reference must not dangle."""
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    assert "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN" in arms


def test_kill_arm_still_refuses_a_redundancy_reason(sandbox):
    """The regression guard the withdrawal depends on.

    If this refusal ever stops firing, the kill option silently becomes available
    again and a sibling clearing the floor could cancel an unrun arm.
    """
    ledger = al.ArmsLedger(cwd=str(sandbox))
    crit = _criterion(ledger)
    with pytest.raises(al.ArmsError, match="REFUSED"):
        ledger.kill_arm(
            crit, "multiresolution-pyramid-4-decades-UNRUN",
            reason="redundant now that the hybrid arm cleared the floor",
            evidence="the hybrid arm measured 0.4254",
        )


def test_kill_arm_still_refuses_without_evidence(sandbox):
    ledger = al.ArmsLedger(cwd=str(sandbox))
    crit = _criterion(ledger)
    with pytest.raises(al.ArmsError, match="evidence"):
        ledger.kill_arm(crit, NEW_ARM, reason="not worth it", evidence="")


def test_killing_every_unrun_arm_would_change_the_verdict(sandbox):
    """Shows the withdrawal is load-bearing, not decorative.

    With evidence and a non-redundancy reason the kills are permitted -- which is
    exactly why a no-evidence kill has to route to user authority rather than being
    something this plan can do on its own.
    """
    ledger = al.ArmsLedger(cwd=str(sandbox))
    crit = _criterion(ledger)
    for name in EXPECTED_PLANNED:
        ledger.kill_arm(crit, name, reason="measured failure on its own run",
                        evidence="synthetic executed evidence (test only)")
    assert ledger.ceiling_verdict(crit).verdict != "blocked:untried-arm"
    # And the real ledger is untouched.
    assert _live_ledger().ceiling_verdict(
        _criterion(_live_ledger())).verdict == "blocked:untried-arm"
