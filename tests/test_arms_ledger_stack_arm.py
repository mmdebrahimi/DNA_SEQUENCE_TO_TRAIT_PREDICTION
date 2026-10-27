"""Guards for the declared orthogonal-modality stack arm.

RE-PINNED TWICE ON 2026-10-27. Second time: the 4-decade arm was RETIRED on its own executed
evidence (`scripts/pyramid_span_feasibility_probe.py` measured the maximum HONEST span at 2.2210
decades against the 4.00 it declares -- 60.1x short at the encoder's 3 bp floor) and the reachable
successor `pyramid-4scale-2.22decades-499bp-fragments` was declared in its place. So the set is now
8 measured / 3 planned / 1 KILLED over 12 arms. That is a legitimate move under R6 -- an arm dying by
its OWN result, not a sibling's -- and the rails follow the ledger rather than freezing a stale set.

First time: the declared stack arm RAN. It was `planned` when these guards were
written and is now `measured` at 0.4404 on its own evidence, so four of them went red
because their SUBJECT legitimately moved -- not because the ledger was tampered with.
Re-pinning to the state that exists is the answer here; the invariants are UNCHANGED and
two are now stronger, asserted over EVERY planned arm rather than over this one.

The EIGHT measured arms are pinned by metric so a future edit cannot quietly restate a
result, and the refusal `kill_arm` makes on a redundancy reason is pinned too -- because
withdrawing the kill option is only safe while that refusal actually fires.

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

# The eight measured arms, pinned by metric. Restating a measured number elsewhere
# would now have to change this line too.
MEASURED = {
    "gc-single-feature": 0.3112,
    "conv-3width-13k-params-150bp": 0.3234,
    "pyramid-4scale-1.70decades-150bp": 0.3167,
    "promoter-calculator-biophysical-346params-fwd_max_tx": 0.3968,
    "promcalc-full-decomposition-ridge-150bp": 0.3925,
    "fragment-rung-321k-rows-size-lever-LB-M9": 0.185,
    "promcalc-plus-conv-rank-average-hybrid-150bp": 0.4254,
    # Ran 2026-10-26: WEAK_STACK, +0.0145 over the 2-channel hybrid at 10/10 seeds,
    # CI95 [0.0108, 0.0172] excluding zero, short only of the INHERITED 0.02 margin.
    NEW_ARM: 0.4404,
}

EXPECTED_PLANNED = {
    "promoteratlas-base-embeddings-frozen-150bp-padded-to-200",
    "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN",
    # The reachable re-scope that replaced multiresolution-pyramid-4-decades-UNRUN.
    "pyramid-4scale-2.22decades-499bp-fragments",
}

#: Retired 2026-10-27 on its own executed evidence, NOT on a sibling's result.
KILLED_ARM = "multiresolution-pyramid-4-decades-UNRUN"
SUCCESSOR = "pyramid-4scale-2.22decades-499bp-fragments"
SPAN_PROBE = "wiki/pyramid_span_feasibility_2026-10-27.json"


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


def test_the_stack_arm_RAN_and_carries_its_own_evidence():
    """Was `planned`; ran 2026-10-26. R6 is satisfied by its OWN executed evidence.

    The arm could not be killed on a sibling's result (see the withdrawal guards below) and
    it was not: it was RUN, and that is the only thing that legitimately moves an arm off
    `planned`. The metric is pinned in MEASURED, so a restatement has to edit that line too.
    """
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    assert NEW_ARM in arms
    arm = arms[NEW_ARM]
    assert arm["status"] == al.STATUS_MEASURED
    assert arm["metric"] == 0.4404
    assert arm["cost_tier"] == "cheap", "running it did not change what it cost"
    assert arm["evidence"].strip(), "a measured arm must carry its evidence"


def test_every_PLANNED_arm_still_has_an_empty_evidence_string():
    """Generalised from the stack arm, which used to be the one carrying this invariant.

    `evidence` means EXECUTED evidence (R6/R7). A substrate fact, a wall classification or a
    rationale is NOT that, so none of them may be parked in this field -- which is why the
    three untried arms' walls live in `wiki/ceiling_arm_wall_classification_2026-10-26.md`
    and not in the ledger. Asserting it over every planned arm is strictly stronger than
    asserting it over one.
    """
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    planned = {n: a for n, a in arms.items() if a["status"] == al.STATUS_PLANNED}
    assert planned, "vacuous: no planned arms left to check"
    for name, arm in planned.items():
        assert arm["evidence"] == "", "%s: a planned arm has no executed evidence" % name
        assert arm["metric"] is None, "%s: a planned arm has no metric" % name


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


def test_the_eight_measured_arms_are_unchanged():
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    measured = {n: a["metric"] for n, a in arms.items()
                if a["status"] == al.STATUS_MEASURED}
    assert measured == MEASURED, "a measured arm's metric moved"
    for name in MEASURED:
        assert arms[name]["evidence"].strip(), f"{name} lost its evidence string"


def test_the_arm_set_is_exactly_eight_measured_and_three_planned():
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    planned = {n for n, a in arms.items() if a["status"] == al.STATUS_PLANNED}
    assert planned == EXPECTED_PLANNED
    killed = {n for n, a in arms.items() if a["status"] == al.STATUS_KILLED}
    assert killed == {KILLED_ARM}, "exactly one arm has ever been retired here"
    # 12 total: 11 plus the declared successor. The arithmetic must reconcile with the whole.
    assert len(arms) == 12
    assert len(EXPECTED_PLANNED) + len(MEASURED) + len(killed) == len(arms)


def test_ceiling_verdict_is_blocked_untried_arm_naming_three():
    ledger = _live_ledger()
    verdict = ledger.ceiling_verdict(_criterion(ledger))
    assert verdict.verdict == "blocked:untried-arm"
    assert set(verdict.planned) == EXPECTED_PLANNED
    assert len(verdict.planned) == 3
    assert not verdict.met


def test_the_retired_arm_cannot_be_READ_as_tested_and_failed():
    """`killed` is the only status the enum offers for "cannot be instantiated", so the REASON
    carries the whole distinction. Without these clauses a later reader concludes the 4-decade
    hypothesis was tested and lost, when the arm was never run at all.
    """
    ledger = _live_ledger()
    arm = ledger._data["criteria"][_criterion(ledger)][KILLED_ARM]
    assert arm["status"] == al.STATUS_KILLED
    assert arm["evidence"] == SPAN_PROBE, "the kill must cite the probe that grounds it"
    reason = arm["reason"]
    low = reason.lower()
    # R7: scope narrow enough that the evidence supports it.
    assert "scope" in low and "this host" in low, "the kill must bound itself to this substrate"
    assert "not a refutation" in low, "it must refuse the tested-and-failed reading outright"
    assert "predictive merit" in low, "it must disclaim saying anything about the arm's merit"
    # R7: the observation that re-opens it.
    assert "revival" in low
    assert "30000" in reason or "30,000" in reason
    # It must name the successor, or the re-scope is lost with the arm.
    assert SUCCESSOR in reason


def test_the_retired_arm_was_NOT_killed_for_redundancy():
    """The exact defect R6 exists to prevent: dying because a sibling cleared the floor."""
    ledger = _live_ledger()
    reason = ledger._data["criteria"][_criterion(ledger)][KILLED_ARM]["reason"]
    assert not al._REDUNDANCY_REASON.search(reason), \
        "the stored reason matches the redundancy pattern kill_arm refuses"
    assert "substrate-gated" in reason.lower()


def test_the_successor_arm_records_WHY_it_is_worth_running():
    """Planned arms carried no reason at all, which is why `blocked:untried-arm` could only ever
    name them. The successor sets `rationale` directly, since `declare_arm` takes no such argument.
    """
    ledger = _live_ledger()
    arm = ledger._data["criteria"][_criterion(ledger)][SUCCESSOR]
    r = arm["rationale"]
    assert "2.22" in r and "1.70" in r, "it must state the span it buys over the measured sibling"
    assert "0.7928" in r, "it must name the FRAGMENT rung's own ceiling"
    assert "0.1313" in r, "it must name the fragment-rung incumbent it has to beat"
    assert "never be read against" in r.lower(), \
        "it must refuse the cross-rung comparison the census ceiling_rule forbids"
    assert "pad/crop" in r, "the variable-length design point must be declared pre-run"
    # And the rationale must NOT have leaked into the field reserved for executed evidence.
    assert arm["evidence"] == ""


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
        ledger.kill_arm(crit, "promoteratlas-base-embeddings-frozen-150bp-padded-to-200",
                        reason="not worth it", evidence="")


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
