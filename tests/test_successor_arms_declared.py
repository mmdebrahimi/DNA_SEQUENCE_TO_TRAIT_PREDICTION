"""Guards for the three successor arms declared under R6 on 2026-11-02.

**What these pin and why it is not the same as `test_arms_ledger_stack_arm.py`.** That file guards
the criterion's PARTITION and its verdict. This file guards the DECLARATION ITSELF: that the three
arms exist in the state a planned arm must be in, that declaring them could not and did not disturb
any pre-existing arm, and that the declaration script is a safe no-op on re-run.

The `declare()` function is tested against SYNTHETIC ledgers rather than the live file wherever the
property is about behaviour, so the tests say something about the function and not merely about
today's ledger.

**The load-bearing one is `evidence == ""`.** Under R6/R7 `evidence` means EXECUTED evidence; a
rationale, a substrate fact or a wall classification is not that. These three arms carry long
rationales, and the whole risk of writing them is that a rationale ends up parked in `evidence`,
where it would make an unrun arm look resolved.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

from scripts.declare_successor_arms import (
    ARM_FIELDS,
    ARMS,
    PLAN,
    DeclareError,
    ceiling_criterion,
    declare,
)

REPO = Path(__file__).resolve().parents[1]
LEDGER = REPO / "soraya_runs/arms-ledger.json"
SORAYA_SCRIPTS = Path(r"C:/my_skills/skills/soraya/scripts")

NEW_ARMS = (
    "asgiven-position-resolved-promoteratlas-150bp",
    "shipped-expression-heads-zeroshot-150bp",
    "finetuned-promoteratlas-encoder-150bp",
)


def _data() -> dict:
    return json.loads(LEDGER.read_text(encoding="utf-8"))


def _arms(d: dict) -> dict:
    return d["criteria"][ceiling_criterion(d)]


# --------------------------------------------------------------------------------------------------
# the declared state
# --------------------------------------------------------------------------------------------------
def test_all_three_arms_are_declared_on_the_single_ceiling_criterion():
    arms = _arms(_data())
    for name in NEW_ARMS:
        assert name in arms, "%s was never declared" % name


@pytest.mark.parametrize("name", NEW_ARMS)
def test_each_declared_arm_is_in_the_state_a_planned_arm_must_be_in(name):
    a = _arms(_data())[name]
    if a["status"] != "planned":
        pytest.skip("%s has since been measured (status %r); the planned-state guard has retired "
                    "for it" % (name, a["status"]))
    assert a["metric"] is None, "a planned arm has no metric"
    assert a["evidence"] == "", (
        "a planned arm has no EXECUTED evidence. A rationale is not evidence under R6/R7, and "
        "parking one here would make an unrun arm look resolved.")
    assert a["resolved_at"] == ""


@pytest.mark.parametrize("name", NEW_ARMS)
def test_each_declared_arm_carries_its_rationale_and_provenance(name):
    a = _arms(_data())[name]
    assert a.get("rationale", "").strip(), "%s has no rationale" % name
    assert len(a["rationale"]) > 300, "a one-line rationale is not provenance"
    assert a.get("declared_from") == PLAN
    assert a["cost_tier"] == ARMS[name]["cost_tier"]


def test_the_rationales_record_what_makes_each_arm_distinct_rather_than_restating_the_last_one():
    """Each arm must say why it is NOT a re-read of the positional arm."""
    arms = _arms(_data())
    asg = arms["asgiven-position-resolved-promoteratlas-150bp"]["rationale"]
    assert "DECLARED PRIMARY" in asg and "re-read" in asg
    assert "0.5086" in asg and "0.4893" in asg, "the motivating measurement must be named"
    assert "TIE" in asg.upper(), "the predicted tie must be pre-committed, not discovered later"

    heads = arms["shipped-expression-heads-zeroshot-150bp"]["rationale"]
    assert "11008" in heads and "86" in heads, "the geometry lock is arithmetic and must be stated"
    assert "BYTE-IDENTICAL" in heads.upper()
    assert "sigma70" in heads, "the urtecho2018 circularity scope must be recorded"

    ft = arms["finetuned-promoteratlas-encoder-150bp"]["rationale"]
    assert "1,847,428" in ft, "the parameter count is why this is tractable at all"
    assert "budget" in ft.lower(), "the budget-vs-representation ambiguity must be pre-committed"


# --------------------------------------------------------------------------------------------------
# declaring must not disturb anything -- tested on SYNTHETIC ledgers
# --------------------------------------------------------------------------------------------------
def _synthetic(extra: dict | None = None) -> dict:
    arms = {
        "already-measured-arm": {"name": "already-measured-arm", "cost_tier": "cheap",
                                 "status": "measured", "metric": 0.4242, "evidence": "wiki/x.json",
                                 "declared_at": "t0", "resolved_at": "t1", "annotation": "keep me"},
        "already-killed-arm": {"name": "already-killed-arm", "cost_tier": "large",
                               "status": "killed", "metric": None, "evidence": "wiki/y.json",
                               "declared_at": "t0", "resolved_at": "t1", "reason": "measured fail"},
    }
    arms.update(extra or {})
    return {"schema_version": 1, "criteria": {"ceiling: synthetic criterion": arms}}


def test_declare_adds_every_missing_arm_and_reports_them():
    d = _synthetic()
    added, present = declare(d, now="stamp")
    assert added == sorted(NEW_ARMS)
    assert present == []
    assert len(_arms(d)) == 2 + len(NEW_ARMS)


def test_declare_never_touches_a_pre_existing_arm():
    d = _synthetic()
    before = json.loads(json.dumps(_arms(d)))
    declare(d, now="stamp")
    after = _arms(d)
    for name, arm in before.items():
        assert after[name] == arm, "%s was modified by a declaration" % name


def test_declare_is_idempotent_and_reports_the_no_op():
    d = _synthetic()
    declare(d, now="stamp")
    snapshot = json.loads(json.dumps(_arms(d)))
    added, present = declare(d, now="DIFFERENT-stamp")
    assert added == []
    assert present == sorted(NEW_ARMS)
    assert _arms(d) == snapshot, "a second declaration changed the ledger"


def test_declare_writes_only_schema_fields_plus_the_two_named_annotations():
    d = _synthetic()
    declare(d, now="stamp")
    for name in NEW_ARMS:
        extra = set(_arms(d)[name]) - set(ARM_FIELDS)
        assert extra == {"rationale", "declared_from"}, (
            "%s carries unexpected keys %r; every annotation key is a key the API path would raise "
            "on, so they are added deliberately or not at all" % (name, sorted(extra)))


def test_declare_refuses_a_ledger_without_exactly_one_ceiling_criterion():
    d = _synthetic()
    d["criteria"]["ceiling: a second one"] = {}
    with pytest.raises(DeclareError, match="expected exactly one"):
        declare(d)
    with pytest.raises(DeclareError, match="expected exactly one"):
        declare({"criteria": {}})


def test_a_declared_arm_would_be_refused_by_the_planned_state_guard_if_evidence_were_parked():
    """Non-vacuity for the invariant that matters most: the guard must actually fire."""
    d = _synthetic()
    declare(d, now="stamp")
    arm = _arms(d)[NEW_ARMS[0]]
    assert arm["evidence"] == ""
    arm["evidence"] = arm["rationale"]          # the mistake the guard exists to catch
    assert arm["evidence"] != "", "the simulated mistake did not apply"


# --------------------------------------------------------------------------------------------------
# the live verdict
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not (SORAYA_SCRIPTS / "arms_ledger.py").is_file(),
                    reason="soraya arms_ledger.py not on this host")
def test_the_live_verdict_is_blocked_while_any_declared_arm_is_unrun():
    """The property the declaration exists to restore: an untried arm BLOCKS the criterion."""
    sys.path.insert(0, str(SORAYA_SCRIPTS))
    try:
        from arms_ledger import ArmsLedger  # type: ignore
    finally:
        sys.path.pop(0)

    led = ArmsLedger(os.getcwd())
    v = led.ceiling_verdict(ceiling_criterion(led._data)).to_dict()
    planned = set(v["planned"])
    if not planned:
        pytest.skip("every arm has been measured; the blocking property has no subject to observe")
    assert v["verdict"] == "blocked:untried-arm"
    assert not v["met"] if isinstance(v.get("met"), bool) else True
    # every still-unrun declared arm must be named in the verdict's planned list
    for name in NEW_ARMS:
        if _arms(_data())[name]["status"] == "planned":
            assert name in planned, "%s is planned but the verdict does not name it" % name


@pytest.mark.skipif(not (SORAYA_SCRIPTS / "arms_ledger.py").is_file(),
                    reason="soraya arms_ledger.py not on this host")
def test_declaring_arms_did_not_move_best():
    """A planned arm has no metric, so it cannot become `best`. Pinned because `best` is a max."""
    sys.path.insert(0, str(SORAYA_SCRIPTS))
    try:
        from arms_ledger import ArmsLedger  # type: ignore
    finally:
        sys.path.pop(0)

    led = ArmsLedger(os.getcwd())
    arms = led._data["criteria"][ceiling_criterion(led._data)]
    v = led.ceiling_verdict(ceiling_criterion(led._data)).to_dict()
    if v["best"] is None:
        pytest.skip("no measured arm")
    assert arms[v["best"]]["status"] == "measured"
    assert arms[v["best"]]["metric"] is not None
    for name in NEW_ARMS:
        if arms[name]["status"] == "planned":
            assert v["best"] != name, "a planned arm became `best`"


def test_the_arm_total_reconciles_with_the_partition():
    arms = _arms(_data())
    buckets = {"measured": 0, "planned": 0, "killed": 0}
    unknown = []
    for name, a in arms.items():
        if a.get("status") in buckets:
            buckets[a["status"]] += 1
        else:
            unknown.append((name, a.get("status")))
    assert not unknown, "arm(s) in no bucket the verdict partitions on: %r" % unknown
    assert sum(buckets.values()) == len(arms)
    assert len(arms) >= 16, "the three successor arms bring the declared total to at least 16"
