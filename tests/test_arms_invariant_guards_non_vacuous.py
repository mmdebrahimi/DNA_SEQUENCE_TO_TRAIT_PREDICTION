"""Proof that the 2026-11-02 invariant conversion did not weaken the arms guards.

Eight guards across three files were converted from EXACT-STATE pins to INVARIANT pins, because the
exact forms went red on every legitimate measurement and every new declaration. The risk of such a
conversion is obvious and is exactly what this file exists to exclude: an invariant can be written
so loosely that it always passes, and a guard that always passes looks identical from the outside to
one that works.

So every converted guard is run here against a DELIBERATELY BROKEN ledger and must FAIL. The real
guard functions are invoked -- not re-implementations of them -- by pointing the modules at a
mutated copy of the committed ledger, or at a stub whose verdict function disagrees with its own
file where the real one cannot be made to.

The second half asserts the converted forms hold in all THREE states this plan passes through: the
`met` state before the three arms were declared, the three-planned state after, and the one-planned
state after some are measured. An invariant that only holds in today's state is just an exact pin
with extra words.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

import tests.test_arms_ledger_stack_arm as g1
import tests.test_glm_tile_arms_declared as g2

REPO = Path(__file__).resolve().parents[1]
LEDGER = REPO / "soraya_runs" / "arms-ledger.json"
POSITIONAL = "positional-non-mean-pooled-promoteratlas-150bp"
SUCCESSORS = (
    "asgiven-position-resolved-promoteratlas-150bp",
    "shipped-expression-heads-zeroshot-150bp",
    "finetuned-promoteratlas-encoder-150bp",
)


def _ledger_at(tmp_path: Path, mutate=None) -> "g1.al.ArmsLedger":
    """A real `ArmsLedger` over a mutated copy of the committed file."""
    dest = tmp_path / "soraya_runs"
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(LEDGER, dest / "arms-ledger.json")
    if mutate is not None:
        d = json.loads((dest / "arms-ledger.json").read_text(encoding="utf-8"))
        crit = next(k for k in d["criteria"] if k.startswith(g1.CRITERION_PREFIX))
        mutate(d["criteria"][crit])
        (dest / "arms-ledger.json").write_text(json.dumps(d, indent=2, sort_keys=True),
                                               encoding="utf-8")
    return g1.al.ArmsLedger(cwd=str(tmp_path))


def _point_g1_at(monkeypatch, ledger):
    monkeypatch.setattr(g1, "_live_ledger", lambda: ledger)


def _arms_of(ledger) -> dict:
    return ledger._data["criteria"][g1._criterion(ledger)]


# --------------------------------------------------------------------------------------------------
# the subset pin still catches a moved published number
# --------------------------------------------------------------------------------------------------
def test_a_historical_metric_moved_by_one_ten_thousandth_is_CAUGHT(tmp_path, monkeypatch):
    """The protection the exact-equality pin really provided, preserved by the subset form."""
    led = _ledger_at(tmp_path, lambda arms: arms["gc-single-feature"].update(metric=0.3113))
    _point_g1_at(monkeypatch, led)
    with pytest.raises(AssertionError, match="published metric moved"):
        g1.test_no_historical_measured_metric_has_moved()


def test_a_historical_arm_leaving_the_measured_set_is_CAUGHT(tmp_path, monkeypatch):
    led = _ledger_at(tmp_path, lambda arms: arms["conv-3width-13k-params-150bp"].update(
        status="planned", metric=None, evidence=""))
    _point_g1_at(monkeypatch, led)
    with pytest.raises(AssertionError, match="no longer measured"):
        g1.test_no_historical_measured_metric_has_moved()


def test_a_measured_arm_losing_its_evidence_is_CAUGHT(tmp_path, monkeypatch):
    led = _ledger_at(tmp_path, lambda arms: arms["gc-single-feature"].update(evidence="   "))
    _point_g1_at(monkeypatch, led)
    with pytest.raises(AssertionError, match="lost its evidence"):
        g1.test_no_historical_measured_metric_has_moved()


def test_a_NEW_measured_arm_with_no_metric_is_CAUGHT(tmp_path, monkeypatch):
    """The subset form permits new measured arms -- but not unresolved ones."""
    led = _ledger_at(tmp_path, lambda arms: arms[SUCCESSORS[0]].update(
        status="measured", metric=None, evidence="wiki/x.json"))
    _point_g1_at(monkeypatch, led)
    with pytest.raises(AssertionError, match="measured with no metric"):
        g1.test_no_historical_measured_metric_has_moved()


# --------------------------------------------------------------------------------------------------
# the partition invariant still catches a bucket that does not reconcile
# --------------------------------------------------------------------------------------------------
def test_an_arm_with_an_off_enum_status_is_CAUGHT(tmp_path, monkeypatch):
    """What reading the three sets could never catch: an arm in none of them."""
    led = _ledger_at(tmp_path, lambda arms: arms[SUCCESSORS[0]].update(status="in-progress"))
    _point_g1_at(monkeypatch, led)
    with pytest.raises(AssertionError, match="no bucket the verdict partitions on"):
        g1.test_the_arm_partition_reconciles_with_the_whole()


def test_a_second_retired_arm_is_CAUGHT(tmp_path, monkeypatch):
    led = _ledger_at(tmp_path, lambda arms: arms[SUCCESSORS[1]].update(
        status="killed", evidence="synthetic", reason="synthetic"))
    _point_g1_at(monkeypatch, led)
    with pytest.raises(AssertionError, match="exactly one arm has ever been retired"):
        g1.test_the_arm_partition_reconciles_with_the_whole()


def test_the_declared_total_falling_below_its_floor_is_CAUGHT(tmp_path, monkeypatch):
    def gut(arms):
        for name in list(arms)[3:]:
            del arms[name]
    led = _ledger_at(tmp_path, gut)
    _point_g1_at(monkeypatch, led)
    with pytest.raises(AssertionError):
        g1.test_the_arm_partition_reconciles_with_the_whole()


# --------------------------------------------------------------------------------------------------
# the biconditional still catches a disagreement in BOTH directions
# --------------------------------------------------------------------------------------------------
class _StubVerdict:
    def __init__(self, verdict, planned, measured, best):
        self.verdict, self.planned, self.measured, self.best = verdict, planned, measured, best
        self.met = verdict == "met"


class _StubLedger:
    """A ledger whose verdict function DISAGREES with its own file.

    The real `ceiling_verdict` is a pure function and cannot be made to return `met` over a planned
    arm, which is precisely why the disagreement the biconditional guards against needs a stub to
    simulate at all.
    """

    def __init__(self, arms: dict, verdict: str, best=None):
        self._data = {"criteria": {g1.CRITERION_PREFIX + " synthetic": arms}}
        planned = sorted(n for n, a in arms.items() if a["status"] == "planned")
        measured = sorted(n for n, a in arms.items() if a["status"] == "measured")
        self._v = _StubVerdict(verdict, planned, measured, best)

    def ceiling_verdict(self, _crit):
        return self._v


def _real_arms() -> dict:
    d = json.loads(LEDGER.read_text(encoding="utf-8"))
    return d["criteria"][next(k for k in d["criteria"] if k.startswith(g1.CRITERION_PREFIX))]


def test_a_met_verdict_over_an_untried_arm_is_CAUGHT(monkeypatch):
    arms = json.loads(json.dumps(_real_arms()))
    assert any(a["status"] == "planned" for a in arms.values()), "fixture needs a planned arm"
    _point_g1_at(monkeypatch, _StubLedger(arms, "met"))
    with pytest.raises(AssertionError, match="must block the criterion"):
        g1.test_the_verdict_is_blocked_IF_AND_ONLY_IF_an_arm_is_untried()


def test_a_blocked_verdict_over_a_fully_measured_set_is_CAUGHT(monkeypatch):
    arms = json.loads(json.dumps(_real_arms()))
    for a in arms.values():
        if a["status"] == "planned":
            a.update(status="measured", metric=0.1, evidence="synthetic")
    _point_g1_at(monkeypatch, _StubLedger(arms, "blocked:untried-arm"))
    with pytest.raises(AssertionError, match="no arm is planned but the verdict reads"):
        g1.test_the_verdict_is_blocked_IF_AND_ONLY_IF_an_arm_is_untried()


def test_a_verdict_whose_planned_list_disagrees_with_the_file_is_CAUGHT(monkeypatch):
    arms = json.loads(json.dumps(_real_arms()))
    stub = _StubLedger(arms, "blocked:untried-arm")
    stub._v.planned = ["an-arm-that-is-not-in-the-file"]
    _point_g1_at(monkeypatch, stub)
    with pytest.raises(AssertionError, match="disagree about what is planned"):
        g1.test_the_verdict_is_blocked_IF_AND_ONLY_IF_an_arm_is_untried()


def test_a_best_that_is_not_the_argmax_is_CAUGHT(monkeypatch):
    arms = json.loads(json.dumps(_real_arms()))
    _point_g1_at(monkeypatch, _StubLedger(arms, "blocked:untried-arm", best="gc-single-feature"))
    with pytest.raises(AssertionError, match="is not the argmax"):
        g1.test_best_is_the_argmax_over_measured_and_a_leader_change_must_be_EXPLICIT()


# --------------------------------------------------------------------------------------------------
# the JSON-level blocking property, both branches
# --------------------------------------------------------------------------------------------------
def test_a_planned_arm_carrying_evidence_is_CAUGHT(monkeypatch):
    arms = json.loads(json.dumps(_real_arms()))
    name = next(n for n, a in arms.items() if a["status"] == "planned")
    arms[name]["evidence"] = "a rationale parked where executed evidence belongs"
    monkeypatch.setattr(g2, "_arms", lambda: arms)
    with pytest.raises(AssertionError, match="planned but carries evidence"):
        g2.test_the_blocking_property_holds_in_the_json_itself_whichever_state_it_is_in()


def test_a_planned_arm_carrying_a_metric_is_CAUGHT(monkeypatch):
    arms = json.loads(json.dumps(_real_arms()))
    name = next(n for n, a in arms.items() if a["status"] == "planned")
    arms[name]["metric"] = 0.42
    monkeypatch.setattr(g2, "_arms", lambda: arms)
    with pytest.raises(AssertionError, match="planned but carries a metric"):
        g2.test_the_blocking_property_holds_in_the_json_itself_whichever_state_it_is_in()


def test_a_SETTLED_criterion_with_an_unresolved_measured_arm_is_CAUGHT(monkeypatch):
    """The branch that only runs once nothing is planned -- proven reachable rather than assumed."""
    arms = json.loads(json.dumps(_real_arms()))
    for n, a in arms.items():
        if a["status"] == "planned":
            a.update(status="measured", metric=None, evidence="synthetic")
    monkeypatch.setattr(g2, "_arms", lambda: arms)
    with pytest.raises(AssertionError, match="measured with no metric"):
        g2.test_the_blocking_property_holds_in_the_json_itself_whichever_state_it_is_in()


# --------------------------------------------------------------------------------------------------
# the converted forms hold in ALL THREE states the plan passes through
# --------------------------------------------------------------------------------------------------
def _state_all_measured() -> dict:
    arms = json.loads(json.dumps(_real_arms()))
    for i, (n, a) in enumerate(arms.items()):
        if a["status"] == "planned":
            a.update(status="measured", metric=0.30 + i * 0.001, evidence="synthetic wiki/x.json")
    return arms


def _state_three_planned() -> dict:
    arms = json.loads(json.dumps(_real_arms()))
    for n in SUCCESSORS:
        arms[n].update(status="planned", metric=None, evidence="")
    return arms


def _state_one_planned() -> dict:
    arms = _state_three_planned()
    for n in SUCCESSORS[1:]:
        arms[n].update(status="measured", metric=0.4, evidence="synthetic wiki/x.json")
    return arms


@pytest.mark.parametrize("name,builder", [
    ("met / nothing planned", _state_all_measured),
    ("three planned", _state_three_planned),
    ("one planned", _state_one_planned),
])
def test_the_json_blocking_invariant_holds_in_every_state(name, builder, monkeypatch):
    monkeypatch.setattr(g2, "_arms", builder)
    g2.test_the_blocking_property_holds_in_the_json_itself_whichever_state_it_is_in()


@pytest.mark.parametrize("name,builder,verdict", [
    ("met / nothing planned", _state_all_measured, "met"),
    ("three planned", _state_three_planned, "blocked:untried-arm"),
    ("one planned", _state_one_planned, "blocked:untried-arm"),
])
def test_the_verdict_biconditional_holds_in_every_state(name, builder, verdict, monkeypatch):
    arms = builder()
    _point_g1_at(monkeypatch, _StubLedger(arms, verdict))
    g1.test_the_verdict_is_blocked_IF_AND_ONLY_IF_an_arm_is_untried()


@pytest.mark.parametrize("name,builder", [
    ("met / nothing planned", _state_all_measured),
    ("three planned", _state_three_planned),
    ("one planned", _state_one_planned),
])
def test_the_partition_invariant_holds_in_every_state(name, builder, monkeypatch, tmp_path):
    arms = builder()

    class _L:
        _data = {"criteria": {g1.CRITERION_PREFIX + " synthetic": arms}}

    _point_g1_at(monkeypatch, _L())
    g1.test_the_arm_partition_reconciles_with_the_whole()
    g1.test_no_historical_measured_metric_has_moved()
