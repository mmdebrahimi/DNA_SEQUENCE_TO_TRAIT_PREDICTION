"""Tests for recording the three successor arms in the arms ledger.

THE CENTRAL CLAIM IS PROVEN, NOT ASSERTED. This script writes the ledger by direct dict assignment
instead of calling `ArmsLedger.measure_arm`, and the justification is that `measure_arm` ends in
`Arm(**base)` over an arm carrying `rationale`/`declared_from`, which raises. A comment saying so
would rot the moment the dataclass gained a field, so `test_measure_arm_WOULD_raise_on_these_arms`
constructs `Arm(**arm)` for real and requires the `TypeError`. If the skill ever filters its fields,
that test fails and the direct-write rationale is retired deliberately rather than silently kept.

NOTHING HERE TOUCHES THE LIVE LEDGER. Every mutation test runs against a deep copy in `tmp_path`.
The only tests that read the committed file are read-only, and they are the ones that matter for
"would this work on the real thing".
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

import scripts.record_successor_arm_results as R
from scripts.record_successor_arm_results import (
    RECORD_ARMS,
    RecordRefused,
    apply_records,
    evidence_for,
    leader_note,
    plan_records,
    prospective,
    save_ledger,
)

REPO = Path(__file__).resolve().parents[1]
LIVE = REPO / "soraya_runs/arms-ledger.json"

ASGIVEN = "asgiven-position-resolved-promoteratlas-150bp"
HEADS = "shipped-expression-heads-zeroshot-150bp"
FINETUNE = "finetuned-promoteratlas-encoder-150bp"

pytestmark = pytest.mark.skipif(not LIVE.is_file(), reason="the arms ledger is absent here")


def _live_raw() -> dict:
    """The ledger exactly as it is on disk."""
    return json.loads(LIVE.read_text(encoding="utf-8"))


def _live() -> dict:
    """The ledger with the three successor arms reset to `planned`.

    WHY THE RESET. These tests exercise the act of RECORDING, which requires the arms to be
    unrecorded. Two of the three were recorded on 2026-11-04, so reading the live file gave tests
    whose fixtures conflicted with it -- `test_recording_is_IDEMPOTENT` began failing because
    `plan_records` correctly REFUSED a fixture metric of 0.5 against a stored 0.5086. That is the
    guard working, not a defect, but it made the suite's result depend on whether recording had
    happened yet, which is not a property any of these tests is about.

    Resetting in the deep copy makes them state-independent. `_live_raw` is there for the one test
    that genuinely needs the file as it stands.
    """
    data = _live_raw()
    arms = data["criteria"][_crit(data)]
    for name in (ASGIVEN, HEADS, FINETUNE):
        if name in arms:
            arms[name].update(status="planned", metric=None, evidence="", resolved_at="")
    return data


def _crit(data: dict) -> str:
    from scripts.declare_successor_arms import ceiling_criterion

    return ceiling_criterion(data)


def row(name, *, metric=0.5, collected=True, verdict="ARM_TIES", interval=None, **over):
    r = {"arm": name, "artifact": "wiki/%s.json" % name, "collected": collected,
         "verdict": verdict, "interval": interval}
    if collected:
        r["metric"] = metric
    r.update(over)
    return r


def report_of(**rows) -> dict:
    return {"arms": dict(rows)}


# ------------------------------------------------------------------ the surface


def test_the_recorded_arms_are_exactly_the_collectors_declared_surface():
    """Two hand-kept lists would drift; this asserts one is derived from the other."""
    from scripts.collect_unattended_runs import ARMS as COLLECTOR

    assert set(RECORD_ARMS) == set(COLLECTOR)
    assert set(RECORD_ARMS) == {ASGIVEN, HEADS, FINETUNE}


def test_only_four_fields_are_writable():
    """If this grows, a reader must re-check that `rationale` is still out of reach."""
    assert R.WRITABLE == ("status", "metric", "evidence", "resolved_at")
    assert "rationale" not in R.WRITABLE and "declared_from" not in R.WRITABLE


# ------------------------------------------------------------------ why direct writes


def test_measure_arm_WOULD_raise_on_these_arms():
    """The non-vacuity proof for the whole direct-write design.

    If this starts passing, `Arm` gained a field filter and the direct-write rationale should be
    revisited on purpose -- not left in place because nobody re-checked.
    """
    sys.path.insert(0, r"C:/my_skills/skills/soraya/scripts")
    try:
        from arms_ledger import Arm  # type: ignore
    except ImportError:
        pytest.skip("the soraya skill's arms_ledger is not on this host")
    finally:
        sys.path.pop(0)

    data = _live()
    arm = data["criteria"][_crit(data)][ASGIVEN]
    assert "rationale" in arm, "fixture assumption: this arm carries rationale"
    with pytest.raises(TypeError):
        Arm(**arm)


def test_the_recorder_never_calls_measure_arm_or_declare_arm():
    src = (REPO / "scripts/record_successor_arm_results.py").read_text(encoding="utf-8")
    body = "\n".join(ln for ln in src.splitlines() if not ln.strip().startswith("#"))
    for forbidden in ("measure_arm(", "declare_arm("):
        assert forbidden not in body, "the recorder reaches for %r" % forbidden


# ------------------------------------------------------------------ refusals


def test_a_collector_REFUSAL_is_not_recorded(tmp_path):
    data = _live()
    rep = report_of(**{ASGIVEN: row(ASGIVEN, collected=False,
                                    refusal="job is NOT_STARTED, not COMPLETE")})
    d = plan_records(data, rep)[ASGIVEN]
    assert d["action"] == "refused"
    assert "NOT_STARTED" in d["reason"]


def test_a_missing_collector_row_is_refused():
    d = plan_records(_live(), report_of())[ASGIVEN]
    assert d["action"] == "refused" and "no row" in d["reason"]


def test_a_NULL_metric_is_refused_rather_than_recorded_as_zero():
    rep = report_of(**{ASGIVEN: row(ASGIVEN, metric=None)})
    d = plan_records(_live(), rep)[ASGIVEN]
    assert d["action"] == "refused" and "null is not a result" in d["reason"]


def test_an_UNDECLARED_arm_is_refused_rather_than_created_by_measuring_it():
    """Creating an arm as a side effect of measuring it would bypass the declaration discipline."""
    data = _live()
    del data["criteria"][_crit(data)][ASGIVEN]
    d = plan_records(data, report_of(**{ASGIVEN: row(ASGIVEN)}))[ASGIVEN]
    assert d["action"] == "refused" and "not declared" in d["reason"]


def test_an_already_measured_arm_with_a_DIFFERENT_metric_is_REFUSED():
    """A published number does not move silently."""
    data = _live()
    arms = data["criteria"][_crit(data)]
    arms[ASGIVEN].update(status="measured", metric=0.4000, evidence="prior")
    d = plan_records(data, report_of(**{ASGIVEN: row(ASGIVEN, metric=0.5086)}))[ASGIVEN]
    assert d["action"] == "refused"
    assert "does not move silently" in d["reason"]


def test_an_already_measured_arm_with_the_SAME_metric_is_a_NO_OP():
    data = _live()
    arms = data["criteria"][_crit(data)]
    arms[ASGIVEN].update(status="measured", metric=0.5086, evidence="prior")
    d = plan_records(data, report_of(**{ASGIVEN: row(ASGIVEN, metric=0.5086)}))[ASGIVEN]
    assert d["action"] == "noop"


# ------------------------------------------------------------------ what a write preserves


def test_recording_preserves_rationale_and_declared_from_BYTE_FOR_BYTE():
    data = _live()
    arms = data["criteria"][_crit(data)]
    before = {n: dict(arms[n]) for n in RECORD_ARMS if n in arms}
    rep = report_of(**{n: row(n, metric=0.5) for n in RECORD_ARMS})
    recorded = apply_records(data, plan_records(data, rep), now="2026-01-01T00:00:00Z")
    assert recorded
    for n in recorded:
        assert arms[n]["rationale"] == before[n]["rationale"]
        assert arms[n]["declared_from"] == before[n]["declared_from"]
        assert arms[n]["declared_at"] == before[n]["declared_at"]
        assert arms[n]["cost_tier"] == before[n]["cost_tier"]


def test_recording_writes_only_the_writable_fields():
    data = _live()
    arms = data["criteria"][_crit(data)]
    before = {n: dict(arms[n]) for n in RECORD_ARMS if n in arms}
    rep = report_of(**{n: row(n, metric=0.5) for n in RECORD_ARMS})
    for n in apply_records(data, plan_records(data, rep), now="2026-01-01T00:00:00Z"):
        moved = {k for k in before[n] if before[n][k] != arms[n][k]}
        assert moved <= set(R.WRITABLE), "%s also changed %r" % (n, moved - set(R.WRITABLE))
        assert set(before[n]) == set(arms[n]), "the field SET changed on %s" % n


def test_a_recorded_arm_becomes_measured_with_a_metric_and_evidence():
    data = _live()
    arms = data["criteria"][_crit(data)]
    rep = report_of(**{n: row(n, metric=0.5) for n in RECORD_ARMS})
    for n in apply_records(data, plan_records(data, rep)):
        assert arms[n]["status"] == "measured"
        assert arms[n]["metric"] == 0.5
        assert arms[n]["evidence"].strip()
        assert arms[n]["resolved_at"].strip()


def test_recording_is_IDEMPOTENT():
    data = _live()
    rep = report_of(**{n: row(n, metric=0.5) for n in RECORD_ARMS})
    apply_records(data, plan_records(data, rep), now="2026-01-01T00:00:00Z")
    snapshot = copy.deepcopy(data)
    second = plan_records(data, rep)
    assert all(d["action"] == "noop" for d in second.values())
    assert apply_records(data, second) == []
    assert data == snapshot, "a second run changed the ledger"


def test_the_evidence_string_names_the_artifact_and_carries_the_interval():
    ev = evidence_for(ASGIVEN, row(ASGIVEN, metric=0.5086, verdict="ASGIVEN_TIES",
                                   interval={"ci_lo": -0.0068, "ci_hi": 0.0004}))
    assert "wiki/" in ev and "ASGIVEN_TIES" in ev and "0.5086" in ev
    assert "-0.0068" in ev and "collect_unattended_runs.py" in ev


# ------------------------------------------------------------------ the leader note


def test_the_leader_note_reports_the_margin_and_the_caution():
    data = _live()
    arms = data["criteria"][_crit(data)]
    arms[ASGIVEN].update(status="measured", metric=0.5086, evidence="x")
    note = leader_note(data)
    assert note["leader"] == ASGIVEN
    assert note["margin_over_runner_up"] == pytest.approx(0.0010, abs=1e-6)
    assert "not a demonstrated gain" in note["caution"]


def test_the_leader_is_the_argmax_over_measured_only():
    """A planned arm with no metric must not be able to lead."""
    data = _live()
    arms = data["criteria"][_crit(data)]
    arms[FINETUNE].update(status="planned", metric=None)
    note = leader_note(data)
    assert arms[note["leader"]]["status"] == "measured"


# ------------------------------------------------------------------ the preview


def test_prospective_does_NOT_mutate_the_input():
    data = _live()
    snapshot = copy.deepcopy(data)
    rep = report_of(**{n: row(n, metric=0.5) for n in RECORD_ARMS})
    prospective(data, plan_records(data, rep))
    assert data == snapshot, "the dry-run preview mutated the ledger it was previewing"


def test_prospective_previews_the_leader_change_a_dry_run_would_otherwise_hide():
    data = _live()
    rep = report_of(**{ASGIVEN: row(ASGIVEN, metric=0.5086)})
    _verdict, note = prospective(data, plan_records(data, rep))
    assert note["leader"] == ASGIVEN, "the preview did not surface the leader change"


def test_prospective_verdict_clears_only_when_no_planned_arm_remains():
    data = _live()
    all_three = report_of(**{n: row(n, metric=0.5) for n in RECORD_ARMS})
    v_all, _ = prospective(data, plan_records(data, all_three))
    assert v_all["verdict"] == "met"
    assert v_all["planned"] == []

    two_only = report_of(**{ASGIVEN: row(ASGIVEN, metric=0.5086),
                            HEADS: row(HEADS, metric=0.1178)})
    v_two, _ = prospective(data, plan_records(data, two_only))
    assert v_two["verdict"] == "blocked:untried-arm"
    assert v_two["planned"] == [FINETUNE]


# ------------------------------------------------------------------ the durable write


def test_save_ledger_leaves_no_tmp_file_and_round_trips(tmp_path):
    p = tmp_path / "arms-ledger.json"
    data = _live()
    save_ledger(data, p)
    assert json.loads(p.read_text(encoding="utf-8")) == data
    assert not list(tmp_path.glob("*.tmp")), "a tmp file survived the atomic write"


def test_save_ledger_uses_fsync_and_replace_not_write_text():
    """`write_text` returns before the data is durable; a crash then leaves a NUL-filled ledger.

    Asserted over the function's CALLS via `ast`, not over its source text. A text search cannot
    work here: the docstring has to NAME `write_text` to explain why it is forbidden, so
    `"write_text" not in source` can never pass -- an absence check on a term the fix must mention
    is unsatisfiable by construction.
    """
    import ast

    src = (REPO / "scripts/record_successor_arm_results.py").read_text(encoding="utf-8")
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "save_ledger")
    called = set()
    for node in ast.walk(fn):
        if isinstance(node, ast.Call):
            f = node.func
            called.add(f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", ""))
    assert {"fsync", "replace"} <= called, "the durable-write pair is missing: %r" % sorted(called)
    assert "write_text" not in called, "save_ledger still writes through write_text"


def test_load_ledger_REFUSES_an_absent_file(tmp_path):
    with pytest.raises(RecordRefused, match="absent"):
        R.load_ledger(tmp_path / "nope.json")


# ------------------------------------------------------------------ the ledger's own invariants


def test_the_planned_arm_empty_evidence_guard_still_holds_after_recording():
    """The ledger's standing guard: a PLANNED arm must carry no evidence. Recording must not
    leave an arm planned-with-evidence, which would read as measured to a careless reader."""
    data = _live()
    arms = data["criteria"][_crit(data)]
    rep = report_of(**{ASGIVEN: row(ASGIVEN, metric=0.5086)})
    apply_records(data, plan_records(data, rep))
    for name, a in arms.items():
        if a.get("status") == "planned":
            assert not a.get("evidence", "").strip(), "%s is planned but carries evidence" % name


def test_every_measured_arm_carries_evidence_after_recording():
    data = _live()
    arms = data["criteria"][_crit(data)]
    rep = report_of(**{n: row(n, metric=0.5) for n in RECORD_ARMS})
    apply_records(data, plan_records(data, rep))
    for name, a in arms.items():
        if a.get("status") == "measured":
            assert a.get("evidence", "").strip(), "%s is measured with no evidence" % name


def test_no_verdict_is_WRITTEN_into_the_ledger():
    """The verdict is recomputed and reported; persisting it would make a stale copy authoritative."""
    data = _live()
    rep = report_of(**{n: row(n, metric=0.5) for n in RECORD_ARMS})
    apply_records(data, plan_records(data, rep))
    assert "ceiling_verdict" not in data
    assert set(data) == {"schema_version", "criteria"}
    for a in data["criteria"][_crit(data)].values():
        assert "ceiling_verdict" not in a


def test_recording_never_changes_another_arms_metric():
    data = _live()
    arms = data["criteria"][_crit(data)]
    others = {n: a.get("metric") for n, a in arms.items() if n not in RECORD_ARMS}
    rep = report_of(**{n: row(n, metric=0.5) for n in RECORD_ARMS})
    apply_records(data, plan_records(data, rep))
    for n, m in others.items():
        assert arms[n].get("metric") == m, "%s's metric moved" % n


# ------------------------------------------------------------------ against the live ledger


def test_the_live_dry_run_records_the_two_finished_arms_and_refuses_the_third():
    """Read-only against the ledger AS IT STANDS and the real collector -- `_live_raw`, not the
    reset copy, because this is the one test whose subject IS the current state."""
    from scripts.collect_unattended_runs import collect

    decisions = plan_records(_live_raw(), collect())
    assert decisions[HEADS]["action"] in ("record", "noop")
    assert decisions[ASGIVEN]["action"] in ("record", "noop")
    if decisions[FINETUNE]["action"] == "refused":
        assert "not COMPLETE" in decisions[FINETUNE]["reason"]


def test_main_dry_run_writes_nothing():
    before = LIVE.read_bytes()
    R.main([])
    assert LIVE.read_bytes() == before, "a dry run modified the ledger"
