"""Tests for the three-arm collector.

MOSTLY REFUSALS, and deliberately so. The collector's job is to hand a metric to the thing that
writes the arms ledger, so its dangerous failure is not an error but a NULL that reads like a
measurement. Every path that cannot produce a trustworthy number must refuse WITH A REASON instead.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pytest

import scripts.collect_unattended_runs as C
from scripts.collect_unattended_runs import CollectRefused

ASGIVEN = "asgiven-position-resolved-promoteratlas-150bp"
HEADS = "shipped-expression-heads-zeroshot-150bp"
FINETUNE = "finetuned-promoteratlas-encoder-150bp"


# --------------------------------------------------------------------- dotted paths


def test_dotted_resolves_a_nested_path():
    assert C.dotted({"a": {"b": {"c": 7}}}, "a.b.c") == 7


def test_dotted_RAISES_rather_than_returning_None():
    """A missing metric must become a refusal with a reason, not a null that reads like a
    measurement. Returning None here is how a null reaches a ledger as a number."""
    with pytest.raises(CollectRefused, match="does not resolve"):
        C.dotted({"a": {"b": 1}}, "a.x.c")


def test_dotted_names_where_it_stopped():
    with pytest.raises(CollectRefused, match="stopped at 'x'"):
        C.dotted({"a": {}}, "a.x")


def test_dotted_refuses_to_walk_through_a_non_dict():
    with pytest.raises(CollectRefused):
        C.dotted({"a": 5}, "a.b")


# --------------------------------------------------------------------- fixtures


def _spec(tmp_path, *, schema="s-v1", job=None, metric_path="m", interval_path="iv"):
    return {"artifact": "art.json", "schema": schema, "metric_path": metric_path,
            "interval_path": interval_path, "job": job,
            "why_no_job": "attended run, no launcher"}


def _write_art(tmp_path, payload):
    (tmp_path / "art.json").write_text(json.dumps(payload), encoding="utf-8")


def _hb(state_dir: Path, job: str, *, pid: int, age=0.0, cadence=60):
    state_dir.mkdir(parents=True, exist_ok=True)
    (state_dir / ("%s.heartbeat.json" % job)).write_text(json.dumps({
        "job": job, "pid": pid, "units_done": 5, "units_total": 10,
        "written_epoch": time.time() - age, "cadence_seconds": cadence}), encoding="utf-8")


def _sent(state_dir: Path, job: str, rc: int):
    state_dir.mkdir(parents=True, exist_ok=True)
    (state_dir / ("%s.sentinel" % job)).write_text("EXIT=%d\r\n" % rc, encoding="utf-8")


def _dead_pid() -> int:
    import subprocess
    import sys as _s
    p = subprocess.Popen([_s.executable, "-c", "pass"],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p.wait()
    return p.pid


GOOD = {"schema": "s-v1", "m": 0.42, "iv": {"ci_lo": -0.01, "ci_hi": 0.05},
        "verdict": "ARM_TIES", "verdict_detail": "x",
        "verdict_is_rederivable_from_this_artifact": True}


# --------------------------------------------------------------------- the happy path


def test_an_attended_arm_with_NO_JOB_collects(tmp_path):
    """Requiring a COMPLETE launcher state would refuse the as-given arm, which was scored
    attended and has no launcher, no heartbeat and no sentinel. That is not a defect."""
    _write_art(tmp_path, GOOD)
    row = C.collect_arm("a", _spec(tmp_path, job=None), repo=tmp_path)
    assert row["collected"] is True
    assert row["metric"] == 0.42
    assert row["job_state"] == "n/a"
    assert "attended" in row["job_detail"]


def test_a_job_backed_arm_collects_when_COMPLETE(tmp_path):
    _write_art(tmp_path, GOOD)
    state = tmp_path / "st"
    _sent(state, "j", 0)
    row = C.collect_arm("a", _spec(tmp_path, job="j"), repo=tmp_path, state_dir=state)
    assert row["collected"] is True and row["job_state"] == "COMPLETE"


def test_the_interval_is_carried_through(tmp_path):
    _write_art(tmp_path, GOOD)
    row = C.collect_arm("a", _spec(tmp_path), repo=tmp_path)
    assert row["interval"] == {"ci_lo": -0.01, "ci_hi": 0.05}


def test_a_missing_interval_is_reported_absent_but_NOT_fatal(tmp_path):
    """Not every arm carries a bootstrap interval; a missing one is information, not a refusal."""
    _write_art(tmp_path, {k: v for k, v in GOOD.items() if k != "iv"})
    row = C.collect_arm("a", _spec(tmp_path), repo=tmp_path)
    assert row["collected"] is True and row["interval"] is None


# --------------------------------------------------------------------- refusals


@pytest.mark.parametrize("state,setup", [
    ("NOT_STARTED", lambda s: None),
    ("FAILED", lambda s: _sent(s, "j", 1)),
    ("REFUSED_LOCK_HELD", lambda s: _sent(s, "j", 9)),
    ("RUNNING", lambda s: _hb(s, "j", pid=os.getpid(), age=1.0)),
    ("STALLED", lambda s: _hb(s, "j", pid=os.getpid(), age=5000, cadence=60)),
    ("DEAD_NO_SENTINEL", lambda s: _hb(s, "j", pid=_dead_pid())),
])
def test_every_NON_COMPLETE_launcher_state_REFUSES(tmp_path, state, setup):
    """Only COMPLETE licenses reading an artifact as the arm's result. A RUNNING job's artifact is
    from a previous run; a DEAD one's is whatever survived."""
    _write_art(tmp_path, GOOD)
    s = tmp_path / "st"
    s.mkdir(parents=True, exist_ok=True)
    setup(s)
    row = C.collect_arm("a", _spec(tmp_path, job="j"), repo=tmp_path, state_dir=s)
    assert row["collected"] is False
    assert row["job_state"] == state
    assert "not COMPLETE" in row["refusal"]
    assert "metric" not in row, "a refused arm must carry NO metric"


def test_an_absent_artifact_refuses(tmp_path):
    row = C.collect_arm("a", _spec(tmp_path), repo=tmp_path)
    assert row["collected"] is False and "absent" in row["refusal"]


def test_malformed_json_refuses(tmp_path):
    (tmp_path / "art.json").write_text('{"schema": "s-v1", "m"', encoding="utf-8")
    row = C.collect_arm("a", _spec(tmp_path), repo=tmp_path)
    assert row["collected"] is False and "not readable JSON" in row["refusal"]


def test_an_UNRECOGNISED_SCHEMA_REFUSES_rather_than_guessing(tmp_path):
    """The load-bearing refusal. A default extractor would read None on an unfamiliar shape and
    report it as a measured absence of signal -- the same trap as mapping NaN to 0.0."""
    _write_art(tmp_path, dict(GOOD, schema="something-else-v9"))
    row = C.collect_arm("a", _spec(tmp_path), repo=tmp_path)
    assert row["collected"] is False
    assert "expected" in row["refusal"]
    assert "guessing an extractor" in row["refusal"]
    assert "metric" not in row


def test_an_unresolvable_metric_path_refuses(tmp_path):
    _write_art(tmp_path, {"schema": "s-v1", "other": 1})
    row = C.collect_arm("a", _spec(tmp_path), repo=tmp_path)
    assert row["collected"] is False and "metric unreadable" in row["refusal"]


def test_a_NULL_metric_refuses(tmp_path):
    """A null metric is a MISSING result, not a result of zero."""
    _write_art(tmp_path, dict(GOOD, m=None))
    row = C.collect_arm("a", _spec(tmp_path), repo=tmp_path)
    assert row["collected"] is False
    assert "not a result" in row["refusal"]


# --------------------------------------------------------------------- the roll-up


def test_all_collected_is_False_while_any_arm_refuses(tmp_path):
    _write_art(tmp_path, GOOD)
    arms = {"a": _spec(tmp_path), "b": dict(_spec(tmp_path), artifact="missing.json")}
    rep = C.collect(repo=tmp_path, arms=arms)
    assert rep["all_collected"] is False
    assert rep["n_collected"] == 1 and rep["n_refused"] == 1


def test_all_collected_is_True_only_when_every_arm_is_in(tmp_path):
    _write_art(tmp_path, GOOD)
    rep = C.collect(repo=tmp_path, arms={"a": _spec(tmp_path)})
    assert rep["all_collected"] is True and rep["n_refused"] == 0


def test_the_report_says_it_records_nothing(tmp_path):
    """The roll-up is read-only; a refused arm must never reach the ledger."""
    _write_art(tmp_path, GOOD)
    rep = C.collect(repo=tmp_path, arms={"a": _spec(tmp_path)})
    assert "READ-ONLY" in rep["what_this_is_not"]
    assert "must not be recorded" in rep["what_this_is_not"]


def test_render_marks_refusals_distinctly(tmp_path):
    _write_art(tmp_path, GOOD)
    arms = {"ok": _spec(tmp_path), "bad": dict(_spec(tmp_path), artifact="missing.json")}
    text = C.render(C.collect(repo=tmp_path, arms=arms))
    assert "OK " in text and "REFUSED" in text


# --------------------------------------------------------------------- the real tree


def test_the_declared_surface_covers_exactly_the_three_successor_arms():
    assert set(C.ARMS) == {ASGIVEN, HEADS, FINETUNE}


def test_only_COMPLETE_licenses_reading():
    """If this ever grows, a reader must re-check which states are being trusted."""
    assert C.TERMINAL_OK == ("COMPLETE",)


def test_the_asgiven_arm_has_no_job_declared():
    """It was scored attended. Declaring a job for it would refuse a committed result."""
    assert C.ARMS[ASGIVEN]["job"] is None
    assert C.ARMS[ASGIVEN]["why_no_job"]


REAL_ASGIVEN = REPO_ASGIVEN = Path(__file__).resolve().parents[1] / \
    "wiki/asgiven_arm_result_2026-11-02.json"


@pytest.mark.skipif(not REAL_ASGIVEN.is_file(), reason="the as-given artifact is absent here")
def test_the_REAL_asgiven_metric_matches_an_INDEPENDENT_read():
    """Real-data check: the collector's extracted metric must equal the artifact field read
    directly, so the declared path is verified against the file rather than trusted."""
    art = json.loads(REAL_ASGIVEN.read_text(encoding="utf-8"))
    expected = art["rung"]["median_spearman"]
    row = C.collect_arm(ASGIVEN, C.ARMS[ASGIVEN])
    assert row["collected"] is True
    assert row["metric"] == expected
    assert row["verdict"] == art["verdict"] == "ASGIVEN_TIES"


def test_the_collector_never_reads_the_heads_SCALAR_BLOCK():
    """The heads block refuses a capped run by design (`load_scalars` raises unless the sidecar
    reads COMPLETE). Reading it directly here would defeat that refusal."""
    src = (Path(__file__).resolve().parents[1]
           / "scripts/collect_unattended_runs.py").read_text(encoding="utf-8")
    for forbidden in ("load_scalars", "head_scalars", "expression_heads_block"):
        assert forbidden not in src, "collector reaches into the scalar block via %r" % forbidden
