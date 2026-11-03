"""Tests for the detached-launch harness.

TWO OF THESE ARE NOT UNIT TESTS AND CANNOT BE. The sentinel-redirect trap is a property of
`cmd.exe`'s parser, and process survival past a closed shell is a property of the OS -- neither
is visible to an in-process assertion. This repo has already been bitten twice by a guard proven
on a path its only caller does not take (`glm_run_lock`'s `acquire`-vs-`run` defect and its
argparse `--` defect, both invisible to in-process tests). So the trap is RUN through real cmd,
and the survival check is opt-in and was run for real once rather than simulated.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import detached_launch as dl  # noqa: E402
from glm_run_lock import HELD_EXIT  # noqa: E402

WINDOWS_ONLY = pytest.mark.skipif(os.name != "nt", reason="cmd.exe parser behaviour")


# --------------------------------------------------------------------- the redirect trap


def _run_sentinel_line(tmp_path: Path, line: str, rc: str = "0") -> bytes:
    """Write a minimal .bat containing exactly one sentinel line, run it, return the bytes."""
    sentinel = tmp_path / "s.sentinel"
    bat = tmp_path / "t.bat"
    bat.write_bytes("\r\n".join([
        "@echo off",
        'set "SENTINEL=%s"' % sentinel,
        'set "RC=%s"' % rc,
        line,
        "",
    ]).encode("ascii"))
    subprocess.run(["cmd", "/c", str(bat)], capture_output=True, timeout=60)
    return sentinel.read_bytes() if sentinel.is_file() else b"__ABSENT__"


@WINDOWS_ONLY
def test_the_natural_ordering_produces_an_EMPTY_sentinel():
    """The trap, measured rather than asserted.

    `echo EXIT=%RC%>>"%SENTINEL%"` -> `echo EXIT=0>>"file"`, and cmd reads the `0` as a FILE
    HANDLE. If this ever stops being empty on some future cmd build, the whole reason the
    generator reverses the order has evaporated and this test is the place that says so.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        got = _run_sentinel_line(Path(d), dl.SENTINEL_LINE_TRAP, rc="0")
    assert got.strip() == b"", (
        "the trap no longer reproduces: cmd wrote %r. The generator's reversed ordering is "
        "justified by this behaviour -- re-derive it before relying on the comment." % got)


@WINDOWS_ONLY
def test_the_reversed_ordering_records_the_exit_code():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        got = _run_sentinel_line(Path(d), dl.SENTINEL_LINE_CORRECT, rc="0")
    assert b"EXIT=0" in got


@WINDOWS_ONLY
def test_the_reversed_ordering_records_a_NONZERO_exit_code_too():
    """A single-digit rc is the trap's trigger; a non-zero one must survive the fix as well."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        got = _run_sentinel_line(Path(d), dl.SENTINEL_LINE_CORRECT, rc="7")
    assert b"EXIT=7" in got


@WINDOWS_ONLY
def test_an_empty_sentinel_reads_as_UNKNOWN_not_as_success(tmp_path):
    """The trap's output must never classify as COMPLETE.

    This is the load-bearing consequence: an empty sentinel parsed as exit 0 converts a loud
    missing exit code into a silent false success, which is strictly worse than the bug.
    """
    p = tmp_path / "x.sentinel"
    p.write_bytes(b"")
    assert dl.read_sentinel(p) is None
    p.write_bytes(b"\r\n  \r\n")
    assert dl.read_sentinel(p) is None


def test_a_garbled_sentinel_reads_as_unknown(tmp_path):
    p = tmp_path / "x.sentinel"
    p.write_text("EXIT=not-a-number", encoding="utf-8")
    assert dl.read_sentinel(p) is None
    p.write_text("something else entirely", encoding="utf-8")
    assert dl.read_sentinel(p) is None


def test_a_sentinel_with_leading_noise_still_parses(tmp_path):
    p = tmp_path / "x.sentinel"
    p.write_text("some banner line\r\nEXIT=3\r\n", encoding="utf-8")
    assert dl.read_sentinel(p) == 3


# --------------------------------------------------------------------- bat generation


def _bat(**kw) -> str:
    return dl.render_bat(kw.pop("job", "demo"),
                         kw.pop("command", ["uv", "run", "python", "scripts/x.py", "--n", "3"]),
                         **kw)


def test_the_generated_bat_wraps_the_payload_in_the_run_lock():
    text = _bat()
    assert "glm_run_lock.py run" in text
    assert "scripts/x.py" in text


def test_the_generated_bat_uses_run_not_acquire():
    """`acquire` from a launcher records a pid that is dead before the next launch looks -- the
    defect glm_run_lock's own docstring records. A generator must not reintroduce it."""
    text = _bat()
    assert "glm_run_lock.py acquire" not in text


def test_the_generated_bat_invokes_the_lock_through_uv_run():
    """Bare `python` picks whichever interpreter is first on PATH. The argparse `--` defect
    reproduces under the venv's 3.11 and not under 3.14, so the launcher must pin the
    interpreter the project actually uses."""
    assert "uv run python scripts/glm_run_lock.py" in _bat()


def test_the_generated_bat_sets_unbuffered_output():
    assert "set PYTHONUNBUFFERED=1" in _bat()


def test_the_generated_bat_deletes_a_stale_sentinel_before_starting():
    assert 'if exist "%SENTINEL%" del "%SENTINEL%"' in _bat()


def test_the_generated_bat_uses_the_reversed_sentinel_ordering():
    text = _bat()
    assert dl.SENTINEL_LINE_CORRECT in text
    assert dl.SENTINEL_LINE_TRAP not in text


def test_the_generated_bat_propagates_the_exit_code():
    text = _bat()
    assert 'set "RC=%ERRORLEVEL%"' in text
    assert "exit /b %RC%" in text


def test_the_generated_bat_anchors_paths_to_its_own_directory():
    """`%~dp0` not a relative path: the launcher is started with -WorkingDirectory set by the
    caller, and a relative log path would land wherever that happened to be."""
    text = _bat()
    assert 'cd /d "%~dp0"' in text
    assert 'set "LOG=%~dp0data\\unattended\\demo.log"' in text
    assert 'set "SENTINEL=%~dp0data\\unattended\\demo.sentinel"' in text


def test_the_generated_bat_sets_requested_env_vars_before_the_payload():
    """`PA_REPO` is read at MODULE IMPORT time by atlas_embed, so it must be in the environment
    before python starts -- a runtime argument could not reach it."""
    text = _bat(env={"PA_REPO": "C:/repo/data/pa_clone_durable"})
    assert 'set "PA_REPO=C:/repo/data/pa_clone_durable"' in text
    lines = text.splitlines()
    assert lines.index('set "PA_REPO=C:/repo/data/pa_clone_durable"') < \
        next(i for i, ln in enumerate(lines) if "glm_run_lock.py run" in ln)


def test_env_vars_are_set_inside_setlocal_so_they_cannot_leak():
    text = _bat(env={"PA_REPO": "C:/x"})
    lines = text.splitlines()
    assert lines.index("setlocal") < lines.index('set "PA_REPO=C:/x"')
    assert lines.index('set "PA_REPO=C:/x"') < lines.index("endlocal")


def test_env_values_are_QUOTED_so_a_path_with_spaces_survives():
    """This host's repo path contains a space. The unquoted `set K=V` form would truncate at it
    and the job would start with a silently wrong root."""
    text = _bat(env={"PA_REPO": "C:/Users/Laptop Portable/x"})
    assert 'set "PA_REPO=C:/Users/Laptop Portable/x"' in text


def test_a_bat_with_no_env_is_unchanged():
    """The env block is omitted entirely when empty, so the generated text for the existing jobs
    does not move."""
    assert _bat(env={}) == _bat()
    assert _bat(env=None) == _bat()


@pytest.mark.parametrize("bad", [{"A B": "x"}, {"A=B": "x"}, {'A"B': "x"}, {"A%B": "x"}])
def test_an_unsettable_env_NAME_is_refused(bad):
    with pytest.raises(dl.DetachedLaunchError, match="not settable"):
        dl.render_bat("demo", ["cmd", "/c", "exit", "0"], env=bad)


@pytest.mark.parametrize("bad", [{"A": 'x"y'}, {"A": "50%done"}])
def test_an_env_VALUE_cmd_would_mangle_is_refused(bad):
    """Refusing beats starting a multi-day job with a silently wrong value."""
    with pytest.raises(dl.DetachedLaunchError, match="mangle"):
        dl.render_bat("demo", ["cmd", "/c", "exit", "0"], env=bad)


def test_the_generated_bat_has_no_caret_continuations():
    """One long command line instead: a trailing space after `^` silently truncates the command
    with no diagnostic, on a generated file nobody proof-reads."""
    for line in _bat().splitlines():
        assert not line.rstrip().endswith("^"), line


def test_the_generated_bat_is_ascii_and_crlf(tmp_path):
    out = dl.write_bat("smoke", ["cmd", "/c", "exit", "0"], repo=tmp_path)
    raw = out.read_bytes()
    raw.decode("ascii")
    assert b"\r\n" in raw
    assert b"\n" not in raw.replace(b"\r\n", b"")


def test_write_bat_creates_the_state_directory(tmp_path):
    dl.write_bat("smoke", ["cmd", "/c", "exit", "0"], repo=tmp_path)
    assert (tmp_path / "data" / "unattended").is_dir()


def test_an_empty_command_is_refused():
    """A launcher with no payload writes EXIT=0 having done nothing, which reads as success."""
    with pytest.raises(dl.DetachedLaunchError):
        dl.render_bat("demo", [])


@pytest.mark.parametrize("bad", ["a b", "a/b", "a\\b", "a:b", 'a"b', "a*b", ""])
def test_an_unsafe_job_name_is_refused(bad):
    with pytest.raises(dl.DetachedLaunchError):
        dl.render_bat(bad, ["cmd", "/c", "exit", "0"])


# --------------------------------------------------------------------- status states


def _write_hb(d: Path, job: str, *, pid: int, age: float = 0.0, done: int = 1,
              total: int = 10, cadence: int = 60) -> None:
    d.mkdir(parents=True, exist_ok=True)
    now = time.time() - age
    (d / ("%s.heartbeat.json" % job)).write_text(json.dumps({
        "job": job, "pid": pid, "units_done": done, "units_total": total,
        "wall_seconds": 100.0, "estimated_finish": "2026-11-04T00:00:00",
        "written_epoch": now, "cadence_seconds": cadence,
    }), encoding="utf-8")


def _write_sentinel(d: Path, job: str, rc: int) -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / ("%s.sentinel" % job)).write_text("EXIT=%d\r\n" % rc, encoding="utf-8")


def _a_definitely_dead_pid() -> int:
    """Spawn a trivial process, wait for it, return its pid.

    A made-up large integer is not reliable -- it can collide with a live process. Reaping a
    real child is the only way to hold a pid that is certainly gone.
    """
    p = subprocess.Popen([sys.executable, "-c", "pass"],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p.wait()
    return p.pid


def test_state_NOT_STARTED_when_nothing_exists(tmp_path):
    st = dl.classify("ghost", state_dir=tmp_path)
    assert st.state == dl.STATE_NOT_STARTED
    assert st.exit_code is None


def test_state_COMPLETE_on_exit_zero(tmp_path):
    _write_sentinel(tmp_path, "j", 0)
    st = dl.classify("j", state_dir=tmp_path)
    assert st.state == dl.STATE_COMPLETE
    assert not st.needs_attention


def test_state_FAILED_on_nonzero_exit(tmp_path):
    _write_sentinel(tmp_path, "j", 1)
    st = dl.classify("j", state_dir=tmp_path)
    assert st.state == dl.STATE_FAILED
    assert st.needs_attention


def test_exit_nine_is_REFUSED_not_FAILED(tmp_path):
    """HELD_EXIT exists precisely so a reader can tell 'refused to start' from 'ran and failed'.
    Folding it into FAILED destroys the only distinction the constant provides."""
    _write_sentinel(tmp_path, "j", HELD_EXIT)
    st = dl.classify("j", state_dir=tmp_path)
    assert st.state == dl.STATE_REFUSED
    assert st.state != dl.STATE_FAILED
    assert "NEVER RAN" in st.reason


def test_state_RUNNING_on_live_pid_and_fresh_heartbeat(tmp_path):
    _write_hb(tmp_path, "j", pid=os.getpid(), age=1.0)
    st = dl.classify("j", state_dir=tmp_path)
    assert st.state == dl.STATE_RUNNING
    assert st.pid_is_alive is True
    assert not st.needs_attention


def test_state_STALLED_on_live_pid_and_stale_heartbeat(tmp_path):
    _write_hb(tmp_path, "j", pid=os.getpid(), age=60 * dl.STALE_MULTIPLE + 120, cadence=60)
    st = dl.classify("j", state_dir=tmp_path)
    assert st.state == dl.STATE_STALLED
    assert st.pid_is_alive is True
    assert st.needs_attention


def test_state_DEAD_NO_SENTINEL_on_dead_pid_with_a_FRESH_heartbeat(tmp_path):
    """The ordering test. A fresh heartbeat does NOT make a dead job look alive -- pid existence
    is the liveness authority, exactly as in glm_run_lock."""
    _write_hb(tmp_path, "j", pid=_a_definitely_dead_pid(), age=0.5)
    st = dl.classify("j", state_dir=tmp_path)
    assert st.state == dl.STATE_DEAD
    assert st.pid_is_alive is False
    assert st.needs_attention


def test_state_DEAD_NO_SENTINEL_on_dead_pid_with_a_stale_heartbeat(tmp_path):
    _write_hb(tmp_path, "j", pid=_a_definitely_dead_pid(), age=99999)
    assert dl.classify("j", state_dir=tmp_path).state == dl.STATE_DEAD


def test_a_sentinel_OUTRANKS_a_live_fresh_heartbeat(tmp_path):
    """Terminal evidence first. A lingering heartbeat beside a finished job must not report
    RUNNING over a COMPLETE."""
    _write_hb(tmp_path, "j", pid=os.getpid(), age=0.1)
    _write_sentinel(tmp_path, "j", 0)
    assert dl.classify("j", state_dir=tmp_path).state == dl.STATE_COMPLETE


def test_an_EMPTY_sentinel_does_not_outrank_a_live_heartbeat(tmp_path):
    """The trap's artifact carries no evidence, so classification must fall through to liveness
    rather than treating the file's existence as terminal."""
    _write_hb(tmp_path, "j", pid=os.getpid(), age=0.1)
    (tmp_path / "j.sentinel").write_bytes(b"")
    assert dl.classify("j", state_dir=tmp_path).state == dl.STATE_RUNNING


def test_a_malformed_heartbeat_reads_as_absent(tmp_path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / "j.heartbeat.json").write_text('{"job": "j", "pi', encoding="utf-8")
    assert dl.classify("j", state_dir=tmp_path).state == dl.STATE_NOT_STARTED


def test_every_classification_returns_a_declared_state(tmp_path):
    """A state string that is not in ALL_STATES would render but mean nothing to a reader."""
    _write_hb(tmp_path, "a", pid=os.getpid())
    _write_sentinel(tmp_path, "b", 2)
    _write_sentinel(tmp_path, "c", 0)
    _write_hb(tmp_path, "d", pid=_a_definitely_dead_pid())
    for st in dl.status_all(tmp_path):
        assert st.state in dl.ALL_STATES


def test_status_all_finds_jobs_known_only_by_sentinel(tmp_path):
    _write_sentinel(tmp_path, "sentinel-only", 0)
    _write_hb(tmp_path, "heartbeat-only", pid=os.getpid())
    assert dl.known_jobs(tmp_path) == ["heartbeat-only", "sentinel-only"]


def test_known_jobs_is_empty_on_a_missing_state_dir(tmp_path):
    assert dl.known_jobs(tmp_path / "nope") == []


def test_render_status_marks_only_the_states_needing_attention(tmp_path):
    _write_sentinel(tmp_path, "good", 0)
    _write_sentinel(tmp_path, "bad", 1)
    text = dl.render_status(dl.status_all(tmp_path))
    bad = [ln for ln in text.splitlines() if ln.startswith("!!")]
    assert len(bad) == 1 and "bad" in bad[0]


# --------------------------------------------------------------------- heartbeat


def test_heartbeat_writes_and_is_readable(tmp_path):
    hb = dl.Heartbeat("j", units_total=4, state_dir=tmp_path, cadence_seconds=3600)
    assert hb.beat(1, force=True)
    rec = dl.read_heartbeat(hb.path)
    assert rec["units_done"] == 1 and rec["units_total"] == 4
    assert rec["pid"] == os.getpid()


def test_heartbeat_respects_its_cadence(tmp_path):
    hb = dl.Heartbeat("j", units_total=4, state_dir=tmp_path, cadence_seconds=3600)
    assert hb.beat(1, force=True)
    assert hb.beat(2) is False, "a sub-cadence beat must not write"
    assert dl.read_heartbeat(hb.path)["units_done"] == 1


def test_heartbeat_force_overrides_the_cadence(tmp_path):
    hb = dl.Heartbeat("j", units_total=4, state_dir=tmp_path, cadence_seconds=3600)
    hb.beat(1, force=True)
    assert hb.beat(4, force=True)
    assert dl.read_heartbeat(hb.path)["units_done"] == 4


def test_heartbeat_leaves_no_tmp_file_behind(tmp_path):
    """Atomic tmp+replace: a leftover .tmp beside the real file would be read by nothing but
    signals a write that never completed."""
    hb = dl.Heartbeat("j", units_total=2, state_dir=tmp_path)
    hb.beat(1, force=True)
    assert list(tmp_path.glob("*.tmp")) == []


def test_heartbeat_carries_extra_fields(tmp_path):
    hb = dl.Heartbeat("j", units_total=2, state_dir=tmp_path)
    hb.beat(1, force=True, seed=3, epoch=7, note="seed 3 epoch 7")
    rec = dl.read_heartbeat(hb.path)
    assert rec["seed"] == 3 and rec["epoch"] == 7


def test_a_zero_denominator_is_refused(tmp_path):
    with pytest.raises(dl.DetachedLaunchError):
        dl.Heartbeat("j", units_total=0, state_dir=tmp_path)


def test_no_finish_estimate_before_any_unit_is_done():
    """None, not a fabricated timestamp. The estimate is the field a returning reader is most
    likely to act on, so an invented one is the costliest possible guess."""
    assert dl.estimate_finish(0, 10, 100.0) is None
    assert dl.estimate_finish(1, 10, 0.0) is None


def test_the_finish_estimate_extrapolates_the_measured_rate():
    now = 1_000_000.0
    # 2 of 10 units in 100 s -> 50 s/unit -> 8 units left -> +400 s.
    got = dl.estimate_finish(2, 10, 100.0, now=now)
    from datetime import datetime
    assert got == datetime.fromtimestamp(now + 400).isoformat(timespec="seconds")


def test_the_finish_estimate_at_completion_is_now():
    from datetime import datetime
    now = 1_000_000.0
    assert dl.estimate_finish(10, 10, 100.0, now=now) == \
        datetime.fromtimestamp(now).isoformat(timespec="seconds")


# --------------------------------------------------------------------- cross-process lifetime


@pytest.mark.skipif(not os.environ.get("DNA_DECODE_RUN_DETACH_SMOKE"),
                    reason="spawns a real detached process; set DNA_DECODE_RUN_DETACH_SMOKE=1")
@WINDOWS_ONLY
def test_a_detached_job_outlives_the_launching_shell():
    """The only check that can establish survival past a closed shell, AND the whole chain.

    IT RUNS IN THE REAL REPO, NOT A TMP DIR, AND THAT IS THE POINT. The first version wrote the
    launcher into `tmp_path` and got EXIT=2: the generated `.bat` resolves
    `scripts/glm_run_lock.py` and `uv run` against its OWN directory, and a tmp dir has neither
    the script nor a pyproject. The sentinel appeared -- detachment worked -- while the payload
    had never run, which is exactly the `REFUSED`-vs-`FAILED` confusion this module exists to
    prevent, arriving from the test harness instead of the job.

    So the launcher is written where launchers actually live. That makes this a check of the
    real chain end to end -- uv resolves, the lock wrapper acquires and releases, the payload
    runs, the reversed sentinel records a real exit code -- rather than of process spawning
    alone. Opt-in because it writes to the repo root and waits on wall-clock; cleaned up in a
    `finally` so a failure cannot leave a stray launcher behind.
    """
    job = "detachsmoke"
    bat = dl.write_bat(job, ["cmd", "/c", "ping", "-n", "6", "127.0.0.1"], repo=dl.REPO)
    state = dl.DEFAULT_STATE_DIR
    try:
        subprocess.run(
            ["powershell", "-NoProfile", "-Command",
             "Start-Process -FilePath '%s' -WorkingDirectory '%s' -WindowStyle Hidden"
             % (bat, dl.REPO)],
            capture_output=True, timeout=60, check=True)

        # The launching shell is already gone by the time Start-Process returns.
        deadline = time.time() + 120
        saw_running = False
        while time.time() < deadline:
            st = dl.classify(job, state_dir=state)
            if st.state == dl.STATE_RUNNING:
                saw_running = True
            if dl.read_sentinel(dl.sentinel_path(job, state)) is not None:
                break
            time.sleep(1.0)

        rc = dl.read_sentinel(dl.sentinel_path(job, state))
        assert rc is not None, (
            "no sentinel after 120 s. Either the job never started or it died without a "
            "trace -- read %s" % dl.log_path(job, state))
        assert rc == 0, (
            "the chain ran but exited %d. 9 means the lock refused (another job holds it); "
            "2 usually means uv or the lock script did not resolve. Log: %s"
            % (rc, dl.log_path(job, state)))
        assert dl.classify(job, state_dir=state).state == dl.STATE_COMPLETE
        # `saw_running` is reported, NOT asserted, and the reason is structural rather than
        # timing: RUNNING requires a heartbeat, and this payload is a bare `ping` that writes
        # none, so mid-flight it classifies NOT_STARTED. Measured False on this host. The
        # RUNNING/STALLED branches are therefore covered only by the synthetic-heartbeat tests
        # above -- a real instrumented payload (the Step 5 runner) is what exercises them in
        # anger. Asserting it here would pin a state this fixture cannot produce.
        print("observed RUNNING mid-flight: %s (expected False: the payload writes no "
              "heartbeat)" % saw_running)
    finally:
        bat.unlink(missing_ok=True)
        for p in (dl.sentinel_path(job, state), dl.heartbeat_path(job, state),
                  dl.log_path(job, state), dl.lock_path(job, state)):
            p.unlink(missing_ok=True)
