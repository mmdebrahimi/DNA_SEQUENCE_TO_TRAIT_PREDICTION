"""The single-job lock for GLM generation runs.

WHY: two concurrent torch processes is this host's documented access-violation hazard, and the
detached runner had nothing preventing a second launch. During the run that motivated this, a
job's liveness was misjudged and relaunched.

The sentinel file cannot answer "is someone running" -- it is written only at EXIT -- so
liveness is by PID existence. A stale lock (PID gone) is RECLAIMED rather than honoured,
because a killed job must not wedge the next launch; that exact wedge cost two restarts.
"""
from __future__ import annotations

import json
import os
import pathlib
import sys

import pytest

from scripts.glm_run_lock import (
    HELD_EXIT,
    GlmLockError,
    acquire,
    holder,
    main,
    pid_alive,
    read_lock,
    release,
    run_guarded,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent


# ------------------------------------------------------------------ pid liveness
def test_our_own_pid_is_alive():
    assert pid_alive(os.getpid()) is True


def test_a_nonsense_pid_is_not_alive():
    assert pid_alive(-1) is False
    assert pid_alive(0) is False


def test_an_almost_certainly_dead_pid_is_not_alive():
    """A high pid no process holds. Skipped rather than guessed if it happens to exist."""
    dead = 999_999
    if pid_alive(dead):
        pytest.skip("pid 999999 happens to exist on this host")
    assert pid_alive(dead) is False


# ------------------------------------------------------------------ acquire / release
def test_acquire_on_a_free_lock_succeeds(tmp_path):
    lock = tmp_path / "x.lock"
    ok, rec = acquire(lock, "first")
    assert ok and rec["pid"] == os.getpid()
    assert lock.is_file()


def test_a_SECOND_acquire_against_a_LIVE_holder_is_REFUSED(tmp_path):
    lock = tmp_path / "x.lock"
    assert acquire(lock, "first")[0] is True
    ok, live = acquire(lock, "second")
    assert ok is False
    assert live["pid"] == os.getpid()


def test_a_STALE_lock_is_RECLAIMED_not_honoured(tmp_path):
    """A killed job must not wedge the next launch -- that wedge cost two restarts."""
    lock = tmp_path / "x.lock"
    dead = 999_999
    if pid_alive(dead):
        pytest.skip("pid 999999 happens to exist on this host")
    lock.write_text(json.dumps({"pid": dead, "label": "killed job"}), encoding="utf-8")
    assert holder(lock) is None
    ok, rec = acquire(lock, "next")
    assert ok is True and rec["pid"] == os.getpid()


def test_a_MALFORMED_lock_is_treated_as_absent_not_as_held(tmp_path):
    """A half-written file is not evidence a job is running; treating it as held wedges."""
    lock = tmp_path / "x.lock"
    lock.write_text("{not json", encoding="utf-8")
    assert read_lock(lock) is None
    assert holder(lock) is None
    assert acquire(lock, "after-crash")[0] is True


def test_a_lock_with_a_NON_INTEGER_pid_is_treated_as_absent(tmp_path):
    lock = tmp_path / "x.lock"
    lock.write_text(json.dumps({"pid": "banana"}), encoding="utf-8")
    assert holder(lock) is None


def test_release_is_idempotent(tmp_path):
    lock = tmp_path / "x.lock"
    acquire(lock, "l")
    assert release(lock) is True
    assert release(lock) is False
    assert not lock.exists()


def test_release_by_a_NON_OWNER_leaves_the_owners_lock_alone(tmp_path):
    lock = tmp_path / "x.lock"
    acquire(lock, "owner")
    assert release(lock, pid=os.getpid() + 1) is False
    assert lock.is_file(), "a non-owner's release must not free someone else's lock"
    assert release(lock, pid=os.getpid()) is True


# ------------------------------------------------------------------ the `run` wrapper
# The shape a launcher must use, for a reason the in-process tests above structurally cannot
# see: `acquire` records the pid of the PROCESS THAT CALLS IT, so a launcher's short-lived
# `acquire` process is already dead when the next launch looks and the lock is reclaimed as
# stale. Measured on the real surface -- with acquire/release the second launch exited 0 and
# its job RAN; with `run` it exits 9 and the job never starts.
def test_run_guarded_propagates_the_childs_exit_code(tmp_path):
    lock = tmp_path / "x.lock"
    assert run_guarded(lock, [sys.executable, "-c", "raise SystemExit(0)"], "ok") == 0
    assert run_guarded(lock, [sys.executable, "-c", "raise SystemExit(3)"], "fail") == 3


def test_run_guarded_RELEASES_even_when_the_child_fails(tmp_path):
    """Released in a finally, so a crashing child does not leave a lock behind."""
    lock = tmp_path / "x.lock"
    run_guarded(lock, [sys.executable, "-c", "raise SystemExit(7)"], "boom")
    assert not lock.exists()


def test_run_guarded_REFUSES_while_a_live_holder_exists(tmp_path):
    """The lock is held by a pid that lives as long as the work it guards."""
    lock = tmp_path / "x.lock"
    acquire(lock, "incumbent")              # this test process stays alive throughout
    rc = run_guarded(lock, [sys.executable, "-c", "print('SHOULD NOT RUN')"], "intruder")
    assert rc == HELD_EXIT == 9


def test_the_held_exit_code_is_distinguishable_from_a_child_failure(tmp_path):
    """9 means refused-to-start; a launcher must tell that from "ran and failed"."""
    assert HELD_EXIT == 9
    lock = tmp_path / "x.lock"
    assert run_guarded(lock, [sys.executable, "-c", "raise SystemExit(1)"], "x") == 1


def test_run_with_no_command_is_REFUSED(tmp_path):
    with pytest.raises(GlmLockError, match="requires a command"):
        run_guarded(tmp_path / "x.lock", [], "empty")


# ------------------------------------------------------------------ CLI contract
def test_cli_acquire_then_acquire_returns_0_then_1(tmp_path):
    lock = str(tmp_path / "x.lock")
    assert main(["acquire", lock, "--label", "a"]) == 0
    assert main(["acquire", lock, "--label", "b"]) == 1


def test_cli_status_reports_free_vs_held(tmp_path):
    lock = str(tmp_path / "x.lock")
    assert main(["status", lock]) == 0
    main(["acquire", lock])
    assert main(["status", lock]) == 1
    main(["release", lock])
    assert main(["status", lock]) == 0


def test_cli_release_always_exits_0(tmp_path):
    assert main(["release", str(tmp_path / "never.lock")]) == 0


# ------------------------------------------------------------------ the .bat wiring
# These read the launcher as TEXT. They cannot prove cmd.exe behaviour -- that is the manual
# check -- but they pin the three things that were actually wrong.
def _bat() -> str:
    p = ROOT / "run_glm_anchor_detached.bat"
    if not p.is_file():
        pytest.skip("detached launcher absent")
    return p.read_text(encoding="utf-8", errors="replace")


def _bat_code() -> str:
    """Executable lines only.

    The launcher documents its traps verbatim in REM blocks, and the first version of these
    guards matched those comments -- a test reading documentation as if it were the
    implementation, which is the same mistake in a different costume.
    """
    keep = [l for l in _bat().splitlines() if not l.strip().upper().startswith("REM")]
    return chr(10).join(keep)


def test_the_launcher_uses_RUN_not_acquire_plus_release():
    code = _bat_code()
    assert "glm_run_lock.py" in code
    assert " run " in code, "the launcher must use the `run` wrapper"
    assert "acquire" not in code, "a launcher must not use bare acquire -- its pid dies at once"
    assert "release" not in code, "release is the wrapper's job, in a finally"


def test_the_launcher_guards_the_REAL_job_not_something_else():
    """Non-vacuity: the wrapper must sit in front of the generation command."""
    code = _bat_code()
    i, j = code.find("glm_run_lock.py"), code.find("glm_generator_falsifier.py")
    assert i != -1 and j != -1 and i < j, "the lock wrapper must precede the guarded command"


def test_the_sentinel_redirect_PRECEDES_the_echo():
    """The measured batch trap.

    `echo EXIT=%RC%>> file` expands to `echo EXIT=0>> file`, and cmd reads a digit immediately
    before `>>` as a FILE HANDLE number -- so it redirected stdin and left the sentinel EMPTY
    while `EXIT=` went to the console. The first detached run's completion was therefore
    undetectable and had to be spotted in `git status`.
    """
    assert '>>"%SENTINEL%" echo EXIT=' in _bat(), "redirect must come before the echo"
    assert "echo EXIT=%RC%>>" not in _bat_code(), "the trailing-redirect form is the trap"
