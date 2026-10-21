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


# ================================================================== the CLI `run` path
# MEASURED 2026-10-21: the launcher's own documented form could not start a job AT ALL.
#
#   $ python scripts/glm_run_lock.py run X.lock --label "glm anchor x" -- python -c "print(1)"
#   glm_run_lock.py: error: unrecognized arguments: -- python -c print(1)   [exit 2]
#
# argparse does not strip the `--` separator when two leading positionals precede a
# `nargs="*"` positional. So the guarded command never ran and the sentinel recorded `EXIT=2`,
# which reads as "the job failed" rather than "the job was never launched".
#
# Every test above is green against that defect, and structurally had to be: they call
# `run_guarded` in-process, or use an action that takes no command, so none of them crossed the
# argparse layer the launcher actually goes through. That is the THIRD instance in this module
# of a guard proven on a path its only caller does not take -- after the dead temperature knob
# and after `acquire`-vs-`run`. Hence these drive `main`, not `run_guarded`.
def _argv_reporter(tmp_path):
    """A child that records BOTH its own argv and the label the lock recorded for it.

    Written as a file rather than `-c`, so the test is not asserting anything about shell
    quoting on this host.
    """
    p = tmp_path / "reporter.py"
    p.write_text(
        "import json, sys\n"
        "lock, out = sys.argv[1], sys.argv[2]\n"
        "rec = json.load(open(lock))\n"
        "json.dump({'argv': sys.argv[3:], 'lock_label': rec.get('label')}, open(out, 'w'))\n",
        encoding="utf-8")
    return str(p)


def test_cli_run_accepts_THE_LAUNCHERS_OWN_FORM_and_the_child_actually_runs(tmp_path):
    """The regression test. Pre-fix this returned 2 with a usage error and ran nothing."""
    lock = str(tmp_path / "x.lock")
    rc = main(["run", lock, "--label", "glm anchor out.json",
               "--", sys.executable, "-c", "raise SystemExit(4)"])
    assert rc == 4, "the child must run and its exit code propagate through the CLI"


def test_the_child_command_is_NOT_swallowed_by_the_LOCKS_OWN_FLAGS(tmp_path):
    """Non-vacuity, and the sharper half: everything after `--` belongs to the CHILD.

    The guarded command carries `--model auto --bins 6` in the real launcher, and a child
    command may legitimately repeat a flag name the lock itself defines. If `--label` after the
    separator reconfigured the LOCK instead of reaching the child, the split would be a
    coincidence rather than a contract -- and a test that only checked the exit code could not
    tell the difference.
    """
    lock, out = str(tmp_path / "x.lock"), str(tmp_path / "seen.json")
    tail = ["--model", "auto", "--bins", "6", "--label", "CHILDS-OWN"]
    rc = main(["run", lock, "--label", "LAUNCHERS-OWN",
               "--", sys.executable, _argv_reporter(tmp_path), lock, out, *tail])
    assert rc == 0
    seen = json.loads(pathlib.Path(out).read_text(encoding="utf-8"))
    assert seen["argv"] == tail, "the post-separator command must reach the child intact"
    assert seen["lock_label"] == "LAUNCHERS-OWN", \
        "a `--label` after `--` must NOT reconfigure the lock"


def test_cli_run_REFUSES_with_9_not_1_while_a_live_holder_exists(tmp_path):
    """The launcher documents 9 as refused-to-start; `acquire` returns 1. Pin the difference.

    Collapsing them would make "another job is running" indistinguishable from "acquire failed"
    at the one layer the .bat can see: %ERRORLEVEL%.
    """
    lock = str(tmp_path / "x.lock")
    acquire(pathlib.Path(lock), "incumbent")        # this test process stays alive throughout
    assert main(["run", lock, "--label", "intruder",
                 "--", sys.executable, "-c", "print('SHOULD NOT RUN')"]) == HELD_EXIT == 9
    assert main(["acquire", lock]) == 1


def test_cli_run_with_NO_command_is_REFUSED_not_silently_SUCCESSFUL(tmp_path):
    """Exiting 0 having run nothing is the failure mode that looks like a completed job."""
    with pytest.raises(GlmLockError, match="requires a command"):
        main(["run", str(tmp_path / "x.lock")])


def _exit6(tmp_path):
    p = tmp_path / "exit6.py"
    p.write_text("raise SystemExit(6)\n", encoding="utf-8")
    return str(p)


def test_the_SEPARATORLESS_form_still_works_for_an_OPTIONLESS_command(tmp_path):
    """Non-vacuity for the split itself: it must not have broken the simple form."""
    lock = str(tmp_path / "x.lock")
    assert main(["run", lock, sys.executable, _exit6(tmp_path)]) == 6


def test_but_WITHOUT_the_separator_an_OPTION_LOOKING_argument_is_still_claimed_by_argparse(
        tmp_path):
    """The limitation, pinned rather than hidden -- and the reason `--` is MANDATORY.

    MEASURED while writing these tests: the first version of the test above used
    `python -c "raise SystemExit(6)"` with no separator and failed with
    `unrecognized arguments: -c raise SystemExit(6)`. A `nargs="*"` positional cannot hold a
    token that looks like a flag, and the real guarded command is full of them
    (`--model auto --bins 6 --compare-to ...`).

    So the separator-less form is a convenience for trivial commands ONLY, and a launcher that
    drops `--` fails LOUDLY (exit 2) rather than silently mis-parsing. Pinning that keeps a
    future author from "simplifying" the .bat by removing the separator.
    """
    lock = str(tmp_path / "x.lock")
    with pytest.raises(SystemExit) as e:
        main(["run", lock, sys.executable, "-c", "raise SystemExit(6)"])
    assert e.value.code == 2, "a dropped separator must fail loudly, not mis-parse"


def test_the_other_actions_are_UNAFFECTED_by_the_split(tmp_path):
    """`acquire`/`status`/`release` take no command; the partition must be a no-op for them."""
    lock = str(tmp_path / "x.lock")
    assert main(["acquire", lock, "--label", "a"]) == 0
    assert main(["status", lock]) == 1
    assert main(["release", lock, "--pid", str(os.getpid())]) == 0
    assert main(["status", lock]) == 0


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


def test_the_launcher_passes_a_BARE_SEPARATOR_before_the_guarded_command():
    """Closes the loop with the CLI tests above.

    Those prove the CLI ACCEPTS `... -- <command>`; this proves the launcher USES that form.
    Either alone is insufficient: the CLI could accept a form nothing sends, and the launcher
    could send a form the CLI rejects -- which is exactly what was true until 2026-10-21.
    """
    code = _bat_code()
    i, j = code.find("glm_run_lock.py"), code.find("glm_generator_falsifier.py")
    between = code[i:j]
    assert any(tok == "--" for tok in between.split()), \
        "the launcher must separate the guarded command with a bare `--`"


def test_the_sentinel_redirect_PRECEDES_the_echo():
    """The measured batch trap.

    `echo EXIT=%RC%>> file` expands to `echo EXIT=0>> file`, and cmd reads a digit immediately
    before `>>` as a FILE HANDLE number -- so it redirected stdin and left the sentinel EMPTY
    while `EXIT=` went to the console. The first detached run's completion was therefore
    undetectable and had to be spotted in `git status`.
    """
    assert '>>"%SENTINEL%" echo EXIT=' in _bat(), "redirect must come before the echo"
    assert "echo EXIT=%RC%>>" not in _bat_code(), "the trailing-redirect form is the trap"


def test_the_launcher_pins_the_INTERPRETER_rather_than_trusting_PATH():
    """The `--`-separator defect is interpreter-dependent, so the launcher must not gamble.

    MEASURED 2026-10-21: the pre-fix code fails under the venv (3.11.17, what `uv run` and the
    test suite use) with `unrecognized arguments: -- ...` and exit 2, and SUCCEEDS under the
    3.14.7 on PATH. Five verification attempts under bare `python` all passed and nearly had the
    defect recorded as unreproducible.

    The fix makes either interpreter work, so this is not what prevents the bug -- it removes
    the launcher's dependence on whichever `python` is first on PATH, which is what made it work
    on this host by luck rather than by design.
    """
    code = _bat_code()
    assert "uv run python" in code and "glm_run_lock.py" in code
    i = code.find("glm_run_lock.py")
    prefix = code[max(0, i - 80):i]
    assert "uv run python" in prefix, (
        "the lock script must be invoked through `uv run python`, not bare `python`")
