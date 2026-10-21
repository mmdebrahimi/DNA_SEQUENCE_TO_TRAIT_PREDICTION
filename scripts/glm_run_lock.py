"""A PID-bearing lock so two GLM generation jobs cannot run at once.

WHY. Two concurrent torch processes is this host's documented access-violation hazard, and the
detached runner had nothing preventing a second launch. During the run that motivated this, a
job's liveness was misjudged once and relaunched; the sentinel file could not help, because it
is written only at EXIT and therefore says nothing at all while a job is running.

So liveness is checked by PID EXISTENCE, not by the sentinel. The two facts a launcher needs are
"is someone running" and "did the last one finish", and they are different questions with
different evidence.

A STALE lock -- one whose PID is gone -- is RECLAIMED rather than honoured. That direction is
deliberate: a killed job must not block the next launch forever, which is exactly the failure
that cost this execution two restarts. The opposite choice (refuse until a human clears it)
trades a crash hazard for a wedge, and the crash hazard is already covered by the PID check.

USE `run`, NOT `acquire` + `release`, FROM A LAUNCHER. This was found by the manual two-launch
check, and the unit tests structurally could not see it: `acquire` records the pid of the
PROCESS THAT CALLS IT, and a launcher's `python glm_run_lock.py acquire` exits immediately. Its
pid is therefore dead by the time a second launch looks, the lock is correctly reclaimed as
stale, and the second job starts. The in-process tests passed because there the caller's pid
really does stay alive for the duration.

`run` holds the lock for the CHILD's lifetime -- it acquires, spawns, waits, releases in a
`finally`, and propagates the child's exit code -- so the recorded pid lives exactly as long as
the work it guards. `acquire`/`release` remain for use INSIDE one long-lived process.

Usage (from a .bat or a shell):
    python scripts/glm_run_lock.py run <lockfile> [--label TEXT] -- <command...>
                                                                       -> child's exit code, or 9 if held
    python scripts/glm_run_lock.py acquire <lockfile> [--label TEXT]   -> 0 ok, 1 held
    python scripts/glm_run_lock.py release <lockfile>                  -> 0 always
    python scripts/glm_run_lock.py status  <lockfile>                  -> 0 free, 1 held
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path


def pid_alive(pid: int) -> bool:
    """True iff a process with this pid currently exists.

    Windows has no `os.kill(pid, 0)` semantics worth relying on, so this uses the one check
    available on both platforms without spawning anything: on Windows, OpenProcess via ctypes;
    elsewhere, signal 0. An UNDECIDABLE answer returns True -- i.e. assume alive -- because the
    fail-safe direction here is to refuse a second launch, not to permit one.
    """
    if pid <= 0:
        return False
    if os.name == "nt":
        try:
            import ctypes
            from ctypes import wintypes

            PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
            k32 = ctypes.WinDLL("kernel32", use_last_error=True)
            k32.OpenProcess.restype = wintypes.HANDLE
            h = k32.OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
            if not h:
                return False
            try:
                code = wintypes.DWORD()
                if k32.GetExitCodeProcess(h, ctypes.byref(code)):
                    STILL_ACTIVE = 259
                    return code.value == STILL_ACTIVE
                return True          # undecidable -> assume alive (fail safe)
            finally:
                k32.CloseHandle(h)
        except Exception:
            return True              # undecidable -> assume alive (fail safe)
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True                  # exists, owned by someone else
    except Exception:
        return True


def read_lock(path: Path) -> dict | None:
    """The lock's contents, or None when absent/unreadable/malformed.

    A malformed lock is treated as ABSENT rather than as held: a half-written file is not
    evidence that a job is running, and treating it as held would wedge the launcher on a
    crash mid-write.
    """
    try:
        if not path.is_file():
            return None
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def holder(path: Path) -> dict | None:
    """The LIVE holder of this lock, or None when free (absent, malformed, or stale)."""
    rec = read_lock(path)
    if not rec:
        return None
    try:
        pid = int(rec.get("pid", -1))
    except (TypeError, ValueError):
        return None
    if not pid_alive(pid):
        return None                  # stale -> reclaimable
    return rec


def acquire(path: Path, label: str = "") -> tuple[bool, dict]:
    """Take the lock unless a live holder has it. Returns (ok, record)."""
    live = holder(path)
    if live is not None:
        return False, live
    path.parent.mkdir(parents=True, exist_ok=True)
    rec = {"pid": os.getpid(), "label": label,
           "acquired": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "argv": " ".join(sys.argv[1:])}
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(rec), encoding="utf-8")
    os.replace(tmp, path)
    return True, rec


def release(path: Path, pid: int | None = None) -> bool:
    """Drop the lock. A release by a NON-owner is a no-op.

    Returns True when a lock was actually removed. Idempotent: releasing a free lock is fine,
    which matters because the launcher releases on every exit path including the ones where it
    never acquired.
    """
    rec = read_lock(path)
    if not rec:
        return False
    if pid is not None and int(rec.get("pid", -1)) != int(pid):
        return False                 # not ours; leave the owner's lock alone
    try:
        path.unlink()
        return True
    except Exception:
        return False


HELD_EXIT = 9
"""Exit code when the lock is held. Distinct from any plausible child exit code, so a launcher
can tell "refused to start" from "ran and failed"."""


def run_guarded(path: Path, command: list[str], label: str = "") -> int:
    """Hold the lock for the CHILD's lifetime, then propagate its exit code.

    This is the correct shape for a launcher. `acquire` + separate `release` is not: a
    launcher's `acquire` process exits at once, so the pid it recorded is dead and the next
    launch reclaims the lock as stale. Measured by the two-launch manual check, which the
    in-process unit tests could not reproduce.

    Releases in a `finally`, so a crashing or killed child does not leave a lock behind -- and
    even if it does, the pid is gone and the next `run` reclaims it.
    """
    if not command:
        raise GlmLockError("run requires a command after `--`")
    ok, live = acquire(path, label)
    if not ok:
        print(f"REFUSING: lock held by live pid={live.get('pid')} since {live.get('acquired')} "
              f"({live.get('label') or 'no label'}). Two concurrent torch jobs is this host's "
              f"access-violation hazard.", file=sys.stderr)
        return HELD_EXIT
    try:
        import subprocess
        return subprocess.call(command)
    finally:
        release(path, os.getpid())


class GlmLockError(ValueError):
    """A lock operation that cannot mean what the caller intends."""


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["run", "acquire", "release", "status"])
    p.add_argument("lockfile")
    p.add_argument("--label", default="")
    p.add_argument("--pid", type=int, default=None,
                   help="release only if the lock is held by this pid")
    p.add_argument("command", nargs="*",
                   help="for `run`: the command to guard (everything after `--`)")
    a = p.parse_args(argv)
    path = Path(a.lockfile)

    if a.action == "run":
        return run_guarded(path, list(a.command), a.label)

    if a.action == "acquire":
        ok, rec = acquire(path, a.label)
        if ok:
            print(f"lock acquired pid={rec['pid']} {a.lockfile}")
            return 0
        print(f"REFUSING: lock held by live pid={rec.get('pid')} "
              f"since {rec.get('acquired')} ({rec.get('label') or 'no label'}). "
              f"Two concurrent torch jobs is this host's access-violation hazard.",
              file=sys.stderr)
        return 1

    if a.action == "release":
        print("released" if release(path, a.pid) else "nothing to release")
        return 0

    live = holder(path)
    if live is None:
        print("free")
        return 0
    print(f"held pid={live.get('pid')} since {live.get('acquired')}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
