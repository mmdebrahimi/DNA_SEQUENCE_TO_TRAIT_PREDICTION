"""Generate, launch and READ BACK long detached jobs on this host.

WHY THIS EXISTS. Two unattended windows are about to be spent -- a ~4.6 h expression-heads
scoring run and a multi-day fine-tune -- and the binding constraint on both is ATTENDANCE, not
compute. Nothing will invoke me while the user is away, so the only thing that makes a 10-day
run worth starting is that its state is READABLE ON RETURN without re-deriving it from a log
tail. That is this module's job: emit the launcher, carry a heartbeat, and classify the result.

THE SIX-DEATH HISTORY IS THE DESIGN INPUT. `run_tile_conv_detached.bat` records six long-run
deaths on this machine: four STATUS_IN_PAGE_ERROR (0xc0000006) with torch mmapping weights off
a USB volume, and TWO THAT LEFT NO EVENT ON ANY OF 118 LOGS and remain unexplained. The one
~12 h survivor kept its data on C:. A reader must therefore be able to tell an unexplained
death apart from a clean finish, because on this host that is not a hypothetical state -- it has
happened twice. `DEAD_NO_SENTINEL` is that state.

SEVEN STATES, NOT FOUR, AND EACH ADDITION CLOSES A CONFLATION THE OTHERS WOULD CREATE.
The plan named four (RUNNING / DEAD_NO_SENTINEL / FAILED / COMPLETE). Three more are not scope
creep; each one exists because folding it into a neighbour would destroy a distinction some
OTHER part of this system was deliberately built to preserve:

  * `REFUSED_LOCK_HELD` (sentinel EXIT=9) is NOT `FAILED`. `glm_run_lock.HELD_EXIT`'s own
    docstring says 9 is "distinct from any plausible child exit code, so a launcher can tell
    'refused to start' from 'ran and failed'." Calling it FAILED would throw away the only
    distinction that constant exists to provide, and would send a reader hunting for a bug in a
    job that never ran a single step.
  * `NOT_STARTED` (no sentinel AND no heartbeat) is NOT `DEAD_NO_SENTINEL`. Conflating
    "never launched" with "died silently" is the exact failure DEAD_NO_SENTINEL was added to
    prevent, one level up: both look like absence, and only one is alarming.
  * `STALLED` (pid alive, heartbeat stale) is NOT `RUNNING`. Over a 10-day window the
    decision-relevant question is not "is the process object alive" but "is it still making
    progress"; a job frozen since hour 3 answers yes to the first and no to the second.

PID EXISTENCE IS THE LIVENESS AUTHORITY, NOT THE HEARTBEAT TIMESTAMP. This mirrors
`glm_run_lock`'s hard-won rule, and the ordering matters: a FRESH heartbeat with a DEAD pid is
`DEAD_NO_SENTINEL`, not RUNNING. The heartbeat is written by the worker and simply stops when
the worker dies, so its last timestamp says when the job was last alive -- never that it is
alive now. A reader that trusted freshness would report a crashed job as healthy for exactly as
long as its cadence allows.

THE HEARTBEAT PID AND THE LOCK PID ARE DIFFERENT PROCESSES, DELIBERATELY. The launcher runs
`glm_run_lock.py run ... -- <payload>`, which records the WRAPPER's pid in the lock and spawns
the payload as a child. The heartbeat is written from inside the payload, so it carries the
WORKER's pid -- the process whose death actually ends the work. Reading liveness off the lock
instead would report the wrapper, which outlives nothing useful.

Usage:
    python scripts/detached_launch.py emit  <job> -- <command...>   # write run_<job>_detached.bat
    python scripts/detached_launch.py status [<job>]                # classify one job or all
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from glm_run_lock import HELD_EXIT, pid_alive  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
DEFAULT_STATE_DIR = REPO / "data" / "unattended"
"""Runtime state lives under the gitignored `data/` tree (`.gitignore:4:/data/*`) and on C:,
which is the one profile a long job on this host has ever SURVIVED with."""

DEFAULT_CADENCE_SECONDS = 60
STALE_MULTIPLE = 5
"""A heartbeat is stale once it is older than this many cadences. Five rather than two because
a CPU-bound epoch boundary on four threads can legitimately straddle a few cadences, and a
false STALLED on a healthy job is the costlier error here -- it invites a kill."""

# --- The sentinel-redirect trap, as two literal lines so a test can RUN both ---------------
#
# `echo EXIT=%RC%>>"%SENTINEL%"` expands to `echo EXIT=0>>"file"`, and cmd reads a digit
# IMMEDIATELY PRECEDING `>>` as a FILE HANDLE number -- it redirects handle 0 (stdin) and leaves
# the sentinel EMPTY. The exit code is then unrecoverable: an empty sentinel is indistinguishable
# from a job that died before writing one. Measured on run_glm_anchor_detached.bat, not theorised.
SENTINEL_LINE_CORRECT = '>>"%SENTINEL%" echo EXIT=%RC%'
SENTINEL_LINE_TRAP = 'echo EXIT=%RC%>>"%SENTINEL%"'

STATE_COMPLETE = "COMPLETE"
STATE_FAILED = "FAILED"
STATE_REFUSED = "REFUSED_LOCK_HELD"
STATE_RUNNING = "RUNNING"
STATE_STALLED = "STALLED"
STATE_DEAD = "DEAD_NO_SENTINEL"
STATE_NOT_STARTED = "NOT_STARTED"

ALL_STATES = (STATE_COMPLETE, STATE_FAILED, STATE_REFUSED, STATE_RUNNING,
              STATE_STALLED, STATE_DEAD, STATE_NOT_STARTED)

ATTENTION_STATES = (STATE_FAILED, STATE_DEAD, STATE_STALLED, STATE_REFUSED)


class DetachedLaunchError(ValueError):
    """A launch or status request that cannot mean what the caller intends."""


def heartbeat_path(job: str, state_dir: Path | None = None) -> Path:
    return (state_dir or DEFAULT_STATE_DIR) / ("%s.heartbeat.json" % job)


def sentinel_path(job: str, state_dir: Path | None = None) -> Path:
    return (state_dir or DEFAULT_STATE_DIR) / ("%s.sentinel" % job)


def log_path(job: str, state_dir: Path | None = None) -> Path:
    return (state_dir or DEFAULT_STATE_DIR) / ("%s.log" % job)


def lock_path(job: str, state_dir: Path | None = None) -> Path:
    return (state_dir or DEFAULT_STATE_DIR) / ("%s.lock" % job)


# --------------------------------------------------------------------------- heartbeat


def estimate_finish(units_done: int, units_total: int, wall_seconds: float,
                    now: float | None = None) -> str | None:
    """ISO finish estimate, or None when no estimate is DERIVABLE.

    Returns None at zero units done or zero elapsed -- there is no rate yet, and inventing one
    would put a confident timestamp in the one field a returning reader is most likely to act
    on. None is the honest answer and reads as such.
    """
    if units_done <= 0 or units_total <= 0 or wall_seconds <= 0:
        return None
    base = now if now is not None else time.time()
    if units_done >= units_total:
        return datetime.fromtimestamp(base).isoformat(timespec="seconds")
    rate = units_done / wall_seconds
    remaining = (units_total - units_done) / rate
    return datetime.fromtimestamp(base + remaining).isoformat(timespec="seconds")


class Heartbeat:
    """A cadence-limited, ATOMICALLY-written progress file for one detached job.

    Written from inside the PAYLOAD process, so `pid` is the worker's -- see the module note on
    why that differs from the lock's pid.

    Atomic tmp+replace, for the same reason `glm_run_lock.read_lock` treats a malformed lock as
    absent: a reader polling a 10-day job will eventually land mid-write, and a half-written JSON
    must not be able to turn a healthy RUNNING into an unreadable state.
    """

    def __init__(self, job: str, units_total: int, *, state_dir: Path | None = None,
                 cadence_seconds: int = DEFAULT_CADENCE_SECONDS, label: str = ""):
        if units_total <= 0:
            raise DetachedLaunchError("units_total must be positive; a progress fraction with "
                                      "no denominator cannot be read as progress")
        self.job = job
        self.units_total = int(units_total)
        self.cadence_seconds = cadence_seconds
        self.label = label
        self.path = heartbeat_path(job, state_dir)
        self.started = time.time()
        self._last_write = 0.0

    def beat(self, units_done: int, *, force: bool = False, **extra) -> bool:
        """Write a heartbeat if the cadence has elapsed. True when a write happened.

        `force=True` is for the first and last beats, which must land regardless of cadence --
        the first because a job with no heartbeat at all reads as NOT_STARTED, the last because
        the final unit is the one a reader most wants confirmed.
        """
        now = time.time()
        if not force and (now - self._last_write) < self.cadence_seconds:
            return False
        wall = now - self.started
        rec = {
            "job": self.job,
            "label": self.label,
            "pid": os.getpid(),
            "units_done": int(units_done),
            "units_total": self.units_total,
            "wall_seconds": round(wall, 1),
            "estimated_finish": estimate_finish(units_done, self.units_total, wall, now=now),
            "written": datetime.fromtimestamp(now).isoformat(timespec="seconds"),
            "written_epoch": now,
            "cadence_seconds": self.cadence_seconds,
        }
        rec.update(extra)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(self.path.suffix + ".tmp")
        tmp.write_text(json.dumps(rec, indent=2), encoding="utf-8")
        os.replace(tmp, self.path)
        self._last_write = now
        return True


def read_heartbeat(path: Path) -> dict | None:
    """The heartbeat record, or None when absent/unreadable/malformed.

    Malformed reads as ABSENT, matching `glm_run_lock.read_lock`: a half-written file is not
    evidence about the job, and the caller's next question (is the pid alive) does not need it.
    """
    try:
        if not path.is_file():
            return None
        rec = json.loads(path.read_text(encoding="utf-8"))
        return rec if isinstance(rec, dict) else None
    except Exception:
        return None


def read_sentinel(path: Path) -> int | None:
    """The recorded exit code, or None when the sentinel is absent OR EMPTY.

    An EMPTY sentinel is the redirect trap's signature and must read as None, not 0. Treating an
    empty file as success would convert the trap from a loud missing exit code into a silent
    false COMPLETE -- which is strictly worse than the bug it comes from.
    """
    try:
        if not path.is_file():
            return None
        text = path.read_text(encoding="utf-8", errors="replace").strip()
        if not text:
            return None
        for tok in text.replace("\r", "\n").split("\n"):
            tok = tok.strip()
            if tok.upper().startswith("EXIT="):
                try:
                    return int(tok.split("=", 1)[1].strip())
                except ValueError:
                    return None
        return None
    except Exception:
        return None


# --------------------------------------------------------------------------- status


@dataclass(frozen=True)
class JobStatus:
    job: str
    state: str
    reason: str
    exit_code: int | None = None
    pid: int | None = None
    pid_is_alive: bool | None = None
    units_done: int | None = None
    units_total: int | None = None
    wall_seconds: float | None = None
    heartbeat_age_seconds: float | None = None
    estimated_finish: str | None = None

    @property
    def needs_attention(self) -> bool:
        return self.state in ATTENTION_STATES

    def as_dict(self) -> dict:
        return asdict(self)


def classify(job: str, *, state_dir: Path | None = None, now: float | None = None) -> JobStatus:
    """Classify one job into exactly one of ALL_STATES.

    ORDER IS THE CONTRACT. The sentinel is consulted FIRST because it is the only terminal
    evidence: once it exists the job is over, and a stale lingering heartbeat must not be able to
    report RUNNING over a finished job. Only in the sentinel's absence does liveness matter, and
    there PID EXISTENCE outranks heartbeat freshness -- so a fresh heartbeat with a dead pid is
    DEAD_NO_SENTINEL.
    """
    now = time.time() if now is None else now
    hb = read_heartbeat(heartbeat_path(job, state_dir))
    rc = read_sentinel(sentinel_path(job, state_dir))

    pid = None
    if hb is not None:
        try:
            pid = int(hb.get("pid", -1))
        except (TypeError, ValueError):
            pid = None

    age = None
    if hb is not None and isinstance(hb.get("written_epoch"), (int, float)):
        age = round(now - float(hb["written_epoch"]), 1)

    common = dict(
        job=job, exit_code=rc, pid=pid,
        units_done=(hb or {}).get("units_done"),
        units_total=(hb or {}).get("units_total"),
        wall_seconds=(hb or {}).get("wall_seconds"),
        estimated_finish=(hb or {}).get("estimated_finish"),
        heartbeat_age_seconds=age,
    )

    # 1. Terminal evidence first.
    if rc is not None:
        if rc == 0:
            return JobStatus(state=STATE_COMPLETE, reason="sentinel recorded EXIT=0", **common)
        if rc == HELD_EXIT:
            return JobStatus(
                state=STATE_REFUSED,
                reason="sentinel recorded EXIT=%d: the run lock was held, so the payload NEVER "
                       "RAN. This is not a failure of the work." % HELD_EXIT, **common)
        return JobStatus(state=STATE_FAILED,
                         reason="sentinel recorded EXIT=%d" % rc, **common)

    # 2. No sentinel. Absence of a heartbeat too means it never got going.
    if hb is None:
        return JobStatus(state=STATE_NOT_STARTED,
                         reason="no sentinel and no heartbeat: this job has not started "
                                "(distinct from having died silently)", **common)

    # 3. Heartbeat present, no sentinel -> liveness decides, by PID.
    alive = pid_alive(pid) if pid else False
    if not alive:
        return JobStatus(
            state=STATE_DEAD, pid_is_alive=False,
            reason="pid %s is gone and no sentinel was written -- the unexplained-death "
                   "signature this host has twice. Last heartbeat was %ss ago."
                   % (pid, age if age is not None else "?"), **common)

    cadence = hb.get("cadence_seconds") or DEFAULT_CADENCE_SECONDS
    if age is not None and age > cadence * STALE_MULTIPLE:
        return JobStatus(
            state=STATE_STALLED, pid_is_alive=True,
            reason="pid %s is alive but its heartbeat is %ss old against a %ss cadence -- "
                   "alive is not the same as progressing" % (pid, age, cadence), **common)
    return JobStatus(state=STATE_RUNNING, pid_is_alive=True,
                     reason="pid %s alive, heartbeat %ss old" % (pid, age), **common)


def known_jobs(state_dir: Path | None = None) -> list[str]:
    d = state_dir or DEFAULT_STATE_DIR
    if not d.is_dir():
        return []
    suffix = ".heartbeat.json"
    jobs = {p.name[: -len(suffix)] for p in d.glob("*" + suffix)}
    jobs |= {p.stem for p in d.glob("*.sentinel")}
    return sorted(jobs)


def status_all(state_dir: Path | None = None, now: float | None = None) -> list[JobStatus]:
    return [classify(j, state_dir=state_dir, now=now) for j in known_jobs(state_dir)]


# --------------------------------------------------------------------------- bat emission


def _win(rel) -> str:
    return str(rel).replace("/", "\\")


def render_bat(job: str, command: list[str], *, rationale: str = "",
               state_rel: str = "data/unattended", env: dict | None = None) -> str:
    """The launcher text for one job, matching the existing convention EXACTLY.

    THE COMMAND IS EMITTED ON ONE LINE, not with `^` continuations as the hand-written launchers
    use. A trailing space after a caret silently ends the continuation and cmd runs a truncated
    command -- a failure mode with no diagnostic, on a generated file nobody will proof-read. One
    long line has no such edge, and the cost is only readability.
    """
    if not job or any(c in job for c in ' \\/:*?"<>|'):
        raise DetachedLaunchError("job name %r must be a bare filename-safe token" % job)
    if not command:
        raise DetachedLaunchError("a launcher with no command would write EXIT=0 having done "
                                  "nothing, which reads as a successful run")

    # Extra environment, set INSIDE `setlocal` so it cannot leak to the launching shell.
    #
    # This exists because `PA_REPO` is read at MODULE IMPORT time by atlas_embed, so it must be
    # in the environment before python starts -- a runtime argument could not reach it. Values are
    # quoted `set "K=V"` form, which is what keeps a path containing spaces intact; this host's
    # repo path has one, so the unquoted form would silently truncate it.
    env = dict(env or {})
    for k, val in env.items():
        if not k or any(c in k for c in ' ="%'):
            raise DetachedLaunchError("env name %r is not settable from a .bat" % k)
        if '"' in str(val) or "%" in str(val):
            raise DetachedLaunchError(
                "env value for %r contains a quote or a %% -- cmd would mangle it and the job "
                "would start with a silently wrong value, which is worse than refusing" % k)

    log = _win("%s/%s.log" % (state_rel, job))
    sentinel = _win("%s/%s.sentinel" % (state_rel, job))
    lock = _win("%s/%s.lock" % (state_rel, job))
    guarded = " ".join(command)

    head = ["REM " + ln for ln in
            (rationale or ("Detached launcher for job `%s`." % job)).splitlines()]
    lines = [
        "@echo off",
        "REM GENERATED by scripts/detached_launch.py -- edit the generator, not this file.",
        *head,
        "REM",
        "REM Launch DETACHED from PowerShell so it outlives the session:",
        "REM   Start-Process -FilePath <repo>\\run_%s_detached.bat "
        "-WorkingDirectory <repo> -WindowStyle Hidden" % job,
        "",
        "setlocal",
        'cd /d "%~dp0"',
        "",
        'set "LOG=%~dp0' + log + '"',
        'set "SENTINEL=%~dp0' + sentinel + '"',
        'set "LOCK=%~dp0' + lock + '"',
        "",
        "REM A sentinel from a PREVIOUS run must go before this one starts: a stale EXIT=0 left",
        "REM beside a live job reads as COMPLETE and would hide the run entirely.",
        'if exist "%SENTINEL%" del "%SENTINEL%"',
        "",
        "REM Unbuffered, or a death loses the last minutes of the log -- which is precisely the",
        "REM window that explains the death.",
        "set PYTHONUNBUFFERED=1",
        *([""] + ["REM Job environment, set inside setlocal so it cannot leak to the caller.",
                  "REM PA_REPO and friends are read at module import, so a runtime argument",
                  "REM could not reach them."]
          + ['set "%s=%s"' % (k, v) for k, v in sorted(env.items())]
          if env else []),
        "",
        'echo ==== start %DATE% %TIME% ==== > "%LOG%"',
        "",
        "REM Guarded by the PID-bearing run lock: a second launch over a long window exits 9",
        "REM rather than running two concurrent torch jobs, which is this host's documented",
        "REM access-violation hazard. `run` (not `acquire`) holds the lock for the CHILD's",
        "REM lifetime -- an `acquire` process exits at once and its pid is reclaimed as stale.",
        'uv run python scripts/glm_run_lock.py run "%LOCK%" --label "' + job + '" -- '
        + guarded + ' >> "%LOG%" 2>&1',
        "",
        'set "RC=%ERRORLEVEL%"',
        "",
        'echo ==== end %DATE% %TIME% rc=%RC% ==== >> "%LOG%"',
        "",
        "REM THE REDIRECT IS WRITTEN BEFORE THE ECHO ON PURPOSE. `echo EXIT=%RC%>>file` expands",
        "REM to `echo EXIT=0>>file`, and cmd reads a digit immediately preceding `>>` as a FILE",
        "REM HANDLE -- it redirects handle 0 (stdin) and leaves the sentinel EMPTY. Measured on",
        "REM run_glm_anchor_detached.bat, not theorised.",
        SENTINEL_LINE_CORRECT,
        "",
        "endlocal",
        "exit /b %RC%",
        "",
    ]
    return "\n".join(lines)


def write_bat(job: str, command: list[str], *, rationale: str = "",
              repo: Path | None = None, state_rel: str = "data/unattended",
              env: dict | None = None) -> Path:
    repo = repo or REPO
    text = render_bat(job, command, rationale=rationale, state_rel=state_rel, env=env)
    out = repo / ("run_%s_detached.bat" % job)
    # CRLF: cmd tolerates LF, but every hand-written launcher here is CRLF and a mixed tree
    # invites a diff nobody can read.
    out.write_bytes(text.replace("\n", "\r\n").encode("ascii"))
    (repo / _win(state_rel)).mkdir(parents=True, exist_ok=True)
    return out


# --------------------------------------------------------------------------- cli


def render_status(rows: list[JobStatus]) -> str:
    if not rows:
        return "no jobs known (no heartbeat or sentinel under the state dir)"
    out = []
    for r in rows:
        mark = "!!" if r.needs_attention else "  "
        prog = ""
        if r.units_done is not None and r.units_total:
            prog = " %s/%s" % (r.units_done, r.units_total)
        out.append("%s %-28s %-18s%s" % (mark, r.job, r.state, prog))
        out.append("     " + r.reason)
        if r.state in (STATE_RUNNING, STATE_STALLED) and r.estimated_finish:
            out.append("     estimated finish: " + r.estimated_finish)
    return "\n".join(out)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    after: list[str] = []
    if "--" in argv:
        i = argv.index("--")
        argv, after = argv[:i], argv[i + 1:]
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("action", choices=["emit", "status"])
    p.add_argument("job", nargs="?")
    p.add_argument("--rationale", default="")
    p.add_argument("--json", action="store_true")
    a = p.parse_args(argv)

    if a.action == "emit":
        if not a.job:
            print("emit requires a job name", file=sys.stderr)
            return 2
        if not after:
            print("emit requires a command after `--`", file=sys.stderr)
            return 2
        out = write_bat(a.job, after, rationale=a.rationale)
        print("wrote %s" % out)
        return 0

    rows = [classify(a.job)] if a.job else status_all()
    if a.json:
        print(json.dumps([r.as_dict() for r in rows], indent=2))
    else:
        print(render_status(rows))
    return 1 if any(r.needs_attention for r in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
