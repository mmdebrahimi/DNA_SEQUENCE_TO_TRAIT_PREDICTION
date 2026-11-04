"""Relaunch the fine-tune arm's detached run if, and only if, it has stopped without finishing.

WHY THIS EXISTS NOW AND NOT BEFORE. `FUTURE_FEATURES.md` named this and deliberately left it
unbuilt: "a task that relaunches a training job unattended is a larger change than the
one-invocation gap it closes, and it should be an explicit decision rather than a convenience." The
decision has now been made explicitly, for a ten-day absence in which the run finishes around day
four. Without it, a single host reboot on day two silently costs the entire window -- and the run's
own design makes the fix cheap, because every epoch writes a checkpoint and a relaunch resumes
bit-for-bit rather than restarting.

IT IS SAFE BECAUSE IT CANNOT STACK A SECOND TORCH JOB. Two concurrent torch jobs are this host's
documented access-violation hazard. Three independent things prevent it: this script refuses while a
LIVE pid holds the run lock (`glm_run_lock.holder` already treats a stale or malformed lock as free,
so a crash does not wedge it); the launcher it invokes is itself lock-guarded and exits 9 rather than
starting a second payload; and the runner carries its own wall-clock ceiling.

THE "DONE" SIGNAL IS THE ARTIFACT, NOT THE SENTINEL, and that distinction is load-bearing. A
sentinel can read EXIT=9 simply because this watchdog raced a healthy run and the launcher refused --
so treating "a sentinel exists" as "the run finished" would stop the watchdog permanently on its own
first collision. The run artifact only appears when `run()` returned a record.

IT IS BOUNDED. A deterministically-failing run would otherwise be relaunched every interval for ten
days. After `MAX_RELAUNCHES` it REFUSES and says so, leaving the diagnosis to a human instead of
burning the window on a loop. Every decision is appended to a log with its reason, so the history is
readable afterwards rather than inferred.

Usage:
    python scripts/watchdog_finetune_arm.py            # report the decision, change nothing
    python scripts/watchdog_finetune_arm.py --apply    # act on it
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
JOB = "finetune_arm"
STATE_DIR = REPO / "data/unattended"
LAUNCHER = REPO / ("run_%s_detached.bat" % JOB)
ARTIFACT = REPO / "wiki/finetune_arm_run_2026-11-13.json"

#: Enough for a plausible number of host reboots across a ten-day window, and few enough that a
#: deterministically-broken run stops rather than looping. Not a tuned number; a bound.
MAX_RELAUNCHES = 12

#: Decisions. Only one of these causes an action.
DONE = "artifact_present"
RUNNING = "running"
EXHAUSTED = "relaunch_budget_exhausted"
RELAUNCH = "relaunch"
NO_LAUNCHER = "launcher_absent"
NOT_OWNER = "not_the_owning_host"

ALL_DECISIONS = (DONE, RUNNING, EXHAUSTED, RELAUNCH, NO_LAUNCHER, NOT_OWNER)

#: AN EXPLICIT OWNERSHIP MARKER, and it exists because of a footgun this script created.
#:
#: The run lives on server2. But `run_finetune_arm_detached.bat` is COMMITTED, so it exists in the
#: laptop's checkout too, and there the cohort cache also resolves -- so a watchdog run on the
#: laptop would find no artifact, no live lock, a present launcher, and cheerfully start a
#: four-thread torch job on the host that hard-powered-off three times in eighteen hours. Every
#: other check in `decide` would have been satisfied.
#:
#: So relaunching requires this marker to exist AND to name this machine. It is deliberately NOT
#: inferred from anything: a path that happens to resolve is not evidence that a host owns a run,
#: which is the same reasoning that stops `external_root.py` adopting a drive by its name.
OWNER_MARKER_NAME = "%s.owner" % JOB


def state_path(state_dir: Path | None = None) -> Path:
    return (state_dir or STATE_DIR) / ("%s.watchdog.json" % JOB)


def read_state(state_dir: Path | None = None) -> dict:
    """The watchdog's own history. A malformed file reads as EMPTY, never as exhausted.

    Reading damage as "budget spent" would silently disable the watchdog, which is the failure
    direction that costs the window -- the same malformed-versus-absent distinction the sidecar
    reader got wrong.
    """
    p = state_path(state_dir)
    try:
        rec = json.loads(p.read_text(encoding="utf-8"))
        if not isinstance(rec, dict):
            return {"relaunches": 0, "history": []}
        rec.setdefault("relaunches", 0)
        rec.setdefault("history", [])
        return rec
    except (OSError, ValueError):
        return {"relaunches": 0, "history": []}


def write_state(rec: dict, state_dir: Path | None = None) -> Path:
    """Atomic durable write: tmp -> flush -> fsync -> replace."""
    final = state_path(state_dir)
    final.parent.mkdir(parents=True, exist_ok=True)
    tmp = final.with_suffix(final.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, indent=1, sort_keys=True))
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, final)
    return final


def owner_path(state_dir: Path | None = None) -> Path:
    return (state_dir or STATE_DIR) / OWNER_MARKER_NAME


def read_owner(state_dir: Path | None = None) -> str | None:
    """The hostname this run belongs to, or None when unclaimed."""
    try:
        txt = owner_path(state_dir).read_text(encoding="utf-8", errors="replace").strip()
    except OSError:
        return None
    return txt or None


def decide(*, artifact: Path, lock_holder: dict | None, relaunches: int,
           launcher: Path, max_relaunches: int = MAX_RELAUNCHES,
           owner: str | None = None, this_host: str = "") -> dict:
    """PURE. What to do, and why. Order matters and is the contract.

    `artifact` is checked FIRST so a finished run can never be relaunched, and the live lock holder
    SECOND so a healthy run is never raced. Both are read-only facts and are reported truthfully on
    ANY host. Ownership is checked only before the branches that would ACT.
    """
    if Path(artifact).is_file():
        return {"decision": DONE,
                "reason": "the run artifact exists, so run() returned a record; nothing to do"}
    if lock_holder is not None:
        return {"decision": RUNNING,
                "reason": "pid %s holds the run lock and is alive (label %r)"
                          % (lock_holder.get("pid"), lock_holder.get("label"))}
    if not Path(launcher).is_file():
        return {"decision": NO_LAUNCHER,
                "reason": "launcher absent at %s; refusing to improvise a command" % launcher}
    if not owner or owner != this_host:
        return {"decision": NOT_OWNER,
                "reason": ("this run is owned by %r and this host is %r, so relaunching here is "
                           "refused. The launcher is committed and its cache paths resolve on more "
                           "than one machine, so every other check would have passed -- a path "
                           "that resolves is not evidence that a host owns a run."
                           % (owner or "(unclaimed)", this_host))}
    if relaunches >= max_relaunches:
        return {"decision": EXHAUSTED,
                "reason": ("already relaunched %d times against a bound of %d. Refusing to keep "
                           "looping: a run that dies this often is failing deterministically and "
                           "needs a human, not another attempt." % (relaunches, max_relaunches))}
    return {"decision": RELAUNCH,
            "reason": ("no artifact and no live lock holder, so the run stopped without "
                       "finishing. Relaunching; it resumes from the last per-epoch checkpoint.")}


def current_decision(*, state_dir: Path | None = None, artifact: Path | None = None,
                     launcher: Path | None = None, max_relaunches: int = MAX_RELAUNCHES) -> dict:
    import platform

    from scripts.glm_run_lock import holder

    sd = state_dir or STATE_DIR
    st = read_state(sd)
    host = platform.node()
    d = decide(artifact=artifact or ARTIFACT,
               lock_holder=holder(sd / ("%s.lock" % JOB)),
               relaunches=int(st.get("relaunches", 0)),
               launcher=launcher or LAUNCHER,
               max_relaunches=max_relaunches,
               owner=read_owner(sd), this_host=host)
    d["relaunches_so_far"] = int(st.get("relaunches", 0))
    d["max_relaunches"] = max_relaunches
    d["this_host"] = host
    d["owner"] = read_owner(sd)
    return d


def launch(launcher: Path | None = None) -> int:
    """Start the launcher DETACHED. Measured on this host: a process started over ssh outlives the
    session, so no scheduled-task indirection is needed for the child itself."""
    lp = Path(launcher or LAUNCHER)
    proc = subprocess.Popen(["cmd", "/c", str(lp)], cwd=str(lp.parent),
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            close_fds=True)
    return proc.pid


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="act on the decision")
    ap.add_argument("--state-dir", default=str(STATE_DIR))
    ap.add_argument("--max-relaunches", type=int, default=MAX_RELAUNCHES)
    ap.add_argument("--claim", action="store_true",
                    help="record THIS host as the one allowed to relaunch this run, then exit")
    a = ap.parse_args(argv)

    sd = Path(a.state_dir)
    if a.claim:
        import platform

        sd.mkdir(parents=True, exist_ok=True)
        owner_path(sd).write_text(platform.node(), encoding="utf-8")
        print("claimed: %s may relaunch job %r" % (platform.node(), JOB))
        return 0
    d = current_decision(state_dir=sd, max_relaunches=a.max_relaunches)
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    print("%s  %s  %s" % (stamp, d["decision"], d["reason"]))

    if not a.apply:
        print("DRY RUN -- nothing started.")
        return 0

    st = read_state(sd)
    entry = {"at": stamp, "decision": d["decision"], "reason": d["reason"]}
    if d["decision"] == RELAUNCH:
        entry["launched_pid"] = launch()
        st["relaunches"] = int(st.get("relaunches", 0)) + 1
        print("relaunched; pid %s; relaunch %d of %d"
              % (entry["launched_pid"], st["relaunches"], a.max_relaunches))
    # Every decision is recorded, not just the actions: a history of "running" entries is how a
    # reader later confirms the watchdog was actually alive the whole time.
    st["history"] = (st.get("history") or [])[-199:] + [entry]
    st["last_decision"] = d["decision"]
    st["last_checked"] = stamp
    write_state(st, sd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
