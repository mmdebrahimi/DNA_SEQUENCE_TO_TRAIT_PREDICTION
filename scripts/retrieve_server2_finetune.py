"""Bring the fine-tune arm's result back from server2, then say what to run next.

WHY THIS EXISTS. The fine-tune run was moved to server2 because the laptop that declared the arm
hard-powered-off three times in eighteen hours under exactly this load. That makes the result
REACHABLE but not LOCAL, and the collector's launcher gate reads the LAPTOP's `data/unattended` --
so without the job state coming back too, a finished run still reads `NOT_STARTED` here and the arm
is refused. Copying the artifact alone would be the subtler failure: the number would be present and
the gate would still refuse it, which looks like the run never happened.

IT VERIFIES BY CONTENT, NEVER BY ARRIVAL. scp exits 0 on this host while transferring nothing when a
path is mis-spelled, so every file is checked after the copy: the run artifact must parse and carry
the expected schema, and the sentinel must be NON-EMPTY -- an empty sentinel is the documented
`echo EXIT=%RC%>>file` redirect trap, and it is indistinguishable from a clean exit if nobody looks.

IT REFUSES TO OVERWRITE A DIFFERENT LOCAL ARTIFACT. If a run artifact is already here and differs,
one of the two is from another run, and silently replacing it would destroy the question of which.

Usage:
    python scripts/retrieve_server2_finetune.py            # dry run: report what is there
    python scripts/retrieve_server2_finetune.py --apply    # copy it back and verify
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

HOST = "server2"
REMOTE_REPO = "D:/finetune_arm_run/repo"
JOB = "finetune_arm"

RUN_ARTIFACT = "wiki/finetune_arm_run_2026-11-13.json"
RUN_SCHEMA = "finetune-arm-run-v1"

#: The job-state files. The sentinel is what turns a finished run into a COMPLETE launcher state,
#: which is what licenses the collector to read the artifact at all.
STATE_FILES = (
    "data/unattended/%s.sentinel" % JOB,
    "data/unattended/%s.heartbeat.json" % JOB,
    "data/unattended/%s.log" % JOB,
)


class RetrieveRefused(RuntimeError):
    """The result cannot be brought back, or cannot be trusted once here."""


def _remote(path: str) -> str:
    """scp source spec. The `/D:/...` form is what works from Git Bash on this host."""
    return "%s:/%s/%s" % (HOST, REMOTE_REPO.replace("D:/", "D:/"), path)


def remote_listing(*, host: str = HOST, remote_repo: str = REMOTE_REPO) -> dict:
    """Which files exist remotely. Read-only, ONE ssh call.

    Deliberately `if exist` rather than a `for` loop over file sizes: a `for` body needs `%%A` in a
    batch file and `%A` under `cmd /c`, and getting that backwards fails with "%%A was unexpected at
    this time" -- which says nothing about the run. Presence plus the sentinel's own text is the
    information that matters anyway, and `verify_local` does the real checking after the copy.
    """
    paths = (RUN_ARTIFACT,) + STATE_FILES
    win = remote_repo.replace("/", "\\")
    out = {}
    for p in paths:
        # ONE call per path, deliberately. Chaining `if exist (..) else (..)` with `&` lets the
        # following command bind INSIDE the else branch, so a row silently vanishes depending on
        # whether the previous file happened to exist -- measured, after one row went missing.
        cmd = 'cmd /c "if exist "%s\\%s" (echo PRESENT) else (echo ABSENT)"' % (
            win, p.replace("/", "\\"))
        try:
            res = subprocess.run(["ssh", "-o", "BatchMode=yes", host, cmd],
                                 capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.SubprocessError) as e:
            raise RetrieveRefused("could not reach %s: %s" % (host, e)) from e
        said = [ln.strip() for ln in (res.stdout or "").splitlines()
                if ln.strip() in ("PRESENT", "ABSENT")]
        if not said:
            raise RetrieveRefused(
                "no answer for %r from %s -- refusing to report it as absent, because a failed "
                "check and a missing file are different facts. stderr=%r"
                % (p, host, (res.stderr or "")[-200:]))
        out[p] = said[-1]
    return out


def copy_back(*, host: str = HOST, remote_repo: str = REMOTE_REPO,
              repo: Path | None = None) -> list:
    """scp each file, then VERIFY it by content. Returns the local paths written."""
    repo = repo or REPO
    written = []
    for rel in (RUN_ARTIFACT,) + STATE_FILES:
        dest = Path(repo) / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        src = "%s:/%s/%s" % (host, remote_repo, rel)
        res = subprocess.run(["scp", "-o", "BatchMode=yes", src, str(dest)],
                             capture_output=True, text=True, timeout=600)
        if res.returncode != 0:
            raise RetrieveRefused("scp of %s failed: %s" % (rel, (res.stderr or "")[-300:]))
        if not dest.is_file():
            raise RetrieveRefused("scp reported success but %s is not here -- an exit code is not "
                                  "evidence of a transfer on this host" % dest)
        written.append(dest)
    return written


def verify_local(*, repo: Path | None = None) -> dict:
    """Check what came back is usable. Raises rather than returning a soft verdict."""
    repo = repo or REPO
    art = Path(repo) / RUN_ARTIFACT
    if not art.is_file():
        raise RetrieveRefused("%s is absent" % RUN_ARTIFACT)
    try:
        run = json.loads(art.read_text(encoding="utf-8"))
    except ValueError as e:
        raise RetrieveRefused("%s is not readable JSON (%s)" % (RUN_ARTIFACT, e)) from e
    if run.get("schema") != RUN_SCHEMA:
        raise RetrieveRefused("run schema is %r, expected %r" % (run.get("schema"), RUN_SCHEMA))

    sentinel = Path(repo) / STATE_FILES[0]
    if not sentinel.is_file():
        raise RetrieveRefused("the sentinel is absent, so the launcher state cannot read COMPLETE")
    text = sentinel.read_text(encoding="utf-8", errors="replace").strip()
    if not text:
        raise RetrieveRefused(
            "the sentinel is EMPTY. That is the documented `echo EXIT=%RC%>>file` trap -- cmd reads "
            "the digit before `>>` as a file handle -- and an empty sentinel is indistinguishable "
            "from a clean exit unless someone looks.")
    return {"schema": run.get("schema"), "status": run.get("status"),
            "n_seeds": len(run.get("per_seed") or {}),
            "horizon_epochs": run.get("horizon_epochs"), "sentinel": text}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="copy the files back")
    ap.add_argument("--host", default=HOST)
    ap.add_argument("--remote-repo", default=REMOTE_REPO)
    ap.add_argument("--force", action="store_true",
                    help="replace a local run artifact that differs from the remote one")
    a = ap.parse_args(argv)

    try:
        listing = remote_listing(host=a.host, remote_repo=a.remote_repo)
    except RetrieveRefused as e:
        print("REFUSED: %s" % e)
        return 2

    print("ON %s:" % a.host)
    for k, v in sorted(listing.items()):
        print("  %-48s %s" % (k, v))
    print()

    if listing.get(RUN_ARTIFACT) in (None, "ABSENT"):
        print("The run has not produced its artifact yet, so there is nothing to score.")
        print("Check progress with:")
        print("  ssh %s \"cmd /c type %s\\data\\unattended\\%s.log\""
              % (a.host, a.remote_repo.replace("/", "\\"), JOB))
        return 1

    local = REPO / RUN_ARTIFACT
    if local.is_file() and not a.force:
        print("A local %s already exists. Not replacing it without --force: if the two differ, one "
              "is from another run and overwriting destroys the question of which." % RUN_ARTIFACT)
        return 1

    if not a.apply:
        print("DRY RUN. Re-run with --apply to copy these back.")
        return 0

    try:
        written = copy_back(host=a.host, remote_repo=a.remote_repo)
        info = verify_local()
    except RetrieveRefused as e:
        print("REFUSED: %s" % e)
        return 2

    for p in written:
        print("  <- %s" % p)
    print()
    print("VERIFIED  schema %s  status %s  seeds %s  horizon %s  sentinel %r"
          % (info["schema"], info["status"], info["n_seeds"], info["horizon_epochs"],
             info["sentinel"]))
    print()
    print("NEXT, in order:")
    print("  uv run python scripts/detached_launch.py status")
    print("  uv run python scripts/finetune_arm_score.py")
    print("  uv run python scripts/collect_unattended_runs.py")
    print("  uv run python scripts/record_successor_arm_results.py --apply")
    print("  uv run python scripts/successor_arms_memo.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
