"""Wait for stage B to finish, validate its table, then run stage C. Unattended.

    python scripts/promcalc_chain_stage_c.py            # the real chain
    python scripts/promcalc_chain_stage_c.py --dry-run  # one pass, decide nothing, run nothing

**Why this exists.** Stage B takes ~12 h and finishes in the middle of the night. Without a chain the
morning state is a finished TABLE; with it, the morning state is a VERDICT. Stage C is cheap once the
table exists, so the only reason it would not already be done is that nobody was awake to start it.

**The gate order is load-bearing and is not negotiable.** validate THEN score, never the reverse:

  1. Stage B complete -- sentinel `EXIT=0`, or `completed_tiles` == 44,106.
  2. No promcalc worker alive. The validator decompresses the whole table, and doing that during a
     live run is the recorded trap that turned a ~1 h job into a 5.0 h projection -- observing the
     sweep slowed the sweep, and the slowed rate was then quoted as the sweep's rate.
  3. `promcalc_validate_output.py` with NO `--expect-partial`. Exit 0 = scoreable. Anything else
     STOPS the chain and stage C is never invoked.
  4. Stage C, the frozen pre-registered comparison. Its own reconcile gate is the validity anchor.

Step 3 is the whole point. `_site_field` defaults on KeyError/AttributeError/TypeError, so a future
field rename yields a ~140 MB all-nan table with 0 errors, full 44,106 coverage and
`sites_per_tile_mean` still 73.0 -- every count reporting success. Scoring that unattended would
produce `INDETERMINATE_RECONCILE_FAILED`, which is indistinguishable from genuine model drift, and
the morning reader would have a plausible wrong diagnosis instead of a defect report.

**What it deliberately does NOT do:** it never kills anything, never retries stage B (the watchdog
owns that, with its own retry cap), never edits the frozen bar, and never scores a partial table. A
stage-B crash leaves the chain waiting, which is correct -- the watchdog relaunches, and the chain
picks up whenever coverage reaches the full count.

Stdlib only, so it runs under either interpreter. Stage C is invoked through `uv run` because it
needs the project env; the chain itself needs nothing.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
IN_DIR = REPO / "data" / "promcalc"
MPRA = REPO / "data" / "mpra"
SENTINEL = IN_DIR / "decompose_sentinel.txt"
LOG = IN_DIR / "chain_stage_c.log"
LOCK = IN_DIR / "chain_stage_c.lock"
VERDICT = IN_DIR / "chain_stage_c_verdict.txt"

EXPECT_TILES = 44106
POLL_S = 300
#: Generous: the sweep's own ETA plus room for a watchdog relaunch. Not a deadline for stage B --
#: just a bound so this process cannot outlive the night it was started for.
MAX_WAIT_S = 14 * 3600
#: After completion, let the parent flush and exit before touching the table.
SETTLE_S = 60


def log(msg: str) -> None:
    line = f"{datetime.now():%Y-%m-%d %H:%M:%S}  {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def completed_tiles() -> int | None:
    """Exact count via stage B's own resume logic. None if it cannot be read."""
    try:
        sys.path.insert(0, str(REPO / "scripts"))
        import promcalc_decompose as pd  # stdlib-only module imports by construction

        return len(pd.completed_tiles(IN_DIR))
    except Exception as exc:  # noqa: BLE001 - a read failure must not end the chain
        log(f"coverage read failed ({type(exc).__name__}: {exc}); will retry")
        return None


def sentinel_exit() -> str | None:
    """`EXIT=<n>` if stage B wrote one. Absent means killed before the exit line."""
    try:
        for ln in reversed(SENTINEL.read_text(encoding="utf-8", errors="replace").splitlines()):
            if ln.strip().startswith("EXIT="):
                return ln.strip()
    except OSError:
        return None
    return None


def workers_alive() -> int:
    """Count live promcalc processes.

    Deliberately over-broad, matching the watchdog's rationale: `uv venv`'s python.exe is a 46 KB
    trampoline, so multiprocessing spawn launches the BASE interpreter and a venv-path grep misses
    every worker. Over-counting only makes the chain wait longer, which is the safe direction.
    """
    ps = (
        "(Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
        "Where-Object { $_.CommandLine -like '*promcalc_decompose*' -or "
        "$_.CommandLine -like '*multiprocessing*' }).Count"
    )
    try:
        out = subprocess.run(
            ["powershell", "-NoProfile", "-Command", ps],
            capture_output=True, text=True, timeout=60,
        ).stdout.strip()
        return int(out or 0)
    except Exception:  # noqa: BLE001
        return -1  # unknown -> treated as "do not proceed"


def stage_b_done() -> tuple[bool, str]:
    sx = sentinel_exit()
    if sx == "EXIT=0":
        return True, "sentinel EXIT=0"
    n = completed_tiles()
    if n is not None and n >= EXPECT_TILES:
        return True, f"coverage {n}/{EXPECT_TILES}"
    return False, f"coverage {n}/{EXPECT_TILES}" + (f", sentinel {sx}" if sx else "")


def run(cmd: list[str], tag: str) -> int:
    log(f"RUN {tag}: {' '.join(cmd)}")
    with LOG.open("a", encoding="utf-8") as fh:
        rc = subprocess.run(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT).returncode
    log(f"{tag} exit={rc}")
    return rc


def finish(verdict: str, detail: str) -> int:
    VERDICT.write_text(
        f"{datetime.now():%Y-%m-%d %H:%M:%S}\n{verdict}\n{detail}\n", encoding="utf-8"
    )
    log(f"=== CHAIN VERDICT: {verdict} === {detail}")
    return 0 if verdict.startswith("SCORED") else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="one pass: report what the chain sees and would do, then exit. Runs "
                         "neither the validator nor stage C.")
    a = ap.parse_args(argv)

    if a.dry_run:
        done, why = stage_b_done()
        log(f"DRY-RUN: stage_b_done={done} ({why}); workers={workers_alive()}")
        return 0

    if LOCK.exists():
        log(f"another chain holds {LOCK.name}; refusing to double-run")
        return 3
    LOCK.write_text(str(os.getpid()), encoding="utf-8")
    try:
        log(f"chain start (poll {POLL_S}s, max wait {MAX_WAIT_S // 3600}h, expect {EXPECT_TILES})")
        deadline = time.time() + MAX_WAIT_S

        while time.time() < deadline:
            done, why = stage_b_done()
            if done:
                log(f"stage B complete: {why}")
                break
            log(f"waiting: {why}")
            time.sleep(POLL_S)
        else:
            return finish("TIMED_OUT", f"stage B not complete within {MAX_WAIT_S // 3600}h")

        # Never decompress the table while the sweep is writing it.
        while time.time() < deadline:
            n = workers_alive()
            if n == 0:
                break
            log(f"{n} promcalc process(es) still alive (or unknown); waiting before validation")
            time.sleep(POLL_S)
        else:
            return finish("TIMED_OUT", "workers never cleared")

        log(f"settling {SETTLE_S}s so the parent can flush and exit")
        time.sleep(SETTLE_S)

        rc = run([sys.executable, "scripts/promcalc_validate_output.py", str(IN_DIR)], "validate")
        if rc != 0:
            return finish(
                "BLOCKED_VALIDATOR",
                f"validator exit={rc} (2=defect, 3=indeterminate). Stage C NOT invoked -- scoring a "
                f"bad table would yield a reconcile failure indistinguishable from model drift. "
                f"See {LOG.name}.",
            )

        rc = run(
            ["uv", "run", "python", "scripts/promcalc_decomposition_arm.py",
             "--in-dir", str(IN_DIR), "--mpra-cache", str(MPRA), "--seeds", "10"],
            "stage-C",
        )
        if rc != 0:
            return finish("STAGE_C_NONZERO",
                          f"stage C exit={rc}. A refusal is a legitimate outcome (incomplete "
                          f"coverage exits 1); read {LOG.name} for which.")
        return finish("SCORED", f"stage C exit=0; verdict is in its own artifact and {LOG.name}.")
    finally:
        LOCK.unlink(missing_ok=True)
        log("lock released")


if __name__ == "__main__":
    sys.exit(main())
