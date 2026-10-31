"""Resumable driver for the 270-seed fragment-null resolution across Kaggle batches.

WHY A DRIVER. The run is ~18 h of T4 across four batches, Kaggle permits exactly TWO concurrent
batch GPU sessions, and a session boundary must not lose compute. This is deliberately NOT a daemon:
every action is idempotent, so `step` can be run at any time, by any session, as often as you like.
It pulls whatever has finished, launches the next batch when a slot is free, and otherwise reports.

    uv run python scripts/fragment_null_batch_driver.py status
    uv run python scripts/fragment_null_batch_driver.py step      # pull done, launch next
    uv run python scripts/fragment_null_batch_driver.py pool      # pooled null over all seeds

THE BATCH PLAN IS NOT RESTATED HERE. It is imported from the kernel that runs it, because a plan
duplicated beside its definition drifts -- and the one thing that must not drift is which seeds have
been drawn. See `scripts/kaggle/fragment_null_batch_kernel.py`.

THE PULL RECIPE IS THE LOAD-BEARING PART. The kernel clones the repo, so its output listing is
~3,600 files and a plain `kernels_output` 429s while paging the listing server-side -- measured
again 2026-10-30, including the trap that `file_pattern` is a REGEX (not a glob) and filters only
what is DOWNLOADED, not what is LISTED, so it does not avoid the 429. It also matched a false
positive inside the cloned repo. The working recipe, verified against the completed 10..49 kernel:
walk `kernels_list_files` at page_size=200, find the page whose listing contains the exact artifact
name (it is the LAST page, 18 of 18, but this searches rather than assuming), then make ONE
`kernels_output` call with that page's request token. Matched files come from signed storage URLs,
not the throttled listing endpoint.

WHAT IT WILL NOT DO. It never pushes to a slug that already holds output (Kaggle keeps one output
per kernel, so re-pushing destroys retrieved artifacts), and it never writes a pooled verdict --
pooling is reported by `pool` for a human to read against
`wiki/fragment_null_270seed_prereg_2026-10-30.md`, whose three outcome branches were fixed before
any of these seeds existed.
"""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KERNEL = ROOT / "scripts" / "kaggle" / "fragment_null_batch_kernel.py"
STAGE = ROOT / "wiki" / "fragment_null_batches"
USER = "emanueleebrahimi"
MAX_CONCURRENT_GPU = 2          # measured: a third push is refused outright
SLUG = "glm-fragment-null-b%d"

#: Seeds already measured BEFORE this campaign, and where they live. Not redrawn.
PRIOR = {
    "wiki/glm_fragment_pyramid_arm_kaggle_full_2026-10-28.json": (0, 10),
    "wiki/fragment_null_ext_summary_2026-10-29.json": (10, 40),
}


def _batches() -> dict:
    """The batch plan, imported from the kernel rather than restated."""
    spec = importlib.util.spec_from_file_location("_fnbk", KERNEL)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return {int(k): tuple(v) for k, v in m.BATCHES.items()}


def _api():
    from kaggle import KaggleApi
    api = KaggleApi()
    api.authenticate()
    return api


#: A kernel that DOES exist and is complete, used only to tell "absent" from "broken auth".
#: Never pushed to -- it holds the retrieved seeds 10..49 output.
KNOWN_GOOD = "glm-fragment-null-ext"

#: What `kernels_status` raises for a kernel that does not exist. Measured 2026-10-30: it is a
#: ValueError reading "Permission 'kernels.get' was denied ... most likely cause is a wrong kernel
#: slug", NOT an HTTPError carrying 404.
_NOT_FOUND = ("kernels.get", "wrong kernel slug", "404")


def _looks_absent(err: Exception) -> bool:
    s = str(err).lower()
    return any(m.lower() in s for m in _NOT_FOUND)


def _api_healthy(api) -> bool:
    """Is the API actually answering, or has the ~12 h OAuth lapsed?

    THIS IS NOT DEFENSIVE PADDING. An expired token raises the SAME not-found-shaped ValueError as a
    kernel that was never created, so without this probe `step` would read every unlaunched batch as
    ABSENT, try to launch one, and fail again in a different place -- reporting "nothing launched
    yet" during an auth outage. The project's recorded discriminator for exactly this is to repeat
    the call on a KNOWN-GOOD slug, so that is what this does.
    """
    try:
        api.kernels_status("%s/%s" % (USER, KNOWN_GOOD))
        return True
    except Exception:
        return False


def _status(api, slug: str, healthy: bool) -> str:
    """Kernel status, or a word describing why there is none. Never raises."""
    try:
        st = api.kernels_status("%s/%s" % (USER, slug))
    except Exception as e:
        if _looks_absent(e):
            # Only absent if the API is demonstrably answering; otherwise this is an auth/API
            # failure wearing the costume of an unlaunched batch.
            return "ABSENT" if healthy else "API_UNHEALTHY"
        return "ERROR:%s" % type(e).__name__
    s = getattr(st, "status", None)
    return str(s).rsplit(".", 1)[-1] if s is not None else "UNKNOWN"


def _healthy_or_explain(api) -> bool:
    if _api_healthy(api):
        return True
    print("REFUSING: the Kaggle API will not answer for a known-good kernel (%s), so an absent\n"
          "  batch cannot be told from an expired token. Kaggle OAuth lapses after ~12 h.\n"
          "  Fix it in YOUR OWN terminal (credentials must not route through a session):\n"
          "    uv run kaggle auth login --force" % KNOWN_GOOD, file=sys.stderr)
    return False


def _artifacts(batch: int) -> list:
    return ["fragment_null_b%d_summary.json" % batch,
            "glm_fragment_pyramid_arm_kaggle_b%d_ext.json" % batch]


def _pull(api, slug: str, names: list, dest: Path) -> list:
    """ONE kernels_output call, aimed at the page that actually lists `names`."""
    dest.mkdir(parents=True, exist_ok=True)
    kid = "%s/%s" % (USER, slug)
    token, hit, found = None, None, False
    while True:
        r = api.kernels_list_files(kid, page_token=token, page_size=200)
        files = [getattr(f, "name", str(f)) for f in (r.files or [])]
        if any(n in files for n in names):
            hit, found = token, True
            break
        nxt = getattr(r, "nextPageToken", None) or getattr(r, "next_page_token", None)
        if not nxt or not files:
            break
        token = nxt
    if not found:
        return []
    pattern = "^(%s)$" % "|".join(re.escape(n) for n in names)
    got, _ = api.kernels_output(kid, path=str(dest), file_pattern=pattern,
                               page_token=hit, page_size=200, quiet=True)
    return [str(g) for g in got]


def _pulled(batch: int):
    p = STAGE / ("fragment_null_b%d_summary.json" % batch)
    return p if p.is_file() else None


def cmd_status() -> int:
    api, batches = _api(), _batches()
    if not _healthy_or_explain(api):
        return 2
    print("270-seed null resolution -- batch status")
    print("  prior seeds (not redrawn): 0..49 from %d artifacts" % len(PRIOR))
    for b, (start, n) in sorted(batches.items()):
        st = _status(api, SLUG % b, True)
        got = "pulled" if _pulled(b) else "-"
        print("  b%-2d seeds %3d..%-3d (n=%d)  %-12s %s" % (b, start, start + n - 1, n, st, got))
    return 0


def cmd_step() -> int:
    api, batches = _api(), _batches()
    if not _healthy_or_explain(api):
        return 2
    running, done, absent = [], [], []
    for b in sorted(batches):
        st = _status(api, SLUG % b, True)
        if st == "RUNNING":
            running.append(b)
        elif st == "COMPLETE":
            done.append(b)
        elif st == "ABSENT":
            absent.append(b)
        else:
            print("  b%d: %s -- needs a human look" % (b, st))

    for b in done:
        if _pulled(b):
            continue
        print("pulling b%d ..." % b, flush=True)
        got = _pull(api, SLUG % b, _artifacts(b), STAGE)
        print("  %d file(s): %s" % (len(got), ", ".join(Path(g).name for g in got) or "NONE"))

    slots = MAX_CONCURRENT_GPU - len(running)
    if not absent:
        print("all %d batches launched; %d running" % (len(batches), len(running)))
        return 0
    if slots <= 0:
        print("%d running, GPU cap is %d -- next batch waits (b%s)"
              % (len(running), MAX_CONCURRENT_GPU, ", b".join(str(b) for b in absent)))
        return 0

    nxt = absent[0]
    print("launching b%d (slot free: %d running) ..." % (nxt, len(running)), flush=True)
    src = KERNEL.read_text(encoding="utf-8")
    cur = re.search(r"^BATCH = (\d+)$", src, re.M)
    assert cur, "BATCH line not found in the kernel"
    if int(cur.group(1)) != nxt:
        with open(KERNEL, "w", encoding="utf-8", newline="\r\n") as fh:
            fh.write(re.sub(r"^BATCH = \d+$", "BATCH = %d" % nxt, src, count=1, flags=re.M))
        print("  set BATCH = %d" % nxt)
    rc = subprocess.run([sys.executable, str(ROOT / "scripts" / "kaggle_push_poll.py"),
                         "push", str(KERNEL), SLUG % nxt, "--gpu"], cwd=str(ROOT)).returncode
    print("  push rc=%d" % rc)
    return rc


def cmd_pool() -> int:
    """Pooled per-seed nulls over every artifact on disk. Reports; decides nothing."""
    import statistics as stats
    seeds = []
    sources = []
    for rel in PRIOR:
        p = ROOT / rel
        if not p.is_file():
            print("  MISSING prior artifact %s -- pool would be incomplete, refusing" % rel)
            return 2
        d = json.loads(p.read_text(encoding="utf-8"))
        null = d.get("ext_null") or d.get("null") or {}
        v = list(null.get("per_seed") or [])
        seeds += v
        sources.append((rel, len(v)))
    batch_files = sorted(STAGE.glob("fragment_null_b*_summary.json")) if STAGE.is_dir() else []
    for p in batch_files:
        d = json.loads(p.read_text(encoding="utf-8"))
        v = list(((d.get("ext_null") or {}).get("per_seed")) or [])
        seeds += v
        sources.append((p.name, len(v)))

    print("pooled null over %d seeds" % len(seeds))
    for s, n in sources:
        print("    %-58s %3d" % (s, n))
    if len(seeds) < 2:
        return 0
    med = stats.median(seeds)
    print("  median %+.4f   sd %.4f   bar |median| < 0.02 -> %s"
          % (med, stats.stdev(seeds), "WITHIN" if abs(med) < 0.02 else "OUTSIDE"))
    print("  NOTE: the decision statistic is a BOOTSTRAP CI on these per-seed values, not this")
    print("        point estimate and not a normal-theory interval -- the two disagreed about")
    print("        zero at n=50. Read wiki/fragment_null_270seed_prereg_2026-10-30.md section 5.")
    return 0


if __name__ == "__main__":
    cmds = {"status": cmd_status, "step": cmd_step, "pool": cmd_pool}
    if len(sys.argv) != 2 or sys.argv[1] not in cmds:
        print("usage: fragment_null_batch_driver.py {%s}" % "|".join(cmds), file=sys.stderr)
        raise SystemExit(2)
    raise SystemExit(cmds[sys.argv[1]]())
