"""Unattended runner that carries the 270-seed campaign to completion while nobody is watching.

WHY. Kaggle permits exactly TWO concurrent batch GPU sessions, so the four batches cannot all be in
flight at once and the campaign otherwise stalls at the cap until a human re-pushes. This fills any
free slot from a queue, every POLL_SECONDS, until the queue is empty and nothing is running.

IT OWNS NO STATE. Every decision is re-derived from Kaggle's own kernel statuses on each pass, and
every action it takes is one `fragment_null_batch_driver` already performs idempotently. So if this
process dies -- session ends, laptop sleeps, token lapses -- NOTHING is lost and nothing is
duplicated: `uv run python scripts/fragment_null_batch_driver.py step` in the morning resumes from
whatever actually happened. That property is the whole design, and it is why this is a convenience
rather than a load-bearing component.

THE QUEUE IS ORDERED SHORTEST-FIRST ON PURPOSE. The attested sampling sweep is ~15 min and
discharges a named blocker (the `--verify-knob` half of G-D's gate); the batches are ~4.5 h each.
Running the sweep in the first freed slot costs the batch schedule ~15 min of latency and buys a
result a day earlier.

WHAT IT WILL NOT DO. It never re-pushes a slug that already exists in any state, because Kaggle
keeps ONE output per kernel and re-pushing destroys retrieved artifacts. It never writes a verdict:
the pooled number and its pre-registered branch are written to a report for a human to read against
`wiki/fragment_null_270seed_prereg_2026-10-30.md`.

    uv run python scripts/fragment_null_overnight.py
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "fragment_null_batch_driver.py"
SWEEP_KERNEL = ROOT / "scripts" / "kaggle" / "glm_sampling_sweep_kernel.py"
SWEEP_SLUG = "glm-sampling-sweep-vk"

#: Family G-D's LoRA fine-tune, authorised by the user 2026-10-31. Queued here rather than pushed
#: directly because both GPU sessions were occupied by the funded null campaign at authorisation
#: time, and pre-empting that campaign to start this would spend the user's own decision against
#: itself. It sits LAST in the queue so any legitimately-ABSENT batch always outranks it.
FINETUNE_KERNEL = ROOT / "scripts" / "kaggle" / "glm_finetune_kernel.py"
FINETUNE_SLUG = "glm-generator-finetune"
#: The FAILURE file is in this list on purpose. The real kernel exits 0 and writes it when the run
#: dies, because Kaggle returns no output at all for an ERROR kernel -- so without it a failed
#: fine-tune is unreadable. It is a DISTINCT filename from the verdict and is never a verdict.
FINETUNE_ARTIFACTS = ["glm_generator_finetune_verdict.json",
                      "arm_a_zeroshot.json", "arm_b_finetuned.json",
                      "glm_finetune_failure.json"]

POLL_SECONDS = 600
MAX_HOURS = 11.5          # a stop, not a schedule: the campaign is ~9 h. This used to be
                          # justified by a ~12 h token lapse; it is not, because the client
                          # is now rebuilt every pass. It is an ATTENTION bound only.


def _driver():
    spec = importlib.util.spec_from_file_location("_fnbd", DRIVER)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


#: Resolved ONCE. An earlier draft re-exec'd the kernel module on every log line, which would have
#: re-run its assertions a few hundred times a night to learn a path that cannot change.
_D = None


def _d():
    global _D
    if _D is None:
        _D = _driver()
    return _D


def log(msg: str) -> None:
    line = "%s  %s" % (datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"), msg)
    print(line, flush=True)
    stage = _d().STAGE
    stage.mkdir(parents=True, exist_ok=True)
    with open(stage / "overnight_log.txt", "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def push(kernel: Path, slug: str) -> int:
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "kaggle_push_poll.py"),
                        "push", str(kernel), slug, "--gpu"],
                       cwd=str(ROOT), capture_output=True, text=True)
    out = ((r.stdout or "") + (r.stderr or "")).strip().splitlines()
    log("    push %s -> rc=%d  %s" % (slug, r.returncode, out[-1] if out else ""))
    return r.returncode


def report(d, statuses: dict) -> None:
    """Rewritten every pass, so it is useful even if this process is killed mid-campaign."""
    d.STAGE.mkdir(parents=True, exist_ok=True)
    pool = subprocess.run([sys.executable, str(DRIVER), "pool"],
                          cwd=str(ROOT), capture_output=True, text=True)
    body = [
        "# Fragment-null 270-seed campaign — overnight report",
        "",
        "<!-- doc-status: LIVE while the campaign runs; rewritten on every poll. The VERDICT READING",
        "     is not here -- it is pre-registered in wiki/fragment_null_270seed_prereg_2026-10-30.md,",
        "     sections 3a and 5. Read those before quoting any number below. -->",
        "",
        "Generated %s by `scripts/fragment_null_overnight.py`." %
        datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
        "",
        "## Kernel state",
        "",
        "| job | seeds | status | pulled |",
        "|---|---|---|---|",
    ]
    for b, (start, n) in sorted(d._batches().items()):
        body.append("| `%s` | %d..%d | %s | %s |"
                    % (d.SLUG % b, start, start + n - 1, statuses.get(d.SLUG % b, "?"),
                       "yes" if d._pulled(b) else "no"))
    body.append("| `%s` | — (attested sampling sweep) | %s | — |"
                % (SWEEP_SLUG, statuses.get(SWEEP_SLUG, "?")))
    body.append("| `%s` | — (G-D LoRA fine-tune) | %s | %s |"
                % (FINETUNE_SLUG, statuses.get(FINETUNE_SLUG, "?"),
                   "yes" if (d.STAGE / FINETUNE_ARTIFACTS[0]).is_file() else "no"))
    body += ["", "## Pooled null", "", "```", (pool.stdout or "").rstrip() or "(no output)", "```",
             "",
             "## The fine-tune (G-D), queued 2026-10-31",
             "",
             "Authorised by the user. It sits LAST in the launch queue, so it takes the first slot",
             "the null campaign frees rather than pre-empting it. Its bar is FROZEN in",
             "`wiki/glm_generator_finetune_acceptance_bar.json` with all four outcomes named --",
             "read that BEFORE reading the verdict, and do not re-specify the 0.6679 threshold,",
             "which is the Markov null's own measured value.",
             "",
             "## If this stalled",
             "",
             "Nothing is lost. Re-run `uv run python scripts/fragment_null_batch_driver.py step`;",
             "it re-derives everything from Kaggle and is safe to repeat. A kernel reading",
             "`API_UNHEALTHY` means the ~12 h Kaggle OAuth lapsed — re-auth in your own terminal",
             "with `uv run kaggle auth login --force`, then step again.",
             ]
    (d.STAGE / "MORNING_REPORT.md").write_text("\n".join(body) + "\n", encoding="utf-8")


def main() -> int:
    d = _d()
    if not d._healthy_or_explain(d._api()):
        return 2

    batches = d._batches()
    t0 = time.time()
    log("overnight runner started; cap %d concurrent GPU, poll %ds, stop after %.1f h"
        % (d.MAX_CONCURRENT_GPU, POLL_SECONDS, MAX_HOURS))

    while True:
        if (time.time() - t0) / 3600.0 > MAX_HOURS:
            log("STOP: %.1f h elapsed. Resume with `fragment_null_batch_driver.py step`." % MAX_HOURS)
            return 0
        # A FRESH client EVERY PASS, and that is the fix for a measured 4-hour outage rather
        # than defensive padding. 2026-10-31: this runner authenticated once at start, went
        # permanently unhealthy ~20 min in, and then logged
        # "API UNHEALTHY (token likely lapsed ~12 h)" on 22 consecutive polls while b3 sat
        # failed and b4 ran -- pulling nothing and relaunching nothing. A FRESH process
        # authenticated fine throughout that window, so the ON-DISK credential was valid the
        # whole time and the HELD CLIENT was the fault: the Kaggle OAuth ACCESS token is far
        # shorter-lived than the ~12 h refresh credential, and nothing here refreshed it. The
        # old message named the wrong cause and sent the operator to re-auth a credential that
        # was not broken. Re-authenticating costs one call per 10 min.
        try:
            api = d._api()
        except Exception as e:
            log("could not build a Kaggle client (%s: %s). Holding; retrying next poll."
                % (type(e).__name__, e))
            time.sleep(POLL_SECONDS)
            continue
        if not d._api_healthy(api):
            log("API UNHEALTHY with a FRESH client, so this is the ON-DISK credential and not a "
                "stale in-process token. Holding; will retry next poll. Fix it in YOUR OWN "
                "terminal (credentials must not route through a session): "
                "`uv run kaggle auth login --force`.")
            time.sleep(POLL_SECONDS)
            continue

        statuses = {d.SLUG % b: d._status(api, d.SLUG % b, True) for b in sorted(batches)}
        statuses[SWEEP_SLUG] = d._status(api, SWEEP_SLUG, True)
        statuses[FINETUNE_SLUG] = d._status(api, FINETUNE_SLUG, True)
        # QUEUED occupies a slot just as RUNNING does -- see driver.OCCUPYING.
        running = [k for k, v in statuses.items() if v in d.OCCUPYING]
        log("statuses: " + "  ".join("%s=%s" % (k.rsplit("-", 1)[-1], v)
                                     for k, v in statuses.items()))

        # The fine-tune's verdict, when it lands. Pulled like a batch artifact so a human reads a
        # file rather than a Kaggle log.
        _ft_done = any((d.STAGE / n).is_file() for n in (FINETUNE_ARTIFACTS[0],
                                                          FINETUNE_ARTIFACTS[-1]))
        if statuses[FINETUNE_SLUG] == "COMPLETE" and not _ft_done:
            log("  pulling the fine-tune verdict ...")
            got = d._pull(api, FINETUNE_SLUG, FINETUNE_ARTIFACTS, d.STAGE)
            log("    %d file(s): %s" % (len(got), ", ".join(Path(g).name for g in got) or "NONE"))

        # Pull anything finished and not yet on disk.
        for b in sorted(batches):
            if statuses[d.SLUG % b] == "COMPLETE" and not d._pulled(b):
                log("  pulling b%d ..." % b)
                got = d._pull(api, d.SLUG % b, d._artifacts(b), d.STAGE)
                log("    %d file(s): %s" % (len(got), ", ".join(Path(g).name for g in got) or "NONE"))

        # Queue, shortest job first. ABSENT is the only launchable state -- anything else exists
        # already, and re-pushing an existing slug destroys its output.
        queue = []
        if statuses[SWEEP_SLUG] == d.LAUNCHABLE:
            queue.append(("sweep", SWEEP_KERNEL, SWEEP_SLUG, None))
        for b in sorted(batches):
            if statuses[d.SLUG % b] == d.LAUNCHABLE:
                queue.append(("batch", d.KERNEL, d.SLUG % b, b))
        # LAST on purpose: the null campaign is the funded work and outranks this.
        if statuses[FINETUNE_SLUG] == d.LAUNCHABLE:
            queue.append(("finetune", FINETUNE_KERNEL, FINETUNE_SLUG, None))

        slots = d.MAX_CONCURRENT_GPU - len(running)
        for kind, kernel, slug, b in queue[:max(0, slots)]:
            if kind == "batch":
                import re
                src = kernel.read_text(encoding="utf-8")
                if int(re.search(r"^BATCH = (\d+)$", src, re.M).group(1)) != b:
                    with open(kernel, "w", encoding="utf-8", newline="\r\n") as fh:
                        fh.write(re.sub(r"^BATCH = \d+$", "BATCH = %d" % b, src, count=1, flags=re.M))
                    log("  set BATCH = %d" % b)
            log("  launching %s (%d slot(s) free)" % (slug, slots))
            push(kernel, slug)
            slots -= 1

        report(d, statuses)

        unfinished = [k for k, v in statuses.items()
                      if v in d.OCCUPYING or v == d.LAUNCHABLE]
        if not unfinished:
            log("CAMPAIGN COMPLETE: nothing running, nothing queued.")
            report(d, statuses)
            log("final pooled result is in %s" % (d.STAGE / "MORNING_REPORT.md"))
            return 0
        time.sleep(POLL_SECONDS)


if __name__ == "__main__":
    raise SystemExit(main())
