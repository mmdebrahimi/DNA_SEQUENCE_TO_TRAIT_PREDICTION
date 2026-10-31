#!/usr/bin/env python
"""Kaggle T4 kernel: the 270-SEED NULL RESOLUTION for the fragment-pyramid arm, in batches.

WHY THIS EXISTS. The arm's validity gate is `|null median| < 0.02`, and it has voided every verdict
so far. Pooled over the 50 seeds drawn to date (0..9 twice-published + 10..49 fresh) the null median
is **+0.0242, sd 0.0654**, with a bootstrap CI of `[-0.0126, +0.0466]` that CONTAINS ZERO -- so the
bar is failing on a POINT ESTIMATE whose own standard error (~0.012) is comparable to the bar
itself. A point-estimate bar narrower than its own SE cannot separate a clean null from a dirty one.

THE USER'S CALL (2026-10-30), taken over two cheaper alternatives: resolve it by SPENDING SEEDS
rather than by re-specifying the bar or accepting no verdict. Target n = 270, which is where
SE(median) ~= 1.2533 * sd / sqrt(n) reaches 0.005 -- a quarter of the bar, enough that the bar
becomes decidable in either direction.

    1.2533 * 0.0654 / sqrt(270) = 0.0050

WHAT WAS PRE-REGISTERED BEFORE THIS RAN. `wiki/fragment_null_270seed_prereg_2026-10-30.md` fixes the
reading of all three possible outcomes BEFORE any of these seeds existed, because the most likely
outcome (CI excluding zero AND above 0.02) CONFIRMS a small real positive bias rather than clearing
the gate. Writing that down first is what stops the result being retro-fitted.

WHY BATCHES. At the MEASURED 4.945 min/seed (197.79 min for 40 seeds on a T4) 220 further seeds is
~18.1 h, against a Kaggle GPU session cap of ~9 h and this kernel's own 8 h projection guard. Four
batches of 55 project to ~4.96 h and land near ~4.5 h actual -- deliberately generous margin,
because the guard exists precisely to avoid losing a near-complete run to a timeout.

POOLING VALIDITY. `PINNED_SHA` is UNCHANGED from the seeds-10..49 run and matches the `head` field
of `wiki/fragment_null_ext_summary_2026-10-29.json` exactly, so all 270 seeds are drawn from
byte-identical code. A bumped pin would silently make the pooled median a mixture of two estimators.

EACH BATCH NEEDS ITS OWN SLUG. Kaggle keeps one output per kernel, so re-pushing an existing slug
DESTROYS its retrieved artifacts -- `glm-fragment-null-ext` in particular still holds the 10..49
result. Set BATCH below, then push to the matching slug:

    uv run python scripts/kaggle_push_poll.py push scripts/kaggle/fragment_null_batch_kernel.py         glm-fragment-null-b1 --gpu

STATUS AND VERDICT. Every batch reports a NON-PREREGISTERED seed range and therefore emits NO
verdict, by design: this resolves the NULL's precision and must not re-adjudicate the
pyramid-vs-kmer4 comparison, whose pre-registered range is seeds 0..9 and whose answer is already
measured.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO_URL = "https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION.git"
#: Pinned so the artifact is traceable to a commit. Bump deliberately, never to "latest".
#: Bumped 2026-10-27 off 74d8af7b, which predates the `learned_contribution` diagnostic -- a run at the
#: old pin could not report how much of the margin over kmer4 the LABELS bought, and would also
#: lose a completed run to a summary bug (the artifact was written AFTER the print).
PINNED_SHA = "74b7a28df74d46974f36c5f55ed601b9cca3b156"

WORK = Path("/kaggle/working")
SRC = WORK / "dna_decode_repo"
CACHE = WORK / "mpra"

SMOKE_LIMIT = 4000
SMOKE_SEEDS = 3
FULL_SEEDS = 10

#: WHICH SLICE THIS PUSH IS. Edit this one line per batch; everything else stays fixed so the
#: batches are byte-comparable apart from the seed range.
BATCH = 1

#: The full 220-seed plan, written out so the file records what it is a part of. Seeds 0..49 are
#: already measured and are NOT redrawn: 0..9 in the two published full runs, 10..49 in
#: wiki/fragment_null_ext_summary_2026-10-29.json. 50..269 are new, giving a pooled n = 270.
BATCHES = {1: (50, 55), 2: (105, 55), 3: (160, 55), 4: (215, 55)}
assert BATCH in BATCHES, "BATCH must be one of %s" % sorted(BATCHES)
assert sum(n for _, n in BATCHES.values()) == 220
assert sorted(s for s, _ in BATCHES.values()) == [50, 105, 160, 215], "batch starts must tile"

SEED_START, EXT_SEEDS = BATCHES[BATCH]

#: MEASURED, not projected -- the 10-seed full run's own wall-clock on a Kaggle T4
#: (wiki/fragment_pyramid_kaggle_summary_2026-10-28.json: full_minutes 54.08).
MEASURED_FULL10_MINUTES = 54.08

#: A SECOND, SMALLER timing point, added 2026-10-28. ONE timing point cannot separate fixed
#: cost from marginal cost, and the projection multiplies whatever it is handed by ~268x
#: (scale 80.3 x 20 fits / 6 fits) -- so ~30 s of one-off setup (import torch, build the model,
#: parse the frames) was projected as +2.2 h of the full run, and the projection STOPPED a run
#: on its own overhead. Two points give a slope, and the slope is the only part that scales.
CAL_LIMIT = 1000
CAL_SEEDS = 1


def sh(cmd, cwd=None, check=True):
    print("$ %s" % (cmd if isinstance(cmd, str) else " ".join(cmd)), flush=True)
    r = subprocess.run(cmd, cwd=cwd, shell=isinstance(cmd, str),
                       capture_output=True, text=True)
    if r.stdout:
        print(r.stdout[-6000:], flush=True)
    if r.stderr:
        print("STDERR:", r.stderr[-6000:], flush=True)
    if check and r.returncode != 0:
        raise SystemExit("command failed with %d" % r.returncode)
    return r


def device_report():
    """Report the real device. A P100 is correct-but-slow here; silence would hide which we got."""
    try:
        import torch
    except ImportError:
        print("!! torch not importable -- this kernel needs the GPU image", flush=True)
        return "cpu"
    print("torch %s  cuda_available=%s" % (torch.__version__, torch.cuda.is_available()),
          flush=True)
    if not torch.cuda.is_available():
        print("!! NO GPU. Set Accelerator to GPU T4 x1. Continuing on CPU would take tens of "
              "hours and is not worth a session.", flush=True)
        return "cpu"
    cc = torch.cuda.get_device_capability(0)
    name = torch.cuda.get_device_name(0)
    print("device: %s  compute capability %d.%d" % (name, cc[0], cc[1]), flush=True)
    if cc[0] < 7:
        print("  !! CC < 7.0 -- the P100 trap: `enable_gpu` alone provisions a P100 (sm_60). "
              "machine_shape must pin NvidiaTeslaT4. This run is fp32 so the result is still "
              "CORRECT, just slower.", flush=True)
    return "cuda"


def clone():
    if SRC.exists():
        print("repo already present", flush=True)
    else:
        sh(["git", "clone", "--quiet", REPO_URL, str(SRC)])
    sh(["git", "checkout", "--quiet", PINNED_SHA], cwd=str(SRC))
    head = sh(["git", "rev-parse", "HEAD"], cwd=str(SRC)).stdout.strip()
    assert head == PINNED_SHA, "checked out %s, expected %s" % (head, PINNED_SHA)
    print("HEAD pinned at %s" % head, flush=True)
    return head


def verify_pinned_tree():
    """Check the tree just checked out carries the fixes this kernel assumes. Costs ~10 ms.

    WHY, measured 2026-10-27: a run of this kernel pinned at 74d8af7b spent 5.1 minutes computing
    real fits and then DISCARDED them -- that arm reads res["headline"] unconditionally, which
    KeyErrors on the PARTIAL_SEED_RANGE result a 3-of-10-seed smoke returns BY DESIGN, and it wrote
    its artifact AFTER the summary print, so the crash took the completed work with it. Nothing
    noticed until the traceback. The pin decides which code runs, so a stale pin must fail in
    seconds, not minutes, and must say WHY rather than surfacing as a KeyError.

    HONEST LIMIT: this proves the fixes are PRESENT in the checked-out source, not that they are
    correct. It is a staleness guard, not a test.
    """
    arm = SRC / "scripts" / "glm_fragment_pyramid_arm.py"
    if not arm.is_file():
        raise SystemExit("REFUSING: %s is absent from the pinned tree -- wrong PINNED_SHA?" % arm)
    s = arm.read_text(encoding="utf-8")

    missing = []
    if "learned_contribution" not in s:
        missing.append(
            "the `learned_contribution` diagnostic (real minus null, paired per seed). Without it "
            "the run cannot say how much of the margin over kmer4 the LABELS bought, which is the "
            "question this run exists to answer.")
    w, pr = s.find("Path(a.out).write_text"), s.find('print("arm: %s" % ARM)')
    if w == -1 or pr == -1 or w > pr:
        missing.append(
            "the write-before-format ordering. At this pin the artifact is written AFTER the "
            "summary print, so a formatting bug discards a completed run -- measured: 5.1 minutes "
            "of fits lost on 2026-10-27.")
    if 'h = res.get("headline")' not in s:
        missing.append(
            "the partial-aware printer. At this pin the summary reads res['headline'] "
            "unconditionally and CRASHES on PARTIAL_SEED_RANGE, which is exactly what a 3-of-10 "
            "seed smoke returns by design.")

    if missing:
        print("REFUSING: the pinned tree is STALE.\n"
              "  PINNED_SHA = %s\n"
              "  Missing from scripts/glm_fragment_pyramid_arm.py at that commit:" % PINNED_SHA,
              flush=True)
        for m in missing:
            print("    * %s" % m, flush=True)
        print("  This almost always means the NOTEBOOK holds an older copy of this kernel than the\n"
              "  repo does. Re-fetch and re-paste the cell:\n"
              "    https://raw.githubusercontent.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION"
              "/main/scripts/kaggle/fragment_pyramid_kernel.py\n"
              "  Refusing now rather than after the 91 MB download and ~5 minutes of fits.",
              flush=True)
        raise SystemExit(3)
    print("pinned tree verified: diagnostic present, write-before-format, partial-aware printer",
          flush=True)


def fetch_data():
    """Download GSE144621 through the PROJECT'S loader, not curl.

    The loader sniffs the delimiter per file (the fragment files are space-delimited while the
    peak-tile file in the same series is tab-delimited, and a global split parses one of them to
    zero rows, which reads as 'no signal' rather than 'no data'). A hand-rolled fetch of this
    project's data wrote two HTML error pages under genuine data filenames earlier today.
    """
    sys.path.insert(0, str(SRC))
    from dna_decode.glm.genomewide import ensure_cached, load_fragments

    CACHE.mkdir(parents=True, exist_ok=True)
    for which in ("LB", "M9"):
        p = ensure_cached(CACHE, which=which, allow_download=True)
        print("cached %s -> %s (%.1f MB)" % (which, p, p.stat().st_size / 1e6), flush=True)

    frags, st = load_fragments(CACHE, which="LB", allow_download=False)
    print("LB parse: n=%d  len %d..%d  distinct_lengths=%d"
          % (len(frags), st.length_min, st.length_max, st.n_distinct_lengths), flush=True)
    # The expected shape, so a wrong or truncated file is caught HERE rather than in a verdict.
    assert (len(frags), st.length_min, st.length_max) == (321123, 48, 499), (
        "expected 321123 fragments at 48..499 bp, got %d at %d..%d -- the data is wrong, refusing "
        "to run the arm on it" % (len(frags), st.length_min, st.length_max))
    print("substrate verified", flush=True)


def run_arm(tag, extra, device):
    out = WORK / ("glm_fragment_pyramid_arm_kaggle_b%d_%s.json" % (BATCH, tag))
    cmd = [sys.executable, "scripts/glm_fragment_pyramid_arm.py",
           "--cache-dir", str(CACHE), "--device", device,
           "--with-encoder", "--out", str(out)] + extra
    t0 = time.time()
    r = sh(cmd, cwd=str(SRC), check=False)
    dt = time.time() - t0
    print("\n[%s] wall-clock %.1f min  rc=%d" % (tag, dt / 60.0, r.returncode), flush=True)
    if not out.is_file():
        print("[%s] NO ARTIFACT written" % tag, flush=True)
        return None, dt
    d = json.loads(out.read_text())
    print("[%s] status=%s verdict=%s" % (tag, d.get("status"), d.get("verdict")), flush=True)
    if d.get("null"):
        print("[%s] null per-seed=%s median=%s sd=%s"
              % (tag, d["null"].get("per_seed"), d["null"].get("median"), d["null"].get("sd")),
              flush=True)
    if d.get("paired"):
        print("[%s] paired=%s" % (tag, d["paired"]), flush=True)
    return d, dt


def on_kaggle() -> bool:
    """Is this actually a Kaggle notebook? Checked because the failure mode is SILENT.

    `Path("/kaggle/working")` is an ABSOLUTE POSIX path, and on Windows an absolute path with no
    drive resolves against the current drive -- so off-platform it becomes `C:\\kaggle\\working` and
    every write SUCCEEDS. Measured 2026-10-27: run on the project laptop this kernel created
    `C:\\kaggle\\working`, cloned the repo into it, downloaded 91 MB of GEO data and printed
    `substrate verified`, looking exactly like a Kaggle run until the arm died on a missing `scipy`.

    It must refuse rather than half-run: a kernel that appears to work off-platform invites the
    conclusion that it HAS been tested on Kaggle when it has not. (The 91 MB download was in fact
    useful that day, but usefulness by accident is not a reason to leave the trap open.)
    """
    return bool(os.environ.get("KAGGLE_KERNEL_RUN_TYPE")) or Path("/kaggle/input").exists()


def main():
    if not on_kaggle() and "--allow-non-kaggle" not in sys.argv:
        print("REFUSING: this is a Kaggle kernel and this host is not Kaggle.\n"
              "  `/kaggle/working` is an absolute POSIX path; off-platform on Windows it silently\n"
              "  resolves against the current drive (C:\\kaggle\\working) and every write succeeds,\n"
              "  so the kernel half-runs and LOOKS like a Kaggle job.\n"
              "  Two further reasons a local run cannot be trusted even when it appears to work:\n"
              "    * it shells out via sys.executable, which here is whatever python invoked it --\n"
              "      not necessarily the project venv that holds torch/scipy;\n"
              "    * there is no GPU contract off-platform, so a 'successful' run may be CPU-only\n"
              "      and silently hours slower.\n"
              "  Run it in a Kaggle notebook (Accelerator GPU T4 x1, Internet ON). To exercise the\n"
              "  arm locally instead, call it directly with the project venv:\n"
              "    uv run python scripts/glm_fragment_pyramid_arm.py --cache-dir <dir> \\\n"
              "        --limit 4000 --seeds 3 --with-encoder\n"
              "  Pass --allow-non-kaggle to override deliberately.", file=sys.stderr)
        return 2

    print("=" * 90, flush=True)
    device = device_report()
    print("=" * 90, flush=True)

    # A WARNING THE CODE THEN OVERRIDES IS NOT A WARNING. device_report() says continuing on CPU
    # "is not worth a session" -- so do not continue. Measured 2026-10-27: a CPU session got as far
    # as a repo clone, a 91 MB download and 5.1 minutes of fits before dying, on a configuration
    # that could never have finished the real run.
    if device != "cuda" and "--allow-cpu" not in sys.argv:
        print("REFUSING: no CUDA device, so this session cannot finish the real run.\n"
              "  Fix it in the notebook: Settings -> Accelerator -> `GPU T4 x1`.\n"
              "  Plain `GPU` is not enough -- it can provision a P100 (sm_60), which is slow but\n"
              "  CORRECT here (this run is fp32 throughout); no accelerator at all is neither.\n"
              "  Refusing BEFORE the clone and the 91 MB GEO download rather than after.\n"
              "  Pass --allow-cpu to exercise the plumbing on CPU deliberately (the smoke runs in\n"
              "  minutes; the full run does not fit a session).", file=sys.stderr)
        return 3

    head = clone()
    verify_pinned_tree()   # seconds, and BEFORE the download -- a stale pin must not cost 5 minutes
    fetch_data()

    print()
    print("=" * 90, flush=True)
    print("SMOKE (plumbing only): --limit %d --seeds %d --seed-start %d"
          % (SMOKE_LIMIT, SMOKE_SEEDS, SEED_START), flush=True)
    print("=" * 90, flush=True)
    smoke, smoke_sec = run_arm("extsmoke",
                               ["--limit", str(SMOKE_LIMIT), "--seeds", str(SMOKE_SEEDS),
                                "--seed-start", str(SEED_START)], device)

    # A healthy smoke here is a PARTIAL/non-preregistered status carrying paired stats + a null. It
    # is NOT a result: at N=4,000 the null's sd is about the size of any delta. It exists only to
    # prove the `--seed-start` path runs, which is the one thing never exercised anywhere.
    ok = smoke is not None and smoke.get("paired") and smoke.get("null")
    if not ok:
        print("\nSMOKE DID NOT REACH A PAIRED RESULT -- stopping before the long run, which is "
              "exactly what the smoke is for. status=%s"
              % (None if smoke is None else smoke.get("status")), flush=True)
        return 1
    print("\nsmoke reached paired stats + null (status=%s, EXPECTED for a shifted seed range)"
          % smoke.get("status"), flush=True)

    # CONFIRM THE SHIFT LANDED rather than trusting the flag. A --seed-start silently ignored would
    # redraw seeds 0..2 -- already measured twice -- and the run would look entirely healthy while
    # adding no information at all. That is the whole premise of this kernel, so it is checked.
    sr = (smoke or {}).get("seed_range") or {}
    if sr.get("start") != SEED_START:
        print("\nREFUSING: the arm reported seed_range start=%r, expected %d. The --seed-start "
              "flag did not take effect, so a long run would silently redraw the SAME seeds that "
              "are already measured and resolve nothing."
              % (sr.get("start"), SEED_START), flush=True)
        return 4
    print("seed shift CONFIRMED: seed_range=%s" % (sr,), flush=True)

    projected_h = MEASURED_FULL10_MINUTES * (EXT_SEEDS / 10.0) / 60.0
    print()
    print("PROJECTED extension: ~%.2f h  (measured anchor %.2f min for 10 seeds at full N, "
          "scaled x%.1f)" % (projected_h, MEASURED_FULL10_MINUTES, EXT_SEEDS / 10.0), flush=True)
    print("   Kaggle GPU sessions cap around 9 h.", flush=True)
    if projected_h > 8.0:
        print("!! Does not comfortably fit one session -- STOPPING rather than losing a "
              "near-complete run to a timeout. Split the seed range across sessions.", flush=True)
        json.dump({"head": head, "decision": "STOPPED_PROJECTION_TOO_LONG",
                   "projected_hours": projected_h},
                  open(WORK / ("fragment_null_b%d_summary.json" % BATCH), "w"), indent=2)
        return 0

    print("\n" + "=" * 90, flush=True)
    print("EXTENSION: --seeds %d --seed-start %d on all 321,123 fragments (seeds %d..%d)"
          % (EXT_SEEDS, SEED_START, SEED_START, SEED_START + EXT_SEEDS - 1), flush=True)
    print("=" * 90, flush=True)
    ext, ext_sec = run_arm("ext", ["--seeds", str(EXT_SEEDS),
                                   "--seed-start", str(SEED_START)], device)

    json.dump({"head": head, "device": device,
               "seed_start": SEED_START, "ext_seeds": EXT_SEEDS,
               "smoke_minutes": round(smoke_sec / 60, 2),
               "projected_hours": projected_h,
               "projection_basis": "measured 10-seed full run %.2f min, scaled x%.1f"
                                   % (MEASURED_FULL10_MINUTES, EXT_SEEDS / 10.0),
               "ext_minutes": None if ext is None else round(ext_sec / 60, 2),
               "ext_status": None if ext is None else ext.get("status"),
               "ext_verdict": None if ext is None else ext.get("verdict"),
               "ext_seed_range": None if ext is None else ext.get("seed_range"),
               "ext_null": None if ext is None else ext.get("null"),
               "ext_paired": None if ext is None else ext.get("paired"),
               "batch": BATCH, "batch_plan": {str(k): v for k, v in BATCHES.items()},
               "pooled_target_n": 270,
               "prereg": "wiki/fragment_null_270seed_prereg_2026-10-30.md",
               "purpose": "batch %d of 4 of the 270-seed NULL resolution; NOT a verdict run"
                          % BATCH},
              open(WORK / ("fragment_null_b%d_summary.json" % BATCH), "w"), indent=2)

    if ext is None:
        print("\nEXTENSION RUN PRODUCED NO ARTIFACT.", flush=True)
        return 1
    print("\nDONE. Artifacts in /kaggle/working: "
          "glm_fragment_pyramid_arm_kaggle_{extsmoke,ext}.json + "
          "fragment_null_ext_summary.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
