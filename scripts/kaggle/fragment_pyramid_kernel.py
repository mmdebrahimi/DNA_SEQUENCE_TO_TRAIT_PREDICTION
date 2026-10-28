#!/usr/bin/env python
"""Kaggle T4 kernel: the fragment-pyramid encoder arm, smoke then full run.

WHY KAGGLE IS THE RIGHT HOME FOR THIS ONE, and why it needs no payload machinery.

The arm trains a 4-scale masked-pooling conv encoder over 321,123 variable-length E. coli MPRA
fragments x 10 seeds, plus a 10-seed shuffled-label null. Measured on the project laptop's CPU a
single fit at N=4,000 takes ~3 minutes, which extrapolates to tens of hours at full N. The laptop's
GPU is a GTX 860M (Maxwell, CC 5.0, 4 GiB) and cannot help. A Kaggle T4 is Turing (CC 7.5, 16 GiB).

**The decisive reason is DATA, not speed.** The substrate is public GEO GSE144621, which the
project's own loader downloads itself -- so this kernel needs nothing uploaded and is completely
independent of the external D: volume, which detached mid-session on the host and is what blocks the
local run. Nothing here depends on that drive, or on any other machine answering.

**It CLONES the repo at a pinned SHA rather than embedding a base64 payload.** The existing
`build_glm_sampling_sweep_kernel.py` embeds sources and says why: *"the repo is private, so the
kernel cannot `git clone`"*. **That premise is no longer true** -- measured 2026-10-27, the GitHub
API reports `private: False` for this repo. Its second reason, keeping the kernel version-exact
against a working tree newer than the published wheel, is real and is satisfied here by pinning the
COMMIT: a clone at a fixed SHA is exactly as reproducible as an embedded tarball and far less
machinery. The kernel PRINTS the SHA it checked out, so any artifact it writes is traceable to a
commit rather than to "whatever main was that day".

SETTINGS THAT MATTER (set these in the notebook, they are not inferable from the code):
  * Accelerator: **GPU T4 x1**. `enable_gpu` alone can provision a P100 (sm_60); the device report
    below warns if CC < 7.0. This run is fp32 either way, so a P100 is slow-but-correct, not wrong.
  * Internet: **ON**. Required for the GitHub clone and the GEO download. A Kaggle account needs
    phone verification before internet can be enabled.
  * Free tier -> no money gate.

ORDER IS SMOKE THEN FULL, deliberately. The arm's own pre-registered validity null caught two
defects in it on the day it was written, and the CORRECTED multi-seed null has never completed a
scored run anywhere -- the host's drive detached first. So this runs a ~10 minute smoke, checks it
produced a verdict at all, and only then starts the long run. Spending hours of GPU on a code path
that has never executed once is the thing being avoided.
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
    out = WORK / ("glm_fragment_pyramid_arm_kaggle_%s.json" % tag)
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
    print("CALIBRATION: --limit %d --seeds %d (the second timing point, for the slope)"
          % (CAL_LIMIT, CAL_SEEDS), flush=True)
    print("=" * 90, flush=True)
    cal, cal_sec = run_arm("cal", ["--limit", str(CAL_LIMIT), "--seeds", str(CAL_SEEDS)],
                           device)
    print("\n" + "=" * 90, flush=True)
    print("SMOKE: --limit %d --seeds %d" % (SMOKE_LIMIT, SMOKE_SEEDS), flush=True)
    print("=" * 90, flush=True)
    smoke, smoke_sec = run_arm("smoke", ["--limit", str(SMOKE_LIMIT), "--seeds", str(SMOKE_SEEDS)],
                               device)

    # A HEALTHY SMOKE IS `PARTIAL_SEED_RANGE`, not `SCORED`. The smoke runs 3 of the 10
    # pre-registered seeds, and the arm deliberately refuses a verdict on a partial range -- judging
    # 3 seeds against a ">= 8 of 10" bar would redefine the acceptance criterion. So the smoke's
    # success criterion is "it got far enough to produce paired stats and a null", which is what
    # exercises every code path. An earlier version of this kernel demanded `SCORED` and would have
    # stopped after every smoke by construction.
    ok = smoke is not None and smoke.get("paired") and smoke.get("null")
    if not ok:
        print("\nSMOKE DID NOT REACH A PAIRED RESULT -- stopping before the long run. That is the "
              "point of the smoke: the corrected multi-seed null had never executed anywhere "
              "before this. status=%s" % (None if smoke is None else smoke.get("status")),
              flush=True)
        return 1
    print("\nsmoke reached paired stats + null (status=%s, which is EXPECTED for a 3-of-10 seed "
          "range)" % smoke.get("status"), flush=True)

    # The smoke's delta is NOT a result and must not be read as one.
    print("\nREAD THE SMOKE CORRECTLY: at N=%d the null's sd is about the same size as any delta "
          "(measured 0.1016 vs 0.1109 on the host), so the smoke tests that the code RUNS and "
          "nothing more. Per-seed nulls swinging to +/-0.13 are EXPECTED -- composition predicts "
          "this rung at spearman(GC, expression) = -0.1169, and a near-initialisation conv net is "
          "close to a random functional of composition." % SMOKE_LIMIT, flush=True)

    # Size the full run from TWO timing points, so fixed cost is not scaled up ~268x.
    #   t(N, S) ~= fixed + rate * N * (2*S)        [real fits + null fits]
    # Two measurements give both unknowns, and only `rate` multiplies up to the full run.
    w_cal = CAL_LIMIT * (2 * CAL_SEEDS)
    w_smoke = SMOKE_LIMIT * (2 * SMOKE_SEEDS)
    w_full = 321123 * (2 * FULL_SEEDS)
    rate = (smoke_sec - cal_sec) / float(w_smoke - w_cal)   # seconds per row-fit
    fixed = cal_sec - rate * w_cal

    naive_per_fit = smoke_sec / max(1, 2 * SMOKE_SEEDS)
    naive_h = naive_per_fit * (321123 / SMOKE_LIMIT) * (2 * FULL_SEEDS) / 3600.0
    inflation = (321123 / SMOKE_LIMIT) * (2 * FULL_SEEDS) / float(2 * SMOKE_SEEDS)

    if rate <= 0:
        # The smaller run came out no faster than the larger one, so the slope is below the
        # timing noise and no marginal rate is recoverable. REPORT that and fall back to the
        # (pessimistic) one-point number rather than projecting from a negative rate.
        projected_h = naive_h
        basis = ("DEGRADED one-point fallback: cal %.1f s was not faster than smoke %.1f s, so "
                 "the slope is below timing noise and no marginal rate is recoverable."
                 % (cal_sec, smoke_sec))
    else:
        projected_h = (fixed + rate * w_full) / 3600.0
        basis = ("two-point: %.1f s one-off + %.4f ms per row-fit over %d row-fits"
                 % (fixed, rate * 1000.0, w_full))

    print()
    print("PROJECTED full run: ~%.1f h  (%s)" % (projected_h, basis), flush=True)
    print("   one-point projection for comparison: ~%.1f h -- it charges the full run %.0fx "
          "for setup that is paid ONCE, which is the defect the second timing point removes."
          % (naive_h, inflation), flush=True)
    print("   Kaggle GPU sessions cap around 9 h.", flush=True)
    if projected_h > 8.0:
        print("!! Projected run does not comfortably fit one session. STOPPING after the smoke. "
              "Re-run with fewer seeds, or split the seed range across sessions, rather than "
              "losing a near-complete run to a timeout.", flush=True)
        json.dump({"head": head, "smoke": smoke, "projected_hours": projected_h,
                   "projection_basis": basis,
                   "projected_hours_one_point": naive_h,
                   "cal_seconds": cal_sec, "smoke_seconds": smoke_sec,
                   "decision": "STOPPED_AFTER_SMOKE_PROJECTION_TOO_LONG"},
                  open(WORK / "fragment_pyramid_kaggle_summary.json", "w"), indent=2)
        return 0

    print("\n" + "=" * 90, flush=True)
    print("FULL: --seeds %d on all 321,123 fragments" % FULL_SEEDS, flush=True)
    print("=" * 90, flush=True)
    full, full_sec = run_arm("full", ["--seeds", str(FULL_SEEDS)], device)

    json.dump({"head": head, "device": device, "smoke_minutes": round(smoke_sec / 60, 2),
               "cal_minutes": round(cal_sec / 60, 2),
               "projected_hours": projected_h, "projection_basis": basis,
               "projected_hours_one_point": naive_h,
               "full_minutes": round(full_sec / 60, 2),
               "smoke_status": smoke.get("status"),
               "full_status": None if full is None else full.get("status"),
               "full_verdict": None if full is None else full.get("verdict"),
               "full_headline": None if full is None else full.get("headline"),
               "full_null": None if full is None else full.get("null")},
              open(WORK / "fragment_pyramid_kaggle_summary.json", "w"), indent=2)

    if full is None:
        print("\nFULL RUN PRODUCED NO ARTIFACT.", flush=True)
        return 1
    print("\nDONE. Artifacts in /kaggle/working: "
          "glm_fragment_pyramid_arm_kaggle_{smoke,full}.json + "
          "fragment_pyramid_kaggle_summary.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
