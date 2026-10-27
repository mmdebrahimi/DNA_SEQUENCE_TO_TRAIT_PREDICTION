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
PINNED_SHA = "74d8af7b4f6d4f9ba9d0d6931277a860f5f472bd"

WORK = Path("/kaggle/working")
SRC = WORK / "dna_decode_repo"
CACHE = WORK / "mpra"

SMOKE_LIMIT = 4000
SMOKE_SEEDS = 3
FULL_SEEDS = 10


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


def main():
    print("=" * 90, flush=True)
    device = device_report()
    print("=" * 90, flush=True)
    head = clone()
    fetch_data()

    print("\n" + "=" * 90, flush=True)
    print("SMOKE: --limit %d --seeds %d" % (SMOKE_LIMIT, SMOKE_SEEDS), flush=True)
    print("=" * 90, flush=True)
    smoke, smoke_min = run_arm("smoke", ["--limit", str(SMOKE_LIMIT), "--seeds", str(SMOKE_SEEDS)],
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

    # Size the full run from the smoke rather than from an estimate.
    per_fit = smoke_min / max(1, 2 * SMOKE_SEEDS)          # real fits + null fits
    scale = 321123 / SMOKE_LIMIT
    projected_h = per_fit * scale * (2 * FULL_SEEDS) / 3600.0
    print("\nPROJECTED full run: ~%.1f h (from %.1f s per fit at N=%d, scaled %.0fx, %d fits). "
          "Kaggle GPU sessions cap around 9 h." % (projected_h, per_fit, SMOKE_LIMIT, scale,
                                                   2 * FULL_SEEDS), flush=True)
    if projected_h > 8.0:
        print("!! Projected run does not comfortably fit one session. STOPPING after the smoke. "
              "Re-run with fewer seeds, or split the seed range across sessions, rather than "
              "losing a near-complete run to a timeout.", flush=True)
        json.dump({"head": head, "smoke": smoke, "projected_hours": projected_h,
                   "decision": "STOPPED_AFTER_SMOKE_PROJECTION_TOO_LONG"},
                  open(WORK / "fragment_pyramid_kaggle_summary.json", "w"), indent=2)
        return 0

    print("\n" + "=" * 90, flush=True)
    print("FULL: --seeds %d on all 321,123 fragments" % FULL_SEEDS, flush=True)
    print("=" * 90, flush=True)
    full, full_min = run_arm("full", ["--seeds", str(FULL_SEEDS)], device)

    json.dump({"head": head, "device": device, "smoke_minutes": round(smoke_min / 60, 2),
               "full_minutes": round(full_min / 60, 2),
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
