#!/usr/bin/env python
"""Build the self-contained Kaggle T4 kernel for the GLM generator SAMPLING SWEEP.

WHY A BUILDER. The sweep needs `scripts/glm_generator_falsifier.py` plus a handful of `dna_decode`
modules, and the repo is private, so the kernel cannot `git clone` and the PyPI wheel does not ship
`scripts/`. Embedding an exact copy of the needed sources keeps the kernel VERSION-EXACT against this
working tree -- which matters here, because the fix that makes the GLM path generate at all
(`generate.py` calling `setup_tokenizer`) is newer than the published wheel.

Emits `scripts/kaggle/glm_sampling_sweep_kernel.py` (payload = base64 tar.gz, ~30 KB, far under the
~1 MB kernel source cap recorded in the project ledger).

    uv run python scripts/kaggle/build_glm_sampling_sweep_kernel.py
    uv run python scripts/kaggle_push_poll.py push \
        scripts/kaggle/glm_sampling_sweep_kernel.py glm-sampling-sweep --gpu

The payload is substituted by TOKEN REPLACEMENT, not str.format: the kernel body is full of f-strings
and escaping every brace for .format() is how this builder broke the first time.
"""
from __future__ import annotations

import base64
import io
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

# The exact import closure of `scripts/glm_generator_falsifier.py`. Kept VERBATIM -- no stubbing of
# __init__ files, so the kernel runs the same code this tree runs.
FILES = [
    "dna_decode/__init__.py",
    "dna_decode/data/__init__.py",
    "dna_decode/data/annotations.py",
    "dna_decode/data/refseq.py",
    "dna_decode/glm/__init__.py",
    "dna_decode/glm/baselines.py",
    "dna_decode/glm/corpus.py",
    "dna_decode/glm/edit.py",
    "dna_decode/glm/falsifier.py",
    "dna_decode/glm/generate.py",
    "scripts/glm_generator_falsifier.py",
]

KERNEL = r'''"""Kaggle T4 kernel: the GLM generator SAMPLING SWEEP (temperature x top_k).

WHY KAGGLE. The GLM is `GenerTeam/GENERator-v2-prokaryote-1.2b-base` (1.2B, fp32 ~4.8 GB). The local
host has no CUDA device at all (Intel UHD 620), so a sweep there runs on 4 CPU threads at ~14 s per
generated sequence. Kaggle's free T4 removes that wall with no money gate.

WHAT THE SWEEP IS FOR. Every recorded result for this generator was produced at DEFAULT sampling BY
OMISSION: `HFGenerator` always accepted `temperature`/`top_k`, but the CLI never passed them. This
kernel varies both and reports whether any setting beats the default.

LOWER distinguishability is BETTER (less distinguishable from real E. coli promoter sequence).

THE SETTINGS ARE RE-DERIVED FROM THE RECORDED ARTIFACT, NOT FROM PROSE. The recorded baseline
(`wiki/glm_generator_falsifier_REAL_2026-10-07.json`) was produced at:
    --n 120 --limit-natural 1200 --bins 6 --modes kmer,both --natural-set promoter --conditional
    --context 300 --length 150
Reproduced EXACTLY on the local CPU for the two baselines that need no GPU:
    kmer:  markov 0.5990  uniform 0.9408      both: markov 0.6679  uniform 0.9641
A sweep run at any other --n / --bins / --limit-natural is NOT comparable to the recorded
0.7237 (kmer) / 0.7706 (both), so those values are pinned below.

THE DEFAULT CELL IS THE SELF-CHECK. The grid includes (temperature=1.0, top_k=0), which IS the
recorded default. If that cell does not land near 0.7237/0.7706 the sweep is not measuring what it
claims and every other cell is suspect. It runs first and is flagged. It is NOT expected to match
bit-for-bit: the recorded run sampled on different hardware, so the RNG stream differs at a fixed seed.

A NOTE ON WHY THIS KERNEL CARRIES ITS OWN `generate.py`. GENERator-v2's remote `modeling_generator.py`
does its bp-level setup by RE-FETCHING its own tokenizer from the Hub inside `from_pretrained`, with no
cache_dir and no local_files_only. When that fetch fails the remote code CATCHES the error, prints a
warning, and returns a model whose `_bp_base_index` buffer was never registered; `generate()` then
raises, and the falsifier's `except Exception` records the generator as UNAVAILABLE with ZERO
sequences -- a complete-looking artifact in which the GLM scored nothing. The embedded `generate.py`
calls `setup_tokenizer` explicitly. The per-cell guard below treats n_candidate == 0 as a FAILED cell
rather than a result.

Everything is printed so the log itself is the audit trail.
"""
import base64
import io
import json
import os
import subprocess
import sys
import tarfile
import time
import traceback

os.environ.setdefault("PYTHONUTF8", "1")
OUT = "/kaggle/working"
REPO = OUT + "/repo"
HF = OUT + "/hf_cache"

# ---- the recorded baseline this sweep is measured against (from the committed artifact) ----
RECORDED = {
    "artifact": "wiki/glm_generator_falsifier_REAL_2026-10-07.json",
    "n": 120, "limit_natural": 1200, "bins": 6, "length": 150,
    "context": 300, "conditional": True, "natural_set": "promoter",
    "hf": {"kmer": 0.7237, "both": 0.7706},
    "markov_null": {"kmer": 0.5990, "both": 0.6679},
    "uniform": {"kmer": 0.9408, "both": 0.9641},
}

# The grid. (1.0, 0) is the recorded default and doubles as the self-check.
GRID = [(1.0, 0), (0.8, 0), (1.2, 0), (1.0, 40), (0.8, 40), (1.2, 40)]

PAYLOAD = "__PAYLOAD__"


def sh(cmd):
    print("$ " + cmd, flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.stdout:
        print(r.stdout[-4000:], flush=True)
    if r.returncode and r.stderr:
        print("STDERR:", r.stderr[-2000:], flush=True)
    return r.returncode


def unpack():
    os.makedirs(REPO, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(base64.b64decode(PAYLOAD.encode())), mode="r:gz") as tf:
        tf.extractall(REPO)
    print("unpacked sources:", flush=True)
    for root, _d, fs in os.walk(REPO):
        for f in sorted(fs):
            p = os.path.join(root, f)
            print("   ", os.path.relpath(p, REPO), os.path.getsize(p), "B", flush=True)


def device_report():
    """Print the GPU. The T4 pin is load-bearing and recorded in the project's lessons."""
    try:
        import torch
        print("torch", torch.__version__, "cuda_available", torch.cuda.is_available(), flush=True)
        if not torch.cuda.is_available():
            print("  !! NO CUDA -- the sweep would run on CPU, ~50x slower. Check the accelerator.",
                  flush=True)
            return "cpu"
        name = torch.cuda.get_device_name(0)
        cc = torch.cuda.get_device_capability(0)
        vram = torch.cuda.get_device_properties(0).total_memory / 2 ** 30
        print("GPU: %s  CC %d.%d  %.1f GiB" % (name, cc[0], cc[1], vram), flush=True)
        if cc[0] < 7:
            print("  !! CC < 7.0 -- this is the P100 trap: `enable_gpu` alone provisions a P100 "
                  "(sm_60) that current Kaggle torch cannot run in fp16. machine_shape must pin "
                  "NvidiaTeslaT4. Continuing in fp32, which is correct but slower.", flush=True)
        return "cuda"
    except Exception:
        traceback.print_exc()
        return "cpu"


def fetch_reference():
    """MG1655 via the repo's own idempotent downloader. The falsifier REFUSES without it."""
    sys.path.insert(0, REPO)
    from dna_decode.data.refseq import download_genome
    dest = REPO + "/data/cache/refseq"
    os.makedirs(dest, exist_ok=True)
    t0 = time.time()
    download_genome("GCF_000005845.2", dest)
    got = dest + "/GCF_000005845.2"
    print("reference fetched in %.0fs:" % (time.time() - t0), flush=True)
    for f in sorted(os.listdir(got)):
        print("   ", f, os.path.getsize(os.path.join(got, f)), "B", flush=True)
    for need in ("genome.fna", "annotations.gff3"):
        if not os.path.exists(os.path.join(got, need)):
            raise RuntimeError("reference incomplete: %s missing" % need)


def load_cli():
    sys.path.insert(0, REPO)
    sys.path.insert(0, REPO + "/scripts")
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "glm_falsifier_cli", REPO + "/scripts/glm_generator_falsifier.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def patch_weight_reuse():
    """Reuse the loaded model across cells. Same weights, so this cannot move a number; it just
    avoids paying the cold load six times. temperature/top_k are read at generate() time from the
    HFGenerator instance, so each cell still samples with its own settings, and
    generate_conditional() re-seeds torch per call, so cells stay independent of ordering."""
    import dna_decode.glm.generate as G
    cache = {}
    orig = G.HFGenerator.load

    def cached_load(self):
        if self._model is not None:
            return self
        key = (self.model_id, self.device, self.dtype, self.cache_dir)
        if key in cache:
            self._tok, self._model, self._dev = cache[key]
            return self
        t0 = time.time()
        orig(self)
        print("  [model loaded in %.0fs, cached for later cells]" % (time.time() - t0), flush=True)
        cache[key] = (self._tok, self._model, self._dev)
        return self

    G.HFGenerator.load = cached_load


def run_cell(cli, temperature, top_k, device):
    """One grid cell, at the RECORDED settings so the number is comparable."""
    out = "%s/cell_t%s_k%s.json" % (OUT, temperature, top_k)
    argv = [
        "--model", "auto", "--device", device, "--conditional",
        "--context", str(RECORDED["context"]),
        "--n", str(RECORDED["n"]),
        "--limit-natural", str(RECORDED["limit_natural"]),
        "--bins", str(RECORDED["bins"]),
        "--length", str(RECORDED["length"]),
        "--modes", "kmer,both",
        "--natural-set", RECORDED["natural_set"],
        "--cache-dir", HF,
        "--temperature", str(temperature),
        "--top-k", str(top_k),
        "--out", out,
    ]
    tag = "   <-- RECORDED DEFAULT (self-check)" if (temperature, top_k) == (1.0, 0) else ""
    print("\n" + "#" * 78, flush=True)
    print("# CELL temperature=%s  top_k=%s%s" % (temperature, top_k, tag), flush=True)
    print("#" * 78, flush=True)
    print("argv: " + " ".join(argv), flush=True)

    t0 = time.time()
    rc = cli.main(argv)
    el = (time.time() - t0) / 60
    print("# cell rc=%s in %.1f min" % (rc, el), flush=True)

    row = {"temperature": temperature, "top_k": top_k, "minutes": round(el, 1)}
    if rc != 0 or not os.path.exists(out):
        row["status"] = "FAILED"
        row["rc"] = rc
        return row
    art = json.load(open(out))
    row["status"] = "ok"
    row["sampling_stamp"] = art.get("sampling")
    row["modes"] = {}
    for r in art["results"]:
        if r["candidate"].startswith("hf:"):
            row["modes"][r["feature_mode"]] = {
                "distinguishability": r["distinguishability_auroc_LOWER_IS_BETTER"],
                "markov_null": r["markov_null_distinguishability"],
                "verdict": r["verdict"],
                "n_candidate": r["n_candidate"],
                "gc_candidate": r["gc_candidate"],
                "gc_natural": r["gc_natural"],
            }
    # A generator that produced ZERO sequences is the silent-failure mode this sweep must never
    # report as a result (see the module docstring).
    for m in row["modes"].values():
        if not m["n_candidate"]:
            row["status"] = "GENERATOR_PRODUCED_NOTHING"
    return row


def fmt(v):
    return ("%.4f" % v) if isinstance(v, (int, float)) else "    --  "


def table(rows):
    base = RECORDED["hf"]
    print("\n" + "=" * 92, flush=True)
    print("SWEEP RESULT  (distinguishability: LOWER = BETTER)", flush=True)
    print("=" * 92, flush=True)
    print("%5s %6s %8s %8s %8s %8s %7s  notes"
          % ("temp", "top_k", "kmer", "both", "mkKmer", "mkBoth", "gcCand"), flush=True)
    for r in rows:
        km = r.get("modes", {}).get("kmer", {})
        bo = r.get("modes", {}).get("both", {})
        note = ""
        if isinstance(km.get("distinguishability"), float) \
                and km["distinguishability"] < base["kmer"]:
            note += " BEATS-kmer"
        if isinstance(bo.get("distinguishability"), float) \
                and bo["distinguishability"] < base["both"]:
            note += " BEATS-both"
        if (r["temperature"], r["top_k"]) == (1.0, 0):
            note += " [DEFAULT/self-check]"
        if r.get("status") != "ok":
            note += " " + str(r.get("status"))
        print("%5.1f %6d %8s %8s %8s %8s %7s %s"
              % (r["temperature"], r["top_k"], fmt(km.get("distinguishability")),
                 fmt(bo.get("distinguishability")), fmt(km.get("markov_null")),
                 fmt(bo.get("markov_null")), fmt(km.get("gc_candidate")), note), flush=True)
    print("\nRECORDED DEFAULT to beat:  kmer %s   both %s" % (base["kmer"], base["both"]),
          flush=True)
    print("RECORDED markov null:      kmer %s   both %s"
          % (RECORDED["markov_null"]["kmer"], RECORDED["markov_null"]["both"]), flush=True)


def main():
    print("=" * 78, flush=True)
    print("GLM GENERATOR SAMPLING SWEEP -- temperature x top_k", flush=True)
    print("LOWER distinguishability is BETTER", flush=True)
    print("=" * 78, flush=True)
    unpack()
    device = device_report()
    sh("nvidia-smi")
    fetch_reference()
    patch_weight_reuse()
    cli = load_cli()

    rows = []
    for (t, k) in GRID:
        try:
            rows.append(run_cell(cli, t, k, device))
        except Exception:
            traceback.print_exc()
            rows.append({"temperature": t, "top_k": k, "status": "EXCEPTION",
                         "traceback": traceback.format_exc()[-1500:]})
        json.dump({"record": "glm-sampling-sweep-v1", "recorded_baseline": RECORDED,
                   "device": device, "grid": rows},
                  open(OUT + "/glm_sampling_sweep.json", "w"), indent=2)
    table(rows)
    print("\nwrote " + OUT + "/glm_sampling_sweep.json", flush=True)
    print("\nSWEEP COMPLETE", flush=True)


if __name__ == "__main__":
    main()
'''


def build() -> Path:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for rel in FILES:
            src = ROOT / rel
            if not src.exists():
                raise SystemExit(f"REFUSE: missing source {src}")
            tf.add(src, arcname=rel)
    payload = base64.b64encode(buf.getvalue()).decode()
    if "__PAYLOAD__" not in KERNEL:
        raise SystemExit("REFUSE: payload token missing from kernel template")
    kernel = KERNEL.replace("__PAYLOAD__", payload)
    out = ROOT / "scripts" / "kaggle" / "glm_sampling_sweep_kernel.py"
    out.write_text(kernel, encoding="utf-8")
    compile(kernel, str(out), "exec")          # syntax-check the emitted kernel, not just the template
    print(f"embedded {len(FILES)} files -> payload {len(payload) / 1024:.1f} KB b64")
    print(f"wrote {out}  ({out.stat().st_size / 1024:.1f} KB total)  [syntax OK]")
    if out.stat().st_size > 900_000:
        print("WARNING: approaching the ~1 MB Kaggle kernel source cap")
    return out


if __name__ == "__main__":
    build()
