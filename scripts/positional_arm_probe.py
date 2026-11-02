"""Validity FIRST: prove the positional arm is runnable and the parity inverse holds, before any bar.

This probe must pass BEFORE `freeze_positional_arm_bar.py` writes a threshold, mirroring how the
sibling arms derived their bars from a feasibility run rather than from judgement. It measures what
the bar then inherits -- the parity deviation, the throughput, the per-rung RAM, and the lag set --
so no figure in the bar is a preference dressed as a measurement.

**It also OWNS the cohort contract, and that is deliberate.** `scripts/promoteratlas_local_reconcile.py`
builds the committed embedding cache from `load_peak_tiles(...)` filtered to
`category == "tile" and len(seq) == 150`, and the cache is indexed by that FILTERED list. Measured
here: 46,713 rows parsed -> 45,201 of category `tile` -> **44,106** after the length filter, which is
exactly the cache's row count. A second copy of that filter in the extract and score runners would
be three chances to drift from the cache by one row and still produce a plausible number, so
`load_arm_cohort()` lives here once and the other two runners import it. `genomewide` stays the
registered owner of GSE144621; this adds no second reader of the file, only one shared reader of the
filter.

**What is measured here rather than assumed:**

- the checkpoint's sha256, because it lives in a VOLATILE temp directory
  (`C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas`) that a reboot or a cleaner can
  remove; a hash lets a replacement be verified identical rather than assumed
- the parameter count against `atlas_embed.EXPECTED_PARAMS`
- the forward contract: `model(X, return_latent=True)[2]` is `(B, 128, L)`
- the parity inverse against the COMMITTED embedding rows for the same tiles
- throughput, extrapolated to the full pass
- free RAM, and the peak design-matrix bytes for every ladder rung against a named headroom factor
- the autocorrelation lag set, chosen by the pre-declared rule below rather than by preference

**The lag-set rule is declared before the numbers, not after.** The plan left dense-versus-log-spaced
open to this probe. The rule is: take the LARGEST candidate lag set whose primary design matrix fits
in RAM under the headroom factor. A larger lag set is a strict superset of the information, so this
only ever gives up resolution to a resource limit; it never trades resolution for a nicer number.

Exits non-zero on ANY failure, because a broken probe that exited 0 would read as a clean
feasibility and the bar would be frozen on top of it.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
import platform
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/positional_arm_feasibility_2026-11-01.json"
EMB_CACHE = REPO / "data/promoteratlas/tile_embeddings_native150_n44106.npy"
SUBSTRATE_CACHE = Path(os.environ.get("MPRA_CACHE", "D:/dna_decode_cache/mpra"))

#: What the committed cache's row count must be. A cohort of a different size is not this arm.
EXPECTED_TILES = 44_106
EXPECTED_CHANNELS = 128
EXPECTED_LENGTH = 150

#: Float64 copies of a design matrix alive at once inside the ridge fit. `StandardScaler.transform`
#: returns a NEW array while the caller still holds the original, so two is the honest count.
RIDGE_COPIES = 2

#: Multiplier applied to the peak requirement before comparing against free RAM. Named, not hidden:
#: this project already has a registered constraint that declared the kmer6 rung infeasible at a
#: 9.8 GB design matrix, so a pre-check is house practice rather than caution.
RAM_HEADROOM = 1.5

#: Candidate autocorrelation lag sets, largest first. The pre-declared rule picks the first that
#: fits; see the module docstring.
LAG_CANDIDATES = {
    "dense_1_20": list(range(1, 21)),
    "log_spaced_9": [1, 2, 3, 4, 5, 7, 10, 14, 20],
    "coarse_5": [2, 5, 10, 15, 20],
}

#: Default parity tolerance. The probe REPORTS the achieved deviation so the bar can freeze a
#: tolerance from the measurement; this is only the pass/fail used by the probe itself.
PARITY_ATOL = 1e-4


class ProbeFailure(RuntimeError):
    """A feasibility precondition did not hold. Never swallowed: the probe exits non-zero."""


# ---------------------------------------------------------------------------------------------------
# PURE helpers -- unit-testable with no model, no D:, no checkpoint
# ---------------------------------------------------------------------------------------------------
def design_matrix_bytes(n_rows: int, n_features: int, *, itemsize: int = 8,
                        copies: int = RIDGE_COPIES) -> int:
    """Peak bytes held at once for a ridge fit over `n_rows x n_features`.

    `copies` is 2 because the scaler's transform allocates a second array while the first is still
    referenced. Train and test are counted together: both are resident during `predict`.
    """
    for name, v in (("n_rows", n_rows), ("n_features", n_features), ("itemsize", itemsize),
                    ("copies", copies)):
        if int(v) < 1:
            raise ValueError("%s must be >= 1, got %r" % (name, v))
    return int(n_rows) * int(n_features) * int(itemsize) * int(copies)


def stored_block_bytes(n_tiles: int, n_features_total: int, *, itemsize: int = 4) -> int:
    """On-disk bytes of the extracted feature block (float32 by default)."""
    if int(n_tiles) < 1 or int(n_features_total) < 1:
        raise ValueError("n_tiles and n_features_total must both be >= 1")
    return int(n_tiles) * int(n_features_total) * int(itemsize)


def fits_in_ram(need_bytes: int, free_bytes: int | None, *, headroom: float = RAM_HEADROOM):
    """True / False / None. **None means UNKNOWN and must never be read as True.**

    A probe that could not measure free RAM has not shown a rung is affordable; the scoring runner
    treats `None` the same as a failed pre-check and records the rung as skipped.
    """
    if headroom < 1.0:
        raise ValueError("headroom must be >= 1.0, got %r" % (headroom,))
    if free_bytes is None:
        return None
    return bool(int(need_bytes) * float(headroom) <= int(free_bytes))


def extrapolate_wall_seconds(n_total: int, n_sample: int, elapsed_seconds: float) -> float:
    """Linear extrapolation of a sampled forward pass to the full cohort."""
    if int(n_sample) < 1:
        raise ValueError("n_sample must be >= 1")
    if float(elapsed_seconds) <= 0:
        raise ValueError("elapsed_seconds must be > 0 to extrapolate a rate")
    return float(elapsed_seconds) * int(n_total) / int(n_sample)


def free_ram_bytes() -> int | None:
    """Available physical memory, or None when it cannot be measured.

    `psutil` is NOT installed here, so this uses the platform call directly: `GlobalMemoryStatusEx`
    on Windows, `/proc/meminfo` elsewhere. Returning None rather than a guess is the point -- see
    `fits_in_ram`.
    """
    try:
        if platform.system() == "Windows":
            class _MemStatusEx(ctypes.Structure):
                _fields_ = [
                    ("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]
            st = _MemStatusEx()
            st.dwLength = ctypes.sizeof(_MemStatusEx)
            if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):
                return None
            return int(st.ullAvailPhys)
        meminfo = Path("/proc/meminfo")
        if meminfo.is_file():
            for line in meminfo.read_text().splitlines():
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
        return None
    except Exception:
        return None


def rung_feature_counts(n_lags: int, *, channels: int = EXPECTED_CHANNELS) -> dict[str, int]:
    """Feature count per ladder rung. Pure arithmetic over the geometry."""
    if int(n_lags) < 1 or int(channels) < 1:
        raise ValueError("n_lags and channels must both be >= 1")
    c = int(channels)
    return {
        "bins10_mean_oriented": c * 10,
        "autocorr": c * int(n_lags),
        "bins10_mean_asgiven": c * 10,
        "bins10_max_oriented": c * 10,
        "bins2_mean_oriented": c * 2,
        "bins5_mean_oriented": c * 5,
    }


def choose_lag_set(candidates: dict[str, list[int]], *, n_rows: int, channels: int,
                   free_bytes: int | None, headroom: float = RAM_HEADROOM):
    """The pre-declared rule: the LARGEST candidate whose design matrix fits under the headroom.

    Returns `(name, lags, decisions)` where `decisions` records every candidate's bytes and verdict,
    so the choice is auditable and a reader can see what was given up. Raises when none fits --
    better than silently scoring a lag set nobody chose.
    """
    decisions = []
    chosen = None
    for name, lags in sorted(candidates.items(), key=lambda kv: -len(kv[1])):
        need = design_matrix_bytes(n_rows, channels * len(lags))
        ok = fits_in_ram(need, free_bytes, headroom=headroom)
        decisions.append({
            "candidate": name, "n_lags": len(lags), "n_features": channels * len(lags),
            "peak_design_bytes": need, "fits_in_ram": ok,
        })
        if ok is True and chosen is None:
            chosen = (name, list(lags))
    if chosen is None:
        raise ProbeFailure(
            "no candidate autocorrelation lag set fits in RAM under headroom %.2f (free=%r). "
            "Refusing to pick one anyway: the arm's primary would then be a rung nobody chose. "
            "Decisions: %s" % (headroom, free_bytes, json.dumps(decisions))
        )
    return chosen[0], chosen[1], decisions


def parity_report(derived, ref, *, atol: float = PARITY_ATOL) -> dict:
    """Compare a DERIVED global mean against the committed embedding rows. Pure, so it is testable.

    Split out of `run` on purpose: the parity gate is the single check that makes every positional
    rung comparable to the incumbent's 0.5076, and a gate whose arithmetic can only be exercised by
    a 686 MB checkpoint and a forward pass is a gate nobody can prove non-vacuous.
    """
    a = np.asarray(derived, dtype=np.float64)
    b = np.asarray(ref, dtype=np.float64)
    if a.shape != b.shape:
        raise ProbeFailure(
            "parity shapes differ: derived %r vs committed %r. A positional pairing is impossible, "
            "and broadcasting them would compare the wrong rows." % (a.shape, b.shape)
        )
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        raise ProbeFailure("parity inputs hold non-finite values; refusing to report a deviation")
    dev = np.abs(a - b)
    max_dev = float(dev.max())
    return {
        "n_rows": int(a.shape[0]),
        "max_abs_deviation": max_dev,
        "mean_abs_deviation": float(dev.mean()),
        "p99_abs_deviation": float(np.percentile(dev, 99)),
        "atol_used": float(atol),
        "passed": bool(max_dev <= float(atol)),
    }


def sha256_of(path: Path, *, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def stride_sample(n_total: int, n_sample: int) -> list[int]:
    """Deterministic evenly-spaced indices. A stride covers the file; a head slice does not."""
    if int(n_sample) < 1 or int(n_total) < 1:
        raise ValueError("n_total and n_sample must both be >= 1")
    k = min(int(n_sample), int(n_total))
    step = int(n_total) / k
    return sorted({min(int(n_total) - 1, int(i * step)) for i in range(k)})


# ---------------------------------------------------------------------------------------------------
# the cohort contract -- ONE definition, imported by the extract and score runners
# ---------------------------------------------------------------------------------------------------
def load_arm_cohort(cache_dir: Path | str = SUBSTRATE_CACHE, *, expect: int = EXPECTED_TILES):
    """The arm's tiles, filtered EXACTLY as the committed embedding cache was built.

    Refuses a count other than `expect`: the cache is indexed positionally, so an off-by-one cohort
    would pair every feature row with the wrong embedding row and still score.
    """
    from dna_decode.glm import genomewide as gw

    tiles_all, st = gw.load_peak_tiles(cache_dir, allow_download=False)
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == EXPECTED_LENGTH]
    if len(tiles) != expect:
        raise ProbeFailure(
            "cohort is %d tiles, the committed cache is indexed for %d. The embedding cache is "
            "positional, so scoring a differently-sized cohort would pair feature rows with the "
            "wrong embedding rows. Parsed %d rows from %s."
            % (len(tiles), expect, st.n_parsed, st.path)
        )
    return tiles, st


# ---------------------------------------------------------------------------------------------------
# the model-touching probe
# ---------------------------------------------------------------------------------------------------
def run(*, n_sample: int, batch: int, device: str, atol: float) -> dict:
    import torch

    from dna_decode.glm.atlas_embed import (
        EXPECTED_PARAMS, NUC, AtlasEmbedError, DEFAULT_PA_REPO, audit_checkpoint, load_atlas,
    )
    from dna_decode.glm.positional_readout import bin_means, global_mean_from_bin_means

    rec: dict = {
        "schema": "positional-arm-feasibility-v1",
        "date": "2026-11-01",
        "arm": "positional-non-mean-pooled-promoteratlas-150bp",
        "host": {"platform": platform.platform(), "python": platform.python_version(),
                 "torch": torch.__version__, "cuda": bool(torch.cuda.is_available()),
                 "threads": torch.get_num_threads()},
    }
    t0 = time.time()

    # --- the checkpoint, which lives in a volatile temp directory -----------------------------
    weights = Path(DEFAULT_PA_REPO) / "trained_weights/base_model/promoteratlas-base.pt"
    if not weights.is_file():
        raise ProbeFailure(
            "checkpoint absent at %s. It lives in a VOLATILE temp directory, so this is a live "
            "fact and not an inherited one; re-clone the MIT-licensed PromoterAtlas repo (or set "
            "PA_REPO) and re-run. No number may be read until it resolves." % weights
        )
    clone_bytes = sum(f.stat().st_size for f in Path(DEFAULT_PA_REPO).rglob("*") if f.is_file())
    rec["checkpoint"] = {
        "path": str(weights),
        "bytes": weights.stat().st_size,
        "mb": round(weights.stat().st_size / (1 << 20), 2),
        "sha256": sha256_of(weights),
        "clone_bytes": clone_bytes,
        "clone_mb": round(clone_bytes / (1 << 20), 0),
        "size_note": ("the BASE CHECKPOINT FILE is ~21 MB; the 686 MB figure in the surrounding "
                      "notes is the whole CLONE, which also ships five expression-prediction heads "
                      "and a segmentation head at ~11 MB each. 1,847,428 params x 4 bytes is only "
                      "7.05 MB, so the 21 MB file carries more than the weights. Recorded because "
                      "a conflated size reads as a checkpoint-integrity fact and is not one."),
        "volatile_location": True,
        "why_hashed": ("a temp-directory checkpoint can be removed by a reboot or a cleaner; the "
                       "hash lets a replacement be verified identical rather than assumed"),
    }
    rec["checkpoint"]["pickle_globals"] = audit_checkpoint(weights)

    try:
        model = load_atlas()
    except AtlasEmbedError as e:
        raise ProbeFailure("load_atlas refused: %s" % e) from e
    n_params = sum(p.numel() for p in model.parameters())
    if n_params != EXPECTED_PARAMS:
        raise ProbeFailure("%d parameters, the arm's bar expects %d" % (n_params, EXPECTED_PARAMS))
    rec["model"] = {"n_params": n_params, "expected_params": EXPECTED_PARAMS,
                    "params_match": True}

    # --- the cohort ----------------------------------------------------------------------------
    tiles, st = load_arm_cohort()
    strands = {}
    for t in tiles:
        strands[t.strand] = strands.get(t.strand, 0) + 1
    rec["substrate"] = {
        "path": st.path, "n_parsed": st.n_parsed, "n_cohort": len(tiles),
        "filter": "category == 'tile' and len(seq) == %d" % EXPECTED_LENGTH,
        "why_this_filter": ("it is the filter scripts/promoteratlas_local_reconcile.py applied "
                            "before building the committed embedding cache, which is indexed "
                            "positionally against the filtered list"),
        "sequence_length": EXPECTED_LENGTH,
        "n_peaks": len({t.peak for t in tiles}),
        "strand_counts": strands,
        "n_duplicate_sequences": len(tiles) - len({t.seq for t in tiles}),
    }

    # --- the forward contract ------------------------------------------------------------------
    probe_seqs = [tiles[0].seq, tiles[1].seq]
    X = torch.zeros((2, 4, EXPECTED_LENGTH))
    for b, s in enumerate(probe_seqs):
        for j, c in enumerate(s):
            k = NUC.get(c)
            if k is not None:
                X[b, k, j] = 1.0
    with torch.no_grad():
        fwd = model(X.to(device), return_latent=True)
    if not isinstance(fwd, tuple) or len(fwd) != 3:
        raise ProbeFailure(
            "model(X, return_latent=True) returned %r of length %r; the reduction seam assumes a "
            "3-tuple whose third element is the latent"
            % (type(fwd).__name__, len(fwd) if hasattr(fwd, "__len__") else "n/a")
        )
    z = fwd[2]
    if tuple(z.shape) != (2, EXPECTED_CHANNELS, EXPECTED_LENGTH):
        raise ProbeFailure(
            "latent shape %r, expected (B, %d, %d). Every reduction in this arm is written against "
            "that geometry." % (tuple(z.shape), EXPECTED_CHANNELS, EXPECTED_LENGTH)
        )
    rec["forward_contract"] = {"returns": "3-tuple", "latent_shape": list(z.shape), "ok": True}

    # --- parity on a deterministic stride sample ------------------------------------------------
    if not EMB_CACHE.is_file():
        raise ProbeFailure("committed embedding cache absent at %s" % EMB_CACHE)
    committed = np.load(EMB_CACHE, mmap_mode="r")
    if committed.shape != (len(tiles), EXPECTED_CHANNELS):
        raise ProbeFailure(
            "embedding cache is %r, cohort is %d tiles; a positional pairing is impossible"
            % ((committed.shape,), len(tiles))
        )

    idx = stride_sample(len(tiles), n_sample)
    seqs = [tiles[i].seq for i in idx]
    derived = np.zeros((len(idx), EXPECTED_CHANNELS), dtype=np.float64)
    t_fwd = time.time()
    with torch.no_grad():
        for s0 in range(0, len(seqs), batch):
            chunk = seqs[s0:s0 + batch]
            Xb = torch.zeros((len(chunk), 4, EXPECTED_LENGTH))
            for b, s in enumerate(chunk):
                for j, c in enumerate(s):
                    k = NUC.get(c)
                    if k is not None:
                        Xb[b, k, j] = 1.0
            zz = model(Xb.to(device), return_latent=True)[2].float().cpu().numpy()
            bm = bin_means(zz, 10)                       # (b, 128, 10) -- what gets STORED
            derived[s0:s0 + len(chunk)] = global_mean_from_bin_means(bm)
    fwd_seconds = time.time() - t_fwd

    ref = np.asarray(committed[idx], dtype=np.float64)
    par = parity_report(derived, ref, atol=atol)
    par["n_sampled"] = par.pop("n_rows")
    par["sampling"] = "deterministic stride over the cohort"
    par["what_it_proves"] = (
        "the global mean DERIVED from the stored 10-bin means reproduces the committed mean-pool "
        "embedding, so every positional rung is measured on the same representation as the 0.5076 "
        "incumbent"
    )
    rec["parity"] = par
    if not par["passed"]:
        raise ProbeFailure(
            "parity deviation %.3e exceeds atol %.3e over %d tiles. The forward pass has drifted "
            "from the one that built the committed cache, so NO number from this arm would be "
            "comparable to 0.5076." % (par["max_abs_deviation"], atol, par["n_sampled"])
        )

    # --- throughput ----------------------------------------------------------------------------
    rate = len(seqs) / fwd_seconds if fwd_seconds else 0.0
    full = extrapolate_wall_seconds(len(tiles), len(seqs), fwd_seconds)
    rec["throughput"] = {
        "n_sampled": len(seqs), "sample_seconds": round(fwd_seconds, 2),
        "seq_per_second": round(rate, 2), "batch": batch, "device": device,
        "full_pass_seconds_extrapolated": round(full, 1),
        "full_pass_minutes_extrapolated": round(full / 60.0, 2),
        "caveat": ("linear extrapolation from a %d-tile sample; contention from another job moves "
                   "it, which is why it is recorded with the host's thread count" % len(seqs)),
    }

    # --- RAM, the binding constraint, and the lag-set choice ------------------------------------
    free = free_ram_bytes()
    # the largest split the arm will ever fit: all tiles resident (train + test together)
    name, lags, lag_decisions = choose_lag_set(
        LAG_CANDIDATES, n_rows=len(tiles), channels=EXPECTED_CHANNELS, free_bytes=free,
    )
    counts = rung_feature_counts(len(lags))
    rungs = []
    for rung, n_feat in sorted(counts.items(), key=lambda kv: -kv[1]):
        need = design_matrix_bytes(len(tiles), n_feat)
        rungs.append({
            "rung": rung, "n_features": n_feat, "peak_design_bytes": need,
            "peak_design_gb": round(need / 1e9, 3),
            "fits_in_ram": fits_in_ram(need, free, headroom=RAM_HEADROOM),
        })
    stored_features = counts["bins10_mean_oriented"] + counts["bins10_max_oriented"] + counts["autocorr"]
    rec["resources"] = {
        "free_ram_bytes": free,
        "free_ram_gb": round(free / 1e9, 2) if free is not None else None,
        "ram_measured_via": "GlobalMemoryStatusEx" if platform.system() == "Windows" else "/proc/meminfo",
        "headroom_factor": RAM_HEADROOM,
        "ridge_float64_copies": RIDGE_COPIES,
        "rungs": rungs,
        "n_rungs_fitting": sum(1 for r in rungs if r["fits_in_ram"] is True),
        "n_rungs_not_fitting": sum(1 for r in rungs if r["fits_in_ram"] is not True),
        "unknown_is_not_ok": ("fits_in_ram=None means UNMEASURED, and the scoring runner treats it "
                              "as a failed pre-check rather than as affordable"),
    }
    rec["lag_set"] = {
        "chosen": name, "lags": lags, "n_lags": len(lags),
        "n_features": EXPECTED_CHANNELS * len(lags),
        "rule": ("the LARGEST candidate whose design matrix fits in RAM under the headroom factor; "
                 "declared before the numbers, and only ever gives up resolution to a resource "
                 "limit rather than trading it for a nicer number"),
        "candidates_considered": lag_decisions,
    }
    rec["stored_block"] = {
        "blocks": ["bins10_mean", "bins10_max", "autocorr"],
        "n_features_total": stored_features,
        "dtype": "float32",
        "bytes": stored_block_bytes(len(tiles), stored_features),
        "gb": round(stored_block_bytes(len(tiles), stored_features) / 1e9, 3),
        "destination": "D:/dna_decode_cache/ per the standing D: rule, not a tight C:",
        "note": ("bins10_mean is stored UNORIENTED; `orient` acts on the stored bin axis, so both "
                 "orientation rungs cost no second forward pass"),
    }
    rec["wall_seconds"] = round(time.time() - t0, 1)
    rec["verdict"] = "FEASIBLE"
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", type=int, default=512, help="tiles for the parity/throughput sample")
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--atol", type=float, default=PARITY_ATOL)
    args = ap.parse_args(argv)

    try:
        rec = run(n_sample=args.sample, batch=args.batch, device=args.device, atol=args.atol)
    except Exception as e:
        print("PROBE FAILED (%s): %s" % (type(e).__name__, e), flush=True)
        return 1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    p, r, res = rec["parity"], rec["throughput"], rec["resources"]
    print("checkpoint   sha256 %s... (%.2f MB file, %.0f MB clone)"
          % (rec["checkpoint"]["sha256"][:16], rec["checkpoint"]["mb"],
             rec["checkpoint"]["clone_mb"]))
    print("params       %d == %d" % (rec["model"]["n_params"], rec["model"]["expected_params"]))
    print("latent       %r" % (rec["forward_contract"]["latent_shape"],))
    print("cohort       %d tiles / %d peaks / %r"
          % (rec["substrate"]["n_cohort"], rec["substrate"]["n_peaks"],
             rec["substrate"]["strand_counts"]))
    print("PARITY       max |dev| %.3e over %d tiles (atol %.0e) -> %s"
          % (p["max_abs_deviation"], p["n_sampled"], p["atol_used"],
             "PASS" if p["passed"] else "FAIL"))
    print("throughput   %.1f seq/s -> full pass ~%.1f min"
          % (r["seq_per_second"], r["full_pass_minutes_extrapolated"]))
    print("lag set      %s = %r (%d features)"
          % (rec["lag_set"]["chosen"], rec["lag_set"]["lags"], rec["lag_set"]["n_features"]))
    print("stored block %d features -> %.2f GB float32"
          % (rec["stored_block"]["n_features_total"], rec["stored_block"]["gb"]))
    print("RAM          %.2f GB free, headroom %.1fx -> %d/%d rungs fit"
          % (res["free_ram_gb"] if res["free_ram_gb"] is not None else float("nan"),
             res["headroom_factor"], res["n_rungs_fitting"],
             res["n_rungs_fitting"] + res["n_rungs_not_fitting"]))
    for rr in res["rungs"]:
        print("   %-22s %5d feat  %6.3f GB  fits=%s"
              % (rr["rung"], rr["n_features"], rr["peak_design_gb"], rr["fits_in_ram"]))
    print("VERDICT      %s  (%.1fs)" % (rec["verdict"], rec["wall_seconds"]))
    print("-> %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
