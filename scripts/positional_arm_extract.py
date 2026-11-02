"""ONE forward pass over the 44,106 tiles, emitting every positional rung's features at once.

The whole ladder is a different REDUCTION over an identical forward pass, so it costs one pass, not
six. This runner does that pass and stores the three reduction blocks; `positional_arm_score.py`
then slices rungs out of the stored block without touching the model again.

**Stored layout, channel-major throughout, 5,120 float32 per tile:**

    [    0 : 1280 ]  bins10_mean   -- 128 channels x 10 bins, UNORIENTED
    [ 1280 : 2560 ]  bins10_max    -- 128 channels x 10 bins, UNORIENTED
    [ 2560 : 5120 ]  autocorr      -- 128 channels x 20 lags

Both bin blocks are stored UNORIENTED on purpose: `orient` acts on the stored bin axis, so the
oriented and as-given rungs both come out of these same bytes with no second pass. `bins2` and
`bins5` are likewise derived from `bins10` at scoring time -- 10 bins nest exactly into 5 and 2, so
re-binning the stored means is arithmetically identical to re-binning the latent, and storing them
separately would be three more chances to disagree.

**The autocorrelation block needs no orientation variant, and that is a fact rather than a
simplification.** An autocorrelation is invariant under reversal of its series: the set of position
pairs separated by `k` is unchanged when the series is read backwards. So the `autocorr` co-primary
is IMMUNE to the arm's unresolved orientation convention -- the one risk the bins rungs carry a
1-bit inference about. That immunity is a genuine argument for having declared it a co-primary, and
a test pins it.

**Crash discipline.** The block is a real `.npy` memmap written in place and the sidecar records how
many rows are done, so a crash costs ONE BATCH rather than the pass. This project has a measured
precedent of a mid-run crash leaving a partial embedding cache that looked complete to its consumer
(`EmbeddingCache.populate(skip_existing=True)` skipping at gene-dataset level), which is why the row
count is recorded explicitly and a resumed run is byte-checked against an uninterrupted one in the
tests rather than assumed equal.

**Refusals.** A ragged cohort is refused, not padded. An existing block whose sidecar geometry
differs is refused, not overwritten. A `--limit` run stamps the sidecar `PARTIAL` and the scoring
runner's loader refuses a PARTIAL block, so a smoke run can never be scored as if it were whole.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np

from scripts.positional_arm_probe import (
    EXPECTED_CHANNELS,
    EXPECTED_LENGTH,
    EXPECTED_TILES,
    ProbeFailure,
    load_arm_cohort,
    sha256_of,
)

REPO = Path(__file__).resolve().parents[1]
FEASIBILITY = REPO / "wiki/positional_arm_feasibility_2026-11-01.json"

#: The standing D: rule: heavy artifacts never land on a chronically-full C:.
DEFAULT_BLOCK_DIR = Path(os.environ.get(
    "POSITIONAL_BLOCK_DIR", "D:/dna_decode_cache/promoteratlas_positional"))

BLOCK_NAME = "positional_block.npy"
SIDECAR_NAME = "positional_block_sidecar.json"

N_BINS = 10
BATCH = 64
DTYPE = "float32"


class ExtractError(RuntimeError):
    """The block cannot be written, resumed, or trusted."""


# ---------------------------------------------------------------------------------------------------
# the stored layout -- ONE definition, imported by the scoring runner
# ---------------------------------------------------------------------------------------------------
def feature_layout(n_lags: int, *, channels: int = EXPECTED_CHANNELS,
                   n_bins: int = N_BINS) -> dict:
    """Byte-free description of the stored column blocks. Pure; the scorer imports this."""
    if int(n_lags) < 1 or int(channels) < 1 or int(n_bins) < 1:
        raise ValueError("n_lags, channels and n_bins must all be >= 1")
    c, b, g = int(channels), int(n_bins), int(n_lags)
    off = 0
    layout = {}
    for name, width in (("bins10_mean", c * b), ("bins10_max", c * b), ("autocorr", c * g)):
        layout[name] = {"start": off, "stop": off + width, "width": width}
        off += width
    layout["_total"] = off
    layout["_channels"] = c
    layout["_n_bins"] = b
    layout["_n_lags"] = g
    return layout


def slice_block(block, layout: dict, name: str) -> np.ndarray:
    """Pull one named block out of the stored matrix. Refuses an unknown name."""
    if name not in layout or name.startswith("_"):
        raise ExtractError(
            "no stored block named %r; stored blocks are %r"
            % (name, [k for k in layout if not k.startswith("_")])
        )
    s = layout[name]
    a = np.asarray(block)
    if a.ndim != 2 or a.shape[1] != layout["_total"]:
        raise ExtractError(
            "block is %r but the layout describes %d columns; refusing to slice a mismatched block"
            % (a.shape, layout["_total"])
        )
    return a[:, s["start"]:s["stop"]]


def as_bins(flat, layout: dict) -> np.ndarray:
    """`(B, C*n_bins)` -> `(B, C, n_bins)`, so `orient` can act on the bin axis."""
    a = np.asarray(flat)
    return a.reshape(a.shape[0], layout["_channels"], layout["_n_bins"])


def rebin(bins, n_target: int) -> np.ndarray:
    """Coarsen `(B, C, n)` bin MEANS to `n_target` bins. Exact when `n_target` divides `n`.

    10 bins nest into 5 and into 2, so the mean of equal groups of equal bins is exactly the mean
    over the underlying positions -- re-binning the stored means equals re-binning the latent, which
    is why `bins2` and `bins5` need no stored block of their own.
    """
    a = np.asarray(bins, dtype=np.float64)
    if a.ndim != 3:
        raise ExtractError("rebin expects (B, C, n_bins), got %r" % (a.shape,))
    n = a.shape[2]
    t = int(n_target)
    if t < 1 or t > n or n % t:
        raise ExtractError(
            "cannot coarsen %d bins to %d: the groups would be unequal, so the result would no "
            "longer equal a mean over positions" % (n, t)
        )
    return a.reshape(a.shape[0], a.shape[1], t, n // t).mean(axis=3)


# ---------------------------------------------------------------------------------------------------
# the encoder seam -- injectable so the runner is testable with a stub
# ---------------------------------------------------------------------------------------------------
def make_atlas_encoder(*, device: str = "cpu"):
    """Return `encode(seqs) -> (b, 128, L) float32`, backed by the frozen PromoterAtlas base."""
    import torch

    from dna_decode.glm.atlas_embed import NUC, load_atlas

    model = load_atlas()

    def encode(seqs: list[str]) -> np.ndarray:
        X = torch.zeros((len(seqs), 4, len(seqs[0])))
        for b, s in enumerate(seqs):
            for j, ch in enumerate(s):
                k = NUC.get(ch)
                if k is not None:
                    X[b, k, j] = 1.0
        with torch.no_grad():
            z = model(X.to(device), return_latent=True)[2]
        return z.float().cpu().numpy()

    return encode


def reduce_batch(z, lags, *, n_bins: int = N_BINS) -> np.ndarray:
    """All three stored blocks for one batch -> `(b, total_features)` float32."""
    from dna_decode.glm.positional_readout import (
        bin_maxes, bin_means, channel_autocorrelation, flatten_bins,
    )

    bm = flatten_bins(bin_means(z, n_bins))
    bx = flatten_bins(bin_maxes(z, n_bins))
    ac = channel_autocorrelation(z, lags)
    return np.concatenate([bm, bx, ac], axis=1).astype(np.float32)


# ---------------------------------------------------------------------------------------------------
# the run
# ---------------------------------------------------------------------------------------------------
def _read_sidecar(path: Path) -> dict | None:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        raise ExtractError(
            "sidecar %s is not parseable JSON (%s). Refusing to guess at the block's provenance; "
            "delete the block and the sidecar together to start clean." % (path, e)
        ) from e


def _geometry(n_rows: int, lags, n_bins: int) -> dict:
    lay = feature_layout(len(lags), n_bins=n_bins)
    return {
        "n_rows": int(n_rows), "n_features": lay["_total"], "dtype": DTYPE,
        "n_bins": int(n_bins), "lags": list(lags), "channels": EXPECTED_CHANNELS,
        "sequence_length": EXPECTED_LENGTH,
    }


def extract(*, encode, seqs: list[str], lags: list[int], block_dir: Path,
            n_bins: int = N_BINS, batch: int = BATCH, partial: bool = False,
            checkpoint_sha: str | None = None, progress=None) -> dict:
    """Write the stored block, resumably. Returns the sidecar it wrote."""
    if not seqs:
        raise ExtractError("extract got zero sequences")
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ExtractError(
            "ragged cohort (%d distinct lengths, %d-%d). This runner REFUSES to pad: a silent pad "
            "is how pretraining encoded a MASKED base, so it is not a neutral operation on this "
            "encoder and would produce a plausible meaningless number."
            % (len(lengths), min(lengths), max(lengths))
        )

    geom = _geometry(len(seqs), lags, n_bins)
    lay = feature_layout(len(lags), n_bins=n_bins)
    block_dir = Path(block_dir)
    block_dir.mkdir(parents=True, exist_ok=True)
    bpath, spath = block_dir / BLOCK_NAME, block_dir / SIDECAR_NAME

    prev = _read_sidecar(spath)
    start = 0
    if prev is not None and bpath.is_file():
        if prev.get("geometry") != geom:
            raise ExtractError(
                "an existing block at %s has geometry %r but this run needs %r. Refusing to "
                "overwrite: a block scored under the wrong geometry would reshape cleanly and mean "
                "nothing. Move or delete the block and its sidecar to re-extract."
                % (bpath, prev.get("geometry"), geom)
            )
        start = int(prev.get("n_rows_done", 0))
        if start >= len(seqs):
            print("block already complete (%d rows); nothing to do" % start, flush=True)
            return prev
        print("resuming at row %d of %d" % (start, len(seqs)), flush=True)

    mode = "r+" if (bpath.is_file() and prev is not None) else "w+"
    if mode == "w+":
        arr = np.lib.format.open_memmap(bpath, mode="w+", dtype=np.float32,
                                        shape=(len(seqs), lay["_total"]))
    else:
        arr = np.lib.format.open_memmap(bpath, mode="r+")
        if arr.shape != (len(seqs), lay["_total"]):
            raise ExtractError(
                "existing block is %r, this run needs %r" % (arr.shape, (len(seqs), lay["_total"]))
            )

    t0 = time.time()
    done = start
    try:
        for s0 in range(start, len(seqs), batch):
            chunk = seqs[s0:s0 + batch]
            z = encode(chunk)
            z = np.asarray(z)
            if z.ndim != 3 or z.shape[0] != len(chunk):
                raise ExtractError(
                    "encoder returned %r for a batch of %d; expected (b, C, L)"
                    % (z.shape, len(chunk))
                )
            arr[s0:s0 + len(chunk)] = reduce_batch(z, lags, n_bins=n_bins)
            done = s0 + len(chunk)
            arr.flush()
            sidecar = {
                "schema": "positional-block-sidecar-v1",
                "date": "2026-11-01",
                "arm": "positional-non-mean-pooled-promoteratlas-150bp",
                "geometry": geom,
                "layout": {k: v for k, v in lay.items() if not k.startswith("_")},
                "n_rows_done": done,
                "complete": done >= len(seqs),
                "status": "PARTIAL" if partial else ("COMPLETE" if done >= len(seqs) else "IN_PROGRESS"),
                "partial_limit_run": bool(partial),
                "checkpoint_sha256": checkpoint_sha,
                "block": str(bpath),
                "wall_seconds": round(time.time() - t0, 1),
                "stored_unoriented": True,
                "why_unoriented": ("orient acts on the stored bin axis, so both orientation rungs "
                                   "come out of these bytes with no second forward pass"),
                "autocorr_orientation_invariant": True,
                "why_invariant": ("the set of position pairs separated by k is unchanged when the "
                                  "series is reversed, so the autocorr co-primary is immune to the "
                                  "arm's unresolved orientation convention"),
            }
            tmp = spath.with_suffix(".json.tmp")
            tmp.write_text(json.dumps(sidecar, indent=2, sort_keys=True), encoding="utf-8")
            tmp.replace(spath)
            if progress is not None:
                progress(done, len(seqs), time.time() - t0)
    finally:
        del arr  # release the memmap handle so a resumed run can reopen it on Windows

    return _read_sidecar(spath)


def load_block(block_dir: Path | str, *, allow_partial: bool = False):
    """Open the stored block read-only. REFUSES a PARTIAL or incomplete block.

    This is the refusal that makes `--limit` safe: a smoke block reshapes and scores exactly like a
    whole one, so only an explicit status check can stop it being read as the arm's result.
    """
    block_dir = Path(block_dir)
    bpath, spath = block_dir / BLOCK_NAME, block_dir / SIDECAR_NAME
    side = _read_sidecar(spath)
    if side is None or not bpath.is_file():
        raise ExtractError(
            "no extracted block at %s (sidecar %s). Run scripts/positional_arm_extract.py first."
            % (bpath, "present" if side else "absent")
        )
    if side.get("partial_limit_run") and not allow_partial:
        raise ExtractError(
            "the block at %s is stamped PARTIAL (a --limit smoke run of %d rows). REFUSING to load "
            "it for scoring: a partial block reshapes and scores exactly like a whole one, so this "
            "status check is the only thing standing between a smoke run and a published number."
            % (bpath, side.get("n_rows_done"))
        )
    if not side.get("complete"):
        raise ExtractError(
            "the block at %s is incomplete (%d of %d rows). An interrupted pass leaves rows at "
            "zero, which would score as a clean null rather than as an error."
            % (bpath, side.get("n_rows_done"), side.get("geometry", {}).get("n_rows"))
        )
    arr = np.load(bpath, mmap_mode="r")
    lay = feature_layout(len(side["geometry"]["lags"]), n_bins=side["geometry"]["n_bins"])
    if arr.shape != (side["geometry"]["n_rows"], lay["_total"]):
        raise ExtractError(
            "block on disk is %r but its sidecar says %r" % (arr.shape,
                                                             (side["geometry"]["n_rows"],
                                                              lay["_total"]))
        )
    return arr, side, lay


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=0,
                    help="smoke run over the first N tiles; stamps the sidecar PARTIAL so the "
                         "scoring loader refuses it")
    ap.add_argument("--block-dir", default=str(DEFAULT_BLOCK_DIR))
    ap.add_argument("--batch", type=int, default=BATCH)
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args(argv)

    if not FEASIBILITY.is_file():
        print("REFUSED: feasibility artifact %s absent; run the probe first (validity before any "
              "extraction)." % FEASIBILITY, flush=True)
        return 2
    feas = json.loads(FEASIBILITY.read_text(encoding="utf-8"))
    if feas.get("verdict") != "FEASIBLE" or not feas["parity"]["passed"]:
        print("REFUSED: feasibility verdict is %r and parity passed=%r; extraction would produce "
              "features that no number could be read from."
              % (feas.get("verdict"), feas["parity"]["passed"]), flush=True)
        return 2

    lags = list(feas["lag_set"]["lags"])
    ck_sha = feas["checkpoint"]["sha256"]

    try:
        tiles, st = load_arm_cohort()
    except ProbeFailure as e:
        print("REFUSED: %s" % e, flush=True)
        return 2

    seqs = [t.seq for t in tiles]
    partial = bool(args.limit)
    if partial:
        seqs = seqs[:args.limit]

    print("cohort %d tiles (%s) | lags %r | block_dir %s"
          % (len(seqs), "PARTIAL smoke" if partial else "full", lags, args.block_dir), flush=True)

    def prog(done, total, el):
        if done % (args.batch * 20) == 0 or done == total:
            rate = done / el if el else 0.0
            print("   %d/%d  %.0fs  %.1f seq/s  eta %.1f min"
                  % (done, total, el, rate, (total - done) / rate / 60 if rate else 0.0),
                  flush=True)

    try:
        encode = make_atlas_encoder(device=args.device)
        side = extract(encode=encode, seqs=seqs, lags=lags, block_dir=Path(args.block_dir),
                       batch=args.batch, partial=partial, checkpoint_sha=ck_sha, progress=prog)
    except Exception as e:
        print("EXTRACT FAILED (%s): %s" % (type(e).__name__, e), flush=True)
        return 1

    g = side["geometry"]
    print()
    print("status      %s (complete=%s)" % (side["status"], side["complete"]))
    print("block       %s" % side["block"])
    print("geometry    %d rows x %d features %s" % (g["n_rows"], g["n_features"], g["dtype"]))
    print("layout      %s" % json.dumps(side["layout"], sort_keys=True))
    print("wall        %.1f s (%.1f min)" % (side["wall_seconds"], side["wall_seconds"] / 60.0))
    if partial:
        print("NOTE        stamped PARTIAL; the scoring loader will refuse this block.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
