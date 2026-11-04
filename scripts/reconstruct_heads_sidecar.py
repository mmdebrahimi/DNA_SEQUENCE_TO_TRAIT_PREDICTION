"""Rebuild the expression-heads sidecar from the block itself, after a crash destroyed it.

WHY THIS EXISTS. On 2026-11-03 the host lost AC power at 18:22:41 (Kernel-Power 41 + Event 6008)
with the heads extraction at 43,328 of 44,106 tiles. The `.npy` block survived intact -- 3.5 MB of
finite float64, cleanly truncated -- but the sidecar came back as 1,193 NUL bytes, because it was
written with a plain `write_text` whose data never reached the platter. `read_sidecar` maps
malformed to None and `extract` reads None as "no prior block", so 3.8 hours of correct compute
would have been recomputed from row zero.

The write path is now atomic and fsynced (`expression_heads_score.write_sidecar`), so this should
not recur. This tool recovers the blocks already damaged, and the host has now died twice in one
day, so it is committed rather than run as a one-off.

IT REFUSES UNLESS THREE INDEPENDENT SOURCES AGREE, and that is the whole design. Writing a sidecar
is asserting how many rows of the block are trustworthy; getting it wrong high means scoring
uninitialised zeros as measurements, which would pass every downstream gate and corrupt the arm's
number silently. So `n_done` must be corroborated by:

  1. THE BLOCK -- the index of the first all-zero row, with every row after it also all-zero. A row
     of ten float64 head outputs that are all exactly 0.0 does not occur in practice; a contiguous
     all-zero tail is the signature of "never written".
  2. THE HEARTBEAT -- written by the worker through a DIFFERENT code path (atomic tmp+replace, which
     is why it survived), recording `units_done` independently of the block.
  3. ARITHMETIC -- the count must be an exact multiple of the batch size, because the sidecar is
     written once per completed batch. A count that is not a batch multiple means the block was
     written by something other than this loop, and the reconstruction's premise is wrong.

Disagreement is a REFUSAL, never a reconciliation. Taking the minimum of two disagreeing sources
would look prudent and would be guessing.

Usage:
    python scripts/reconstruct_heads_sidecar.py [--block-dir DIR] [--job JOB] [--apply]
Dry-run by default: it prints what it would write and exits non-zero until `--apply`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from scripts.expression_heads_score import (
    BATCH,
    BLOCK_DIR,
    BLOCK_NAME,
    REDUCTIONS,
    ROUTE_TRUE,
    SIDECAR_NAME,
    column_layout,
    read_sidecar,
    write_sidecar,
)

REPO = Path(__file__).resolve().parents[1]


class ReconstructRefused(ValueError):
    """The three sources do not agree, so no sidecar can be written honestly."""


def filled_rows_from_block(path: Path) -> dict:
    """Where the block stops, read from the block alone.

    Returns the first all-zero row index plus the evidence a caller needs to judge it: whether the
    tail is contiguously zero, and whether the filled region is finite.
    """
    arr = np.load(path, mmap_mode="r")
    if arr.ndim != 2:
        raise ReconstructRefused("block is %dD, expected 2D" % arr.ndim)
    nonzero = np.asarray((arr != 0).any(axis=1))
    n_rows = int(arr.shape[0])
    first_zero = n_rows if nonzero.all() else int(np.argmin(nonzero))
    tail_all_zero = bool((~nonzero[first_zero:]).all()) if first_zero < n_rows else True
    filled_finite = bool(np.isfinite(np.asarray(arr[:first_zero])).all()) if first_zero else True
    return {
        "n_rows": n_rows,
        "n_cols": int(arr.shape[1]),
        "first_all_zero_row": first_zero,
        "n_nonzero_rows": int(nonzero.sum()),
        "tail_is_contiguously_zero": tail_all_zero,
        "filled_region_is_finite": filled_finite,
    }


def corroborate(block: dict, heartbeat_units: int | None, *, batch: int = BATCH) -> dict:
    """Agree-or-refuse across the three sources. Pure, so the refusals are testable."""
    n_done = block["first_all_zero_row"]

    if not block["tail_is_contiguously_zero"]:
        raise ReconstructRefused(
            "the block has a non-zero row AFTER its first all-zero row, so the all-zero row is a "
            "gap rather than the end of the written region. This block was not produced by a "
            "simple forward pass and its stopping point cannot be inferred.")
    if not block["filled_region_is_finite"]:
        raise ReconstructRefused(
            "the filled region contains non-finite values, so the surviving rows are not "
            "trustworthy as measurements regardless of how many there are.")
    if n_done <= 0:
        raise ReconstructRefused("the block holds no written rows; there is nothing to resume from")
    if n_done % batch != 0 and n_done != block["n_rows"]:
        raise ReconstructRefused(
            "n_done %d is not a multiple of the batch size %d and is not the full row count. The "
            "sidecar is written once per completed batch, so a non-multiple means this block was "
            "written by a different loop and this tool's premise does not hold."
            % (n_done, batch))
    if heartbeat_units is None:
        raise ReconstructRefused(
            "no heartbeat count to corroborate against. The block alone is ONE source; writing a "
            "sidecar off it would assert how many rows are trustworthy on unconfirmed evidence.")
    if int(heartbeat_units) != n_done:
        raise ReconstructRefused(
            "the block stops at %d but the heartbeat recorded %d. These are independent writers, "
            "so a disagreement means one of them is wrong and this tool cannot tell which. "
            "Refusing rather than taking the minimum, which would be a guess wearing the costume "
            "of caution." % (n_done, int(heartbeat_units)))

    return {
        "n_done": n_done,
        "n_rows": block["n_rows"],
        "status": "COMPLETE" if n_done >= block["n_rows"] else "PARTIAL",
        "agreed_by": ["block.first_all_zero_row", "heartbeat.units_done",
                      "batch_multiple(%d)" % batch],
        "batch": batch,
    }


def build_payload(agreed: dict, head_names: list, *, block: dict) -> dict:
    """The sidecar `extract` will read. Layout is DERIVED from the head names, never copied."""
    layout = column_layout(head_names)
    if layout["_total"] != block["n_cols"]:
        raise ReconstructRefused(
            "the layout derived from %d head names has %d columns but the block has %d. A sidecar "
            "whose layout disagrees with its block maps every head to the wrong column."
            % (len(head_names), layout["_total"], block["n_cols"]))
    return {
        "schema": "expression-heads-scalars-v1",
        "status": agreed["status"],
        "n_rows": agreed["n_rows"],
        "n_done": agreed["n_done"],
        "layout": layout,
        "heads": sorted(head_names),
        "reductions": list(REDUCTIONS),
        "route": ROUTE_TRUE,
        "batch": agreed["batch"],
        "reconstructed": True,
        "reconstruction_note": (
            "REBUILT by scripts/reconstruct_heads_sidecar.py after the original was destroyed by "
            "an unclean shutdown (Kernel-Power 41, 2026-11-03 18:22:41). n_done is corroborated by "
            + ", ".join(agreed["agreed_by"]) + ". The block's own data is untouched."),
    }


def head_names_from_repo() -> list:
    """The five shipped head names, read from the durable clone rather than hardcoded."""
    from dna_decode.glm import atlas_heads as AH
    from dna_decode.glm.atlas_embed import DEFAULT_PA_REPO

    base = AH.load_base_state_dict()
    heads = AH.load_heads(Path(DEFAULT_PA_REPO) / "trained_weights/expression_prediction",
                          base_state_dict=base)
    return [h.name for h in heads]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--block-dir", default=str(BLOCK_DIR))
    ap.add_argument("--job", default="expression_heads",
                    help="heartbeat job name to corroborate against")
    ap.add_argument("--apply", action="store_true",
                    help="write the sidecar; without it this is a dry run")
    a = ap.parse_args(argv)

    block_dir = Path(a.block_dir)
    block_path = block_dir / BLOCK_NAME
    if not block_path.is_file():
        print("REFUSED: no block at %s" % block_path)
        return 2

    existing = read_sidecar(block_dir)
    if existing is not None:
        print("REFUSED: %s already reads as valid JSON (status=%r, n_done=%r). This tool only "
              "repairs a DESTROYED sidecar; overwriting a readable one would discard the "
              "authority it carries."
              % (block_dir / SIDECAR_NAME, existing.get("status"), existing.get("n_done")))
        return 3

    from scripts.detached_launch import heartbeat_path, read_heartbeat

    hb = read_heartbeat(heartbeat_path(a.job))
    hb_units = hb.get("units_done") if hb else None

    block = filled_rows_from_block(block_path)
    print("BLOCK      %d rows x %d cols | first all-zero row %d | nonzero rows %d"
          % (block["n_rows"], block["n_cols"], block["first_all_zero_row"],
             block["n_nonzero_rows"]))
    print("           tail contiguously zero: %s | filled region finite: %s"
          % (block["tail_is_contiguously_zero"], block["filled_region_is_finite"]))
    print("HEARTBEAT  units_done %r (job %r)" % (hb_units, a.job))

    try:
        agreed = corroborate(block, hb_units)
    except ReconstructRefused as e:
        print()
        print("REFUSED: %s" % e)
        return 4

    print("AGREED     n_done %d / %d  -> status %s  (%s)"
          % (agreed["n_done"], agreed["n_rows"], agreed["status"],
             ", ".join(agreed["agreed_by"])))

    try:
        names = head_names_from_repo()
    except Exception as e:                                     # noqa: BLE001
        print("REFUSED: cannot read the head set (%s: %s), so the column layout cannot be derived"
              % (type(e).__name__, e))
        return 5
    print("HEADS      %d: %s" % (len(names), ", ".join(sorted(names))))

    try:
        payload = build_payload(agreed, names, block=block)
    except ReconstructRefused as e:
        print("REFUSED: %s" % e)
        return 6

    if not a.apply:
        print()
        print("DRY RUN -- would write %s" % (block_dir / SIDECAR_NAME))
        print("re-run with --apply to write it")
        return 1

    out = write_sidecar(block_dir, payload)
    back = read_sidecar(block_dir)
    if back is None or back.get("n_done") != agreed["n_done"]:
        print("REFUSED: wrote %s but it did not read back correctly" % out)
        return 7
    print()
    print("WROTE %s" % out)
    print("  status=%s n_done=%d/%d -- extract will resume at row %d"
          % (back["status"], back["n_done"], back["n_rows"], back["n_done"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
