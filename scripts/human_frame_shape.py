"""The human lentiMPRA label's SHAPE, and whether a human reproducibility curve exists at all.

WHAT THIS IS FOR
----------------
The eukaryotic arm reports a conv-over-GC gain of +0.2471 (HepG2) and +0.3016 (K562) against
+0.0122 on bacterial tiles, and the README has called that "roughly 20-25x more on human
regulatory DNA". The decisive question for whether that ratio is a FRAME artifact is whether the
human label is detection-dominated the way the bacterial one is.

TWO SEPARATE THINGS, NEVER MERGED INTO ONE HEADLINE
---------------------------------------------------
(a) `shape_summary` on each cell's activity and on the bacterial expression -- modal-bin share
    and distinct-value count only. These are the project's shipped `assay_degeneracy`
    quantities and are not determined by the design.
(b) `reproducibility_curve` on the human label IF AND ONLY IF the deposit carries per-replicate
    columns. `dna_decode.glm.lentimpra.LentiRecord` carries `name`, `seq`, `activity`, `fold`
    and NOTHING replicate-level, so this script scans the raw table headers rather than trusting
    the loader's record shape, and records `human_reproducibility_curve: null` WITH THE COLUMNS
    IT FOUND when none qualify.

There is no fallback, no single-measurement proxy, and no inherited bacterial curve. Both
eukaryotic artifacts currently carry the BACTERIAL 0.941 ceiling as context, and that
substitution must not be repeated here: a human ceiling is either derived from human replicates
or reported as `null` with a reason.

This script does NOT invent a human `active` flag. lentiMPRA publishes no such call, so the
comparison is on distributional shape, never on a threshold this plan chose.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.glm.active_frame import (  # noqa: E402
    reproducibility_curve,
    shape_summary,
    spearman_brown_ceiling,
)

DEFAULT_DEST = "data/mpra_human"
CELLS = ("HepG2", "K562", "WTC11")
BACTERIAL_CACHE = "data/mpra"

#: Header patterns that would indicate PER-REPLICATE data rather than a summary.
#: `fold` is excluded deliberately: in this deposit it is the authors' published 10-fold
#: cross-validation assignment, not a measurement replicate, and treating it as one would
#: manufacture a ceiling out of a CV index.
REPLICATE_PATTERNS = (
    re.compile(r"(?:^|[_.])rep(?:licate)?[_.]?\d+$", re.I),
    re.compile(r"(?:^|[_.])r\d+$", re.I),
    re.compile(r"rna[_.].*[_.]\d+$", re.I),
    re.compile(r"dna[_.].*[_.]\d+$", re.I),
)
NOT_A_REPLICATE = re.compile(r"fold", re.I)

EXIT_OK = 0
EXIT_FETCH = 3
EXIT_MD5 = 4
EXIT_JOIN_RATE = 5
EXIT_BACTERIAL_ABSENT = 6


class HumanShapeError(RuntimeError):
    """A refusal; never a partial measurement."""


# --------------------------------------------------------------------------------------
# the fetch, delegated to the existing MD5-verifying fetcher
# --------------------------------------------------------------------------------------

def fetch_argv(dest: str) -> list[str]:
    """The exact argv handed to the fetcher. Pinned by test so it cannot drift.

    `--select tables` is the ~74 MB four-table subset. `--select all` is 7.68 GB and `C:` has
    33.8 GB free at 93.4% full, so the flag is load-bearing rather than cosmetic. `--dest`
    exists precisely so the detached `D:/dna_decode_cache/mpra_human` default is bypassed.
    """
    return [sys.executable, "scripts/fetch_lentimpra.py",
            "--select", "tables", "--dest", dest, "--yes"]


def run_fetch(dest: str, *, timeout: int = 1800) -> dict:
    """Invoke the fetcher. Its MD5 verification is the integrity guard and is NOT reimplemented."""
    argv = fetch_argv(dest)
    t0 = time.time()
    proc = subprocess.run(argv, cwd=REPO, capture_output=True, text=True, timeout=timeout)
    blob = (proc.stdout or "") + (proc.stderr or "")
    md5_failed = bool(re.search(r"md5|digest|checksum", blob, re.I)) and (
        bool(re.search(r"mismatch|fail|differ|bad", blob, re.I)))
    return {
        "argv": argv,
        "returncode": proc.returncode,
        "ok": proc.returncode == 0,
        "md5_failure_detected": md5_failed,
        "wall_seconds": round(time.time() - t0, 1),
        "tail": blob.strip().splitlines()[-12:],
    }


# --------------------------------------------------------------------------------------
# does the deposit carry per-replicate columns at all?
# --------------------------------------------------------------------------------------

def scan_replicate_columns(data_dir: str | Path) -> dict:
    """Enumerate every table header in the deposit and classify it.

    Scans the RAW headers rather than the loader's record shape, because `LentiRecord` carries
    no replicate fields by construction and absence there says nothing about the deposit.
    """
    root = Path(data_dir)
    if not root.exists():
        return {"data_dir": str(root), "exists": False, "files": [],
                "replicate_like_columns": [], "has_replicates": False,
                "reason": "the deposit directory does not exist"}
    files = []
    rep_cols: list[str] = []
    excluded_as_fold: list[str] = []
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        entry = {"file": str(p.relative_to(root)), "bytes": p.stat().st_size, "columns": None,
                 "note": None}
        if p.suffix.lower() in (".txt", ".tsv", ".csv"):
            try:
                head = p.open("r", errors="replace").readline().rstrip("\n")
                sep = "\t" if "\t" in head else ","
                cols = [c.strip().strip('"') for c in head.split(sep)]
                entry["columns"] = cols[:60]
                for c in cols:
                    looks_replicate = any(pat.search(c) for pat in REPLICATE_PATTERNS)
                    if NOT_A_REPLICATE.search(c):
                        # record ONLY where the exclusion actually changed the outcome: a
                        # filter that removes nothing is not a control, and this project has a
                        # recorded instance of exactly that shipping as a safeguard
                        if looks_replicate:
                            excluded_as_fold.append(f"{entry['file']}::{c}")
                        continue
                    if looks_replicate:
                        rep_cols.append(f"{entry['file']}::{c}")
            except Exception as e:  # noqa: BLE001
                entry["note"] = f"unreadable header: {e!r}"
        else:
            entry["note"] = "not a delimited text table; header not scanned"
        files.append(entry)
    return {
        "data_dir": str(root),
        "exists": True,
        "n_files": len(files),
        "files": files,
        "replicate_like_columns": rep_cols,
        "has_replicates": bool(rep_cols),
        "patterns_used": [p.pattern for p in REPLICATE_PATTERNS],
        "excluded_pattern": NOT_A_REPLICATE.pattern,
        "excluded_as_fold": excluded_as_fold,
        "n_excluded_as_fold": len(excluded_as_fold),
        "fold_exclusion_was_inert": not excluded_as_fold,
        "what_inert_means_here": (
            "the fold exclusion removed NOTHING on this deposit, so it is a guard against a "
            "column shape that did not occur rather than an active control. Reported because a "
            "filter that removes nothing is not a control, and presenting it as one would be "
            "the vacuous-control pattern this project has shipped once already"
        ) if not excluded_as_fold else
        "the fold exclusion actively removed the listed columns from the replicate set",
        "why_fold_is_excluded": "in this deposit `fold` is the authors' published 10-fold "
                                "cross-validation assignment, not a measurement replicate; "
                                "treating it as one would manufacture a ceiling from a CV index",
        "reason": None if rep_cols else
                  "no header matched a per-replicate pattern, so no human reproducibility "
                  "curve and no human ceiling are derivable from this deposit",
    }


# --------------------------------------------------------------------------------------
# per-cell shape
# --------------------------------------------------------------------------------------

def cell_shape(cell: str, data_dir: str | Path, *, min_join_rate: float = 0.95) -> dict:
    """Join rate plus distributional shape for one cell line, or a recorded refusal."""
    from dna_decode.glm.lentimpra import LentiMpraDataError, load_cell

    try:
        recs, stats = load_cell(cell, data_dir, min_join_rate=min_join_rate)
    except LentiMpraDataError as e:
        return {"cell": cell, "loaded": False, "refusal": str(e),
                "refusal_is_expected_for_WTC11": cell == "WTC11",
                "note": "recorded, not worked around"}
    except Exception as e:  # noqa: BLE001
        return {"cell": cell, "loaded": False, "refusal": repr(e), "unexpected": True}
    acts = [float(r.activity) for r in recs]
    sd = {k: v for k, v in vars(stats).items()} if hasattr(stats, "__dict__") else {}
    return {
        "cell": cell,
        "loaded": True,
        "n_records": len(recs),
        "join_rate": sd.get("join_rate"),
        "loader_stats": {k: v for k, v in sd.items() if isinstance(v, (int, float, str, bool))},
        "shape": shape_summary(acts),
        "activity_is": "log2(RNA/DNA); goes negative, so it is NOT a positive-only intensity",
        "no_active_flag_invented": True,
    }


def bacterial_shape(cache: str = BACTERIAL_CACHE) -> dict:
    """Bacterial label shape on BOTH scales, because the raw one is not comparable.

    `modal_bin_share` uses equal-WIDTH bins, so it is not scale-invariant. `expn_med` is a raw
    ratio with a heavy right tail while human activity is already `log2(RNA/DNA)`, so comparing
    their raw modal-bin shares compares a skew against a log. Both scales are reported and the
    comparable one is named, rather than leaving a reader to pick.

    `log1p` is not an arbitrary choice: it is the transform `scripts/glm_encoder_gate.ridge_gc`
    and the conv actually fit against, so it is the scale on which the models see this label.
    """
    import math

    from dna_decode.glm.genomewide import load_peak_tiles

    loaded, _st = load_peak_tiles(cache, allow_download=False)
    tiles = [t for t in loaded if t.category == "tile" and len(t.seq) == 150]
    raw = [float(t.expression) for t in tiles]
    logged = [math.log1p(x) for x in raw]
    s_raw, s_log = shape_summary(raw), shape_summary(logged)
    return {
        "n_tiles": len(tiles),
        "shape": s_raw,
        "shape_on_model_scale": s_log,
        "model_scale_is": "log1p(expn_med) -- the target ridge_gc and the conv are fit against",
        "comparability": {
            "compare_this_row_to_the_human_cells": "shape_on_model_scale",
            "why_not_the_raw_row": "modal_bin_share uses equal-WIDTH bins and so is NOT "
                                   "scale-invariant. expn_med is a raw ratio on [0, 158.3] with "
                                   "median 0.4338; human activity is log2(RNA/DNA) on about "
                                   "[-3, +7]. Comparing those raw modal-bin shares compares a "
                                   "heavy right skew against an already-logged variable",
            "distinct_value_share_is_scale_invariant": True,
            "note_on_distinct_value_share": "a monotone transform preserves distinctness, so "
                                            "that measure is already comparable across the two "
                                            "and log1p does not move it",
        },
    }


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dest", default=DEFAULT_DEST)
    ap.add_argument("--bacterial-cache", default=BACTERIAL_CACHE)
    ap.add_argument("--min-join-rate", type=float, default=0.95)
    ap.add_argument("--skip-fetch", action="store_true",
                    help="assume the deposit is already at --dest")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    t0 = time.time()
    fetch = {"skipped": True} if args.skip_fetch else run_fetch(args.dest)
    if not args.skip_fetch:
        print(f"fetch returncode {fetch['returncode']} in {fetch['wall_seconds']}s")
        for line in fetch["tail"]:
            print(f"   | {line}")
        if not fetch["ok"]:
            code = EXIT_MD5 if fetch["md5_failure_detected"] else EXIT_FETCH
            label = ("INDETERMINATE_MD5_VERIFICATION_FAILED" if code == EXIT_MD5
                     else "INDETERMINATE_FETCH_FAILED")
            print(f"{label}: the deposit was not obtained, so NO shape is reported")
            Path(args.out or f"wiki/human_frame_shape_{date.today()}.json").write_text(
                json.dumps({"schema": "human-frame-shape-v1", "date": str(date.today()),
                            "status": label, "fetch": fetch,
                            "human_reproducibility_curve": None,
                            "human_ceiling": None,
                            "partial_measurement_withheld": True}, indent=2),
                encoding="utf-8")
            return code

    scan = scan_replicate_columns(args.dest)
    print(f"deposit scan: {scan.get('n_files', 0)} file(s), "
          f"{len(scan['replicate_like_columns'])} replicate-like column(s), "
          f"has_replicates={scan['has_replicates']}")

    cells = {c: cell_shape(c, args.dest, min_join_rate=args.min_join_rate) for c in CELLS}
    for c, r in cells.items():
        if r["loaded"]:
            print(f"  {c:<7s} n={r['n_records']:<7d} join_rate={r['join_rate']} "
                  f"modal_bin_share={r['shape']['modal_bin_share']} "
                  f"distinct={r['shape']['n_distinct_values']}")
        else:
            print(f"  {c:<7s} REFUSED: {str(r['refusal'])[:110]}")

    usable = [c for c, r in cells.items() if r["loaded"]]
    if not usable:
        print("INDETERMINATE_JOIN_RATE_BELOW_THRESHOLD: no cell line loaded, so NO shape is "
              "reported")
        Path(args.out or f"wiki/human_frame_shape_{date.today()}.json").write_text(
            json.dumps({"schema": "human-frame-shape-v1", "date": str(date.today()),
                        "status": "INDETERMINATE_JOIN_RATE_BELOW_THRESHOLD",
                        "fetch": fetch, "deposit_scan": scan, "cells": cells,
                        "human_reproducibility_curve": None, "human_ceiling": None,
                        "partial_measurement_withheld": True}, indent=2),
            encoding="utf-8")
        return EXIT_JOIN_RATE

    try:
        bact = bacterial_shape(args.bacterial_cache)
    except Exception as e:  # noqa: BLE001
        print(f"INDETERMINATE_BACTERIAL_SUBSTRATE_ABSENT: {e}")
        return EXIT_BACTERIAL_ABSENT

    out = {
        "schema": "human-frame-shape-v1",
        "date": str(date.today()),
        "status": "MEASURED",
        "fetch": fetch,
        "deposit_scan": scan,
        "cells": cells,
        "bacterial_comparison": bact,
        "human_reproducibility_curve": None,
        "human_reproducibility_curve_reason": scan["reason"],
        "human_ceiling": None,
        "human_ceiling_reason": (
            None if scan["has_replicates"] else
            "no per-replicate columns were found in the deposit, so no Spearman-Brown step-up "
            "is derivable. Reported as null rather than filled with the BACTERIAL 0.941, which "
            "both wiki/eukaryotic_arm_result_*.json currently carry as context"
        ),
        "the_two_halves_are_reported_separately": {
            "a_distributional_shape": "modal-bin share and distinct-value count, per cell and "
                                      "for the bacterial label, directly comparable",
            "b_reproducibility_curve": "available only with per-replicate columns; see "
                                       "human_reproducibility_curve_reason",
            "why_separate": "merging them into one headline is how the retracted "
                            "single-threshold statistic travelled; a shape comparison and a "
                            "reliability comparison answer different questions",
        },
        "no_active_flag_invented": True,
        "wall_seconds": round(time.time() - t0, 1),
    }
    if scan["has_replicates"]:
        out["human_reproducibility_curve_reason"] = (
            "per-replicate columns WERE found: " + ", ".join(scan["replicate_like_columns"][:8])
            + ". Computing the curve needs a per-cell join of those columns, which this script "
              "does not yet implement -- reported as a named gap rather than approximated"
        )

    dest = Path(args.out) if args.out else Path(f"wiki/human_frame_shape_{date.today()}.json")
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print("\n" + "=" * 78)
    print("SHAPE SIDE BY SIDE   modal-bin share / distinct-value share / n")
    print(f"  {'bacterial RAW expn_med':<26s} {bact['shape']['modal_bin_share']:.4f}  "
          f"{bact['shape']['distinct_value_share']:.4f}  {bact['n_tiles']}"
          f"   <- NOT comparable (raw ratio, heavy right tail)")
    sl = bact["shape_on_model_scale"]
    print(f"  {'bacterial log1p (model)':<26s} {sl['modal_bin_share']:.4f}  "
          f"{sl['distinct_value_share']:.4f}  {bact['n_tiles']}"
          f"   <- COMPARE THIS ROW")
    for c in usable:
        s = cells[c]["shape"]
        print(f"  {'human ' + c + ' log2(RNA/DNA)':<26s} {s['modal_bin_share']:.4f}  "
              f"{s['distinct_value_share']:.4f}  {s['n']}")
    print("=" * 78)
    print(f"human reproducibility curve: {out['human_reproducibility_curve']} "
          f"({out['human_reproducibility_curve_reason']})")
    print(f"human ceiling: {out['human_ceiling']}")
    print(f"wrote {dest}")
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
