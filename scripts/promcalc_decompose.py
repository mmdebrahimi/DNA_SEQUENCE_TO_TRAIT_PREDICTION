"""Stage B: run the Promoter Calculator over exported tiles and KEEP THE WHOLE DECOMPOSITION.

**RUN THIS WITH THE ISOLATED INTERPRETER, NOT `uv run`:**

    D:/tools/promcalc_env/Scripts/python.exe scripts/promcalc_decompose.py --workers 7

`promoter_calculator` pulls sklearn 1.9.1 against this project's pinned 1.8.0, so installing it in the
project venv would move sklearn under all 6,409 tests. This module therefore imports NOTHING from
`dna_decode` — it reads the TSV stage A wrote and writes a CSV stage C reads. Keep it that way.

**Why it exists at all.** The 2026-10-12 measurement reduced each tile to ONE scalar (`max_tx`) and
scored 0.3968, beating every learned arm. But the model returns ~71 candidate forward promoters per
150 bp tile, each decomposed into nine component free energies (UP, -35, spacer, -10, extended -10,
discriminator, ITR, binding, total) plus an initiation rate — roughly 710 numbers per tile, of which
the winning arm used one. That sweep's own script was never committed and its per-site values were
discarded, so the information could not be revisited without re-running. This script re-runs it once
and PERSISTS the per-site table, so every future mechanism arm is a cheap fit rather than another
74-minute sweep.

**Forward strand only, on measured grounds.** The published sweep found forward-only (0.3968) beats
best-of-both-strands (0.3800) at full n: the assay has a defined orientation, so reverse-strand
promoters should not drive the reporter and measurably add noise. Keeping only forward also halves the
table. A future arm wanting the reverse strand must re-sweep; that is a stated limit, not an oversight.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

DEFAULT_DIR = Path("D:/dna_decode_cache/promcalc")

#: The nine component energies plus the rate, in a FIXED order so the CSV header is stable.
DG_FIELDS = ("dG_total", "dG_10", "dG_35", "dG_ext10", "dG_spacer", "dG_disc", "dG_ITR",
             "dG_UP", "dG_bind")
SITE_FIELDS = ("tile_id", "tss", *DG_FIELDS, "Tx_rate", "hex35_start", "hex10_start", "spacer_len")


def _int_start(v) -> int:
    """`hex35_position` arrives as the string '[25, 31]'. Take the start, or -1 if unparseable."""
    try:
        return int(str(v).strip().lstrip("[").split(",")[0])
    except Exception:
        return -1


def _load_calculator():
    """Resolve the `Promoter_Calculator` CLASS across promotercalculator versions.

    The install name is `promotercalculator`; the import name is `promoter_calculator`. On top of that
    the CLASS moved relative to the package root: the version used for the 2026-10-12 run re-exported
    it from the top level, and **1.2.4 does not** -- its `__init__.py` exports only the lowercase
    `promoter_calculator` FUNCTION from `wrapper`, so the recorded `from promoter_calculator import
    Promoter_Calculator` raises ImportError on every tile. Measured 2026-10-22: that failure produced
    40 errored tiles and zero rows in 24 s.

    The class itself is UNCHANGED -- `wrapper.py` in 1.2.4 constructs it with the same
    `Promoter_Calculator(threads=, verbosity=)` and reads `output()['Forward_Predictions_per_TSS']`
    exactly as this worker does. So this is an import-PATH fix, not a model swap, which is what lets
    the frozen reconcile anchors in stage C stay the right acceptance test.

    Top level is tried FIRST so an older install keeps its recorded import path.
    """
    try:
        from promoter_calculator import Promoter_Calculator  # the recorded path
        return Promoter_Calculator
    except ImportError:
        from promoter_calculator.promoter_calculator import Promoter_Calculator  # >=1.2.4
        return Promoter_Calculator


def _site_field(d, name: str, default=None):
    """One per-TSS field, for EITHER of the two shapes this value has had.

    In the recorded run each per-TSS value was a plain dict and this worker called `d.get(...)`. In
    1.2.4 it is a `PromoCalcResults` object that defines `__getitem__`, `__setitem__` and `keys()` but
    **NOT `get`** -- so `d.get(...)` raises AttributeError per site and the tile is counted as failed.
    Subscript works on both shapes; a missing key is the only thing that needs the default.
    """
    try:
        return d[name]
    except (KeyError, AttributeError, TypeError):
        return default


def decompose_one(job: tuple[int, str]) -> tuple[int, list[list], str]:
    """Return `(tile_id, per-site rows, error)`. An error yields NO rows and is counted, never silently
    dropped: a tile that failed is not a tile with no promoters, and the two must stay distinguishable."""
    tile_id, seq = job
    try:
        Promoter_Calculator = _load_calculator()
        c = Promoter_Calculator()
        c.run(seq, TSS_range=[0, len(seq)])
        fwd = c.output().get("Forward_Predictions_per_TSS") or {}
        rows = []
        for tss, d in fwd.items():
            spacer = _site_field(d, "spacer") or ""
            rows.append([
                tile_id, int(tss),
                *[float(_site_field(d, f, "nan")) for f in DG_FIELDS],
                float(_site_field(d, "Tx_rate", "nan")),
                _int_start(_site_field(d, "hex35_position")),
                _int_start(_site_field(d, "hex10_position")),
                len(spacer),
            ])
        return tile_id, rows, ""
    except Exception as e:  # noqa: BLE001 - the error text is the artifact
        return tile_id, [], f"{type(e).__name__}: {e}"


def completed_tiles(in_dir: Path, *, smoke: bool = False) -> set[int]:
    """tile_ids already decomposed, across the main file AND every resume shard.

    **Why shards exist.** The first full run was killed by session teardown at 88.4% (38,974 of 44,106
    tiles, ~5 h of compute) and a plain re-run would have started from zero, because a single gzip
    stream cannot be appended to safely. Each resume writes its own `decomposition.part-NN.csv.gz`
    and stage C reads them all, so an interrupted sweep costs only the tiles still missing.

    **THE LAST TILE OF A TRUNCATED SHARD IS DROPPED, and that is the whole correctness argument.**
    A tile contributes ~73 rows, so a writer killed mid-tile leaves that tile with SOME of its rows -
    and "has at least one row" would read as "complete". Measured on the real interruption: tile 38,973
    is present, and the cut deflate stream replayed its sliding window into **2,199 copies of one
    mangled row** (a float emitted as `999999999993236`). Keeping it would have scored one tile on a
    partial, partly-garbage site list while every count said 44,106/44,106.

    This is the HDF5 embedding-cache trap on a different substrate: `populate(skip_existing=True)`
    skipped at gene level, so a mid-flush crash left a partial gene set that looked like a complete
    strain to the consumer. Same shape, same fix - verify at the unit the consumer aggregates over.
    """
    done: set[int] = set()
    for p in shard_paths(in_dir, smoke=smoke):
        seen: set[int] = set()
        cut = False
        try:
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
                for r in csv.DictReader(fh):
                    try:
                        seen.add(int(r["tile_id"]))
                    except (TypeError, ValueError):
                        continue
        except (EOFError, OSError):
            cut = True  # interrupted writer
        if cut and seen:
            seen.discard(max(seen))  # its row set may be partial; re-sweep it
        done |= seen
    return done


def shard_paths(in_dir: Path, *, smoke: bool = False) -> list[Path]:
    """The shards for one mode, EXPLICITLY - never a bare `decomposition*` glob.

    `decomposition*.csv.gz` also matches `decomposition_smoke.csv.gz`, so a full run silently merged
    the 20-tile smoke output into its own substrate. Harmless here only because the smoke tiles are a
    prefix of the real ones; the fix is to name the two sets rather than rely on that luck.
    """
    if smoke:
        return sorted(in_dir.glob("decomposition_smoke*.csv.gz"))
    out = [in_dir / "decomposition.csv.gz"]
    out += sorted(in_dir.glob("decomposition.part-*.csv.gz"))
    return [p for p in out if p.exists()]


def run(in_dir: Path | str = DEFAULT_DIR, *, workers: int = 7, smoke: bool = False,
        limit: int | None = None, resume: bool = False) -> dict:
    in_dir = Path(in_dir)
    tsv = in_dir / ("tiles_smoke.tsv" if smoke else "tiles.tsv")
    if not tsv.exists():
        raise SystemExit(f"{tsv} missing - run scripts/promcalc_export_tiles.py first (project env)")

    jobs: list[tuple[int, str]] = []
    with tsv.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            jobs.append((int(r["tile_id"]), r["seq"]))
    if limit:
        jobs = jobs[:limit]

    n_total = len(jobs)
    already: set[int] = set()
    if resume:
        already = completed_tiles(in_dir, smoke=smoke)
        jobs = [j for j in jobs if j[0] not in already]
        print(f"resume: {len(already)} tiles already done, {len(jobs)} remaining of {n_total}",
              flush=True)
        if not jobs:
            print("nothing to resume - every tile is already decomposed", flush=True)

    base = "decomposition_smoke" if smoke else "decomposition"
    if resume:
        n = 1
        while (in_dir / f"{base}.part-{n:02d}.csv.gz").exists():
            n += 1
        out = in_dir / f"{base}.part-{n:02d}.csv.gz"
    else:
        out = in_dir / f"{base}.csv.gz"
    t0 = time.time()
    n_ok = n_err = n_rows = 0
    errors: list[str] = []

    with gzip.open(out, "wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(SITE_FIELDS)
        if workers <= 1:
            it = (decompose_one(j) for j in jobs)
        else:
            pool = Pool(processes=workers)
            # chunksize must not exceed jobs/workers, or workers sit idle and the run is serial in
            # disguise. MEASURED: at chunksize=16 a 20-tile smoke put 16 tiles on one worker and
            # nothing on five of seven, took 0.95 min, and extrapolated to 2,095 min for the full
            # sweep -- 28x the published 74. The per-tile cost is 0.882 s (constructor only 0.035 s,
            # 3.8%), i.e. ~93 min at 7 workers. A small smoke measures pool spawn and load imbalance,
            # NOT throughput; size the chunk to the job count and project from the per-tile cost.
            it = pool.imap_unordered(decompose_one, jobs,
                                     chunksize=max(1, min(8, len(jobs) // max(workers, 1) or 1)))
        for tile_id, rows, err in it:
            if err:
                n_err += 1
                if len(errors) < 20:
                    errors.append(f"tile {tile_id}: {err}")
                continue
            n_ok += 1
            n_rows += len(rows)
            w.writerows(rows)
            if n_ok % 2000 == 0:
                el = time.time() - t0
                print(f"  {n_ok}/{len(jobs)} tiles | {n_rows} sites | {el/60:.1f} min "
                      f"| eta {(el/max(n_ok,1))*(len(jobs)-n_ok)/60:.1f} min", flush=True)
        if workers > 1:
            pool.close()
            pool.join()

    meta = {
        "attempted": len(jobs), "ok": n_ok, "err": n_err,
        "resumed": bool(resume), "already_done_before_this_run": len(already),
        "tiles_in_tsv": n_total,
        # the COVERAGE that matters is over the whole substrate, not over this resume slice
        "total_complete": len(already) + n_ok,
        "substrate_coverage": round((len(already) + n_ok) / n_total, 6) if n_total else None,
        "n_site_rows": n_rows,
        "sites_per_tile_mean": round(n_rows / n_ok, 2) if n_ok else None,
        "minutes": round((time.time() - t0) / 60, 2),
        "coverage": round(n_ok / len(jobs), 6) if jobs else None,
        "strand": "forward_only",
        "fields": list(SITE_FIELDS),
        "out": str(out),
        "errors_sample": errors,
        "workers": workers,
        "smoke": smoke,
    }
    (in_dir / ("decompose_smoke.json" if smoke else "decompose.json")).write_text(
        json.dumps(meta, indent=1), encoding="utf-8")
    return meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", default=str(DEFAULT_DIR))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--smoke", action="store_true", help="read tiles_smoke.tsv, write *_smoke.*")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--resume", action="store_true",
                    help="skip tiles already present in any decomposition*.csv.gz and write a new shard")
    a = ap.parse_args()
    m = run(a.in_dir, workers=a.workers, smoke=a.smoke, limit=a.limit, resume=a.resume)
    print(json.dumps(m, indent=1))
    # A sweep that lost tiles, or that did not cover the substrate, is not a usable arm substrate.
    if m["err"]:
        return 1
    return 0 if (m["substrate_coverage"] or 0) >= 1.0 else 2


if __name__ == "__main__":
    sys.exit(main())
