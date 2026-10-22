"""Validate stage B's decomposition table BEFORE stage C scores it.

    python scripts/promcalc_validate_output.py data/promcalc
    python scripts/promcalc_validate_output.py data/promcalc --expect-partial   # mid-run check

Exit 0 = scoreable · 2 = a defect that must not reach stage C · 3 = indeterminate.

**THE DEFECT THIS EXISTS FOR: a table of NaN that looks perfect.** `decompose_one` reads each per-site
field through `_site_field`, which returns its default on KeyError / AttributeError / TypeError. That is
correct for the two 2026-10-22 version breaks, but it means a FUTURE rename of a field produces a
~140 MB table where every value is `nan`, with **0 errors reported, full 44,106 coverage, and
`sites_per_tile_mean` still 73.0**. Every count says success. Only stage C's frozen reconcile gate would
catch it, ~12 h later, and it would surface as `INDETERMINATE_RECONCILE_FAILED` -- indistinguishable
from genuine model drift, which is the one diagnosis you most want to keep clean.

**Why this is a separate validator and NOT a counter inside stage B.** Three reasons, in order:
  1. Stage B is the LIVE module during a ~12 h run, and a watchdog relaunch picks up new code. A guard
     that protects future runs is not worth editing the running one for.
  2. An in-loop counter cannot protect data already written. This can be run on any existing table.
  3. Stage C must not be handed a bad table in the first place, so the check belongs at the boundary
     between the two stages, which is a file.

**Run it IDLE.** It decompresses the whole table. Doing that during a live run is the recorded trap that
turned a ~1 h job into a 5.0 h projection -- observing the sweep slowed the sweep, and the slowed rate
was then quoted as the sweep's rate. Mid-run use is supported via `--expect-partial`, which downgrades
the coverage shortfall from a defect to an expected state so the NaN and shape checks can still be read;
it still decompresses, so use it once, not in a loop.

Stdlib only, so it runs under either interpreter.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import Counter
from pathlib import Path

EXPECT_TILES = 44106

#: The nine component energies plus the rate -- every field whose value is a float and could silently
#: become nan. `tile_id` / `tss` / the two start positions / `spacer_len` are ints and excluded.
FLOAT_FIELDS = ("dG_total", "dG_10", "dG_35", "dG_ext10", "dG_spacer", "dG_disc", "dG_ITR",
                "dG_UP", "dG_bind", "Tx_rate")

EXPECTED_HEADER = ("tile_id", "tss", "dG_total", "dG_10", "dG_35", "dG_ext10", "dG_spacer",
                   "dG_disc", "dG_ITR", "dG_UP", "dG_bind", "Tx_rate", "hex35_start",
                   "hex10_start", "spacer_len")

#: A field whose values are MORE than this fraction non-finite is treated as a defect, not noise.
#: Set deliberately low: the measured rate on a healthy table is 0.0000, so anything approaching this
#: is categorically different from "a few pathological tiles".
NAN_FRACTION_LIMIT = 0.01


def shard_paths(in_dir: Path) -> list[Path]:
    """The full-run shards, explicitly -- never a bare `decomposition*` glob.

    That glob also matches `decomposition_smoke*.csv.gz`, and merging a 40-tile smoke into a real
    substrate is a recorded near-miss in this chain. Mirrors `promcalc_decompose.shard_paths`.
    """
    out = [in_dir / "decomposition.csv.gz"]
    out += sorted(in_dir.glob("decomposition.part-*.csv.gz"))
    return [p for p in out if p.exists()]


def scan(in_dir: Path) -> dict:
    """Read every shard once and return the facts. Tolerates a truncated final shard."""
    paths = shard_paths(in_dir)
    if not paths:
        return {"error": f"no decomposition shards in {in_dir}"}

    headers: list[tuple[str, ...]] = []
    nonfinite = Counter()
    present = Counter()
    rows_per_tile = Counter()
    truncated: list[str] = []
    unparsed = 0
    total_rows = 0

    for p in paths:
        try:
            with gzip.open(p, "rt", encoding="utf-8", errors="replace", newline="") as fh:
                rd = csv.DictReader(fh)
                headers.append(tuple(rd.fieldnames or ()))
                for row in rd:
                    total_rows += 1
                    try:
                        rows_per_tile[int(row["tile_id"])] += 1
                    except (TypeError, ValueError, KeyError):
                        unparsed += 1
                        continue
                    for f in FLOAT_FIELDS:
                        raw = row.get(f)
                        present[f] += 1
                        try:
                            if not math.isfinite(float(raw)):
                                nonfinite[f] += 1
                        except (TypeError, ValueError):
                            nonfinite[f] += 1
        except (EOFError, OSError):
            truncated.append(p.name)

    return {
        "shards": [p.name for p in paths],
        "truncated_shards": truncated,
        "headers": headers,
        "n_rows": total_rows,
        "n_tiles": len(rows_per_tile),
        "unparsed_tile_ids": unparsed,
        "nonfinite": dict(nonfinite),
        "present": dict(present),
        "rows_per_tile_hist": Counter(rows_per_tile.values()).most_common(5),
    }


def judge(facts: dict, *, expect_partial: bool) -> tuple[str, list[str], list[str]]:
    """`(verdict, defects, notes)`. Pure, so the thresholds are testable without a table."""
    defects: list[str] = []
    notes: list[str] = []

    if facts.get("error"):
        return "INDETERMINATE", [facts["error"]], notes

    # --- header ---------------------------------------------------------------------------------
    for h in facts["headers"]:
        if tuple(h) != EXPECTED_HEADER:
            defects.append(f"header mismatch: {h} != {EXPECTED_HEADER}")
            break

    # --- THE defect this exists for -------------------------------------------------------------
    worst = 0.0
    for f in FLOAT_FIELDS:
        n = facts["present"].get(f, 0)
        if not n:
            defects.append(f"field {f} never appeared in any row")
            continue
        frac = facts["nonfinite"].get(f, 0) / n
        worst = max(worst, frac)
        if frac > NAN_FRACTION_LIMIT:
            defects.append(
                f"field {f} is {frac:.1%} non-finite ({facts['nonfinite'].get(f,0)}/{n}), over the "
                f"{NAN_FRACTION_LIMIT:.0%} limit. A field-rename break produces exactly this while "
                f"every count still reports success."
            )
    notes.append(f"worst non-finite fraction across the {len(FLOAT_FIELDS)} float fields: {worst:.4%}")

    # --- shape ----------------------------------------------------------------------------------
    if facts["n_tiles"]:
        spt = facts["n_rows"] / facts["n_tiles"]
        notes.append(f"sites/tile mean {spt:.2f} (published ~73)")
        if not 60.0 <= spt <= 85.0:
            defects.append(f"sites/tile mean {spt:.2f} is outside 60-85; the model's output shape moved")
    if facts["unparsed_tile_ids"]:
        notes.append(f"{facts['unparsed_tile_ids']} row(s) had an unparseable tile_id")

    # --- coverage -------------------------------------------------------------------------------
    cov = facts["n_tiles"] / EXPECT_TILES
    notes.append(f"coverage {cov:.4f} ({facts['n_tiles']}/{EXPECT_TILES})")
    if facts["n_tiles"] < EXPECT_TILES:
        msg = f"coverage short by {EXPECT_TILES - facts['n_tiles']} tiles"
        if expect_partial:
            notes.append(msg + " (expected: --expect-partial)")
        else:
            defects.append(msg + "; stage C would score a partial substrate")

    if facts["truncated_shards"]:
        notes.append(
            f"truncated shard(s) {facts['truncated_shards']} - normal for a live or killed run; "
            f"stage B's resume drops such a shard's LAST tile because a tile killed mid-write holds "
            f"only some of its ~73 rows"
        )

    return ("DEFECT" if defects else "SCOREABLE"), defects, notes


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("in_dir")
    ap.add_argument("--expect-partial", action="store_true",
                    help="mid-run check: a coverage shortfall is a NOTE, not a defect. Still "
                         "decompresses the whole table, so use it once rather than in a loop.")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)

    facts = scan(Path(a.in_dir))
    verdict, defects, notes = judge(facts, expect_partial=a.expect_partial)

    if a.json:
        print(json.dumps({"verdict": verdict, "defects": defects, "notes": notes, "facts": facts},
                         indent=1, default=str))
    else:
        print(f"verdict: {verdict}")
        for n in notes:
            print(f"  note   : {n}")
        for d in defects:
            print(f"  DEFECT : {d}")

    return {"SCOREABLE": 0, "DEFECT": 2, "INDETERMINATE": 3}[verdict]


if __name__ == "__main__":
    sys.exit(main())
