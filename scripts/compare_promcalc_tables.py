#!/usr/bin/env python
"""Compare two promcalc decomposition tables cell-by-cell, keyed on (tile_id, tss).

WHY THIS EXISTS. server1 recomputed the 44,106-tile decomposition while an equivalent table already
sat on server2 (truncated, then rebuilt from its own part-files). Keeping the recompute was a
deliberate call: the server2 copy is 16 days old and its provenance cannot be fully established, and
this project has ALREADY measured that the fragment-pyramid arm is bit-identical across hosts at
3 seeds / N=4,000 and NOT identical at N=321,123 x 10 seeds (one null seed moved 0.0328, larger than
the entire 0.02 validity bar, while the median moved only 0.0041). So "two hosts agree" is an
empirical question here, not an assumption.

WHAT IT REFUSES TO DO. It does not declare agreement from row counts. Counts are exactly what the
`decompose.json` manifest got right while the file it described was truncated -- it reported
coverage 1.0 / 44,106 tiles / sites_per_tile 73.0 / zero errors on a gzip that cannot be read to the
end. Agreement here means VALUES agree, per cell, within a stated tolerance.

    uv run python scripts/compare_promcalc_tables.py A.csv.gz B.csv.gz [--tol 1e-9]

Exit 0 = tables agree within tolerance. 1 = they differ (the report says where and by how much).
2 = a table could not be read to the end (truncation), which is NOT a disagreement and is reported
as its own state rather than folded into one.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys


KEY = ("tile_id", "tss")


def load(path):
    """Stream a table into {(tile_id, tss): [floats]}. Returns (rows, truncated, header)."""
    out, truncated = {}, False
    with gzip.open(path, "rt", newline="", encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh)
        header = next(r)
        ki = [header.index(k) for k in KEY]
        num = [i for i, h in enumerate(header) if i not in ki]
        while True:
            try:
                row = next(r)
            except StopIteration:
                break
            except (EOFError, OSError, gzip.BadGzipFile):
                truncated = True
                break
            if len(row) < len(header):
                continue
            key = tuple(row[i] for i in ki)
            vals = []
            for i in num:
                try:
                    vals.append(float(row[i]))
                except ValueError:
                    vals.append(math.nan)
            out[key] = vals
    return out, truncated, header


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--tol", type=float, default=1e-9,
                    help="absolute tolerance per cell (default 1e-9: bit-identical in practice)")
    ap.add_argument("--max-report", type=int, default=10)
    a = ap.parse_args(argv)

    A, a_trunc, ah = load(a.a)
    B, b_trunc, bh = load(a.b)
    print("A: %s  rows=%d truncated=%s" % (a.a, len(A), a_trunc))
    print("B: %s  rows=%d truncated=%s" % (a.b, len(B), b_trunc))

    if ah != bh:
        print("HEADER MISMATCH -- not comparable.")
        print("  A:", ah)
        print("  B:", bh)
        return 1
    if a_trunc or b_trunc:
        print("REFUSING A VERDICT: at least one table is TRUNCATED, so any agreement would be "
              "scoped to the readable prefix and would read as a clean comparison. Repair first.")
        return 2

    only_a, only_b = set(A) - set(B), set(B) - set(A)
    shared = set(A) & set(B)
    print("shared keys=%d  only_in_A=%d  only_in_B=%d" % (len(shared), len(only_a), len(only_b)))

    fields = [h for h in ah if h not in KEY]
    worst = {f: (0.0, None) for f in fields}
    n_diff, examples = 0, []
    for k in shared:
        va, vb = A[k], B[k]
        row_bad = False
        for j, f in enumerate(fields):
            x, y = va[j], vb[j]
            if math.isnan(x) and math.isnan(y):
                continue
            d = abs(x - y) if not (math.isnan(x) or math.isnan(y)) else float("inf")
            if d > worst[f][0]:
                worst[f] = (d, k)
            if d > a.tol:
                row_bad = True
        if row_bad:
            n_diff += 1
            if len(examples) < a.max_report:
                examples.append(k)

    print()
    print("%-12s %14s   %s" % ("field", "max_abs_delta", "at (tile_id, tss)"))
    for f in fields:
        d, k = worst[f]
        print("%-12s %14.6g   %s" % (f, d, k))
    print()
    print("rows differing beyond tol=%g: %d of %d shared" % (a.tol, n_diff, len(shared)))

    if only_a or only_b or n_diff:
        print("VERDICT: TABLES DIFFER")
        if examples:
            print("  first differing keys:", examples)
        return 1
    print("VERDICT: TABLES AGREE within %g on every shared cell, with identical key sets" % a.tol)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
