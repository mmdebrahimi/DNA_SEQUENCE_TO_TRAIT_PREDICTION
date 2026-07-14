"""Demo: the generalized position-novelty self-awareness flag across cells (2026-07-13).

Applies the cell-agnostic flag (dna_decode.eval.position_novelty) to a genotype against EACH named
cell's committed catalog, showing it generalizes beyond HIV NNRTI (where it was validated,
FLAG_RECOVERS_BLINDSPOT lift 3.98) to SARS-CoV-2 Mpro + fungal ERG11. The flag says "the catalog call
may be incomplete here" — never R/S. Read-only over committed catalogs; frozen surface untouched.

  uv run python scripts/position_novelty_demo.py --cell hiv-nnrti-rt --observed K103R,V179D
  uv run python scripts/position_novelty_demo.py --all   # summarize every known cell's catalog coverage
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.eval import position_novelty as PN  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cell", default=None, help="a known cell (hiv-nnrti-rt / sarscov2-mpro / fungal-<drug>-erg11)")
    ap.add_argument("--observed", default="", help="comma-separated substitutions, e.g. K103R,V179D")
    ap.add_argument("--all", action="store_true", help="summarize catalog coverage for every known cell")
    a = ap.parse_args(argv)

    if a.all:
        print(f"{'cell':<28} {'n_drms':>7} {'n_positions':>12}")
        for cell in PN.KNOWN_CELLS:
            try:
                drms = PN.catalog_drms_for(cell)
                print(f"{cell:<28} {len(drms):>7} {len(PN.catalog_positions(drms)):>12}")
            except Exception as e:  # noqa: BLE001
                print(f"{cell:<28} {'ERR':>7}  {e}")
        return 0

    if not a.cell:
        ap.error("pass --cell (with --observed) or --all")
    observed = [s for s in a.observed.split(",") if s.strip()]
    res = PN.flag_for_cell(observed, a.cell)
    out = {"cell": a.cell, "observed": observed, **res.as_dict(),
           "interpretation": ("position_novel=True means the genotype carries a substitution at a "
                              "catalogued resistance residue that is NOT itself catalogued -> the "
                              "susceptible-by-absence call is LEAST trustworthy here. NOT an R/S call.")}
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
