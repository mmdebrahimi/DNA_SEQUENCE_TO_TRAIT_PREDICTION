"""CLI + artifact writer for the colour-cell substrate screen.

The PURE LOGIC lives in `dna_decode/pigment/substrate_screen.py` -- see that module for what the screen
does and why. It was moved in-package because `dna_decode/data/colour_cell_freeze.py` consumes
`verdicts()` as its single derivation, and an in-package module cannot import from `scripts/`.

This file keeps the parts that are genuinely CLI: the self-check anchor, the catalog-gap record, the
table rendering, and the JSON artifact.

SELF-CHECK. The classifier is a text heuristic, so it is anchored against the DOG cell -- the one case
with measured ground truth. `--self-check` asserts the classifier reproduces what the dog CATALOG
records (not what the literature knows; conflating the two is what made the first version fail).

Run: uv run python scripts/colour_cell_substrate_screen.py [--self-check]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.pigment.substrate_screen import (  # noqa: E402
    classify_variant, collect, snv_panel_scorable, summarise, verdicts,
)

__all__ = ["classify_variant", "collect", "snv_panel_scorable", "summarise", "verdicts",
           "self_check", "main"]

OUT = ROOT / "wiki" / "colour_cell_substrate_screen_2026-08-26.json"

# The DOG cell is the ONE case with measured ground truth, so it anchors the classifier. What is asserted
# here is what the classifier must read OUT OF THE CATALOG TEXT -- not what the literature knows.
_DOG_TRUTH = {"K": "INDEL",          # catalog records CBD103 c.67_69delGGT (3 bp deletion)
              "A": "UNRECORDED",     # catalog records NO causal variant -- see _CATALOG_GAPS below
              "B": "SNV",            # catalog records TYRP1 bs p.Gln331Ter
              "D": "SNV",            # catalog records MLPH c.-22G>A
              "E": "SNV"}            # catalog records MC1R p.Arg306Ter

# Where the MEASURED artifact knows more than the catalog records. This is the honest gap: the screen can
# only classify what is written down, so a locus here is invisible to it even though the substrate answer
# is known. Found by running this self-check -- the classifier said UNRECORDED for dog A and it was RIGHT;
# the initial expectation (STRUCTURAL) had encoded the literature rather than the catalog.
_CATALOG_GAPS = {
    "dog/A": ("measured artifact records ASIP A^y/a^t as a SINE insertion + coding change (STRUCTURAL, "
              "absent from a biallelic-SNV panel); the catalog `source` names only the locus and papers"),
}


def self_check(data: dict) -> list[str]:
    """Does the classifier reproduce what the DOG catalog RECORDS? Returns failures.

    Anchoring on a real cell keeps a text heuristic from silently agreeing with itself. It deliberately
    asserts what the CATALOG says, not what the literature says -- conflating the two is what made the
    first version of this check fail.
    """
    got = {r["locus"]: r["variant_class"] for r in data.get("dog", [])}
    fails = []
    for locus, expect in _DOG_TRUTH.items():
        actual = got.get(locus)
        if actual is None:
            fails.append(f"dog locus {locus} absent from the screen (present: {sorted(got)})")
        elif actual != expect:
            fails.append(f"dog {locus}: classifier said {actual}, dog catalog records {expect}")
    return fails


def main() -> int:
    data = collect()
    fails = self_check(data)

    if "--self-check" in sys.argv:
        for f in fails:
            print("FAIL:", f)
        print("self-check:", "PASS" if not fails else f"{len(fails)} mismatch(es)")
        return 1 if fails else 0

    report = {"_schema": "colour-cell-substrate-screen-v1",
              "ground_truth": "wiki/dog_coat_darwins_ark_measured_2026-07-30.md",
              "self_check_failures": fails,
              "catalog_gaps_vs_measured_artifact": _CATALOG_GAPS,
              "honest_scope": (
                  "Derived from the causal variant each locus records in its OWN catalog `source` string. "
                  "UNRECORDED means the catalog does not write the causal variant down -- that cell cannot "
                  "be screened at all, which is a finding about the catalog, NOT evidence about the "
                  "substrate. A blocked locus means a biallelic-SNV panel cannot REPRESENT the variant; it "
                  "does not by itself predict a cell would fail, and no cell here has been measured except "
                  "dog."),
              "cells": {}}

    print(f"{'cell':12s} {'loci':>5s} {'SNV':>4s} {'INDEL':>6s} {'STRUCT':>7s} {'UNREC':>6s} {'blocked':>8s}  verdict")
    tot = {"SNV": 0, "INDEL": 0, "STRUCTURAL": 0, "UNRECORDED": 0}
    for species, rows in sorted(data.items()):
        s = summarise(rows)
        report["cells"][species] = {**s, "loci": rows}
        for k in tot:
            tot[k] += sum(1 for r in rows if r["variant_class"] == k)
        print(f"{species:12s} {s['n_loci']:5d} {s['n_snv']:4d} {s['n_indel']:6d} {s['n_structural']:7d} "
              f"{s['n_unrecorded']:6d} {s['n_snv_panel_blocked']:8d}  {s['verdict']}")

    n_cells = len(data)
    n_all = sum(tot.values())
    report["totals"] = {"n_cells": n_cells, "n_loci": n_all, **{f"n_{k.lower()}": v for k, v in tot.items()},
                        "n_snv_panel_blocked": tot["INDEL"] + tot["STRUCTURAL"]}
    print(f"\n{n_cells} colour cells / {n_all} loci: SNV {tot['SNV']} | INDEL {tot['INDEL']} | "
          f"STRUCTURAL {tot['STRUCTURAL']} | UNRECORDED {tot['UNRECORDED']}")
    print(f"self-check vs the dog measured artifact: {'PASS' if not fails else str(fails)}")

    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
