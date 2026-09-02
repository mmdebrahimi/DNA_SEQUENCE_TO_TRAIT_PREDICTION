"""Screen PEAR's G6 gate on the REAL extracted fitness values, and re-run the full ten-gate screen.

G6 was the one gate PEAR left open. The memo could not close it because nobody had read a single PEAR
fitness value -- the repo ships no plain text, only .RData workspaces holding ggplot objects. Those are
now extracted (scripts/pear_extract_fitness.R, outputs on D:), so the gate can be MEASURED instead of
declared open.

G6's L4 form is assay DEGENERACY, and the bar is the shipped one from the forward/inverse cell
(`assay_degeneracy`: mode-share > 25% or fewer than 20 distinct levels). That gate exists because a
censored assay does not fail loudly -- CcdB, 79.3% tied at its ceiling, posted the sweep's BEST number.

Offline. Reads only the extracted TSVs on D:; writes wiki/pear_g6_screen_<date>.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.eval.rejection_gates import L4_FORWARD_CONTINUOUS, screen_candidate  # noqa: E402

EXTRACTED = Path("D:/dna_decode_cache/pear/extracted")


def _load(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _floats(rows: list[dict], col: str, drop_gt_mean: bool = False,
            drop_wt: bool = False) -> list[float]:
    """Keep only rows that are a MEASURED VARIANT.

    Two kinds of row are not one, and both would corrupt a degeneracy screen:

    - `gt == "mean"` is a per-position AGGREGATE.
    - `gt == isWt` is the WILD-TYPE base at that position. Its relative growth is 1.0 BY CONSTRUCTION
      (the normalizer measured against itself), so it is one identical value repeated once per position.
      Left in, it alone puts the mode share at 792/2957 = 26.8% and trips the >25% degeneracy bar --
      a censored-assay verdict manufactured entirely by the baseline. Verified: all 792 rows carrying
      the modal value are WT rows, and ZERO non-WT rows carry it.

    Same class of defect as the NNRTI `L234L` self-to-self entries: a non-variant admitted as a variant.
    """
    out = []
    for r in rows:
        if drop_gt_mean and r.get("gt") == "mean":
            continue
        if drop_wt and r.get("gt") is not None and r.get("gt") == r.get("isWt"):
            continue
        v = r.get(col, "")
        if v in ("", "NA", "NaN", None):
            continue
        try:
            out.append(float(v))
        except ValueError:
            continue
    return out


def assay_degeneracy(values: list[float]) -> dict:
    """Imported behaviourally from scripts/forward_inverse_roundtrip.py -- same bars, same verdict."""
    from forward_inverse_roundtrip import assay_degeneracy as shipped  # type: ignore
    return shipped(values)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--extracted", type=Path, default=EXTRACTED)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))   # so the shipped gate imports

    if not a.extracted.is_dir():
        print(f"extracted dir not found: {a.extracted}\nRun scripts/pear_extract_fitness.R first.")
        return 2

    # (label, file, column, drop the per-position 'mean' pseudo-rows)
    targets = [
        # (label, file, column, drop 'mean' aggregates, drop wild-type baseline rows)
        ("figure2A_ceftazidime_per_nt", "Figure.2A__data.tsv", "effect_size", True, True),
        ("figure2B_cefotaxime_per_nt", "Figure.2B__data.tsv", "effect_size", True, True),
        ("figure3A_cefotaxime_per_variant", "Figure3.A__data.tsv", "CTX", False, False),
        ("figure3A_ceftazidime_per_variant", "Figure3.A__data.tsv", "CAZ", False, False),
        # kept deliberately: the SAME table with the WT baseline left in, to show what it does
        ("figure2A_ceftazidime_WT_INCLUDED_diagnostic", "Figure.2A__data.tsv", "effect_size", True, False),
    ]

    measured = {}
    for label, fname, col, drop_mean, drop_wt in targets:
        p = a.extracted / fname
        if not p.is_file():
            measured[label] = {"error": f"missing {fname}"}
            continue
        vals = _floats(_load(p), col, drop_gt_mean=drop_mean, drop_wt=drop_wt)
        deg = assay_degeneracy(vals)
        measured[label] = {"n": deg["n"], "n_distinct_values": deg["n_distinct_values"],
                           "mode_value": deg["mode_value"], "mode_share": deg["mode_share"],
                           "degenerate": deg["degenerate"]}
        print(f"{label:36} n={deg['n']:5d}  distinct={deg['n_distinct_values']:5d}  "
              f"mode_share={deg['mode_share']:.4f}  degenerate={deg['degenerate']}")

    # The gate is screened on the substrate a forward-cell comparison would actually use:
    # the per-variant table (Figure3.A), cefotaxime.
    primary = measured.get("figure3A_cefotaxime_per_variant", {})

    evidence = {
        "label_provenance_evidence": "relative growth of constructed blaCTX-M-14 variants under "
                                     "cefotaxime/ceftazidime selection, read out by barcode sequencing "
                                     "(Zhang 2022 MBE). No genomic tool produced the label.",
        "label_is_measured": True,
        "label_semantics_evidence": "relative growth is an assay reading, not a collection context",
        "label_is_assay_reading": True,
        "variation_is_constructed": True,
        "genotype_defined_by_construction": True,
        "loci_without_recorded_variant_fraction": 0.0,
        "off_panel_variant_fraction": 0.0,
        "mode_share": primary.get("mode_share"),
        "n_distinct_values": float(primary.get("n_distinct_values", 0)) or None,
    }
    res = screen_candidate("PEAR (Zhang 2022 blaCTX-M-14 DMS)", L4_FORWARD_CONTINUOUS, evidence)

    out = {"schema": "pear-g6-screen-v1",
           "substrate": "ggplot $data slots extracted from the authors' .RData via R 4.6.1",
           "measured": measured,
           "g6_screened_on": "figure3A_cefotaxime_per_variant",
           "screen": res.as_dict(),
           "honest_limits": [
               "These are the AGGREGATED per-variant effect sizes the authors plotted, NOT the ~23,000 "
               "raw barcoded strains. U5 is answered: the plot objects carry the full per-variant scan "
               "for the positions shown, not the raw library.",
               "Figure.2A/2B are NUCLEOTIDE-level (gt in A/C/G/T over 792 positions = 3,168 real "
               "substitutions per drug, plus 789 per-position means which are excluded here).",
               "Figure3.A is per-variant in C648T notation (2,114 variants) with CTX and CAZ columns -- "
               "the shape our genome-edit forward path consumes.",
               "A passing G6 bounds the ASSAY. It is not a build recommendation and says nothing about "
               "whether the forward cell will transfer to this protein.",
           ]}
    print(f"\nten-gate verdict: {res.verdict} -- {res.reason}")

    dest = a.out or (Path(__file__).resolve().parents[1] / "wiki" / "pear_g6_screen_2026-09-01.json")
    dest.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
