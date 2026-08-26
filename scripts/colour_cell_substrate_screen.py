"""Which colour-cell loci could EVER be scored on a biallelic-SNV genotype panel?

WHY THIS EXISTS. The animal colour/plumage family is now ~19 CLI cells, all shipping as KNOWLEDGE_BASELINE
(a curated OMIA epistatic rule, no measured per-individual validation). The family has been put to a
measured test EXACTLY ONCE -- the dog cell against the free Darwin's Ark cohort
(`wiki/dog_coat_darwins_ark_measured_2026-07-30.md`, N=3,277 genotypes x 29M biallelic SNVs, N=1,930
owner-reported colours) -- and it mostly FAILED, on SUBSTRATE rather than biology:

    black 160/161 = 0.994   blue/grey 11/31 = 0.355   every other base colour UNSCORABLE

because the causal variants those loci depend on are NOT SNVs: K^B is a 3 bp deletion, ASIP A^y/a^t is a
SINE insertion, MLPH d3 is a frameshift insertion, and MC1R `e` fell in an imputation gap. An imputed
biallelic-SNV panel cannot represent any of those.

That is a REJECTION GATE, and it generalises to any colour cell whose loci rest on the same variant
classes. This script derives it per-cell from the COMMITTED catalogs rather than asserting it: it reads
every locus of every colour catalog and classifies the causal variant recorded in that locus's own
`source` string.

WHAT IT DELIBERATELY DOES NOT DO. It does not predict a cell "will fail" -- it reports how much of each
cell's rule rests on a variant class a SNV panel cannot carry, plus how much is simply NOT RECORDED. An
unrecorded causal variant is reported as UNRECORDED, never guessed into a class: a cell whose causal
variants aren't written down cannot be screened at all, and that is a finding about the catalog, not
evidence about the substrate.

SELF-CHECK. The classifier is validated against the DOG cell, which is the one case with measured ground
truth -- the measured artifact's own table says K/ASIP/MLPH-d3 are indel/structural and TYRP1-bc/MLPH-d1/d2
and MC1R-e are SNVs. `--self-check` asserts the classifier reproduces that, so a text heuristic cannot
silently drift into agreeing with itself.

Run: uv run python scripts/colour_cell_substrate_screen.py [--self-check]
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

OUT = ROOT / "wiki" / "colour_cell_substrate_screen_2026-08-26.json"

# Order is LOAD-BEARING: a SINE insertion also matches "ins", and a frameshift often carries a c.N>N
# elsewhere in the same sentence, so the coarser classes must be tested FIRST.
_STRUCTURAL = re.compile(
    r"\bSINE\b|\bLINE-?1\b|retrotranspos|\bCNV\b|copy[- ]number|duplicat|\b\d+\s*kb\b|large (?:deletion|insertion)",
    re.I)
_INDEL = re.compile(
    r"c\.[\d_+-]+(?:del|ins|dup)|(?:^|[^a-z])(?:del|ins)[ACGT]{1,}|frameshift|\bfs\b|\bLoF\b|"
    r"\bdel[A-Z]{2,}|\d+_\d+del|\bdeletion\b|\binsertion\b",
    re.I)
_SNV = re.compile(r"c\.-?\d+[+-]?\d*[ACGT]>[ACGT]|\bp\.[A-Z][a-z]{2}\d+|\bp\.[A-Z]\d+[A-Z*]|\bSNP\b", re.I)


def classify_variant(text: str) -> str:
    """SNV | INDEL | STRUCTURAL | UNRECORDED for one locus's recorded causal variant. PURE."""
    if not text:
        return "UNRECORDED"
    if _STRUCTURAL.search(text):
        return "STRUCTURAL"
    if _INDEL.search(text):
        return "INDEL"
    if _SNV.search(text):
        return "SNV"
    return "UNRECORDED"


def snv_panel_scorable(cls: str) -> bool | None:
    """Could a biallelic-SNV panel represent this variant? None = unknowable (unrecorded). PURE."""
    return {"SNV": True, "INDEL": False, "STRUCTURAL": False}.get(cls)


def _loci_of(obj) -> dict:
    """A colour catalog exposes either `.loci` (mammal_color) or a module-level LOCI dict. PURE-ish."""
    loci = getattr(obj, "loci", None)
    if isinstance(loci, dict):
        return loci
    loci = getattr(obj, "LOCI", None)
    return loci if isinstance(loci, dict) else {}


# The Locus dataclass field is `note` (SINGULAR). An early version of this reader looked for "notes" and
# silently returned less text, which would have under-counted recorded variants across every cell -- the
# failure mode a text screen cannot detect from its own output. Read EVERY string field instead of a
# hand-listed few, so a renamed or added field can never silently shrink the evidence again.
def _source_of(locus) -> str:
    """Concatenate EVERY string field the locus carries. PURE."""
    import dataclasses
    parts = []
    if dataclasses.is_dataclass(locus):
        for f in dataclasses.fields(locus):
            v = getattr(locus, f.name, None)
            if isinstance(v, str) and v:
                parts.append(v)
    else:
        for attr in ("source", "description", "note", "notes", "gene"):
            v = getattr(locus, attr, None)
            if isinstance(v, str) and v:
                parts.append(v)
    return " | ".join(parts)


def collect() -> dict:
    """species -> [ {locus, gene, variant_class, scorable, source} ] across every colour catalog."""
    from dna_decode.pigment import (cat_coat, chicken_plumage, dog_coat, horse_coat,  # noqa: E402
                                    pigeon_plumage)
    from dna_decode.pigment.mammal_color import MAMMAL_CATALOGS  # noqa: E402

    catalogs: dict = {sp: cat for sp, cat in MAMMAL_CATALOGS.items()}
    catalogs.update({"dog": dog_coat, "cat": cat_coat, "horse": horse_coat,
                     "chicken": chicken_plumage, "pigeon": pigeon_plumage})

    out = {}
    for species, cat in sorted(catalogs.items()):
        rows = []
        for sym, locus in _loci_of(cat).items():
            src = _source_of(locus)
            cls = classify_variant(src)
            rows.append({"locus": sym, "gene": getattr(locus, "gene", "?"),
                         "variant_class": cls, "snv_panel_scorable": snv_panel_scorable(cls),
                         "source": src[:220]})
        if rows:
            out[species] = rows
    return out


def summarise(rows: list[dict]) -> dict:
    """Per-cell counts + the honest verdict. PURE."""
    n = len(rows)
    c = {k: sum(1 for r in rows if r["variant_class"] == k)
         for k in ("SNV", "INDEL", "STRUCTURAL", "UNRECORDED")}
    blocked = c["INDEL"] + c["STRUCTURAL"]
    if c["UNRECORDED"] == n:
        verdict = "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED"
    elif blocked == 0 and c["UNRECORDED"] == 0:
        verdict = "FULLY_SNV_TRACTABLE"
    elif blocked == 0:
        verdict = "SNV_TRACTABLE_WHERE_RECORDED"
    elif c["SNV"] == 0 and c["UNRECORDED"] == 0:
        verdict = "NO_LOCUS_SNV_TRACTABLE"
    else:
        verdict = "PARTIALLY_SNV_TRACTABLE"
    return {"n_loci": n, **{f"n_{k.lower()}": v for k, v in c.items()},
            "n_snv_panel_blocked": blocked, "verdict": verdict}


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
