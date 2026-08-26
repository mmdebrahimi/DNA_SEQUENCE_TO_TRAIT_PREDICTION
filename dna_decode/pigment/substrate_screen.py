"""Which colour-cell loci could EVER be scored on a biallelic-SNV genotype panel?

WHY THIS EXISTS. The animal colour/plumage family is 19 CLI cells, all shipping as KNOWLEDGE_BASELINE
(a curated OMIA epistatic rule, no measured per-individual validation). The family has been put to a
measured test EXACTLY ONCE -- the dog cell against the free Darwin's Ark cohort
(`wiki/dog_coat_darwins_ark_measured_2026-07-30.md`, N=3,277 genotypes x 29M biallelic SNVs, N=1,930
owner-reported colours) -- and it mostly FAILED, on SUBSTRATE rather than biology:

    black 160/161 = 0.994   blue/grey 11/31 = 0.355   every other base colour UNSCORABLE

because the causal variants those loci depend on are NOT SNVs: K^B is a 3 bp deletion, ASIP A^y/a^t is a
SINE insertion, MLPH d3 is a frameshift insertion, and MC1R `e` fell in an imputation gap. An imputed
biallelic-SNV panel cannot represent any of those.

That is a REJECTION GATE (map gates G9/G10), and it generalises to any colour cell whose loci rest on the
same variant classes. This module derives it per-cell from the COMMITTED catalogs rather than asserting
it: it reads every locus of every colour catalog and classifies the causal variant recorded in that
locus's own provenance text.

WHY IT LIVES IN THE PACKAGE (not in scripts/). `dna_decode/data/colour_cell_freeze.py` consumes
`verdicts()` as its single derivation, and an in-package module cannot import from `scripts/`. The script
`scripts/colour_cell_substrate_screen.py` remains the CLI + artifact writer over this logic.

WHAT IT DELIBERATELY DOES NOT DO. It does not predict a cell "will fail" -- it reports how much of each
cell's rule rests on a variant class a SNV panel cannot carry, plus how much is simply NOT RECORDED. An
unrecorded causal variant is reported as UNRECORDED, never guessed into a class: a cell whose causal
variants aren't written down cannot be screened at all, and that is a finding about the catalog, not
evidence about the substrate.
"""
from __future__ import annotations

import dataclasses
import re

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

# The screen keys cells by SPECIES; the CLI routes them by TRAIT. 17 of 19 are `<species>color`; dog and
# chicken are the exceptions. Derived checks cross this map, so a new colour cell that forgets an entry
# surfaces as a missing route rather than silently dropping out of the screen.
SPECIES_TO_TRAIT: dict[str, str] = {
    "dog": "coatcolor",
    "chicken": "plumage",
}

VERDICTS = (
    "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED",
    "NO_LOCUS_SNV_TRACTABLE",
    "PARTIALLY_SNV_TRACTABLE",
    "SNV_TRACTABLE_WHERE_RECORDED",
    "FULLY_SNV_TRACTABLE",
)


def trait_for_species(species: str) -> str:
    """Screen species key -> CLI trait name. PURE."""
    return SPECIES_TO_TRAIT.get(species, f"{species}color")


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
    """species -> [ {locus, gene, variant_class, snv_panel_scorable, source} ] across every catalog."""
    from dna_decode.pigment import (cat_coat, chicken_plumage, dog_coat, horse_coat,
                                    pigeon_plumage)
    from dna_decode.pigment.mammal_color import MAMMAL_CATALOGS

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


def verdicts(data: dict | None = None) -> dict[str, str]:
    """CLI trait name -> screen verdict. THE single derivation every consumer reads. PURE given `data`.

    Keyed by TRAIT (not species) so the freeze module and the contract guards can join against
    `cli.TRAITS` / `cell_registry` routes without re-deriving the naming exception for dog and chicken.
    """
    data = collect() if data is None else data
    return {trait_for_species(sp): summarise(rows)["verdict"] for sp, rows in data.items()}
