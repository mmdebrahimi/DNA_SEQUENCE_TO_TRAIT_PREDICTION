"""Clade and regime specializations — the second layer of the constraint hierarchy.

TAXON AXIS (molecular machinery). Shipped here are the genetic-code specializations, each `sourced` with
real NCBI provenance fetched 2026-10-01. A deeper bucket overrides a shallower one, so
`bacteria.mycoplasma` resolves to table 4 while `bacteria` resolves to table 11.

WHAT IS DELIBERATELY NOT SHIPPED, and why that is the gate working. Several further clade properties
belong on this axis — splicing in eukaryotes, polycistronic operons and horizontal gene transfer in
bacteria, ploidy — and all are well established. But the registry requires a `sourced` constraint to carry
`source_url` + `verbatim_quote`, and supplying those from recall is precisely the fabrication hazard this
project refuses (a volunteered ECOFF value was rejected for that reason). They are declared in
`DEFERRED_CLADE_PROPERTIES` with the evidence each would need, and are NOT registered. The gate refused
its own author's under-sourced entries, which is the behaviour it exists for.

REGIME AXIS (study design) — and a correction to the plan's shape. The plan said "regime specializations
keyed by the imported regime keys", implemented as constraints. That is the wrong object: a regime is a
property of the STUDY, not of the sequence, so it can never REFUTE a sequence. Registering
always-inapplicable constraints would have manufactured exactly the inert filters `report.py` exists to
catch. Instead `regime_scope()` exposes the regime's recorded verdict as advisory scoping metadata, read
live from `eval.regime`. The registry's cross-axis conflict machinery is unchanged and still exercised by
the registry's own tests.

NO "HIGHER LIFE FORM" AXIS. Clade nuances SWAP rather than accumulate — bacteria have HGT and operons,
eukaryotes have splicing and chromatin. Depth in the bucket tree does not mean "more constrained".
"""
from __future__ import annotations

from typing import Any

from dna_decode.constraints.codon_tables import TABLE_IDS, table_for, table_meta
from dna_decode.constraints.registry import (
    INAPPLICABLE,
    SATISFIED,
    SOURCED,
    Constraint,
    ConstraintVerdict,
    register,
)

GENETIC_CODE_KEY = "genetic_code_table"

#: Clades that get an explicit genetic-code specialization. Each must exist in TABLE_IDS and in the
#: registry's TAXON_BUCKETS. Restricted to buckets whose table actually differs from, or is explicitly
#: confirmed equal to, the standard code — all four are NCBI-sourced.
_CODE_CLADES: tuple[str, ...] = (
    "bacteria",                  # table 11, internal assignments equal to standard
    "bacteria.mycoplasma",       # table 4, UGA -> Trp
    "bacteria.spiroplasma",      # table 4
    "mitochondrion.vertebrate",  # table 2, AGA/AGG -> Ter, AUA -> Met, UGA -> Trp
    "eukaryote.ciliate",         # table 6, UAA/UAG -> Gln
)


def _make_code_check(clade: str):
    """A specialization that asserts WHICH code applies; it reports the table rather than refuting."""
    def check(ctx: dict[str, Any]) -> ConstraintVerdict:
        meta = table_meta(clade)
        if not ctx.get("cds"):
            return ConstraintVerdict(
                INAPPLICABLE, f"no cds; clade {clade} would use NCBI table {meta.table_id}")
        return ConstraintVerdict(
            SATISFIED, f"clade {clade} uses NCBI table {meta.table_id} ({meta.name})")
    return check


def genetic_code_specializations() -> tuple[Constraint, ...]:
    """Build (do not register) the sourced genetic-code specializations."""
    out = []
    for clade in _CODE_CLADES:
        assert clade in TABLE_IDS, clade
        meta = table_meta(clade)
        out.append(Constraint(
            key=GENETIC_CODE_KEY,
            scope=f"taxon:{clade}",
            kind=SOURCED,
            check=_make_code_check(clade),
            source_url=meta.source_url,
            verbatim_quote=meta.verbatim_quote,
            note=f"NCBI translation table {meta.table_id}: {meta.differences_from_standard}",
            payload={"table_id": meta.table_id, "clade": clade,
                     "differs_from_standard": sorted(
                         c for c, aa in table_for(clade).items()
                         if aa != table_for("standard")[c])},
        ))
    return tuple(out)


def register_specializations() -> tuple[Constraint, ...]:
    out = []
    for c in genetic_code_specializations():
        try:
            out.append(register(c))
        except Exception:
            out.append(c)
    return tuple(out)


#: Clade properties that BELONG on the taxon axis but are not shipped, because the registry requires
#: provenance for a `sourced` fact and inventing a citation is the hazard the gate exists to stop.
DEFERRED_CLADE_PROPERTIES: dict[str, dict[str, str]] = {
    "splicing_applies": {
        "buckets": "eukaryote (and descendants); absent in bacteria/archaea",
        "why_it_matters": "changes what a 'coding change' means and adds splice-site constraints",
        "evidence_needed": "a citable statement of intron presence/absence per clade",
    },
    "operon_polycistronic": {
        "buckets": "bacteria, archaea",
        "why_it_matters": "a frameshift can affect downstream genes in the same transcript",
        "evidence_needed": "a citable statement of polycistronic transcription per clade",
    },
    "hgt_acquisition_possible": {
        "buckets": "bacteria",
        "why_it_matters": "whole-gene acquisition is why a distributed mobile-element resistance "
                          "mechanism is bacteria-specific (the measured tet failure mode)",
        "evidence_needed": "a citable statement on HGT/plasmid-borne gene acquisition",
    },
    "ploidy": {
        "buckets": "eukaryote (diploid+), vs haploid bacteria/archaea",
        "why_it_matters": "dominance logic and heterozygote interpretation",
        "evidence_needed": "a citable per-clade ploidy statement",
    },
}


def regime_scope(regime_key: str) -> dict[str, Any]:
    """Advisory scoping metadata for a regime, read live from `eval.regime`.

    NOT a constraint: a regime describes the study, so it cannot refute a sequence. It tells a caller
    whether the question is well-posed in that regime — e.g. `natural_organism_zeroshot` is a recorded
    CLOSED_NEGATIVE, so a proposal landing there should be redirected rather than scored.
    """
    from dna_decode.eval.regime import REGIMES
    by_key = {r.key: r for r in REGIMES}
    if regime_key not in by_key:
        raise KeyError(f"unknown regime {regime_key!r}; known: {sorted(by_key)}")
    r = by_key[regime_key]
    return {"key": r.key, "population": r.population, "endpoint": r.endpoint, "method": r.method,
            "verdict": r.verdict, "artifact": r.artifact,
            "is_closed": r.verdict == "CLOSED_NEGATIVE",
            "advisory": ("a regime scopes whether a question is well-posed; it never refutes a "
                         "sequence, which is why it is not registered as a constraint")}


def regime_table() -> dict[str, dict[str, Any]]:
    """Every regime's scoping record, so a caller can see the whole axis at once."""
    from dna_decode.eval.regime import REGIMES
    return {r.key: regime_scope(r.key) for r in REGIMES}
