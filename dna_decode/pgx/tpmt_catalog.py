"""TPMT curated catalog — COMPOUND star alleles (*3A = *3B + *3C) + CPIC thiopurine phenotype (v0).

Fifth CYP-cluster-independent PGx gene (chr6p22.3). Thiopurine (azathioprine/6-MP) dosing. TPMT is the
project's first true COMPOUND-allele cell: *3A is DEFINED by two SNPs on the same haplotype (*3B + *3C),
and each SNP ALONE is a different allele — so a single-SNP tag would mis-call. Resolved via
`compound_caller.assemble_compound_diplotype` (per-haplotype component-set -> star).

PROVENANCE (grounded, NO fabrication):
  * Defining variants + GRCh38 coords VERIFIED via Ensembl REST; TPMT is on the chr6 MINUS strand, so the
    cDNA change is the reverse-complement of the forward-strand genomic ref>alt recorded here. ALT
    orientation EMPIRICALLY CONFIRMED on 1000G (freqs: *3B 0.014, *3C 0.040 — TPMT*3 is rare, *3C > *3B):
      *3B  rs1800460  c.460G>A (p.Ala154Thr)  chr6:18138997 C>T
      *3C  rs1142345  c.719A>G (p.Tyr240Cys)  chr6:18130687 T>C
      *3A  = BOTH rs1800460 + rs1142345 in cis (the common European no-function haplotype)
  * Allele FUNCTION + phenotype: CPIC thiopurine guideline (Relling 2019) — *1 normal; *2/*3A/*3B/*3C
    no function. Diplotype -> Normal/Intermediate/Poor Metabolizer.

HONESTY TIER: star-allele CALLING is INDEPENDENTLY validatable vs the GeT-RM CDC multi-lab consensus
(the consolidated 363-sample PGx table). Validated 85/85 core-comparable on 1000G-overlap
(truth in *1/*3A/*3B/*3C), and the compound *3A path IS exercised (6 *3A + 8 *3C truth samples).
caller_is_independent_baseline=True. Phenotype FAITHFUL-TO-CPIC. Reference tool: PharmCAT.

SCOPE (v0): CORE *3A/*3B/*3C + *1. Rarer no-function alleles (*2/*8/*16/*46...) are non-core -> mis-called
*1 (documented residual; scored separately). NOT a clinical tool.
"""
from __future__ import annotations

from dna_decode.pgx.compound_caller import CompoundAllele
from dna_decode.pgx.cyp2c19_catalog import DefiningVariant

GENE = "TPMT"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

# Component SNPs (the `.star` field carries the COMPONENT TAG used by the compound rules, not a star).
# GRCh38 forward-strand ref>alt, Ensembl-verified + AF-confirmed.
COMPONENTS: list[DefiningVariant] = [
    DefiningVariant("c460T", "rs1800460", "6", 18138997, "C", "T", "c.460G>A (*3B)"),
    DefiningVariant("c719C", "rs1142345", "6", 18130687, "T", "C", "c.719A>G (*3C)"),
]

# Per-haplotype compound rules (most-specific first is handled by the resolver):
#   both component SNPs in cis -> *3A ; c460T alone -> *3B ; c719C alone -> *3C ; neither -> *1
COMPOUND_RULES: list[CompoundAllele] = [
    CompoundAllele("*3A", frozenset({"c460T", "c719C"})),
    CompoundAllele("*3B", frozenset({"c460T"})),
    CompoundAllele("*3C", frozenset({"c719C"})),
]

# `CORE_DEFINING` kept as the component list so the report/registry wiring that reads a defining list works.
CORE_DEFINING = COMPONENTS
SENTINELS: list = []

ALLELE_FUNCTION: dict[str, str] = {
    "*1": "normal",
    "*3A": "none",
    "*3B": "none",
    "*3C": "none",
}

_PHENOTYPE_BY_FUNCTION_PAIR: dict[tuple[str, str], str] = {
    ("normal", "normal"): "Normal Metabolizer",
    ("none", "normal"): "Intermediate Metabolizer",
    ("none", "none"): "Poor Metabolizer",
}

PHENOTYPE_ABBREV: dict[str, str] = {
    "Normal Metabolizer": "NM",
    "Intermediate Metabolizer": "IM",
    "Poor Metabolizer": "PM",
    "Indeterminate": "IND",
}

UNDETECTABLE = sorted({
    "non_core_star_allele",            # *2/*8/*16/*46... beyond the *3-family core -> mis-called *1
    "cnv_or_gene_deletion",
    "novel_uncatalogued_variant",
    "phase_ambiguity_unphased_input",  # >=2 unphased het components -> cis(*3A) vs trans(*3B/*3C) ambiguous
})


def diplotype_phenotype(allele1: str, allele2: str) -> str:
    """CPIC TPMT metabolizer phenotype. 'Indeterminate' if either allele unknown. Order-independent."""
    f1 = ALLELE_FUNCTION.get(allele1, "unknown")
    f2 = ALLELE_FUNCTION.get(allele2, "unknown")
    if "unknown" in (f1, f2):
        return "Indeterminate"
    return _PHENOTYPE_BY_FUNCTION_PAIR[tuple(sorted((f1, f2)))]
