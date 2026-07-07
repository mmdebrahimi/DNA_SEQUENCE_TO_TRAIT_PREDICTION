"""CYP2C8 curated catalog — defining variants + honest CALLING-only design (v0, SNP-defined core).

Third PGx CYP2C-cluster gene (joins CYP2C9 + CYP2C19 on chr10). CYP2C8 is the one remaining SNP-defined
gene in the GeT-RM / ursaPGx benchmark truth set (`CYP2C8_getrm_ngs`), so its star-allele CALLING is
INDEPENDENTLY validatable vs the GeT-RM consensus — the same tier as the C9/C19 cells.

PROVENANCE (grounded, NO fabrication, NO build-on-unverified):
  * Defining variants + GRCh38 coordinates: VERIFIED via the Ensembl REST variation API (forward-strand
    ref/alt at each rsID; CYP2C8 is on the chr10 MINUS strand, so the coding cDNA change is the
    reverse-complement of the forward-strand genomic ref>alt recorded here):
      *2  rs11572103  c.805A>T (p.Ile269Phe)   chr10:95058349 T>A  (fwd; coding A>T)
      *3  rs11572080  c.416G>A (p.Arg139Lys)   chr10:95067273 C>T  (fwd; coding G>A)  [primary *3 tag SNP]
      *4  rs1058930   c.792C>G (p.Ile264Met)   chr10:95058362 G>C  (fwd; coding C>G)
    *3 is a two-variant haplotype: rs11572080 (R139K) + rs10509681 (c.1196A>G, K399R, chr10:95038992 T>C)
    in near-complete linkage (Speed/Dahl 2009; the ScienceDirect CYP2C8 review). v0 tags *3 on the single
    functional variant rs11572080 (R139K); rs10509681 is the co-defining LD partner (NOT separately scanned
    in the single-SNP-per-star framework). The GeT-RM concordance test falsifies this tag choice.
  * Star->rsID->protein mapping: the CYP2C8 pharmacogenetics review (ScienceDirect CYP2C8 overview) +
    Rodriguez-Antona 2008 / the *3-*4 phenotype literature; cross-checked against dbSNP.

HONESTY TIER (load-bearing — DIFFERENT from C9/C19):
  * Star-allele CALLING is INDEPENDENTLY validatable vs the GeT-RM consensus (a free, measured/consensus
    label). caller_is_independent_baseline=True for the calling step.
  * There is NO CPIC standardized metabolizer-phenotype (PM/IM/NM) table for CYP2C8. CYP2C8 allele
    function is SUBSTRATE-DEPENDENT (e.g. *3 shows DECREASED metabolism of paclitaxel / R-ibuprofen /
    arachidonic acid but INCREASED metabolism of pioglitazone / rosiglitazone / repaglinide). So this cell
    reports the star-allele DIPLOTYPE + a per-allele PharmVar-clinical FUNCTION annotation, NOT a
    metabolizer phenotype. HAS_CPIC_PHENOTYPE = False (unlike CYP2C9's activity-score / CYP2C19's
    function-pair). Do NOT emit or infer a PM/IM/NM call for CYP2C8.

SCOPE (v0 — honest): CORE SNP-defined *2/*3/*4 + *1 (reference). CYP2C8's common functional alleles ARE
*2/*3/*4, so core coverage is good; rarer non-core CYP2C8 alleles would be mis-called *1 (a documented
residual — no sentinel layer in v0, mirroring the C9 v0 arc). NOT a clinical tool.
"""
from __future__ import annotations

from dna_decode.pgx.cyp2c19_catalog import DefiningVariant, SentinelVariant  # reuse the dataclasses

GENE = "CYP2C8"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

# CYP2C8 has NO CPIC metabolizer-phenotype system (substrate-dependent function). Callers/consumers MUST
# NOT derive a PM/IM/NM phenotype from a CYP2C8 diplotype.
HAS_CPIC_PHENOTYPE = False

# Core SNP-defined CYP2C8 alleles (v0). GRCh38 forward-strand ref>alt, VERIFIED via Ensembl REST.
CORE_DEFINING: list[DefiningVariant] = [
    DefiningVariant("*2", "rs11572103", "10", 95058349, "T", "A", "c.805A>T"),
    DefiningVariant("*3", "rs11572080", "10", 95067273, "C", "T", "c.416G>A"),
    DefiningVariant("*4", "rs1058930", "10", 95058362, "G", "C", "c.792C>G"),
]

# No sentinel layer in v0 (CYP2C8's common alleles are all in the core set). Documented residual.
SENTINELS: list[SentinelVariant] = []

# PharmVar/literature CLINICAL function per allele (substrate-dependent; NOT a CPIC activity value).
ALLELE_FUNCTION: dict[str, str] = {
    "*1": "normal",                        # reference, normal function
    "*2": "decreased",                     # decreased (paclitaxel / arachidonic acid); Sub-Saharan African-common
    "*3": "decreased_or_substrate_dependent",  # decreased paclitaxel/R-ibuprofen; INCREASED pioglitazone/repaglinide
    "*4": "decreased",                     # decreased; European-common
}


def diplotype_function(allele1: str, allele2: str) -> str:
    """Honest CYP2C8 diplotype output: the per-allele FUNCTION pair, NOT a CPIC metabolizer phenotype.

    Returns a stable descriptor string (order-independent). 'Indeterminate' if either allele is unknown.
    CYP2C8 function is substrate-dependent — this is a function annotation, not a PM/IM/NM call.
    """
    f1 = ALLELE_FUNCTION.get(allele1, "unknown")
    f2 = ALLELE_FUNCTION.get(allele2, "unknown")
    if "unknown" in (f1, f2):
        return "Indeterminate"
    a1, a2 = sorted((allele1, allele2))
    if f1 == f2 == "normal":
        return "Normal function (both alleles; substrate-dependent metabolism)"
    return (f"{a1} + {a2} — carries a non-normal-function CYP2C8 allele "
            "(substrate-dependent; NO CPIC metabolizer phenotype)")


# Backward-compat alias: the generic caller/runner expect a `phenotype_fn` seam. For CYP2C8 this seam
# returns the FUNCTION descriptor above (never a CPIC phenotype). Named `diplotype_phenotype` only so the
# generic `call_diplotype(..., phenotype_fn=...)` wiring is uniform across genes.
diplotype_phenotype = diplotype_function

# Mechanisms a CORE-SNP CYP2C8 scan cannot see (the honest blind-spot list).
UNDETECTABLE = sorted({
    "non_core_star_allele",            # rare CYP2C8 alleles beyond *2/*3/*4 -> mis-called *1
    "cnv_or_gene_deletion",
    "novel_uncatalogued_variant",
    "phase_ambiguity_unphased_input",
})
