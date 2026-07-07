"""UGT1A1 curated catalog — SNP-defined core + CPIC activity-score phenotype (v0, tag-SNP surface).

UGT1A1 drives irinotecan toxicity + atazanavir hyperbilirubinemia (and Gilbert syndrome). CPIC (Gammal 2016)
assigns a metabolizer phenotype from the diplotype's allele functions — the activity-score shape again.

*** THE STRUCTURAL WALL (load-bearing honesty) ***
The MAJOR reduced-function allele UGT1A1*28 is a TA-DINUCLEOTIDE REPEAT in the promoter (TA6 = *1 normal,
TA7 = *28 reduced) — a SHORT TANDEM REPEAT, NOT a SNP. A short-read SNP VCF CANNOT reliably call the TA
repeat length (the CYP2D6-structural situation in a different gene). So v0 uses the validated TAG-SNP
approach labs actually deploy: **rs887829 (*80, C>T) is in strong LD with the TA7 (*28) repeat** (EUR r^2 ~
0.9+; EUR AF ~30% == the *28 frequency), so a *80 genotype is the SNP-callable PROXY for *28 reduced-function
status. This is a WRAPPER-OVER-LD-TAG cell (same posture as the project's HLA tag-SNP cells) — NOT a direct
repeat call. Named residual: the tag is imperfect (non-EUR LD breakdown + recombinants), and the TRUE *28
repeat length is unresolved without a repeat-aware caller.

PROVENANCE (grounded, no fabrication):
  * Coords VERIFIED via Ensembl GRCh38 REST 2026-07-07 (chr2 plus strand):
      *80  rs887829   promoter C>T (*28-tag)  chr2:233759924 C>T  -- DECREASED (LD-tag for *28), activity 0.5
      *6   rs4148323  c.211G>A p.G71R          chr2:233760498 G>A  -- DECREASED (EAS-common),     activity 0.5
  * Activity values + score->phenotype: CPIC UGT1A1 (Gammal 2016) — *1=1.0, decreased=0.5; AS 2.0 = Normal
    Metabolizer, AS 1.5 = Intermediate, AS 1.0 = Poor (irinotecan toxicity risk / Gilbert).

SCOPE (v0 -- honest): the SNP-callable actionable alleles — *80 (rs887829, the *28 LD-tag) + *6 (rs4148323)
+ *1. Direct *28/*37/*36 TA-repeat lengths are NOT called (structural wall). NOT a clinical tool.
"""
from __future__ import annotations

from dna_decode.pgx.cyp2c19_catalog import DefiningVariant, SentinelVariant  # reuse the dataclasses

GENE = "UGT1A1"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

CORE_DEFINING: list[DefiningVariant] = [
    DefiningVariant("*80", "rs887829", "2", 233759924, "C", "T", "promoter C>T (*28 LD-tag; decreased)"),
    DefiningVariant("*6", "rs4148323", "2", 233760498, "G", "A", "c.211G>A p.G71R (decreased; EAS)"),
]

SENTINELS: list[SentinelVariant] = []

ACTIVITY_VALUE: dict[str, float] = {"*1": 1.0, "*80": 0.5, "*6": 0.5}

PHENOTYPE_ABBREV = {"Normal Metabolizer": "NM", "Intermediate Metabolizer": "IM",
                    "Poor Metabolizer": "PM", "Indeterminate": "IND"}

# The load-bearing structural caveat every consumer must surface.
STRUCTURAL_CAVEAT = (
    "UGT1A1*28 is a promoter TA-repeat (STR), NOT a SNP — unresolvable from a short-read SNP VCF. v0 uses "
    "rs887829 (*80) as the validated LD-tag PROXY for *28 reduced-function status (EUR r^2 ~0.9+). This is a "
    "tag-SNP wrapper (like the HLA cells), NOT a direct repeat call; the TRUE *28 length is UNASSESSED.")

UNDETECTABLE = sorted({
    "ugt1a1_star28_ta_repeat_length",   # the promoter TA-repeat itself (structural wall)
    "ugt1a1_star37_star36_ta_repeat",
    "non_core_uncertain_function_star_allele",
    "novel_uncatalogued_variant",
})


def activity_score(allele1: str, allele2: str) -> float | None:
    """UGT1A1 gene activity score = sum of allele activity values. None if any allele unknown."""
    v1 = ACTIVITY_VALUE.get(allele1)
    v2 = ACTIVITY_VALUE.get(allele2)
    if v1 is None or v2 is None:
        return None
    return v1 + v2


def diplotype_phenotype(allele1: str, allele2: str) -> str:
    """CPIC UGT1A1 metabolizer phenotype from the activity score (Gammal 2016).

    AS 2.0 -> Normal; 1.5 -> Intermediate; 1.0 -> Poor (irinotecan toxicity / Gilbert). Tag-SNP based."""
    score = activity_score(allele1, allele2)
    if score is None:
        return "Indeterminate"
    if score == 2.0:
        return "Normal Metabolizer"
    if score == 1.5:
        return "Intermediate Metabolizer"
    if score == 1.0:
        return "Poor Metabolizer"
    return "Indeterminate"
