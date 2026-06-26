"""CYP2C9 curated catalog — defining variants + CPIC ACTIVITY-SCORE phenotype (v0, SNP-defined core).

Second PGx gene (warfarin pair with VKORC1). Unlike CYP2C19 (function-pair), CYP2C9 uses CPIC's
ACTIVITY-SCORE method: each allele carries an activity VALUE, the diplotype score = sum, and the score
maps to a metabolizer phenotype.

PROVENANCE (grounded, no fabrication):
  * Defining variants + GRCh38 coords: PharmVar/CPIC + dbSNP (verified NC_000010.11):
      *2  rs1799853  (c.430C>T, p.R144C)   chr10:94942290 C>T  -- DECREASED function, activity 0.5
      *3  rs1057910  (c.1075A>C, p.I359L)  chr10:94981296 A>C  -- NO function, activity 0.0
  * Activity values + score->phenotype thresholds: CPIC CYP2C9 allele-functionality + the warfarin
    guideline (Johnson 2017) -- *1=1.0 normal, *2=0.5 decreased, *3=0.0 no function; AS 2.0 = Normal
    Metabolizer, AS 1.0-1.5 = Intermediate, AS 0.0-0.5 = Poor.

SCOPE (v0 -- honest): CORE SNP-defined *2/*3 + *1 (reference). The GeT-RM consensus for CYP2C9 also
contains non-core alleles (*5/*6/*8/*9/*11/*61); v0 has NO sentinel layer for CYP2C9 yet, so a non-core
allele is mis-called *1 (a documented residual surfaced by the GeT-RM validation -- the sentinel layer is
a v0.1 follow-up, mirroring the CYP2C19 arc). NOT a clinical tool.
"""
from __future__ import annotations

from dna_decode.pgx.cyp2c19_catalog import DefiningVariant  # reuse the dataclass

GENE = "CYP2C9"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

CORE_DEFINING: list[DefiningVariant] = [
    DefiningVariant("*2", "rs1799853", "10", 94942290, "C", "T", "c.430C>T"),
    DefiningVariant("*3", "rs1057910", "10", 94981296, "A", "C", "c.1075A>C"),
]

# CPIC activity VALUE per allele (sum -> activity score).
ACTIVITY_VALUE: dict[str, float] = {"*1": 1.0, "*2": 0.5, "*3": 0.0}

PHENOTYPE_ABBREV = {"Normal Metabolizer": "NM", "Intermediate Metabolizer": "IM",
                    "Poor Metabolizer": "PM", "Indeterminate": "IND"}

SENTINELS: list = []   # v0: none (non-core *5/*6/*8/*9/*11 -> mis-called *1; sentinel layer = v0.1)

UNDETECTABLE = sorted({
    "non_core_star_allele",       # *5/*6/*8/*9/*11/*61 etc. -> mis-called *1 (no sentinel in v0)
    "novel_uncatalogued_variant",
    "cnv_or_gene_deletion",
})


def activity_score(allele1: str, allele2: str) -> float | None:
    """Diplotype activity score = sum of allele activity values. None if any allele unknown."""
    v1 = ACTIVITY_VALUE.get(allele1)
    v2 = ACTIVITY_VALUE.get(allele2)
    if v1 is None or v2 is None:
        return None
    return v1 + v2


def diplotype_phenotype(allele1: str, allele2: str) -> str:
    """CPIC CYP2C9 metabolizer phenotype from the activity score (Johnson 2017 thresholds)."""
    score = activity_score(allele1, allele2)
    if score is None:
        return "Indeterminate"
    if score == 2.0:
        return "Normal Metabolizer"
    if score in (1.0, 1.5):
        return "Intermediate Metabolizer"
    if score in (0.0, 0.5):
        return "Poor Metabolizer"
    return "Indeterminate"
