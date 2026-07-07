"""NUDT15 curated catalog — defining variant + CPIC ACTIVITY-SCORE phenotype (v0, SNP-defined core).

NUDT15 is the second thiopurine-toxicity pharmacogene (pairs with TPMT): NUDT15-deficient patients suffer
severe myelosuppression from standard azathioprine / mercaptopurine / thioguanine doses. CPIC (Relling 2019)
assigns a metabolizer phenotype from the diplotype's allele functions — the same activity-score shape as
TPMT/CYP2C9/DPYD. The dominant no-function allele is *3 (rs116855232, c.415C>T p.R139C), common in East
Asian / Hispanic populations (~10% EAS) and the highest-VOI actionable NUDT15 variant.

PROVENANCE (grounded, no fabrication):
  * Defining variant + GRCh38 coord: PharmVar/CPIC + dbSNP, VERIFIED via Ensembl GRCh38 REST 2026-07-07:
      *3  rs116855232  c.415C>T p.R139C  chr13:48045719 C>T  -- NO function, activity 0.0
  * Activity values + score->phenotype: CPIC NUDT15 guideline (Relling 2019) — *1=1.0 normal, *3=0.0 no
    function; AS 2.0 = Normal Metabolizer, AS 1.0 = Intermediate (reduce thiopurine), AS 0.0 = Poor (drastically
    reduce / avoid).

SCOPE (v0 -- honest): the dominant actionable no-function haplotype *3 + *1 (reference). NUDT15*2 SHARES
rs116855232 with *3 (it is *3 + an additional c.52G>A), so the single-SNP proxy calls *2 as *3 — a documented
residual, exactly like TPMT*3A needing its second SNP; the phenotype (no-function) is IDENTICAL for *2/*3 so
the CPIC metabolizer call is unaffected. Rarer *4/*5/*6/*9 (uncertain/no-function) are NOT covered -> called
*1. NOT a clinical tool.
"""
from __future__ import annotations

from dna_decode.pgx.cyp2c19_catalog import DefiningVariant, SentinelVariant  # reuse the dataclasses

GENE = "NUDT15"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

CORE_DEFINING: list[DefiningVariant] = [
    DefiningVariant("*3", "rs116855232", "13", 48045719, "C", "T", "c.415C>T p.R139C (no function)"),
]

SENTINELS: list[SentinelVariant] = []

ACTIVITY_VALUE: dict[str, float] = {"*1": 1.0, "*3": 0.0}

PHENOTYPE_ABBREV = {"Normal Metabolizer": "NM", "Intermediate Metabolizer": "IM",
                    "Poor Metabolizer": "PM", "Indeterminate": "IND"}

UNDETECTABLE = sorted({
    "nudt15_star2_vs_star3_second_snp",   # *2 shares rs116855232 -> called *3 (same no-function phenotype)
    "non_core_uncertain_function_star_allele",
    "novel_uncatalogued_variant",
    "cnv_or_gene_deletion",
})


def activity_score(allele1: str, allele2: str) -> float | None:
    """NUDT15 gene activity score = sum of allele activity values. None if any allele unknown."""
    v1 = ACTIVITY_VALUE.get(allele1)
    v2 = ACTIVITY_VALUE.get(allele2)
    if v1 is None or v2 is None:
        return None
    return v1 + v2


def diplotype_phenotype(allele1: str, allele2: str) -> str:
    """CPIC NUDT15 thiopurine metabolizer phenotype from the activity score (Relling 2019).

    AS 2.0 -> Normal; 1.0 -> Intermediate (reduce thiopurine); 0.0 -> Poor (drastically reduce / avoid)."""
    score = activity_score(allele1, allele2)
    if score is None:
        return "Indeterminate"
    if score == 2.0:
        return "Normal Metabolizer"
    if score == 1.0:
        return "Intermediate Metabolizer"
    if score == 0.0:
        return "Poor Metabolizer"
    return "Indeterminate"
