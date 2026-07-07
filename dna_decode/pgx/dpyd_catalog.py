"""DPYD curated catalog — defining variants + CPIC ACTIVITY-SCORE phenotype (v0, SNP-defined core).

DPYD (dihydropyrimidine dehydrogenase) is THE fluoropyrimidine-toxicity pharmacogene: DPD-deficient
patients suffer severe/fatal toxicity from standard 5-fluorouracil / capecitabine doses. CPIC gives a hard
dosing guideline keyed on a GENE ACTIVITY SCORE — the same activity-score method as CYP2C9/CYP2D6: each
allele carries an activity VALUE, the diplotype score = sum, and the score maps to a metabolizer phenotype.

PROVENANCE (grounded, no fabrication):
  * Defining variants + GRCh38 coords: PharmVar/CPIC + dbSNP, coords VERIFIED via Ensembl GRCh38 REST
    2026-07-07 (DPYD is chr1 MINUS strand, so genomic ref>alt = reverse-complement of the coding cDNA):
      *2A       rs3918290   c.1905+1G>A splice   chr1:97450058 C>T   -- NO function,       activity 0.0
      *13       rs55886062  c.1679T>G p.I560S    chr1:97515787 A>C   -- NO function,       activity 0.0
      c.2846A>T rs67376798  c.2846A>T p.D949V    chr1:97082391 T>A   -- DECREASED function, activity 0.5
      HapB3     rs75017182  c.1129-5923C>G       chr1:97579893 G>C   -- DECREASED function, activity 0.5
  * Activity values + score->phenotype thresholds: CPIC DPYD guideline (Amstutz 2018, updated) — the four
    clinically-actionable variants above; AS 2.0 = Normal Metabolizer, AS 1.0-1.5 = Intermediate (reduce
    fluoropyrimidine dose ~50%), AS 0.0-0.5 = Poor (avoid / drastically reduce).

SCOPE (v0 -- honest): the FOUR CPIC-actionable haplotypes (the "big four" DPD-deficiency variants) + *1
(reference). DPYD has 20+ rarer star alleles of uncertain/normal function; CPIC does NOT dose on those, and
v0 has NO sentinel layer, so any non-core variant is called *1 (Normal) — CPIC's own posture (only the
actionable four change dosing). HapB3 is tagged here by its CAUSAL intronic SNP rs75017182 (c.1129-5923C>G),
which a WGS VCF calls directly. NOT a clinical tool.
"""
from __future__ import annotations

from dna_decode.pgx.cyp2c19_catalog import DefiningVariant, SentinelVariant  # reuse the dataclasses

GENE = "DPYD"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

# The four CPIC-actionable DPD-deficiency haplotypes, each defined by its key SNP (GRCh38, minus strand).
CORE_DEFINING: list[DefiningVariant] = [
    DefiningVariant("*2A", "rs3918290", "1", 97450058, "C", "T", "c.1905+1G>A splice (no function)"),
    DefiningVariant("*13", "rs55886062", "1", 97515787, "A", "C", "c.1679T>G p.I560S (no function)"),
    DefiningVariant("c.2846A>T", "rs67376798", "1", 97082391, "T", "A", "c.2846A>T p.D949V (decreased)"),
    DefiningVariant("HapB3", "rs75017182", "1", 97579893, "G", "C", "c.1129-5923C>G intronic causal (decreased)"),
]

# v0: NO sentinel layer — CPIC doses ONLY on the four core haplotypes; rarer DPYD alleles are
# uncertain/normal-function and non-actionable, so calling them *1 (Normal) matches CPIC's own posture.
SENTINELS: list[SentinelVariant] = []

# CPIC activity VALUE per allele (sum -> gene activity score).
ACTIVITY_VALUE: dict[str, float] = {"*1": 1.0, "*2A": 0.0, "*13": 0.0, "c.2846A>T": 0.5, "HapB3": 0.5}

PHENOTYPE_ABBREV = {"Normal Metabolizer": "NM", "Intermediate Metabolizer": "IM",
                    "Poor Metabolizer": "PM", "Indeterminate": "IND"}

UNDETECTABLE = sorted({
    "non_core_uncertain_function_star_allele",  # 20+ rarer DPYD alleles CPIC does not dose on -> called *1
    "novel_uncatalogued_variant",
    "cnv_or_gene_deletion",
})


def activity_score(allele1: str, allele2: str) -> float | None:
    """DPYD gene activity score = sum of allele activity values. None if any allele unknown."""
    v1 = ACTIVITY_VALUE.get(allele1)
    v2 = ACTIVITY_VALUE.get(allele2)
    if v1 is None or v2 is None:
        return None
    return v1 + v2


def diplotype_phenotype(allele1: str, allele2: str) -> str:
    """CPIC DPYD metabolizer phenotype from the gene activity score (Amstutz 2018 thresholds).

    AS 2.0 -> Normal; 1.0-1.5 -> Intermediate (reduce fluoropyrimidine ~50%); 0.0-0.5 -> Poor (avoid)."""
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
