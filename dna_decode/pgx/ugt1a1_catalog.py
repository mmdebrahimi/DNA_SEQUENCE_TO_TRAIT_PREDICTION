"""UGT1A1 curated catalog — SNP-defined core + CPIC activity-score phenotype (v0, tag-SNP surface).

UGT1A1 drives irinotecan toxicity + atazanavir hyperbilirubinemia (and Gilbert syndrome). CPIC (Gammal 2016)
assigns a metabolizer phenotype from the diplotype's allele functions — the activity-score shape again.

*** THE STRUCTURAL WALL — MEASURED 2026-09-22, AND IT IS PANEL-DEPENDENT, NOT ABSOLUTE ***
The MAJOR reduced-function allele UGT1A1*28 is a TA-DINUCLEOTIDE REPEAT in the promoter (TA6 = *1 normal,
TA7 = *28 reduced) — a SHORT TANDEM REPEAT, NOT a SNP. So v0 uses the TAG-SNP approach labs deploy:
**rs887829 (*80, C>T) as the SNP-callable proxy** for *28 reduced-function status. This is a
WRAPPER-OVER-LD-TAG cell (same posture as the project's HLA tag-SNP cells).

**This file previously asserted that a short-read SNP VCF "CANNOT reliably call the TA repeat length" and
that the true *28 length is "UNASSESSED". Both were over-claims, corrected after measurement**
(`wiki/ugt1a1_getrm_concordance_2026-09-22.{md,json}`, bar frozen beforehand). Two findings on 65 1000G
samples carrying CDC GeT-RM multi-lab consensus diplotypes (32 *28 carriers):

  1. **The repeat IS directly callable FROM THIS PANEL.** The 1000G 30x phased panel represents it at
     chr2:233760233 as indel alleles — C>CAT (+1 TA = TA7 = *28), C>CATAT (+2 TA = TA8 = *37),
     CAT>C (-1 TA = TA5 = *36). A direct call scored sens 1.00 / spec 1.00 / PPV 1.00 on *28 carriership
     and 64/65 exact zygosity — it STRICTLY DOMINATES the tag. **SCOPE: this shows the repeat is callable
     from the 1000G 30x phased panel, NOT from an arbitrary short-read SNP VCF.** A VCF that does not
     genotype the indel still cannot resolve it; the wall is panel-dependent, and the direct path is NOT
     wired into this catalog (v0 remains the tag).
  2. **The tag WORKS, and it tags the CLASS not the allele.** Contrary to the HLA-tag precedent, sens for
     *28 carriership was 1.00 overall and 1.00 in each of EUR (n=26) / AFR (n=21) / EAS (n=16) — the
     feared non-EUR LD breakdown did NOT appear at this n. Its 3 apparent false positives are all
     haplotypes carrying *37 (TA8), i.e. **rs887829-T marks the reduced-function promoter class (TA7 +
     TA8), not *28 specifically.** Scored against reduced FUNCTION — the axis irinotecan dosing actually
     turns on — the tag is sens 1.00 / spec 1.00 / PPV 1.00 (post-hoc framing, labelled as such).

Named residual (unchanged in kind, now sized): one sample (NA20509) is *28/*28 by GeT-RM but 1 copy by BOTH
independent assays here; that discrepancy is unadjudicated.

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
    "UGT1A1*28 is a promoter TA-repeat (STR), NOT a SNP. v0 calls it via rs887829 (*80) as an LD-tag PROXY "
    "-- a tag-SNP wrapper, NOT a direct repeat call. MEASURED vs GeT-RM CDC consensus on 65 1000G samples "
    "(2026-09-22): the tag reached sens 1.00 for *28 carriership (1.00 in each of EUR/AFR/EAS) and "
    "64/65 exact zygosity. It marks the reduced-function promoter CLASS -- *28 (TA7) AND *37 (TA8) -- not "
    "*28 specifically, so a *37 carrier is reported as reduced-function (clinically right, but it is not "
    "a *28-specific call). Do NOT read a *80 result as a resolved TA-repeat LENGTH.")

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
