"""CYP3A5 curated catalog — defining variants + CPIC expressor/non-expressor phenotype (v0, SNP+indel core).

Fourth PGx gene, first OUTSIDE the chr10 CYP2C cluster (CYP3A5 is chr7q22.1). Highest-value tacrolimus
(transplant) pharmacogene. Validated against the REAL GeT-RM CDC reference-material consensus (the
free multi-lab truth beyond the ursaPGx 4-gene set — see wiki/pgx_free_label_sources_research_2026-07-05).

PROVENANCE (grounded, NO fabrication):
  * Defining variants + GRCh38 coords VERIFIED via Ensembl REST (forward-strand ref>alt; CYP3A5 is on the
    chr7 MINUS strand, so the coding cDNA change is the reverse-complement):
      *3  rs776746     splice (intron-3 6986A>G, cryptic splice)  chr7:99672916 T>C  -- NO function
      *6  rs10264272   c.624G>A (exon-7 skip splice defect)       chr7:99665212 C>T  -- NO function
      *7  rs41303343   c.1035_1036insT (frameshift, 346 Ter)      chr7:99652770 T>TA -- NO function (indel)
    ref/alt orientation EMPIRICALLY CONFIRMED on the 1000G panel: *3 ALT (C) global freq 0.615 (CYP3A5*3
    is the common non-functional allele, esp. in Europeans → high), *6 ALT (T) freq 0.047 (African-specific
    → low). A reversed ref/alt would invert those frequencies.
  * Allele FUNCTION + phenotype terms: CPIC CYP3A5/tacrolimus guideline (Birdwell 2015, updated) — *1
    normal function (EXPRESSOR); *3/*6/*7 no function; the diplotype maps to Normal/Intermediate/Poor
    Metabolizer (expressor / intermediate / non-expressor).

HONESTY TIER:
  * Star-allele CALLING is INDEPENDENTLY validatable vs the GeT-RM CDC multi-lab consensus (a free,
    wet-lab-derived truth set). caller_is_independent_baseline=True for calling. Validation is UNDERPOWERED
    (only ~8-9 GeT-RM CYP3A5 samples overlap 1000G) — reported honestly, never as a large-N number.
  * Metabolizer PHENOTYPE is FAITHFUL-TO-CPIC (assigned from the diplotype via CPIC, not a measured
    tacrolimus-PK phenotype). phenotype_is_faithful_to_cpic=True. Reference tool: PharmCAT.

SCOPE (v0 — honest): CORE *3/*6/*7 + *1 (the CPIC-relevant no-function alleles; *3/*6/*7 cover the common
functional variation, *6/*7 being the African no-function alleles a *3-only caller would mis-read as *1).
Rarer CYP3A5 alleles (*8/*9 etc.) mis-called *1 (documented residual; no sentinel layer). NOT a clinical tool.
"""
from __future__ import annotations

from dna_decode.pgx.cyp2c19_catalog import DefiningVariant, SentinelVariant  # reuse the dataclasses

GENE = "CYP3A5"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

# Core defining variants (v0). GRCh38 forward-strand ref>alt, VERIFIED via Ensembl REST + AF-confirmed.
# *7 is an INSERTION (VCF left-aligned to the preceding base: pos 99652770 T>TA); the caller matches
# d.alt="TA" against the VCF ALT natively (alts.index) — no single-base assumption.
CORE_DEFINING: list[DefiningVariant] = [
    DefiningVariant("*3", "rs776746", "7", 99672916, "T", "C", "6986A>G (splice)"),
    DefiningVariant("*6", "rs10264272", "7", 99665212, "C", "T", "c.624G>A (splice)"),
    DefiningVariant("*7", "rs41303343", "7", 99652770, "T", "TA", "c.1035_1036insT (frameshift)"),
]

# No sentinel layer in v0 (the *3/*6/*7 no-function alleles ARE the core set). Documented residual.
# Non-core-allele SENTINELS: a proven non-core CYP3A5 allele the *1/*3/*6/*7 core proxy cannot resolve
# WITHHOLDS the phenotype rather than a silent *1 (expressor) mis-call -- the wrong direction for a non-
# expressor allele (tacrolimus dosing). Each row's defining SNP is sourced VERBATIM from the PharmCAT/CLINPGX
# CYP3A5_translation.json (GRCh38.p14, chr7); each (rsid -> GRCh38 pos, ref, alt) machine-verified vs Ensembl
# (4/4 OK, 2026-07-28). None of *8/*9/*10/*11's sites are defined by a core allele (*3/*6/*7) -> no false-
# withhold; accounted_by_core stays None. UNDERPOWERED VALIDATION: the GeT-RM 1000G overlap has 0 non-core
# CYP3A5 carriers (n=9), so concordance shows 0 withheld -- the sentinels are correct+safe but unexercised on
# this cohort (a powered cohort would validate the withhold; core concordance is UNCHANGED = no regression).
SENTINELS: list[SentinelVariant] = [
    SentinelVariant("rs55817950", "7", 99676198, "G", "A", "*8", "CYP3A5*8 non-core"),
    SentinelVariant("rs28383479", "7", 99660516, "C", "T", "*9", "CYP3A5*9 non-core"),
    SentinelVariant("rs150999943", "7", 99676142, "C", "A", "*10", "CYP3A5*10 non-core"),
    SentinelVariant("rs762725013", "7", 99652773, "G", "C", "*11", "CYP3A5*11 non-core"),
]

# CPIC allele FUNCTION: *1 EXPRESSOR (normal), *3/*6/*7 non-functional.
ALLELE_FUNCTION: dict[str, str] = {
    "*1": "normal",   # expressor
    "*3": "none",
    "*6": "none",
    "*7": "none",
}

# CPIC CYP3A5 function-pair -> metabolizer phenotype (expressor terminology in parens).
_PHENOTYPE_BY_FUNCTION_PAIR: dict[tuple[str, str], str] = {
    ("normal", "normal"): "Normal Metabolizer",        # *1/*1 — expressor
    ("none", "normal"): "Intermediate Metabolizer",    # *1/*3, *1/*6, *1/*7 — intermediate
    ("none", "none"): "Poor Metabolizer",              # *3/*3, *6/*7, ... — non-expressor
}

PHENOTYPE_ABBREV: dict[str, str] = {
    "Normal Metabolizer": "NM",           # CYP3A5 expressor
    "Intermediate Metabolizer": "IM",
    "Poor Metabolizer": "PM",             # CYP3A5 non-expressor
    "Indeterminate": "IND",
}

UNDETECTABLE = sorted({
    "non_core_star_allele",            # rare CYP3A5 alleles beyond *3/*6/*7 -> mis-called *1
    "cnv_or_gene_deletion",
    "novel_uncatalogued_variant",
    "phase_ambiguity_unphased_input",
})


def diplotype_phenotype(allele1: str, allele2: str) -> str:
    """CPIC CYP3A5 metabolizer phenotype for a diplotype. 'Indeterminate' if either allele unknown.
    Order-independent. *1/*1 = Normal (expressor); *1/no-function = Intermediate; no-function/no-function
    = Poor (non-expressor)."""
    f1 = ALLELE_FUNCTION.get(allele1, "unknown")
    f2 = ALLELE_FUNCTION.get(allele2, "unknown")
    if "unknown" in (f1, f2):
        return "Indeterminate"
    return _PHENOTYPE_BY_FUNCTION_PAIR[tuple(sorted((f1, f2)))]
