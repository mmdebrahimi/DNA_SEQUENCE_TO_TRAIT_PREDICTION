"""CYP2C19 curated catalog — defining variants + CPIC diplotype->phenotype (v0, SNP-defined core).

PROVENANCE (no fabrication, no build-on-unverified):
  * Defining variants + GRCh38 coordinates: PharmVar CYP2C19 allele definitions (pharmvar.org) +
    dbSNP, cross-checked against Botton et al. 2021 "PharmVar GeneFocus: CYP2C19" (Clin Pharmacol Ther)
    and Gaedigk et al. 2022 (J Mol Diagn, the GeT-RM CYP2C characterization). Positions GRCh38.13 chr10:
      *2  rs4244285  (c.681G>A)   chr10:94781859 G>A  -- splice defect, NO function
      *3  rs4986893  (c.636G>A)   chr10:94780653 G>A  -- p.W212X premature stop, NO function
      *17 rs12248560 (c.-806C>T)  chr10:94761900 C>T  -- promoter GATA site, INCREASED function
  * Allele FUNCTION assignments + diplotype->phenotype: CPIC standardized terms (Caudle et al. 2020,
    Genet Med) + the CPIC/PharmGKB CYP2C19 allele-functionality + diplotype-phenotype tables
    (cpicpgx.org / pharmgkb.org/page/cyp2c19RefMaterials).

HONESTY TIER (load-bearing — two distinct claims, never conflated):
  * Star-allele CALLING is INDEPENDENTLY validatable vs the GeT-RM consensus panel (a free, measured/
    consensus label) -> the genuine number (the P3 cohort run). caller_is_independent_baseline=True for
    the calling step.
  * Metabolizer PHENOTYPE (PM/IM/NM/RM/UM) is FAITHFUL-TO-CPIC -- ASSIGNED from the diplotype via CPIC's
    table, NOT a measured probe-drug PK phenotype (those are not free at scale). Like the serotype/ktype
    cells. caller_is_independent_baseline=False for the phenotype step. The reference tool is PharmCAT.

SCOPE (v0 — honest):
  * Core SNP-defined alleles only: *2, *3, *17, and *1 (= the reference functional allele when no core
    defining variant is present). Many more CYP2C19 star alleles exist (*4-*35, ...); a sample carrying a
    non-core defining variant is NOT *1 -- v0 cannot see it (a documented blind spot, surfaced as a flag).
  * The diplotype is the deterministic combination of the called alleles; CPIC's *2/*17 = Intermediate
    Metabolizer is a PROVISIONAL CPIC assignment (increased *17 does not fully compensate no-function *2).
"""
from __future__ import annotations

from dataclasses import dataclass

GENE = "CYP2C19"
ASSEMBLY = "GRCh38"


@dataclass(frozen=True)
class DefiningVariant:
    star: str        # star allele, e.g. "*2"
    rsid: str        # dbSNP id
    chrom: str       # GRCh38, normalized (no "chr" prefix), e.g. "10"
    pos: int         # GRCh38 1-based position
    ref: str
    alt: str
    cdna: str        # HGVS c. nomenclature


# Core SNP-defined CYP2C19 alleles (v0). One defining SNP each, relative to the reference.
CORE_DEFINING: list[DefiningVariant] = [
    DefiningVariant("*2", "rs4244285", "10", 94781859, "G", "A", "c.681G>A"),
    DefiningVariant("*3", "rs4986893", "10", 94780653, "G", "A", "c.636G>A"),
    DefiningVariant("*17", "rs12248560", "10", 94761900, "C", "T", "c.-806C>T"),
]

# The reference allele when no core defining variant is present on a haplotype.
# v0 CAVEAT: *38 is the true variant-free reference allele; *1 carries rs3758581 (c.991A>G, chr10:94842866
# GRCh38, verified at dbSNP NC_000010.11:g.94842866A>C). Both *1 and *38 are NORMAL function, so the
# phenotype is identical -- v0 reports the normal-function reference allele as *1 (the clinical convention)
# and records this caveat. Distinguishing *1 from *38 requires rs3758581 (a v0.1 refinement;
# phenotype-irrelevant).
REFERENCE_ALLELE = "*1"

# CPIC standardized allele FUNCTION (Caudle 2020).
ALLELE_FUNCTION: dict[str, str] = {
    "*1": "normal",       # normal function (reference)
    "*2": "none",         # no function
    "*3": "none",         # no function
    "*17": "increased",   # increased function
}

# SENTINEL non-core variants (v0.1) -- NOT called as star alleles, but their presence proves a NON-CORE
# allele the single-SNP core proxy CANNOT resolve, so the caller WITHHOLDS the phenotype rather than
# emit a confident mis-call. Targeted at the two demonstrated aliasing failures (NOT a full AMP Tier-2
# screen): *4 (aliases *17 via the shared rs12248560) and *35 (aliases *1 via rs12769205). Coordinates
# GRCh38, confirmed in PharmCAT fixtures + the 1000 Genomes panel.
@dataclass(frozen=True)
class SentinelVariant:
    rsid: str
    chrom: str
    pos: int
    ref: str
    alt: str
    implies: str        # the non-core allele family this ALT signals
    note: str


SENTINELS: list[SentinelVariant] = [
    SentinelVariant("rs28399504", "10", 94762706, "A", "G", "*4",
                    "*4 initiation-codon SNP; *4b ALSO carries the *17 SNP rs12248560, so a core '*17' "
                    "call on this haplotype is actually *4b (reduced/no function, not increased)."),
    SentinelVariant("rs12769205", "10", 94775367, "A", "G", "*35",
                    "shared by *2 (with rs4244285) and *35 (alone); an rs12769205 copy in EXCESS of the "
                    "rs4244285 (*2) copies marks a *35 haplotype the core proxy mis-calls as *1."),
]

# CPIC CYP2C19 function-pair -> metabolizer phenotype (keyed by the SORTED function pair).
# Grounded vs the CPIC diplotype-phenotype table (Caudle 2020 standardized terms).
_PHENOTYPE_BY_FUNCTION_PAIR: dict[tuple[str, str], str] = {
    ("normal", "normal"): "Normal Metabolizer",          # *1/*1
    ("increased", "normal"): "Rapid Metabolizer",        # *1/*17
    ("increased", "increased"): "Ultrarapid Metabolizer",  # *17/*17
    ("none", "normal"): "Intermediate Metabolizer",      # *1/*2, *1/*3
    ("increased", "none"): "Intermediate Metabolizer",   # *2/*17, *3/*17 (provisional IM per CPIC)
    ("none", "none"): "Poor Metabolizer",                # *2/*2, *2/*3, *3/*3
}

PHENOTYPE_ABBREV: dict[str, str] = {
    "Poor Metabolizer": "PM",
    "Intermediate Metabolizer": "IM",
    "Normal Metabolizer": "NM",
    "Rapid Metabolizer": "RM",
    "Ultrarapid Metabolizer": "UM",
    "Indeterminate": "IND",
}

# Mechanisms / alleles a CORE-SNP CYP2C19 scan CANNOT see -- a "*1" (no-core-variant) call cannot rule
# these out. The PGx analogue of the bacterial efflux / HIV minor-mutation blind spot.
UNDETECTABLE = sorted({
    "non_core_star_allele",            # *4-*35 etc. not in the v0 core SNP set -> mis-called as *1
    "cnv_or_gene_deletion",            # CYP2C19 structural variants (rare)
    "novel_uncatalogued_variant",
    "phase_ambiguity_unphased_input",  # >=2 het sites without phasing -> cis/trans ambiguity
})


def diplotype_phenotype(allele1: str, allele2: str) -> str:
    """CPIC metabolizer phenotype for a diplotype. 'Indeterminate' if either allele's function is unknown
    (e.g. a compound/non-core haplotype). Order-independent."""
    f1 = ALLELE_FUNCTION.get(allele1, "unknown")
    f2 = ALLELE_FUNCTION.get(allele2, "unknown")
    if "unknown" in (f1, f2):
        return "Indeterminate"
    return _PHENOTYPE_BY_FUNCTION_PAIR[tuple(sorted((f1, f2)))]
