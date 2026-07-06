"""CYP2D6 curated catalog — SNP-defined star alleles + CPIC ACTIVITY-SCORE phenotype (v0).

The last major pharmacogene, and the project's hardest PGx cell — CYP2D6 is the canonical
"structural-variant" gene (gene deletions *5, duplications *xN, CYP2D6-CYP2D7 hybrids *13/*36/*68).

R2 REFRAME (load-bearing — derived, not asserted): a full CYP2D6 "structural caller" is INFEASIBLE from a
phased SNP VCF (structural detection needs read-depth / breakpoint analysis on a BAM/CRAM — the
Cyrius/Aldy/StellarPGx approach). This cell is therefore a CYP2D6 **SNP-defined star-allele** caller,
validated on the SNP-decodable GeT-RM subset. Structural alleles are NOT VCF-decodable and are NOT
"withheld" (the SNP proxy cannot even see them) -> they are EXCLUDED from the scored denominator and may be
SILENTLY MIS-CALLED; every record carries `cnv_hybrid_unassessed=true`. Never imply full CYP2D6 typing.

PROVENANCE (grounded, NO fabrication, NO build-on-unverified):
  * Defining variants + GRCh38 coords VERIFIED via the NCBI variation service (refsnp SPDI on NC_000022.11)
    AND empirically confirmed against the fetched 1000G 30x panel record + ALT-orientation (CYP2D6 is on the
    chr22 MINUS strand, so the forward-strand genomic ref>alt is the reverse-complement of the coding cDNA
    change). AF on 1000G (global, n=6404 haplotypes) confirms each ALT orientation:
      *4  rs3892097  c.1846G>A splice   chr22:42128945 C>T   AF 0.094  (common EUR no-function)
      *10 rs1065852  c.100C>T  P34S     chr22:42130692 G>A   AF 0.228  (common; ALSO on the *4 background)
      *2  rs16947    c.2851C>T R296C    chr22:42127941 G>A   AF 0.366  (the *2 background SNP)
      486 rs1135840  c.4181G>C S486T    chr22:42126611 C>G   AF 0.594  (background co-marker, majority allele)
      *17 rs28371706 c.1023C>T T107I    chr22:42129770 G>A   AF 0.064  (African-common)
      *41 rs28371725 c.2988G>A splice   chr22:42127803 C>T   AF 0.062
      *29 rs59421388 c.1659G>A V136M    chr22:42127608 C>T   AF 0.031  (African)
      *35 rs769258   c.31G>A   V11M     chr22:42130761 C>T   AF 0.017  (normal function, *2-like)
      *3  rs35742686 c.2549delA fs      chr22:42128241 CT>C  AF 0.005  (1000G left-anchored indel)
      *6  rs5030655  c.1707delT fs      chr22:42129083 CA>C  AF 0.005  (1000G left-anchored indel)
      *9  rs5030656  c.2615_2617del     chr22:42128173 CCTT>C AF 0.007 (K281del; 1000G left-anchored indel)
  * Activity VALUES + score->phenotype thresholds: CPIC "Standardizing CYP2D6 Genotype to Phenotype
    Translation" consensus (Caudle 2020, Clin Transl Sci) — *1/*2/*35 = 1.0 (normal); *9/*17/*29 = 0.5,
    *10/*41 = 0.25 (decreased); *3/*4/*6 = 0.0 (no function). AS 0 = PM, 0<AS<=1.0 = IM,
    1.25<=AS<=2.25 = NM (UM > 2.25 requires a duplication, which v0 cannot call).

CALLING MODEL (why NOT the subset-largest compound resolver): CYP2D6 star alleles sit on a shared SNP
BACKGROUND (*2's 2851/486 is carried by *4/*17/*29/*35/*41; *4 carries *10's 100C>T). A subset-largest
resolver over-matches on that background. This cell uses a PRIORITY-ORDERED per-haplotype resolver instead:
the most-specific defining SNP wins (1846 *4 before 100 *10; every allele-specific SNP before the 2851 *2
background), so a shared background SNP never mis-calls the more-specific allele. See `cyp2d6_caller.py`.

SCOPE (v0 — honest): CORE SNP-defined {*2,*3,*4,*6,*9,*10,*17,*29,*35,*41} + *1. Non-core SNP alleles
(*14/*15/*21/*40/*46) are mis-called (documented residual; no sentinel layer in v0, mirroring the C8/C9/
C3A5/TPMT v0 arc). Structural alleles (*5/*13/*36/*61/*63/*68/*xN) are BAM-required and out of scope
(cnv_hybrid_unassessed). NOT a clinical tool.
"""
from __future__ import annotations

from dna_decode.pgx.cyp2c19_catalog import DefiningVariant  # reuse the dataclass

GENE = "CYP2D6"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

# Component SNPs (the `.star` field carries the COMPONENT TAG used by the priority resolver, not a star).
# GRCh38 forward-strand ref>alt, NCBI-verified + AF-confirmed on 1000G. Indels use the EXACT 1000G
# left-anchored representation (REF/ALT verified in the fetched panel record).
COMPONENTS: list[DefiningVariant] = [
    DefiningVariant("1846A",    "rs3892097",  "22", 42128945, "C",    "T", "c.1846G>A splice (*4)"),
    DefiningVariant("100T",     "rs1065852",  "22", 42130692, "G",    "A", "c.100C>T P34S (*10; also *4 bkg)"),
    DefiningVariant("2851T",    "rs16947",    "22", 42127941, "G",    "A", "c.2851C>T R296C (*2 background)"),
    DefiningVariant("486T",     "rs1135840",  "22", 42126611, "C",    "G", "c.4181G>C S486T (background co-marker)"),
    DefiningVariant("1023T",    "rs28371706", "22", 42129770, "G",    "A", "c.1023C>T T107I (*17)"),
    DefiningVariant("2988A",    "rs28371725", "22", 42127803, "C",    "T", "c.2988G>A splice (*41)"),
    DefiningVariant("1659A",    "rs59421388", "22", 42127608, "C",    "T", "c.1659G>A V136M (*29)"),
    DefiningVariant("31A",      "rs769258",   "22", 42130761, "C",    "T", "c.31G>A V11M (*35)"),
    DefiningVariant("2549delA", "rs35742686", "22", 42128241, "CT",   "C", "c.2549delA frameshift (*3)"),
    DefiningVariant("1707delT", "rs5030655",  "22", 42129083, "CA",   "C", "c.1707delT frameshift (*6)"),
    DefiningVariant("2615del",  "rs5030656",  "22", 42128173, "CCTT", "C", "c.2615_2617delAAG K281del (*9)"),
]

# `CORE_DEFINING` kept as the component list so report/registry wiring that reads a defining list works.
CORE_DEFINING = COMPONENTS

# PRIORITY-ORDERED star resolution: (star, defining_component_tag). The FIRST tag present on a haplotype
# wins. Ordered most-specific -> background so a shared background SNP never mis-calls a more-specific
# allele: 1846(*4) before 100(*10) (every *4 haplotype carries 100C>T); every allele-specific SNP before
# the 2851(*2) background (every *17/*41/*29/*35 haplotype carries 2851). 486T is a background co-marker
# ONLY -> never a standalone star determinant.
STAR_PRIORITY: list[tuple[str, str]] = [
    ("*4",  "1846A"),
    ("*3",  "2549delA"),
    ("*6",  "1707delT"),
    ("*9",  "2615del"),
    ("*17", "1023T"),
    ("*41", "2988A"),
    ("*29", "1659A"),
    ("*35", "31A"),
    ("*10", "100T"),
    ("*2",  "2851T"),
]

# Background tags (not allele-specific): present on many stars, never decide the star alone.
BACKGROUND_TAGS: frozenset[str] = frozenset({"2851T", "486T", "100T"})
# Allele-specific tags: >=2 of these on ONE haplotype is a data anomaly (flagged, not silently resolved).
SPECIFIC_TAGS: frozenset[str] = frozenset({"1846A", "1023T", "2988A", "1659A", "31A",
                                           "2549delA", "1707delT", "2615del"})

# CPIC standardized activity VALUE per allele (Caudle 2020). Diplotype activity score = sum.
ACTIVITY_VALUE: dict[str, float] = {
    "*1": 1.0, "*2": 1.0, "*35": 1.0,   # normal function
    "*9": 0.5, "*17": 0.5, "*29": 0.5,  # decreased function
    "*10": 0.25, "*41": 0.25,           # decreased function (CPIC-downgraded)
    "*3": 0.0, "*4": 0.0, "*6": 0.0,    # no function
}

PHENOTYPE_ABBREV: dict[str, str] = {
    "Ultrarapid Metabolizer": "UM",
    "Normal Metabolizer": "NM",
    "Intermediate Metabolizer": "IM",
    "Poor Metabolizer": "PM",
    "Indeterminate": "IND",
}

# Alleles / mechanisms a CORE-SNP CYP2D6 scan CANNOT see (the honest blind-spot list).
UNDETECTABLE = sorted({
    "structural_gene_deletion",         # *5 — needs read-depth (BAM); may be SILENTLY mis-called
    "structural_gene_duplication",      # *xN / *2x2 — needs copy-number (BAM)
    "cyp2d6_cyp2d7_hybrid",             # *13 / *36 / *61 / *63 / *68 — needs breakpoint analysis (BAM)
    "non_core_snp_star_allele",         # *14/*15/*21/*40/*46 beyond the core set -> mis-called
    "novel_uncatalogued_variant",
    "phase_ambiguity_unphased_input",   # >=2 unphased het specific sites -> cis/trans ambiguous
})


def activity_score(allele1: str, allele2: str) -> float | None:
    """Diplotype activity score = sum of allele activity values. None if any allele unknown."""
    v1 = ACTIVITY_VALUE.get(allele1)
    v2 = ACTIVITY_VALUE.get(allele2)
    if v1 is None or v2 is None:
        return None
    return round(v1 + v2, 3)


def diplotype_phenotype(allele1: str, allele2: str) -> str:
    """CPIC CYP2D6 metabolizer phenotype from the activity score (Caudle 2020 consensus thresholds).
    Order-independent. 'Indeterminate' if either allele unknown. UM (AS>2.25) requires a duplication,
    which the v0 SNP caller cannot call -> never emitted here."""
    score = activity_score(allele1, allele2)
    if score is None:
        return "Indeterminate"
    if score == 0.0:
        return "Poor Metabolizer"
    if score <= 1.0:
        return "Intermediate Metabolizer"
    if score <= 2.25:
        return "Normal Metabolizer"
    return "Ultrarapid Metabolizer"
