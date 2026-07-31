"""Dog (Canis lupus familiaris) body-SIZE decoder — the curated-catalog paradigm applied to a QUANTITATIVE
visible trait. Where coat colour is a fixed-order EPISTASIS across allele CALLS (dog_coat.py), body size is
the opposite shape: a small ADDITIVE polygenic score over a handful of large-effect SNP dosages.

THE BIOLOGY (the rule this encodes — Sutter 2007 / Rimbault 2013 / Plassais 2019, curated in OMIA 001968 +
002524): a *handful* of loci explain roughly HALF the across-breed variance in adult skeletal size. Each
carries a "big" and a "small" allele acting additively (more big-alleles -> larger dog). The dominant loci:

    IGF1  (CFA15) insulin-like growth factor 1 — the single largest size determinant; intron-2 SNP (in
                  complete LD with a SINE insertion) marks the small-size haplotype. Sutter 2007 Science.
    HMGA2 (CFA10) high-mobility-group AT-hook 2 — 5'UTR SNP, second-largest effect. Rimbault 2013.
    STC2  (CFA4)  stanniocalcin 2 — a SNP ~20 kb downstream lowers size. Rimbault 2013 / Plassais 2019.
    GHR   (CFA4)  growth-hormone receptor — non-synonymous exon SNPs in the extracellular domain. Rimbault 2013.

    (Further secondary loci — IGF1R/CFA3, SMAD2/CFA7 — are NOT in v0; they add smaller increments.)

PINNING PROVENANCE (verified, NOT from memory — the anti-fabrication rail that caught the coat MC1R mispin):
each canFam4 panel variant below was pinned by OMIA/literature canFam3.1 coord -> UCSC canFam3ToCanFam4
liftover (pyliftover) -> exact match in the Darwin's Ark canFam4 `.bim` -> FUNCTIONALLY validated: its
dosage tracks owner-reported height (Q121 z-score) on N=3276 dogs with the sign the literature predicts.
The measured single-SNP r and the combined 4-locus polygenic-score r are recorded (see the wiki artifact).

WHY THIS SUBSTRATE WORKS (where coat colour FAILED): the coat causal variants are indels / SVs / imputation
gaps (CBD103 3bp-del, ASIP SINE, MLPH frameshift) — structurally ABSENT from a biallelic-SNV panel. The
body-size causal variants are SNPs, present in the 29M-SNV imputed panel at their exact lifted positions.
So the SAME free Dryad substrate that could only validate BLACK for coat validates the FULL body-size score.

HONEST SCOPE (load-bearing):
  - This is a RELATIVE size score (a polygenic rank), NOT a calibrated absolute-height predictor. Q121 is a
    covariate-adjusted quantile-normalised z-score, so the validated claim is "more big-alleles -> taller
    RANK", not "this dog is 20 inches". An absolute-height model needs a raw-inches label + breed covariate.
  - v0 = the FOUR dominant loci (IGF1/HMGA2/STC2/GHR). Secondary size loci (IGF1R, SMAD2) + the chondro-
    dysplasia leg-length retrogene (FGF4, a structural insertion absent from a SNV panel) are OUT.
  - Input is per-locus BIG-ALLELE DOSAGE (0/1/2), the natural PLINK/VCF panel shape — NOT breeder allele
    symbols (that is coat's input shape). A genome/VCF caller for these 4 canFam4 coords is a v0.1 follow-on.
  - Faithful-to-literature: applies published loci + measured directions; it is not a new GWAS.

Pure-python, wheel-only, offline, deterministic. Regime-A curated catalog, NOT a learned embedding.
Scope: benign visible-trait genetics of a companion animal — NOT any human/forensic application.
The frozen decoder surface (amr_rules / calibrated_amr_rules / forward) is untouched (imports nothing from it).
"""
from __future__ import annotations

from dataclasses import dataclass, field

RULES_VERSION = "dog-body-size-v0.1.0"


@dataclass(frozen=True)
class SizeLocus:
    gene: str
    canfam4_variant: str        # panel ID chrN:pos:ref:alt on UU_Cfam_GSD_1.0 (canFam4)
    canfam3_coord: str          # literature canFam3.1 coord that was lifted
    big_allele: str             # the allele whose dose INCREASES size (larger dog)
    small_allele: str           # the reduced-size allele
    functional_r: float         # measured r(big-allele dose, Q121 height z) on Darwin's Ark N=3276
    source: str
    note: str = ""


# ---- curated + PINNED + functionally-validated catalog (Darwin's Ark canFam4, 2026-07-30) ----------------
# functional_r values are the MEASURED single-SNP correlations with owner-reported height (Q121, N=3276).
SIZE_LOCI: dict[str, SizeLocus] = {
    "IGF1": SizeLocus(
        "IGF1", "chr15:41513523:G:A", "chr15:41221438 (CanFam3.1; Hoopes/Rimbault 2013 fine-map)",
        big_allele="G", small_allele="A", functional_r=0.505,
        source="OMIA 002524-9615 (IGF1-AS body size); Sutter 2007 Science; Rimbault 2013 Genome Res",
        note="intron-2 SNP in complete LD with the SINEC_Cf small-size insertion (the SINE itself is a "
             "structural variant absent from a SNV panel; the SNP tags it). Largest single size locus."),
    "HMGA2": SizeLocus(
        "HMGA2", "chr10:8703415:G:A", "chr10:8348804 (CanFam3.1; Rimbault 2013)",
        big_allele="G", small_allele="A", functional_r=0.542,
        source="OMIA 001968-9615 (Height, HMGA2-associated body-size variation); Rimbault 2013 Genome Res",
        note="5'UTR SNP; ancestral G -> larger, derived A -> smaller. Second-largest effect; strongest "
             "measured single-SNP height correlation in this cohort."),
    "STC2": SizeLocus(
        "STC2", "chr4:40070215:T:A", "STC2 gene chr4:39,151,951-39,165,514 (CanFam3.1; NCBI gene 489112)",
        big_allele="T", small_allele="A", functional_r=0.369,
        source="Rimbault 2013 Genome Res; Plassais 2019 Nat Commun (STC2 CFA4 body-size locus)",
        note="best height-correlated panel SNP ~17 kb downstream of the STC2 gene body — matches the "
             "literature 'SNP ~20 kb downstream of STC2'."),
    "GHR": SizeLocus(
        "GHR", "chr4:67710295:C:T", "GHR gene chr4:66,705,544-66,845,096 (CanFam3.1; NCBI gene 403721)",
        big_allele="C", small_allele="T", functional_r=0.299,
        source="Rimbault 2013 Genome Res (GHR non-synonymous exon-5 SNPs, extracellular domain)",
        note="best height-correlated panel SNP within the GHR gene body; growth-hormone-receptor locus."),
}

# Measured combined 4-locus polygenic-score correlation with Q121 height on Darwin's Ark (N=3276):
POLYGENIC_SCORE_R = 0.619          # r; R^2 = 0.383 (~38% of cross-breed height variance from 4 SNPs)
VALIDATION_N = 3276
VALIDATION_COHORT = "Darwin's Ark (Dryad doi:10.5061/dryad.83bk3jb4r), canFam4 gp-0.70 biallelic imputed"

# ---- validated SINGLE-SNP morphology loci (beyond the polygenic size score) ----------------------------
# These are separate visible morphology traits each driven by ONE known SNP (not the additive size score).
# EAR is the only one validated on Darwin's Ark so far; pinned + functionally validated exactly like the
# size loci (canFam4 SNP -> dosage tracks the owner-reported morphology ordinal). CLEANLY resolved from
# body size: the ear lead (chr10:8,612,500) sits 91 kb from the HMGA2 size SNP (chr10:8,703,415) and the
# ear phenotype's correlation with the SIZE SNP is near-zero/opposite — the exact MSRB3-vs-HMGA2 confound
# Morrill 2022 had to untangle in the diverse cohort. FGF5 (coat length) + KRT71 (curl) were scanned but
# showed NO strong signal on any of the 9 Darwin's Ark morphology questions (those coat-texture traits are
# not among the measured morph Qs); the 4 covariate-adjusted "rerun" morph traits likewise did not map to a
# classic single SNP (max |r|~0.17) — likely SV-caused or different traits (see the validation artifact).


@dataclass(frozen=True)
class MorphLocus:
    trait: str                  # visible morphology trait
    gene: str
    canfam4_variant: str        # panel ID chrN:pos:ref:alt (canFam4)
    high_allele: str            # allele whose dose moves the ordinal UP (e.g. toward erect ears)
    functional_r: float         # measured r(high-allele dose, owner-reported ordinal) on Darwin's Ark
    darwins_ark_question: str   # the survey Q this validated against (identity inferred functionally)
    source: str
    note: str = ""


MORPH_LOCI: dict[str, MorphLocus] = {
    "EAR": MorphLocus(
        "ear_type_erect_vs_drop", "MSRB3", "chr10:8612500:A:G", high_allele="A",
        functional_r=0.543, darwins_ark_question="Q125",
        source="OMIA 000319-9615 (Ears, folded/drop vs prick); Boyko 2010; Vaysse 2011; ear lead "
               "chr10:8,612,500 (canFam4) — Sci Rep 2025 s41598-025-33036-0",
        note="intergenic 3' of MSRB3 / 5' of HMGA2; the exact published canFam4 ear lead variant is in-panel "
             "and its dosage tracks Q125 monotonically (dose 0 -1.13 -> dose 2 +0.43), r=+0.543. Resolved "
             "from body size: r with the HMGA2 size SNP is only -0.13 (opposite sign)."),
}

# secondary / structural size loci intentionally NOT in v0 (naming them keeps the abstention honest)
UNMODELLED_SIZE_LOCI = {
    "IGF1R": "CFA3 splice variant — secondary size locus, smaller increment (not in v0)",
    "SMAD2": "CFA7 ~9.9 kb downstream deletion — structural, secondary (not in v0)",
    "FGF4": "chondrodysplasia leg-length retrogene insertion (dachshund/corgi) — a STRUCTURAL variant "
            "absent from a SNV panel; changes leg length not overall skeletal size axis",
}

UNSEEN_MECHANISMS = (
    "a structural size variant not tagged by these 4 SNPs (FGF4 leg-length retrogene, SMAD2 downstream "
    "deletion) — this reads SNP dosages, so an SV with no LD-tag SNP is invisible",
    "breed-specific modifier loci outside the 4 catalogued genes — v0 gives a relative rank, not a breed call",
    "absolute height in inches — Q121 is a covariate-adjusted z-score, so v0 ranks size, it does not dose it",
    "a dog whose size is set by a locus at very low panel frequency (dropped in imputation)",
)


class SizeInputError(ValueError):
    """Malformed dosage input (never a silent wrong call)."""


@dataclass
class BodySizeCall:
    size_rank: str                  # relative tier: "toy/small" | "below-average" | "average" | "large/giant"
    polygenic_score: int            # sum of big-allele dosages across the loci scored (0..2*n_loci)
    max_score: int                  # 2 * number of loci actually scored
    n_loci_scored: int
    confidence: str                 # "high" | "medium" | "low"
    per_locus: dict                 # locus -> big-allele dose used, for audit
    abstains_on: list[str]
    rule: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "organism": "Canis_lupus_familiaris", "trait": "body_size",
            "regime": "A_curated_catalog_additive_polygenic", "rule": self.rule,
            "size_rank": self.size_rank, "polygenic_score": self.polygenic_score,
            "max_score": self.max_score, "n_loci_scored": self.n_loci_scored,
            "confidence": self.confidence, "per_locus": self.per_locus,
            "abstains_on": self.abstains_on, "notes": self.notes,
            "measure": "RELATIVE size rank (polygenic), NOT calibrated absolute height",
            "validation": {"cohort": VALIDATION_COHORT, "n": VALIDATION_N,
                           "polygenic_score_r": POLYGENIC_SCORE_R},
            "scope_limit": ("v0: 4 dominant loci (IGF1/HMGA2/STC2/GHR); relative rank not absolute inches; "
                            "reads big-allele SNP dosages not breeder allele symbols"),
            "undetectable_mechanisms": list(UNSEEN_MECHANISMS),
        }


def _size_rank(score: int, max_score: int) -> str:
    """Map a polygenic score to a RELATIVE size tier (validated gradient on Darwin's Ark height z).
    This is a rank within the cross-breed distribution, NOT a calibrated absolute height."""
    if max_score <= 0:
        return "undetermined"
    frac = score / max_score                       # 0 (all small) .. 1 (all big)
    if frac <= 0.375:
        return "toy/small"
    if frac <= 0.625:
        return "below-average"
    if frac <= 0.8125:
        return "average"
    return "large/giant"


def polygenic_size_score(dosages: dict[str, int], present_unmodelled: list[str] | None = None) -> BodySizeCall:
    """Deterministic RELATIVE body-size call from big-allele dosages at the catalogued loci.

    `dosages`: {locus -> big-allele copies 0/1/2} for any subset of IGF1/HMGA2/STC2/GHR. Missing loci are
    simply not scored (max_score shrinks); at least one is required. `present_unmodelled`: optional list of
    v0-unmodelled size loci declared present (e.g. ['FGF4']) -> the affected axis is added to abstains_on.
    """
    rule = f"dog_body_size_additive_v0 ({RULES_VERSION})"
    notes: list[str] = []
    per: dict = {}
    score = 0
    n = 0
    for loc, dose in dosages.items():
        L = loc.strip().upper()
        if L not in SIZE_LOCI:
            if L in UNMODELLED_SIZE_LOCI:
                raise SizeInputError(
                    f"locus {loc!r} is a v0-unmodelled size locus ({UNMODELLED_SIZE_LOCI[L]}); pass it via "
                    f"present_unmodelled=[...] so the axis ABSTAINS instead of skewing the additive score")
            raise SizeInputError(f"unknown size locus {loc!r}; v0 loci: {list(SIZE_LOCI)}")
        if dose not in (0, 1, 2):
            raise SizeInputError(f"{L} dosage {dose!r} must be 0, 1 or 2 (big-allele copies)")
        per[L] = dose
        score += dose
        n += 1
    if n == 0:
        raise SizeInputError("at least one catalogued locus dosage is required (IGF1/HMGA2/STC2/GHR)")

    abstains: list[str] = []
    for pl in (present_unmodelled or []):
        P = pl.strip().upper()
        if P in UNMODELLED_SIZE_LOCI:
            abstains.append(f"{P}: {UNMODELLED_SIZE_LOCI[P]}")
        elif P not in SIZE_LOCI:
            raise SizeInputError(f"unknown present locus {pl!r}")

    max_score = 2 * n
    rank = _size_rank(score, max_score)
    # confidence: high with both dominant loci (IGF1+HMGA2) present, else capped medium/low
    dominant_present = sum(1 for k in ("IGF1", "HMGA2") if k in per)
    if dominant_present == 2 and n >= 3:
        conf = "high"
    elif dominant_present >= 1:
        conf = "medium"
    else:
        conf = "low"
        notes.append("neither dominant locus (IGF1/HMGA2) scored -> low confidence relative rank")
    if abstains:
        conf = "medium" if conf == "high" else conf
    if n < len(SIZE_LOCI):
        notes.append(f"scored {n}/{len(SIZE_LOCI)} catalogued loci; max_score={max_score} (absent loci not penalised)")
    notes.append(f"additive polygenic score {score}/{max_score} -> RELATIVE tier '{rank}' "
                 f"(validated gradient r={POLYGENIC_SCORE_R} vs height on N={VALIDATION_N}; rank not inches)")
    return BodySizeCall(rank, score, max_score, n, conf, per, abstains, rule, notes)


def call_ear(dose: int) -> dict:
    """RELATIVE ear-morphology call from the MSRB3 ear-locus high-allele dosage (0/1/2).

    The SIGNAL is validated (dosage tracks the owner-reported ear ordinal Q125 at r=+0.543 on Darwin's Ark,
    N=2834, cleanly resolved from body size). The erect/drop NAMING follows the MSRB3 literature (Boyko 2010:
    two copies of the erect-associated allele -> erect/prick ears; zero -> drop/folded) — it is NOT
    independently label-confirmed here (the dataset has no codebook for the Q-number), so confidence is
    capped MEDIUM and the polarity caveat ships in the output. A future label-confirmed polarity is a
    one-line flip, not a re-derivation.
    """
    if dose not in (0, 1, 2):
        raise SizeInputError(f"EAR dosage {dose!r} must be 0, 1 or 2 (high-allele copies)")
    tier = {2: "erect (prick)", 1: "semi-erect/intermediate", 0: "drop (folded)"}[dose]
    return {
        "trait": "ear_type", "gene": "MSRB3", "locus": MORPH_LOCI["EAR"].canfam4_variant,
        "high_allele_dose": dose, "ear_type": tier, "confidence": "medium",
        "functional_r": MORPH_LOCI["EAR"].functional_r,
        "polarity_caveat": ("signal validated (r=+0.543 vs Q125); erect/drop NAMING is MSRB3-literature-anchored "
                            "(Boyko 2010), not independently label-confirmed on this codebook-less cohort"),
        "source": MORPH_LOCI["EAR"].source,
    }


def reference_integrity_ok() -> bool:
    """Catalog/rule contract guard (offline — no genotype data). Pins: all 4 loci present with well-formed
    canFam4 panel IDs whose big/small alleles are the ref/alt of that ID, positive measured directions, and
    a MONOTONIC additive rule (all-small -> toy/small, all-big -> large/giant, heterozygous -> middle)."""
    if set(SIZE_LOCI) != {"IGF1", "HMGA2", "STC2", "GHR"}:
        return False
    for L in SIZE_LOCI.values():
        parts = L.canfam4_variant.split(":")
        if len(parts) != 4 or not parts[0].startswith("chr"):
            return False
        ref, alt = parts[2], parts[3]
        if {L.big_allele, L.small_allele} != {ref, alt}:   # alleles must be the panel ref/alt
            return False
        if not (0.0 < L.functional_r < 1.0):
            return False
    # monotonic additive rule across the full 4-locus panel
    all_small = polygenic_size_score({k: 0 for k in SIZE_LOCI})
    all_het = polygenic_size_score({k: 1 for k in SIZE_LOCI})
    all_big = polygenic_size_score({k: 2 for k in SIZE_LOCI})
    # MORPH_LOCI (single-SNP morphology) well-formedness: canFam4 id whose high-allele is the ref/alt, r in range
    for M in MORPH_LOCI.values():
        parts = M.canfam4_variant.split(":")
        if len(parts) != 4 or not parts[0].startswith("chr") or M.high_allele not in (parts[2], parts[3]):
            return False
        if not (0.0 < M.functional_r < 1.0):
            return False
    return (all_small.polygenic_score == 0 and all_small.size_rank == "toy/small"
            and all_big.polygenic_score == 8 and all_big.size_rank == "large/giant"
            and all_het.polygenic_score == 4
            and all_small.polygenic_score < all_het.polygenic_score < all_big.polygenic_score
            and POLYGENIC_SCORE_R > max(L.functional_r for L in SIZE_LOCI.values())  # combined beats best single
            and "EAR" in MORPH_LOCI)
