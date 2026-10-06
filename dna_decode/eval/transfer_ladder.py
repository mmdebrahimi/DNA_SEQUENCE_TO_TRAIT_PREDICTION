"""The cross-organism transfer-ladder design, frozen — coverage_lift primary, thresholds parameterised.

THE QUESTION. The project's north star ends at *climb to higher life forms*, and the entire empirical
basis for that step is ONE measurement: the E. coli-tuned conserved-core essentiality decoder
(`dna_decode/essentiality/core_decoder.py`) applied UNCHANGED to human. One point at the maximum
available phylogenetic distance cannot separate

    (A) continuous decay with phylogenetic distance  from  (B) a cliff at the eukaryote boundary

and the two imply different work. A ladder of rungs at increasing distance can.

WHY coverage_lift AND NOT AUROC -- measured, not stylistic. The decay is not a discrimination failure,
it is SILENCE: **83.1% of human essential genes score EXACTLY zero** (E. coli 62.1%), and in human
**11 of 13 core patterns have P(hit | non-essential) = 0.0000** -- the patterns never misfire, they
never fire. AUROC over a distribution that is 83% tied is mostly a tie-mass statistic. Worse, its two
components move in OPPOSITE directions across the two rungs we already have, so AUROC hides both:

    rung      coverage(ess)   coverage(non)   coverage_lift   precision where it fires
    E. coli      0.3789          0.0664          0.3125              0.3684
    human        0.1689          0.0056          0.1633              0.9583

`coverage_lift` = coverage(essential) - coverage(non-essential). BOTH terms are conditioned on the
class, so the statistic is base-rate robust and therefore genuinely cross-rung comparable -- which
AUROC is NOT here (spectrum bias: BAGEL CEGv2/NEGv1 are two curated extremes at base rate 0.431 while
E. coli is genome-wide at 0.0928) and sens / spec / precision are not either. It also has an obvious
permutation null: shuffle the labels and the expected lift is 0.

THE EUKARYOTE TIE, found by doing the derivation rather than assuming it. Ordering is DERIVED from
NCBI Taxonomy lineages (fetched + name-verified, never asserted from memory), and shared-rank depth
from E. coli comes out: P. aeruginosa 5, S. aureus 2, **S. cerevisiae 1, H. sapiens 1**. The two
eukaryotes are EQUIDISTANT from E. coli by this measure, so the ladder has FOUR distance classes, not
five. That is not a defect to paper over -- it makes the eukaryote pair a natural CONSISTENCY CHECK:
both (A) and (B) predict they should behave SIMILARLY (equidistant under A, both past the cliff under
B), so a large yeast-vs-human gap would mean neither account is doing the work and something
organism-specific dominates. `classify_ladder` therefore reports the pair's spread rather than
pretending to an ordering the lineages do not support.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# ---------------------------------------------------------------------------------------------------
# THE BAR, written before any new rung is fetched. Mirrors scripts/transfer_benchmark_falsifier.py's
# PREREGISTERED block. Thresholds are deliberately NOT here -- see classify_ladder.
# ---------------------------------------------------------------------------------------------------
SMOOTH_DECAY_WITH_DISTANCE = "SMOOTH_DECAY_WITH_DISTANCE"
DOMAIN_CLIFF = "DOMAIN_CLIFF"
NO_DECAY_DISTANCE_INSENSITIVE = "NO_DECAY_DISTANCE_INSENSITIVE"
INDETERMINATE_INSUFFICIENT_RUNGS = "INDETERMINATE_INSUFFICIENT_RUNGS"
INDETERMINATE_NULL_NOT_CLEARED = "INDETERMINATE_NULL_NOT_CLEARED"
INDETERMINATE_BAR_SENSITIVE = "INDETERMINATE_BAR_SENSITIVE"

VERDICTS = (SMOOTH_DECAY_WITH_DISTANCE, DOMAIN_CLIFF, NO_DECAY_DISTANCE_INSENSITIVE,
            INDETERMINATE_INSUFFICIENT_RUNGS, INDETERMINATE_NULL_NOT_CLEARED,
            INDETERMINATE_BAR_SENSITIVE)

PRIMARY_TEST = "bacterial_technology_matched_subladder"

PREREGISTERED = {
    "question": "is the conserved-core decoder's cross-organism decay (A) continuous with phylogenetic "
                "distance or (B) a cliff at the eukaryote boundary?",
    "primary_metric": "coverage_lift = coverage(essential) - coverage(non-essential); class-conditioned "
                      "so base-rate robust and cross-rung comparable",
    "secondary_metric": "AUROC -- reported but NOT compared across rungs (different sampling frames; "
                        "83% tie mass in human makes it a tie-mass statistic)",
    "primary_test": PRIMARY_TEST,
    "primary_test_reason": "the three bacterial rungs are all transposon-insertion screens, so distance "
                           "varies while label TECHNOLOGY is held constant",
    "secondary_test": "full 5-rung ladder -- technology-confounded AND sampling-frame-confounded",
    "verdicts": list(VERDICTS),
    "falsified_if": "coverage_lift does not clear its own label-permutation null on the rungs that score",
    "eukaryote_tie": "S. cerevisiae and H. sapiens share equal lineage depth with E. coli (1), so they "
                     "are one distance class; their SPREAD is a consistency check, not a rung ordering",
    "decoder_must_be_unchanged": "a per-organism synonym layer would destroy the transfer condition; "
                                 "sha256(core_decoder.py) is pinned by test",
    "frozen_at": "plan v2 Step 1, before any new rung was fetched",
}

# The one authoritative source for every lineage below. Verified 2026-10-05: each taxid was fetched and
# its <ScientificName> checked against the organism it is claimed to be, so these strings are SOURCED,
# not remembered. Mirrors the provenance gate in dna_decode/data/ecoff_catalog.py.
TAXONOMY_SOURCE = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                   "?db=taxonomy&id={taxid}&retmode=xml")

LINEAGE_SEP = "; "


class RungProvenanceError(ValueError):
    """A rung without a sourced lineage. Raised at construction, never defaulted."""


@dataclass(frozen=True)
class Rung:
    key: str
    taxid: int
    scientific_name: str
    label_source_id: str
    technology: str          # transposon_insertion | deletion_collection | crispr_ko
    lineage: str             # verbatim NCBI Taxonomy <Lineage>
    source_url: str
    annotation_kind: str
    key_space: str

    def __post_init__(self):
        if not self.lineage or not self.lineage.strip():
            raise RungProvenanceError(f"{self.key}: lineage is required and must be sourced")
        if not self.source_url or not self.source_url.startswith("http"):
            raise RungProvenanceError(f"{self.key}: source_url is required (got {self.source_url!r})")
        if self.technology not in ("transposon_insertion", "deletion_collection", "crispr_ko"):
            raise RungProvenanceError(f"{self.key}: unknown label technology {self.technology!r}")

    def ranks(self) -> tuple[str, ...]:
        return tuple(r.strip() for r in self.lineage.split(LINEAGE_SEP) if r.strip())


def _u(taxid: int) -> str:
    return TAXONOMY_SOURCE.format(taxid=taxid)


RUNGS: tuple[Rung, ...] = (
    Rung("ecoli", 562, "Escherichia coli", "goodall_tradis", "transposon_insertion",
         "cellular organisms; Bacteria; Pseudomonadati; Pseudomonadota; Gammaproteobacteria; "
         "Enterobacterales; Enterobacteriaceae; Escherichia",
         _u(562), "ncbi_feature_table", "gene_symbol"),
    Rung("paeruginosa", 287, "Pseudomonas aeruginosa", "plos_gold_pa14", "transposon_insertion",
         "cellular organisms; Bacteria; Pseudomonadati; Pseudomonadota; Gammaproteobacteria; "
         "Pseudomonadales; Pseudomonadaceae; Pseudomonas; Pseudomonas aeruginosa group",
         _u(287), "ncbi_feature_table", "pa14_locus_tag"),
    Rung("saureus", 1280, "Staphylococcus aureus", "ntml_nebraska", "transposon_insertion",
         "cellular organisms; Bacteria; Bacillati; Bacillota; Bacilli; Caryophanales; "
         "Staphylococcaceae; Staphylococcus",
         _u(1280), "ncbi_feature_table", "sausa300_locus_tag"),
    Rung("scerevisiae", 4932, "Saccharomyces cerevisiae", "sgd_phenotype", "deletion_collection",
         "cellular organisms; Eukaryota; Opisthokonta; Fungi; Dikarya; Ascomycota; saccharomyceta; "
         "Saccharomycotina; Saccharomycetes; Saccharomycetales; Saccharomycetaceae; Saccharomyces",
         _u(4932), "sgd_features", "systematic_orf"),
    # label_source_id names the ESSENTIAL-class source; the human rung's negative class comes from a
    # SECOND file (NEGv1). The full per-rung plan lives in
    # dna_decode/essentiality/label_sources.py::RUNG_LABEL_PLAN, because the rungs genuinely differ
    # (one-file-two-columns / two-files / essential-plus-complement) and one id cannot express that.
    Rung("human", 9606, "Homo sapiens", "bagel_ceg", "crispr_ko",
         "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; "
         "Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; "
         "Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; "
         "Euarchontoglires; Primates; Haplorrhini; Simiiformes; Catarrhini; Hominoidea; Hominidae; "
         "Homininae; Homo",
         _u(9606), "ncbi_gene_info", "entrez_gene_id"),
)

RUNGS_BY_KEY = {r.key: r for r in RUNGS}
REFERENCE_RUNG = "ecoli"          # the organism the decoder was tuned on
BACTERIAL_RUNGS = ("ecoli", "paeruginosa", "saureus")


# ---------------------------------------------------------------------------------------------------
# Pure metric. The whole reason this module exists.
# ---------------------------------------------------------------------------------------------------
def coverage(scores) -> float:
    """Fraction of genes the decoder scores NON-ZERO. 'Did it fire at all', not 'how high'."""
    scores = list(scores)
    if not scores:
        raise ValueError("coverage() on an empty score list -- an empty class is not 0.0 coverage, "
                         "it is an absent measurement")
    return sum(1 for s in scores if s != 0) / len(scores)


def coverage_lift(scores_essential, scores_nonessential) -> float:
    """coverage(essential) - coverage(non-essential).

    BASE-RATE ROBUST BY CONSTRUCTION, which is the entire point: each term is conditioned on its own
    class, so duplicating or thinning either class leaves the statistic unchanged. That is what makes
    it comparable across rungs whose prevalences differ by ~5x (E. coli 0.0928 genome-wide vs BAGEL's
    curated-extremes 0.431), where AUROC is not comparable and sens/spec/precision are not either.
    """
    return coverage(scores_essential) - coverage(scores_nonessential)


def shared_rank_depth(lineage_a: str, lineage_b: str) -> int:
    """Number of shared LEADING taxonomic ranks. DERIVES the ladder order instead of asserting it.

    Not an intuition check: E. coli and P. aeruginosa share 5 ranks (both Gammaproteobacteria) while
    S. aureus shares only 2 (different phylum), so P. aeruginosa is the NEARER rung even though
    S. aureus is the more familiar organism. And S. cerevisiae and H. sapiens both share exactly 1
    with E. coli -- they are one distance class, not two rungs.
    """
    a = [r.strip() for r in lineage_a.split(LINEAGE_SEP) if r.strip()]
    b = [r.strip() for r in lineage_b.split(LINEAGE_SEP) if r.strip()]
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def depth_from_reference(rung: Rung, reference: Rung | None = None) -> int:
    """Shared-rank depth from the tuning organism. The reference's depth is its own rank count."""
    reference = reference or RUNGS_BY_KEY[REFERENCE_RUNG]
    if rung.key == reference.key:
        return len(reference.ranks())
    return shared_rank_depth(rung.lineage, reference.lineage)


def distance_classes() -> list[tuple[int, list[str]]]:
    """Rungs grouped by shared-rank depth, NEAREST first. Groups of >1 are one distance class."""
    by_depth: dict[int, list[str]] = {}
    for r in RUNGS:
        by_depth.setdefault(depth_from_reference(r), []).append(r.key)
    return [(d, sorted(by_depth[d])) for d in sorted(by_depth, reverse=True)]


# ---------------------------------------------------------------------------------------------------
# Verdict. Every threshold is a PARAMETER -- the runner freezes them AFTER the distribution analysis,
# per the ratified dna-identify pattern (this repo has two recorded cases where the pre-registered bar
# was itself the error). Nothing here defaults a threshold.
# ---------------------------------------------------------------------------------------------------
@dataclass
class LadderVerdict:
    verdict: str
    reason: str
    n_scored: int
    detail: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"verdict": self.verdict, "reason": self.reason, "n_scored": self.n_scored,
                "detail": self.detail}


def classify_ladder(rungs, *, cliff_drop: float, plateau_tol: float, min_rungs: int) -> LadderVerdict:
    """Pure verdict over scored rung records.

    Each record: {key, depth, coverage_lift, null_p95, coverage_lift_adjusted (optional), scored: bool}.
    Ordered checks, and POWERING IS CHECKED FIRST -- a run with too few rungs necessarily has no decay
    to find, so testing decay before powering would publish an unpowered run as a finding.
    """
    scored = [r for r in rungs if r.get("scored")]
    if len(scored) < min_rungs:
        return LadderVerdict(INDETERMINATE_INSUFFICIENT_RUNGS,
                             f"{len(scored)} scored rung(s) < min_rungs={min_rungs}", len(scored))

    not_cleared = [r["key"] for r in scored if r["coverage_lift"] <= r["null_p95"]]
    if not_cleared:
        return LadderVerdict(INDETERMINATE_NULL_NOT_CLEARED,
                             f"coverage_lift inside its own permutation p95 for: {not_cleared}",
                             len(scored), {"not_cleared": not_cleared})

    # Bar sensitivity: if subtracting the phrasing floor moves a rung into a different band relative to
    # the nearest rung, the instrument is not deciding -- say so rather than pick the flattering number.
    def band(v):
        return "high" if v >= ref_lift - plateau_tol else ("low" if v <= ref_lift - cliff_drop else "mid")

    by_depth = sorted(scored, key=lambda r: -r["depth"])
    ref_lift = by_depth[0]["coverage_lift"]
    disagree = []
    for r in by_depth:
        adj = r.get("coverage_lift_adjusted")
        if adj is not None and band(r["coverage_lift"]) != band(adj):
            disagree.append(r["key"])
    if disagree:
        return LadderVerdict(INDETERMINATE_BAR_SENSITIVE,
                             f"raw and phrasing-adjusted lifts fall in different bands for: {disagree}",
                             len(scored), {"bar_sensitive": disagree})

    lifts = [r["coverage_lift"] for r in by_depth]
    spread = max(lifts) - min(lifts)
    if spread <= plateau_tol:
        return LadderVerdict(NO_DECAY_DISTANCE_INSENSITIVE,
                             f"all scored rungs within plateau_tol={plateau_tol} (spread {spread:.4f})",
                             len(scored), {"spread": spread})

    # DOMAIN_CLIFF: the bacterial rungs sit together and the eukaryote class drops away from them.
    bact = [r for r in by_depth if r["key"] in BACTERIAL_RUNGS]
    euk = [r for r in by_depth if r["key"] not in BACTERIAL_RUNGS]
    if len(bact) >= 2 and euk:
        b_lifts = [r["coverage_lift"] for r in bact]
        if max(b_lifts) - min(b_lifts) <= plateau_tol and min(b_lifts) - max(
                r["coverage_lift"] for r in euk) >= cliff_drop:
            return LadderVerdict(DOMAIN_CLIFF,
                                 f"bacterial rungs within {plateau_tol} and eukaryote class "
                                 f">={cliff_drop} below them", len(scored),
                                 {"bacterial_spread": max(b_lifts) - min(b_lifts)})

    return LadderVerdict(SMOOTH_DECAY_WITH_DISTANCE,
                         f"spread {spread:.4f} exceeds plateau_tol with no bacterial/eukaryote cliff",
                         len(scored), {"spread": spread, "lifts_by_depth": lifts})
