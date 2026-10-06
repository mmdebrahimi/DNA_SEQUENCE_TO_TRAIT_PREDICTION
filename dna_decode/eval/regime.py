"""Which genotype->phenotype regime is a proposal in, and does that regime have a measured positive?

WHY THIS IS CODE AND NOT PROSE. The regime boundary has been mis-stated three separate times, always
the same way: compressing a SCOPED negative into a general one ("organism-level g->p is a closed
negative"), which hid a live direction each time. Prose went stale; a function with a test does not.
The compression is now impossible to make accidentally, because the scope is a parameter.

THE DISCRIMINATING VARIABLE IS POPULATION DESIGN, NOT ORGANISM COMPLEXITY. Constructed variation
randomises ancestry by construction; a natural population cannot. The yeast segregant cross decoded
12/12 quantitative traits at r 0.46-0.80 -- a clean organism-level positive -- while zero-shot
embeddings on natural populations are 0-for-5 de-confounded. "Eukaryotes are too complex" is the
wrong reading of the same data.

THE NEGATIVE IS ZERO-SHOT-SCOPED, and that scope is load-bearing. The shipped architecture is a
deterministic catalogue plus a SUPERVISED complement; refusing a supervised natural-population
proposal as "closed" would re-commit the exact over-generalisation this module exists to prevent. A
supervised proposal gets REQUIRES_DECONFOUNDING -- a condition to meet, not a wall.

Pure + dependency-free. Every cited number is traceable to a committed artifact; `scripts/regime_map.py`
REFUSES to emit a regime whose artifact is missing.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# --- the axes that were MEASURED to matter (organism complexity is deliberately not one) ----------
POPULATIONS = ("constructed", "natural")
ENDPOINTS = ("molecular", "organism", "organism_condition_switch")
METHODS = ("zero_shot", "supervised", "deterministic_catalog")

# --- WHICH SPLIT produced each number: an UNORDERED fact, plus ONE narrow ordered field -----------
#
# WHY TWO FIELDS AND NOT ONE LADDER. The first draft put every held-out unit on a single
# weakest-to-strongest scale (position < study < protein < organism < clade) so that "stronger than"
# was computable. That is wrong: a cross-PROTEIN molecular split and a leave-one-STUDY viral split are
# not weaker or stronger than each other, they answer DIFFERENT questions. Ranking them would let this
# module announce that one piece of evidence beats another when no such relation was measured --
# manufacturing exactly the confident-but-unmeasured claim it exists to prevent.
#
# So `split_unit` is an unordered SET and nothing compares its members. Only `organism_transfer` is
# ordered, and it answers one question: does this row carry held-out-ORGANISM evidence?
SPLIT_UNITS = frozenset({"position", "study", "protein", "organism", "clade", "condition", "none"})

ORGANISM_TRANSFER_UNMEASURED = "unmeasured"
ORGANISM_TRANSFER_LEVELS = (ORGANISM_TRANSFER_UNMEASURED, "held_out_organism", "held_out_clade")

CLOSED_NEGATIVE = "CLOSED_NEGATIVE"
WORKS = "WORKS"
OPEN = "OPEN"
LOSES_TO_CATALOG = "LOSES_TO_CATALOG"
REQUIRES_DECONFOUNDING = "REQUIRES_DECONFOUNDING"


@dataclass(frozen=True)
class Regime:
    key: str
    population: str
    endpoint: str
    method: str
    verdict: str
    evidence: str
    artifact: str
    note: str = ""
    # KEYWORD-ONLY and REQUIRED. `note` already carries a default, so a plain required field appended
    # after it is a TypeError; kw_only sidesteps that AND forces every row to name its split rather
    # than inheriting one by position. An empty frozenset is a real statement -- "no split was
    # performed" -- which is why there is no default to fall back on.
    split_unit: frozenset[str] = field(kw_only=True)
    # DEFAULTED, deliberately: `unmeasured` is simultaneously the fail-closed value and the true value
    # for every row today, so a default costs nothing and breaks no constructor.
    organism_transfer: str = field(kw_only=True, default=ORGANISM_TRANSFER_UNMEASURED)

    def as_dict(self) -> dict:
        return {"key": self.key, "population": self.population, "endpoint": self.endpoint,
                "method": self.method, "verdict": self.verdict, "evidence": self.evidence,
                "artifact": self.artifact, "note": self.note,
                # sorted list, not the frozenset: this dict is json.dumps'd by scripts/regime_map.py
                "split_unit": sorted(self.split_unit),
                "organism_transfer": self.organism_transfer}


REGIMES: tuple[Regime, ...] = (
    Regime("natural_organism_zeroshot", "natural", "organism", "zero_shot", CLOSED_NEGATIVE,
           "0-for-5 de-confounded; Arabidopsis FT10 embedding within-group r2 negative (-0.13) vs "
           "structure-only spearman 0.48 -- the embedding learned POPULATION STRUCTURE",
           "wiki/organism_gp_regime_correction_2026-08-29.md",
           "Do NOT scale this on a bigger GPU. A negative de-confounded metric is a signal-vs-structure "
           "problem, not a window-budget one.",
           # `none`, not `clade`: the within-group r2 GROUPS by kinship/PC to de-confound, but the
           # held-out unit inside each group is an INDIVIDUAL. A de-confounding group is not a
           # held-out split, and recording it as one would be a category error.
           split_unit=frozenset({"none"})),
    Regime("natural_organism_supervised", "natural", "organism", "supervised", REQUIRES_DECONFOUNDING,
           "not closed -- the 0-for-5 is ZERO-SHOT-only; a supervised complement is the shipped "
           "architecture. Condition: report WITHIN-GROUP performance against each group's own null",
           "wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md",
           "Pooled accuracy is dominated by any grouping variable the genotype tracks (clone, ancestry, "
           "submitter). Pooled numbers here are uninformative, not encouraging.",
           # no number at all on this row -- it carries conditions, not a measurement
           split_unit=frozenset({"none"})),
    Regime("constructed_molecular", "constructed", "molecular", "supervised", WORKS,
           "TEM-1 genome-edit path, Spearman 0.761 vs measured ampicillin fitness; externally "
           "replicated on a SECOND beta-lactamase (CTX-M-14/cefotaxime, independent lab) at 0.352",
           "wiki/pear_forward_replication_2026-09-02.md",
           "The one working molecular regime -- but the MAGNITUDE is protein-specific. Measured range "
           "0.35-0.76 across two beta-lactamases; do not quote 0.761 as the path's general strength. "
           "Direction holds in both (ESM2 beats BLOSUM62: 0.352 vs 0.198 on CTX-M-14). Lift comes from "
           "ORTHOGONAL MODALITIES, not scale: ESM2+GEMME+ProSST beats ESM2 on 90.5% of proteins paired, "
           "while 650M > 3B > 15B. A DAMAGE predictor cannot score a GAIN-of-function axis -- CTX-M-14 "
           "on ceftazidime is 0.078 for exactly that reason. **THE CITED EVIDENCE IS ZERO-SHOT** (ESM2 masked-marginal scoring); nothing in it trains on measured phenotype, so do NOT read this row as 'the supervised path is achieved'. Measured 2026-09-29 on PEAR CTX-M-14 with a held-out-POSITION split: a SUPERVISED alphabet head on top of ESM2 gives NO reliable lift (unstable across strides 2/3/4/5, position-clustered bootstrap 95% CI [-0.022,+0.099] spans zero), and a POSITION-FREE learned alphabet LOSES to BLOSUM62 on all four splits (0.13-0.27 vs 0.23-0.31) even trained on that protein's own 1,017 variants. See wiki/pear_genotype_alphabet_2026-09-29.md.",
           # TEM-1 -> CTX-M-14, independent lab: a cross-PROTEIN split. NOT comparable to a
           # leave-one-study split -- that incomparability is why SPLIT_UNITS is unordered.
           split_unit=frozenset({"protein"})),
    Regime("constructed_molecular_zeroshot", "constructed", "molecular", "zero_shot", WORKS,
           "the SAME evidence as constructed_molecular, correctly attributed: ESM2-650M "
           "masked-marginal ZERO-SHOT scoring, TEM-1 0.761 / CTX-M-14 0.352, and 0.31-0.41 "
           "across four held-out-position splits on CTX-M-14",
           "wiki/pear_genotype_alphabet_2026-09-29.md",
           "Declared explicitly because the default verdict for this cell was OPEN with the "
           "condition 'measure a de-confounded baseline first' -- a condition the 0.761/0.352 "
           "numbers ALREADY satisfy. A pretrained protein LM scoring variants zero-shot IS the "
           "working method here, and it BEAT every learned-alphabet variant tested.",
           split_unit=frozenset({"position"})),
    Regime("constructed_organism_per_condition", "constructed", "organism", "supervised", WORKS,
           "FBA iML1515 conditional essentiality, MCC 0.70-0.74 across four media",
           "wiki/organism_gp_regime_correction_2026-08-29.md",
           "Mechanistic, not learned. The organism-level positive that the 'too complex' reading misses.",
           # one organism (iML1515), four media. Nothing organism-level was held out.
           split_unit=frozenset({"none"})),
    Regime("constructed_organism_condition_switch", "constructed", "organism_condition_switch",
           "supervised", OPEN,
           "within-gene AUROC 0.73/0.81/0.71 on three axes (all p<=0.001), but the model emits ONE "
           "identical ratio for 61-76% of genes -- silent, not wrong",
           "wiki/fba_within_gene_ranking_2026-08-29.md",
           "The one genuinely OPEN cell. The readout lever is closed (+1.8pp oracle ceiling on the "
           "best-measured axis); the measured bottleneck is condition coverage in the expression data.",
           # `condition` is a SPLIT_UNITS member precisely so this row needs no ladder position --
           # it was the open question the single-ladder draft could not place.
           split_unit=frozenset({"condition"})),
    Regime("curated_catalog_exists", "natural", "molecular", "zero_shot", LOSES_TO_CATALOG,
           "HIV NNRTI: curated catalog AUC 0.926-0.962 vs ESM2 0.454 -- BELOW CHANCE",
           "wiki/hiv_esm_vs_catalog_2026-07-09.md",
           "Antagonistic endpoints INVERT a plausibility scorer: resistance is reached via chemically "
           "CONSERVATIVE substitutions at averagely-conserved sites, so likelihood calls them benign.",
           split_unit=frozenset({"none"})),
    # Added 2026-09-30. The cell above was the ZERO-SHOT half of (natural, molecular); the SUPERVISED half
    # was never measured, and the map's own warning -- that a zero-shot negative says nothing about
    # supervised -- applied to the map itself. Measured on the IDENTICAL isolate set with the IDENTICAL
    # pre-registered bar the zero-shot arm failed.
    Regime("natural_molecular_supervised_blindspot", "natural", "molecular", "supervised", WORKS,
           "HIV catalog blind spot, now on THREE genes: supervised genotype-token model 0.8102 (RT/EFV) "
           "and 0.8923 (integrase/RAL) leave-one-STUDY-out, both PASSING the bar zero-shot ESM2 failed "
           "at 0.4485; generalizes across 7 of 8 RT drugs. Protease is UNSCOREABLE -- its own catalog "
           "calls 74.2% positive, leaving 3 R of 910 in its blind spot",
           "wiki/glm_alphabet_headroom_2026-09-30.md",
           "SCOPED TO THE BLIND SPOT, not to replacing the catalog: the subset is catalog-NEGATIVE by "
           "construction, so the model cannot be rediscovering the catalog, and it still does NOT beat "
           "the catalog's 0.926 full-cohort AUC. The LINEAR one-hot model was called a FLOOR on the "
           "supervised family until 2026-09-30; it is now MEASURED to be the CEILING for this feature "
           "space -- pairwise interactions and a non-linear learner both LOSE on both genes (4 of 4 "
           "negative, 2 with a CI entirely below zero; wiki/hiv_context_vs_linear_floor_2026-09-30.md), "
           "corroborating the independent 2026-07-11 epistasis negative. That is NOT 'resistance is "
           "additive' (capacity may be unaffordable at this sample size) and does NOT test attention "
           "over RAW SEQUENCE, which remains unmeasured. All three genes are ONE virus. ALSO CLOSED 2026-10-01: the ORTHOGONAL-MODALITY lever has no room either -- one-hot's blindness to zero-training-signal columns is real (coefficients measured at exactly 0.0) but ABSORBED by per-isolate redundancy (median 1 blind column of 12), so the model scores slightly BETTER where it is blind (0.8229 vs 0.8163); see wiki/hiv_onehot_unseen_headroom_2026-10-01.md. A modality adds where the incumbent representation is STARVED of signal, not where it is merely blind in principle -- which is why the forward cell's single-variant modality lift does not transfer here.",
           # leave-one-STUDY-out across three genes of ONE virus. The one-virus limit this row's own
           # note already states in prose is now machine-readable: organism_transfer stays `unmeasured`.
           split_unit=frozenset({"study"})),
    # THE FIRST ROW IN THIS TABLE TO CARRY HELD-OUT-ORGANISM EVIDENCE (added 2026-10-06). Verdict is set
    # FROM the measurement, not chosen: the cross-organism number EXISTS and clears its own null by ~4x,
    # so this is not CLOSED_NEGATIVE; but the ladder built to decide decay-with-distance vs
    # eukaryote-cliff returned INDETERMINATE_INSUFFICIENT_RUNGS (2 scored of 5), so it is not WORKS
    # either. OPEN with a named condition is the honest cell.
    Regime("curated_catalog_cross_organism", "constructed", "organism", "curated_catalog", OPEN,
           "the E. coli-tuned conserved-core essentiality decoder applied UNCHANGED to human: "
           "coverage_lift 0.1633 vs its own label-permutation null MAX 0.0446 (~4x), against 0.3125 "
           "on the tuning organism -- so function-catalogue transfer across the bacteria/eukaryote "
           "boundary is REAL but roughly halves. The 5-rung ladder meant to separate smooth decay from "
           "a domain cliff returned INDETERMINATE_INSUFFICIENT_RUNGS: only 2 rungs scored, the other 3 "
           "are WALL_SCHEMA_UNVERIFIED (no label file on this host)",
           "wiki/essentiality_transfer_ladder_2026-10-06.md",
           "READ coverage_lift, NOT AUROC. The decay is SILENCE, not error: 83.1% of human essentials "
           "score EXACTLY zero and 11 of 13 core patterns have P(hit|non-essential)=0.0000, so AUROC is "
           "largely a tie-mass statistic and its two components move in OPPOSITE directions across the "
           "two rungs (coverage 0.3789->0.1689 while precision-where-it-fires 0.3684->0.9583). Cross-rung "
           "AUROC is ALSO a different sampling frame (BAGEL two curated extremes at base rate 0.431 vs "
           "genome-wide 0.0928) and is refused. ~10% of the human miss is a catalogue PHRASING gap, not "
           "phylogeny (floor 57/566 vs host-specific 141/566; wiki/essentiality_missed_vocabulary_"
           "2026-10-05.md), so the phrasing-adjusted lift ships beside the raw one. CONDITION to close "
           "this cell: score a third rung -- and note S. cerevisiae and H. sapiens are EQUIDISTANT from "
           "E. coli by shared lineage depth, so they are one distance class, not two rungs.",
           # held out the ORGANISM: tuned on E. coli, applied unchanged to human with nothing re-fit.
           split_unit=frozenset({"organism"}),
           organism_transfer="held_out_organism"),
)

# Where a learned layer is pointed. The catalog-beats-learning result is about REPLACING a catalog; it was
# being applied to every learned proposal, which over-refused the one shape that is now measured to work.
TARGETS = ("replace", "blind_spot_complement")

_BY_KEY = {r.key: r for r in REGIMES}


def split_units_for(key: str) -> frozenset[str]:
    """Which unit(s) this regime's number held out. UNORDERED -- nothing compares the members."""
    if key not in _BY_KEY:
        raise KeyError(f"unknown regime {key!r}; known: {sorted(_BY_KEY)}")
    return _BY_KEY[key].split_unit


def organism_transfer_is_unmeasured() -> tuple[str, ...]:
    """The regimes carrying NO held-out-organism evidence.

    Today this is every row of THIS TABLE, which is the headline F1 exists to make machine-readable.
    Callers should treat a SHRINKING return value as the event worth noticing, not the current length.

    SCOPE CORRECTION 2026-10-05 -- the earlier wording here said "the project has never measured a
    cross-organism number", and that is FALSE as a project-level claim. One genuine cross-kingdom
    transfer measurement exists and predates this table: `wiki/essentiality_e4_transfer_2026-07-28.json`
    applies the E. coli-tuned conserved-core essentiality decoder UNCHANGED to human (BAGEL CEGv2 vs
    NEGv1, 681 essential / 899 non) and reports AUROC 0.5805 against a 0.5 null, with spec 0.998 and
    sens 0.157 -- high-precision, very-low-recall. Its mechanism is measured, not guessed: the UNIVERSAL
    core (ribosome / tRNA-synthetase / translation / polymerase) transfers cross-kingdom, while the
    human-specific core is missed entirely (proteasome 0/53, spliceosome 0/49).

    So the honest claim is "no regime IN THIS TABLE carries held-out-organism evidence", NOT "the project
    has never measured one". The essentiality cell's transfer number is simply not represented as a
    regime row. DO NOT read a 0.911 human AUROC anywhere as transfer -- that figure
    (`essentiality_e3_human_2026-07-28.json`) is a WITHIN-human 5-fold CV on human labels, trained and
    tested on human, and misreading it as cross-organism transfer is a live trap.
    """
    return tuple(r.key for r in REGIMES
                 if r.organism_transfer == ORGANISM_TRANSFER_UNMEASURED)


def _transfer_rank(level: str) -> int:
    """Position within ORGANISM_TRANSFER_LEVELS. Defined ONLY for that narrow field -- there is
    deliberately no equivalent for SPLIT_UNITS, because its members are incomparable."""
    return ORGANISM_TRANSFER_LEVELS.index(level)


def _transfer_gap_conditions(regime: "Regime", claim: str | None) -> list[str]:
    """AUGMENT-ONLY: a condition when the caller claims more organism-transfer evidence than the
    matched regime carries. Returns [] when there is no claim or the claim is already supported."""
    if claim is None or _transfer_rank(claim) <= _transfer_rank(regime.organism_transfer):
        return []
    # DERIVED, never restated. This used to end with the literal "no regime in this map carries one
    # today", which became FALSE the moment `curated_catalog_cross_organism` landed (2026-10-06). A
    # restated fact goes stale silently -- the same failure mode as the retired-lock filename this repo
    # already records -- so the sentence is computed from the live table instead.
    carriers = sorted(r.key for r in REGIMES
                      if _transfer_rank(r.organism_transfer) >= _transfer_rank("held_out_organism"))
    if carriers:
        where = ("regime(s) that DO carry held-out-organism evidence: %s -- cite one of those, or "
                 "measure a held-out-ORGANISM number for this one" % ", ".join(carriers))
    else:
        where = ("measure a held-out-ORGANISM number before claiming organism transfer; no regime in "
                 "this map carries one today")
    return [f"this proposal claims {claim!r} but regime {regime.key!r} carries "
            f"organism_transfer={regime.organism_transfer!r} (split unit(s): "
            f"{sorted(regime.split_unit)}) -- the claim is NOT supported by the cited evidence",
            where]


@dataclass
class ScreenResult:
    regime: str | None
    verdict: str
    reason: str
    evidence: str = ""
    artifact: str = ""
    conditions: list = field(default_factory=list)

    @property
    def refused(self) -> bool:
        return self.verdict == CLOSED_NEGATIVE

    def as_dict(self) -> dict:
        return {"regime": self.regime, "verdict": self.verdict, "reason": self.reason,
                "evidence": self.evidence, "artifact": self.artifact,
                "conditions": list(self.conditions), "refused": self.refused}


def classify_regime(population: str, endpoint: str, method: str) -> Regime | None:
    """Exact (population, endpoint, method) match, or None. Pure."""
    p, e, m = (str(x).strip().lower() for x in (population, endpoint, method))
    for r in REGIMES:
        if (r.population, r.endpoint, r.method) == (p, e, m):
            return r
    return None


def screen_proposal(population: str, endpoint: str, method: str,
                    curated_catalog_exists: bool = False,
                    target: str = "replace",
                    claims_organism_transfer: str | None = None) -> ScreenResult:
    """Screen a learned-decoder proposal against the measured regime map.

    A CLOSED_NEGATIVE verdict is a refusal: that exact regime has been tested and failed under
    de-confounding, and more scale does not address it. Every other verdict names conditions rather
    than blocking -- the boundary exists to stop the ONE repeated mistake, not to forbid learning.

    `claims_organism_transfer` is AUGMENT-ONLY. When a caller claims more organism-transfer evidence
    than the matched regime actually carries, a CONDITION is appended naming the gap. It never changes
    `verdict`, `reason`, `evidence` or `artifact` on any path -- the same discipline every other
    disclosure layer in this repo follows, and it is proved by a with/without diff rather than
    promised.
    """
    p, e, m = (str(x).strip().lower() for x in (population, endpoint, method))
    if p not in POPULATIONS:
        return ScreenResult(None, "UNKNOWN", f"population must be one of {POPULATIONS}; got {p!r}")
    if e not in ENDPOINTS:
        return ScreenResult(None, "UNKNOWN", f"endpoint must be one of {ENDPOINTS}; got {e!r}")
    if m not in METHODS:
        return ScreenResult(None, "UNKNOWN", f"method must be one of {METHODS}; got {m!r}")

    if target not in TARGETS:
        return ScreenResult(None, "UNKNOWN", f"target must be one of {TARGETS}; got {target!r}")

    # Validated BEFORE the catalog short-circuit so a bad value cannot slip through on either path.
    if claims_organism_transfer is not None and claims_organism_transfer not in ORGANISM_TRANSFER_LEVELS:
        return ScreenResult(None, "UNKNOWN",
                            f"claims_organism_transfer must be one of {ORGANISM_TRANSFER_LEVELS}; "
                            f"got {claims_organism_transfer!r}")

    # A curated catalog beats a learned scorer wherever one exists -- checked BEFORE the regime match,
    # because it is the strongest measured result and it inverts (ESM 0.454, below chance).
    #
    # SCOPED 2026-09-30: this applies to REPLACING the catalog. Applying it to a blind-spot COMPLEMENT
    # over-refused the one shape now measured to work (supervised tokens, 0.8142 leave-study-out on the
    # catalog-negative subset) -- the same over-compression this module exists to prevent, committed by
    # this module. A complement proposal falls through to the regime match instead of being refused here.
    if curated_catalog_exists and m != "deterministic_catalog" and target == "replace":
        r = _BY_KEY["curated_catalog_exists"]
        return ScreenResult(r.key, LOSES_TO_CATALOG,
                            "a curated catalog exists for this endpoint, and a learned scorer has been "
                            "measured to LOSE to it here",
                            r.evidence, r.artifact,
                            ["beat the curated catalog on held-out data before proposing to replace it",
                             "if the endpoint is antagonistic (drug resistance), expect BELOW-chance"]
                            + _transfer_gap_conditions(r, claims_organism_transfer))

    r = classify_regime(p, e, m)
    if r is None:
        return ScreenResult(None, OPEN,
                            f"no measured result for ({p}, {e}, {m}) -- unscreened, which is not the "
                            "same as promising", conditions=["measure a de-confounded baseline first"])
    conds = []
    if r.verdict == REQUIRES_DECONFOUNDING:
        conds = ["report WITHIN-GROUP performance, not pooled",
                 "compare each group against its OWN null",
                 "mark single-class groups unscorable rather than scoring them"]
    elif r.verdict == CLOSED_NEGATIVE:
        conds = ["do NOT re-run this at larger scale — the failure is signal-vs-structure",
                 "changing population DESIGN (constructed variation) moves it to a different regime"]
    conds = conds + _transfer_gap_conditions(r, claims_organism_transfer)
    return ScreenResult(r.key, r.verdict, r.note or r.evidence, r.evidence, r.artifact, conds)
