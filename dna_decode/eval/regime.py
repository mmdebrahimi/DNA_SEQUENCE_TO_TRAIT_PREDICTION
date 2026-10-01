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

    def as_dict(self) -> dict:
        return {"key": self.key, "population": self.population, "endpoint": self.endpoint,
                "method": self.method, "verdict": self.verdict, "evidence": self.evidence,
                "artifact": self.artifact, "note": self.note}


REGIMES: tuple[Regime, ...] = (
    Regime("natural_organism_zeroshot", "natural", "organism", "zero_shot", CLOSED_NEGATIVE,
           "0-for-5 de-confounded; Arabidopsis FT10 embedding within-group r2 negative (-0.13) vs "
           "structure-only spearman 0.48 -- the embedding learned POPULATION STRUCTURE",
           "wiki/organism_gp_regime_correction_2026-08-29.md",
           "Do NOT scale this on a bigger GPU. A negative de-confounded metric is a signal-vs-structure "
           "problem, not a window-budget one."),
    Regime("natural_organism_supervised", "natural", "organism", "supervised", REQUIRES_DECONFOUNDING,
           "not closed -- the 0-for-5 is ZERO-SHOT-only; a supervised complement is the shipped "
           "architecture. Condition: report WITHIN-GROUP performance against each group's own null",
           "wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md",
           "Pooled accuracy is dominated by any grouping variable the genotype tracks (clone, ancestry, "
           "submitter). Pooled numbers here are uninformative, not encouraging."),
    Regime("constructed_molecular", "constructed", "molecular", "supervised", WORKS,
           "TEM-1 genome-edit path, Spearman 0.761 vs measured ampicillin fitness; externally "
           "replicated on a SECOND beta-lactamase (CTX-M-14/cefotaxime, independent lab) at 0.352",
           "wiki/pear_forward_replication_2026-09-02.md",
           "The one working molecular regime -- but the MAGNITUDE is protein-specific. Measured range "
           "0.35-0.76 across two beta-lactamases; do not quote 0.761 as the path's general strength. "
           "Direction holds in both (ESM2 beats BLOSUM62: 0.352 vs 0.198 on CTX-M-14). Lift comes from "
           "ORTHOGONAL MODALITIES, not scale: ESM2+GEMME+ProSST beats ESM2 on 90.5% of proteins paired, "
           "while 650M > 3B > 15B. A DAMAGE predictor cannot score a GAIN-of-function axis -- CTX-M-14 "
           "on ceftazidime is 0.078 for exactly that reason. **THE CITED EVIDENCE IS ZERO-SHOT** (ESM2 masked-marginal scoring); nothing in it trains on measured phenotype, so do NOT read this row as 'the supervised path is achieved'. Measured 2026-09-29 on PEAR CTX-M-14 with a held-out-POSITION split: a SUPERVISED alphabet head on top of ESM2 gives NO reliable lift (unstable across strides 2/3/4/5, position-clustered bootstrap 95% CI [-0.022,+0.099] spans zero), and a POSITION-FREE learned alphabet LOSES to BLOSUM62 on all four splits (0.13-0.27 vs 0.23-0.31) even trained on that protein's own 1,017 variants. See wiki/pear_genotype_alphabet_2026-09-29.md."),
    Regime("constructed_molecular_zeroshot", "constructed", "molecular", "zero_shot", WORKS,
           "the SAME evidence as constructed_molecular, correctly attributed: ESM2-650M "
           "masked-marginal ZERO-SHOT scoring, TEM-1 0.761 / CTX-M-14 0.352, and 0.31-0.41 "
           "across four held-out-position splits on CTX-M-14",
           "wiki/pear_genotype_alphabet_2026-09-29.md",
           "Declared explicitly because the default verdict for this cell was OPEN with the "
           "condition 'measure a de-confounded baseline first' -- a condition the 0.761/0.352 "
           "numbers ALREADY satisfy. A pretrained protein LM scoring variants zero-shot IS the "
           "working method here, and it BEAT every learned-alphabet variant tested."),
    Regime("constructed_organism_per_condition", "constructed", "organism", "supervised", WORKS,
           "FBA iML1515 conditional essentiality, MCC 0.70-0.74 across four media",
           "wiki/organism_gp_regime_correction_2026-08-29.md",
           "Mechanistic, not learned. The organism-level positive that the 'too complex' reading misses."),
    Regime("constructed_organism_condition_switch", "constructed", "organism_condition_switch",
           "supervised", OPEN,
           "within-gene AUROC 0.73/0.81/0.71 on three axes (all p<=0.001), but the model emits ONE "
           "identical ratio for 61-76% of genes -- silent, not wrong",
           "wiki/fba_within_gene_ranking_2026-08-29.md",
           "The one genuinely OPEN cell. The readout lever is closed (+1.8pp oracle ceiling on the "
           "best-measured axis); the measured bottleneck is condition coverage in the expression data."),
    Regime("curated_catalog_exists", "natural", "molecular", "zero_shot", LOSES_TO_CATALOG,
           "HIV NNRTI: curated catalog AUC 0.926-0.962 vs ESM2 0.454 -- BELOW CHANCE",
           "wiki/hiv_esm_vs_catalog_2026-07-09.md",
           "Antagonistic endpoints INVERT a plausibility scorer: resistance is reached via chemically "
           "CONSERVATIVE substitutions at averagely-conserved sites, so likelihood calls them benign."),
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
           "over RAW SEQUENCE, which remains unmeasured. All three genes are ONE virus. ALSO CLOSED 2026-10-01: the ORTHOGONAL-MODALITY lever has no room either -- one-hot's blindness to zero-training-signal columns is real (coefficients measured at exactly 0.0) but ABSORBED by per-isolate redundancy (median 1 blind column of 12), so the model scores slightly BETTER where it is blind (0.8229 vs 0.8163); see wiki/hiv_onehot_unseen_headroom_2026-10-01.md. A modality adds where the incumbent representation is STARVED of signal, not where it is merely blind in principle -- which is why the forward cell's single-variant modality lift does not transfer here."),
)

# Where a learned layer is pointed. The catalog-beats-learning result is about REPLACING a catalog; it was
# being applied to every learned proposal, which over-refused the one shape that is now measured to work.
TARGETS = ("replace", "blind_spot_complement")

_BY_KEY = {r.key: r for r in REGIMES}


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
                    target: str = "replace") -> ScreenResult:
    """Screen a learned-decoder proposal against the measured regime map.

    A CLOSED_NEGATIVE verdict is a refusal: that exact regime has been tested and failed under
    de-confounding, and more scale does not address it. Every other verdict names conditions rather
    than blocking -- the boundary exists to stop the ONE repeated mistake, not to forbid learning.
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
                             "if the endpoint is antagonistic (drug resistance), expect BELOW-chance"])

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
    return ScreenResult(r.key, r.verdict, r.note or r.evidence, r.evidence, r.artifact, conds)
