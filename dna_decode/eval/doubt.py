"""L2 DOUBT — "this call may be incomplete, and here is why". Never a competing call.

WHY. The catalog's measured failure mode is COMPLETENESS, not accuracy, and it has been found twice
with the same shape (gentamicin `rmt*` unrepresentable by a `Subclass == GENTAMICIN` rule; HIV NNRTI
drivers outside the 8 catalogued positions). Both were invisible until an independent label set
arrived. This module carries the signals that say so, in the record, at call time.

THE LOAD-BEARING CONSTRAINT. A doubt signal may QUALIFY a call, or EXPLAIN itself. It may never
overrule L1 and never emit a resistance prediction of its own. That is not a convention here -- it is
enforced: `DoubtBlock.as_dict()` runs `assert_no_call` on its own output and raises rather than emit.
The constraint is what keeps L2 out of the learned-predictor regime that is 0-for-5 de-confounded.

WHY A SURPRISE MEASURE, NOT A COUNT. The raw signature "many R carriers, no S carriers" is what
identified `rmt`, but it fires on families far too small to mean anything: the full-index run flags
ciprofloxacin `qnrA1` at 4R/0S, where the cohort's own base S-rate is 0.583 and P(all four R) = 0.030
by chance -- against ~125 uncounted families tested for that drug. Reporting that beside `rmtE1`
(36R/0S, p = 4e-12) as though they were the same kind of finding would make the doubt layer noise.
`completeness_tier` applies the family-wise correction, so STRONG means survived it.

Pure + dependency-free. No model, no network, no structures, no I/O.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import comb

# A doubt block describes UNCERTAINTY about a call. These are the shapes a CALL takes; if one appears
# inside a doubt block, the layer has crossed the line it exists to hold.
_CALL_VALUES = {"R", "S", "I", "RESISTANT", "SUSCEPTIBLE", "INTERMEDIATE", "NONSUSCEPTIBLE"}
_CALL_KEYS = {"prediction", "predicted", "call", "phenotype", "resistance", "susceptibility",
              "predicted_phenotype", "resistance_call"}

STRONG, WEAK, NONE = "strong", "weak", "none"
FAMILYWISE_ALPHA = 0.05


def binom_lower_tail(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p). Exact, pure, no scipy.

    Used one-sided: how surprising is it to see AT MOST `k` susceptible carriers among `n` labelled
    carriers of a determinant family, if that family were labelled at the cohort's own base rate?
    """
    if n <= 0:
        return 1.0
    p = min(max(p, 0.0), 1.0)
    k = min(max(k, 0), n)
    return sum(comb(n, i) * (p ** i) * ((1.0 - p) ** (n - i)) for i in range(k + 1))


def completeness_surprise(carriers_labelled_r: int, carriers_labelled_s: int,
                          base_s_rate: float) -> float | None:
    """P(zero susceptible carriers among n), or None when the purity signature is absent.

    ENRICHMENT IS THE WRONG NULL, and that is measured, not stylistic. A lower-tail binomial on the
    observed susceptible count calls gentamicin `aph(6)-Id` (62R/28S) STRONG at p ~ 5e-5, because 28
    susceptible carriers really are fewer than a 0.517 base rate predicts. But `aph(6)-Id` is a
    CORRECT exclusion -- a streptomycin determinant that travels with gentamicin resistance by
    linkage. EVERY co-occurring determinant is R-enriched, so an enrichment null floods the layer
    with correct exclusions, exactly the way raw volume ranking buried the answer in step 1.

    The `rmt` signature is PURITY: zero susceptible carriers among a well-powered set. A single
    susceptible carrier is positive evidence that the exclusion is deliberate, so it ENDS the signal
    rather than weakening it. Returns None for both "no labelled carriers" (unassessable) and "has a
    susceptible carrier" (signature absent) -- `completeness_signal` distinguishes them in prose.
    """
    n = carriers_labelled_r + carriers_labelled_s
    if n <= 0 or carriers_labelled_s > 0:
        return None
    return (1.0 - min(max(base_s_rate, 0.0), 1.0)) ** n


def completeness_tier(p: float | None, n_families_tested: int) -> str:
    """STRONG only after the family-wise correction; WEAK is nominal-only; None -> NONE.

    `n_families_tested` is the count of uncounted determinant families screened for THIS drug. A
    screen that probes ~125 families will turn up nominally-significant purity by chance, so an
    uncorrected threshold would manufacture doubt signals rather than find them.
    """
    if p is None:
        return NONE
    if n_families_tested > 0 and p <= FAMILYWISE_ALPHA / n_families_tested:
        return STRONG
    return WEAK if p <= FAMILYWISE_ALPHA else NONE


@dataclass(frozen=True)
class DoubtSignal:
    """One reason a call may be incomplete. Carries evidence and a tier -- never a prediction."""
    kind: str                 # "determinant_completeness" | "position_novelty"
    tier: str                 # strong | weak | none
    reason: str
    evidence: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"kind": self.kind, "tier": self.tier, "reason": self.reason,
                "evidence": dict(self.evidence)}


@dataclass
class DoubtBlock:
    """The `doubt` field of a decoder record. Qualifies the call; never replaces it."""
    signals: list = field(default_factory=list)

    @property
    def any_doubt(self) -> bool:
        return any(s.tier in (STRONG, WEAK) for s in self.signals)

    @property
    def max_tier(self) -> str:
        for t in (STRONG, WEAK):
            if any(s.tier == t for s in self.signals):
                return t
        return NONE

    def as_dict(self) -> dict:
        """Emit the block, REFUSING to emit anything call-shaped.

        The check runs here rather than only in tests: a doubt block that carried a prediction would
        be a product-surface falsehood, and the cheapest place to make that impossible is the one
        function every consumer goes through.
        """
        out = {
            "schema": "decoder-doubt-v1",
            "contract": ("A doubt signal qualifies the call and explains itself. It NEVER overrules "
                         "the call and NEVER emits a resistance prediction of its own."),
            "any_doubt": self.any_doubt,
            "max_tier": self.max_tier,
            "signals": [s.as_dict() for s in self.signals],
        }
        assert_no_call(out)
        return out


def doubt_one_line(block: dict) -> str | None:
    """Compact human-readable doubt line, or None when silence is the honest answer.

    WHY SILENCE IS ONLY EVER HONEST IN ONE CASE. A machine-readable block a human never sees is not
    a disclosure -- the target-site CLI carried the block in JSON for a day while its human output
    showed only a STATIC blind-spot list ("an S call can't rule out an uncatalogued substitution"),
    never that this genotype actually HAS one. So a line prints whenever the block has something to
    say, and is omitted ONLY for the assessed-and-quiet case, which is the single situation where "we
    checked and found nothing" is what silence would truthfully mean.
    """
    if not block or not block.get("signals"):
        return None
    sigs = block["signals"]
    tier = block.get("max_tier", NONE)
    # Report the reason belonging to the signal that ACTUALLY fired, not signals[0]. A block now
    # carries more than one signal (position-novelty + target-site completeness, added 2026-09-02),
    # and hardcoding the first paired a STRONG tier with the other signal's "found nothing" prose --
    # a self-contradicting line, and the worst possible failure for a human-facing disclosure.
    sig = next((s for s in sigs if s.get("tier") == tier), sigs[0])
    ev = sig.get("evidence") or {}
    if tier in (STRONG, WEAK):
        return f"DOUBT [{tier}]: {sig['reason']}"
    if ev.get("applicable") is False:
        return ("doubt: n/a -- this catalog is position-based, so the completeness flag could never "
                "fire here (NOT an absence of doubt)")
    if ev.get("assessed") is False:
        return ("doubt: NOT ASSESSED -- this input path does not surface the observed substitutions, "
                "so the completeness flag could not be evaluated (NOT a clean result)")
    return None                     # assessed, nothing found: the one honest silence


def assert_no_call(obj, _path: str = "doubt") -> None:
    """Raise if anything call-shaped appears anywhere in a doubt block. Recursive, fail-closed."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = str(k).strip().lower()
            if key in _CALL_KEYS:
                raise ValueError(f"doubt block carries a call-shaped key at {_path}.{k!r} -- L2 may "
                                 "qualify a call, never emit one")
            assert_no_call(v, f"{_path}.{k}")
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            assert_no_call(v, f"{_path}[{i}]")
    elif isinstance(obj, str):
        if obj.strip().upper() in _CALL_VALUES:
            raise ValueError(f"doubt block carries a call-shaped value {obj!r} at {_path} -- L2 may "
                             "qualify a call, never emit one")


# --- the two shipped signal kinds -----------------------------------------------------------------

def completeness_signal(symbol: str, subclass: str, carriers_labelled_r: int,
                        carriers_labelled_s: int, base_s_rate: float,
                        n_families_tested: int) -> DoubtSignal:
    """Doubt from a determinant family the deployed rule cannot represent (the `rmt` shape)."""
    p = completeness_surprise(carriers_labelled_r, carriers_labelled_s, base_s_rate)
    tier = completeness_tier(p, n_families_tested)
    n = carriers_labelled_r + carriers_labelled_s
    if n == 0:
        reason = (f"determinant {symbol!r} ({subclass}) is present but the deployed rule cannot "
                  "represent it, and no labelled carrier exists to assess it -- unassessable, not clean")
    elif carriers_labelled_s > 0:
        reason = (f"determinant {symbol!r} ({subclass}) is not represented by the deployed rule, but "
                  f"{carriers_labelled_s} of {n} labelled carriers are susceptible -- positive "
                  "evidence that the exclusion is deliberate, not a gap")
    else:
        reason = (f"determinant {symbol!r} ({subclass}) is present but the deployed rule cannot "
                  f"represent it; {carriers_labelled_r}/{n} labelled carriers are resistant against a "
                  f"cohort base susceptible-rate of {base_s_rate:.3f}")
    return DoubtSignal(kind="determinant_completeness", tier=tier, reason=reason,
                       evidence={"symbol": symbol, "subclass": subclass,
                                 "carriers_labelled_r": carriers_labelled_r,
                                 "carriers_labelled_s": carriers_labelled_s,
                                 "base_s_rate": round(base_s_rate, 4),
                                 "purity_surprise_p": p,
                                 "n_families_tested": n_families_tested,
                                 "familywise_alpha": FAMILYWISE_ALPHA})


# ONLY a MUTANT-LEVEL catalog can carry a position-novelty signal, and that is load-bearing. For a
# POSITION-based catalog (HIV NRTI / PI / INSTI) every substitution at a catalogued position is
# already called, so "a novel substitution at a catalogued position" is empty BY CONSTRUCTION. The
# flag there would be permanently silent -- and a permanently-silent signal reads as "no doubt", which
# is the exact failure this layer exists to prevent. Those cells report not-applicable instead.
_MUTANT_LEVEL_CELLS = {
    "efavirenz": "hiv-nnrti-rt", "nevirapine": "hiv-nnrti-rt", "etravirine": "hiv-nnrti-rt",
    "rilpivirine": "hiv-nnrti-rt", "doravirine": "hiv-nnrti-rt",
    "nirmatrelvir": "sarscov2-mpro", "ensitrelvir": "sarscov2-mpro", "lufotrelvir": "sarscov2-mpro",
    "fluconazole": "fungal-fluconazole-erg11", "voriconazole": "fungal-voriconazole-erg11",
}
_CELL_GENE = {"hiv-nnrti-rt": "RT", "sarscov2-mpro": "Mpro",
              "fungal-fluconazole-erg11": "ERG11", "fungal-voriconazole-erg11": "ERG11"}


def doubt_cell_for(drug: str) -> str | None:
    """The position-novelty cell for a drug, or None when the flag cannot apply. Pure."""
    return _MUTANT_LEVEL_CELLS.get(str(drug).strip().lower())


def target_site_doubt(drug: str, observed_by_gene: dict | None) -> DoubtBlock:
    """The doubt block for a target-site call. Distinguishes THREE states, never collapsing them.

    not-applicable (position-based catalog) / not-assessable (this path did not surface the observed
    substitutions) / assessed. Reporting "no doubt" for either of the first two would be a false
    clean bill of health -- exactly the failure mode the layer exists to surface.
    """
    cell = doubt_cell_for(drug)
    if cell is None:
        return DoubtBlock([DoubtSignal(
            kind="position_novelty", tier=NONE,
            reason=(f"position-novelty does not apply to {drug}: its catalog is position-based or "
                    "unregistered, so every substitution at a catalogued position is already called "
                    "and the flag could never fire -- this is NOT an absence of doubt"),
            evidence={"cell": None, "applicable": False})])
    if observed_by_gene is None:
        return DoubtBlock([DoubtSignal(
            kind="position_novelty", tier=NONE,
            reason=("observed substitutions were not surfaced on this input path, so NEITHER the "
                    "position-novelty flag nor the catalog-completeness screen could be evaluated "
                    "-- not assessable, NOT clean"),
            evidence={"cell": cell, "applicable": True, "assessed": False})])
    subs = observed_by_gene.get(_CELL_GENE.get(cell, ""), set()) or set()
    sig = position_novelty_signal(sorted(subs), cell)
    # SECOND, COMPLEMENTARY signal (2026-09-02). position-novelty fires only at CATALOGUED positions,
    # so it is silent on a gap at a position the catalog does not carry -- verified: V179F returns
    # position_novel=False. Appended, never merged: the two answer different questions and collapsing
    # them would hide whichever fired.
    return DoubtBlock([sig, target_site_completeness_signal(sorted(subs), cell)])


def target_site_completeness_signal(observed_substitutions, cell: str) -> DoubtSignal:
    """Doubt from a substitution at a position the target-site catalog does NOT carry (the V179F shape).

    Mirrors the AMR arm's `completeness_signal` vocabulary -- the 2026-09-02 probe measured that the
    purity signature is well-formed here, so this is one vocabulary rather than two. Three states, never
    collapsed: not-measured for this cell / measured-and-no-hit / hit.
    """
    from ..data.target_site_completeness import (completeness_units_for, is_measured,
                                                  matching_units)

    if not is_measured(cell):
        return DoubtSignal(
            kind="target_site_completeness", tier=NONE,
            reason=(f"catalog-completeness has NOT been measured for {cell} -- no free isolate-level "
                    "phenotype source, or measured and found underpowered. This is an absence of "
                    "measurement, NOT an absence of doubt"),
            evidence={"cell": cell, "measured": False})

    hits = matching_units(observed_substitutions, cell)
    if not hits:
        return DoubtSignal(
            kind="target_site_completeness", tier=NONE,
            reason=(f"no observed substitution matches a measured completeness gap in {cell}; the "
                    "screen ran and found nothing, which is a result rather than an assumption"),
            evidence={"cell": cell, "measured": True,
                      "n_known_gaps": len(completeness_units_for(cell))})

    sub, st = hits[0]
    n = st["carriers_labelled_r"] + st["carriers_labelled_s"]
    tier = completeness_tier(st["purity_surprise_p"], st["n_units_tested"])
    reason = (f"substitution {sub} sits at a position the {cell} catalog does not carry, and "
              f"{st['carriers_labelled_r']}/{n} of its labelled carriers are resistant against a base "
              f"susceptible-rate of {st['base_s_rate']:.3f} -- a measured catalog-completeness gap, so a "
              "susceptible call here is the least trustworthy kind")
    return DoubtSignal(kind="target_site_completeness", tier=tier, reason=reason,
                       evidence={"cell": cell, "measured": True, "substitution": sub,
                                 "all_matches": [h[0] for h in hits], **{
                                     k: st[k] for k in ("carriers_labelled_r", "carriers_labelled_s",
                                                        "purity_surprise_p", "n_units_tested",
                                                        "base_s_rate", "scored_on", "artifact")}})


def position_novelty_signal(observed_substitutions, cell: str) -> DoubtSignal:
    """Doubt from a novel substitution at a catalogued target-site position (the HIV shape).

    Delegates to the shipped `position_novelty` flag rather than restating its logic. Measured
    incumbent: median sensitivity 0.604 on the EFV catalog-negative blind spot, lift 4.69, no tools.
    """
    from .position_novelty import flag_for_cell

    res = flag_for_cell(observed_substitutions, cell)
    tier = WEAK if res.position_novel else NONE
    if res.position_novel:
        reason = (f"substitution(s) {', '.join(res.novel_substitutions)} sit at catalogued "
                  f"{cell} resistance positions but are not themselves catalogued -- a "
                  "susceptible-by-absence call is least trustworthy here")
    else:
        reason = f"no uncatalogued substitution at any of the {res.n_catalog_positions} catalogued " \
                 f"{cell} positions"
    return DoubtSignal(kind="position_novelty", tier=tier, reason=reason,
                       evidence={"cell": cell, **res.as_dict()})
