"""Growth-COUPLED strain design — the DESIGN direction of the FBA cell.

The rest of `fba/` answers "given an edit, what happens?" (edit -> growth / essentiality). This module
answers the inverse, which is what strain engineering actually needs:

    given a PRODUCT, which gene knockouts make producing it NECESSARY for the cell to grow?

That property is **growth coupling**, and it is the thing that makes an engineered strain stable: if the
product can only be made when the cell grows, then selection for growth selects for production, instead
of fighting it. The test is a two-sided LP at a fixed growth floor:

    fix   biomass >= growth_frac * THIS STRAIN'S OWN maximum growth   (near-optimal, e.g. 90%)
    then  MINIMISE the product exchange   -> `min_flux`
    and   MAXIMISE the product exchange   -> `max_flux`

`min_flux > 0` means **every** flux distribution that grows at that rate also secretes the product —
production is obligatory, not merely allowed. `max_flux > 0` with `min_flux == 0` means the cell *can*
make it and *can* avoid it; that is the un-engineered case and is NOT a design.

Honest scope, stated once and stamped into every record:
  * This is a stoichiometric prediction. A coupled design is a **hypothesis for the bench**, never a
    validated strain. FBA sees mass balance and reaction knockouts; it does not see regulation, enzyme
    kinetics, toxicity, burden, or whether the knockout strain is actually constructible.
  * Coupling is evaluated at the model's default medium/exchange bounds unless the caller changes them.

`cobra.design`'s OptKnock was removed from cobrapy (absent in 0.31), so the search here is an explicit
bounded enumeration rather than a MILP: exhaustive single knockouts, then pairs and triples over the most
promising singles. Slower asymptotically but transparent, and it never truncates without saying so.

**Three settings are load-bearing, each verified against a known design on 2026-08-07.** Get any one wrong
and the search silently returns zero designs while appearing to work:
  1. **REACTION-level knockouts** (default). GPR isozymes mean a single gene deletion often leaves the
     reaction fully open (`apply_knockout`).
  2. **A NEAR-OPTIMAL growth floor, relative to each strain's OWN maximum** (default 0.9). At a low floor
     the cell has enough slack to avoid producing, so `min_flux` collapses to ~0 for everything
     (`evaluate_knockouts`).
  3. **Depth 3.** The classic fermentation-route designs need three knockouts (`find_coupled_designs`).
End-to-end gate: with all three, the search recovers the OptKnock-lineage anaerobic succinate design
`PFL + LDH_D + ALCD2x` at a guaranteed 9.26 mmol/gDW/h against a wild-type floor of 0.047
(`tests/test_fba_design.py::test_reproduces_the_literature_anaerobic_succinate_design`).
"""
from __future__ import annotations

from dataclasses import dataclass, field

# A flux below this is numerically zero (LP tolerances put "no flux" at ~1e-9, not exactly 0).
FLUX_TOL = 1e-6

OBLIGATORY = "OBLIGATORY"   # min_flux > tol -> growth REQUIRES production (a growth-coupled design)
POSSIBLE = "POSSIBLE"       # max_flux > tol but min_flux ~ 0 -> can produce, can also avoid it
INFEASIBLE = "INFEASIBLE"   # max_flux ~ 0 -> cannot produce at this growth floor at all


def classify_coupling(min_flux: float | None, max_flux: float | None, tol: float = FLUX_TOL) -> str:
    """PURE: (min,max) product flux at a growth floor -> the coupling verdict.

    A None/NaN bound means the LP was infeasible at that growth floor, which is not "no production" —
    it is "cannot grow", so it cannot be a design either. Treated as INFEASIBLE.
    """
    if min_flux is None or max_flux is None:
        return INFEASIBLE
    if min_flux != min_flux or max_flux != max_flux:  # NaN
        return INFEASIBLE
    if max_flux <= tol:
        return INFEASIBLE
    return OBLIGATORY if min_flux > tol else POSSIBLE


@dataclass
class Design:
    """One candidate strain design: a knockout set + what it does to product flux and growth."""
    knockouts: tuple[str, ...]
    min_product: float
    max_product: float
    growth: float
    coupling: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "knockouts": list(self.knockouts),
            "n_knockouts": len(self.knockouts),
            "min_product_flux": round(self.min_product, 6),
            "max_product_flux": round(self.max_product, 6),
            "growth_per_h": round(self.growth, 6),
            "coupling": self.coupling,
            "notes": self.notes,
        }


def improves_on_baseline(design_min: float, baseline_min: float, tol: float = FLUX_TOL) -> bool:
    """PURE: does this knockout set actually RAISE the guaranteed product floor above wild type?

    Load-bearing, and learned the hard way (2026-08-07): anaerobic *E. coli* already secretes a little
    succinate obligatorily, so the wild type is OBLIGATORY before any edit. Counting "is OBLIGATORY" as
    "is a design" then reported **2096 of 2096 evaluated knockouts as growth-coupled designs** — every
    one of them merely inheriting the baseline. A design is only a design if it beats the baseline.
    """
    return design_min > baseline_min + tol


def rank_designs(designs: list[Design]) -> list[Design]:
    """PURE: best design first.

    Ranked by GUARANTEED product (`min_product`) before anything else — a design's whole value is the
    floor it enforces, not the ceiling it permits. Ties break toward higher growth (a coupled strain
    that barely grows is not useful), then toward FEWER knockouts (each one is real bench work).
    """
    return sorted(designs, key=lambda d: (-d.min_product, -d.growth, len(d.knockouts)))


def biomass_reaction(model):
    """The model's objective (biomass) reaction. Raises if the objective is not a single reaction."""
    from cobra.util.solver import linear_reaction_coefficients  # noqa: PLC0415 (lazy heavy import)

    coeffs = linear_reaction_coefficients(model)
    rxns = [r for r, c in coeffs.items() if c != 0]
    if len(rxns) != 1:
        raise ValueError(
            f"expected exactly one objective reaction, found {len(rxns)}: {[r.id for r in rxns]}. "
            "Set a single biomass objective before running a design search."
        )
    return rxns[0]


def resolve_target(model, target: str):
    """Resolve a product target to a reaction. Accepts a reaction id or a bare metabolite-ish name.

    Tries, in order: exact reaction id -> `EX_<target>_e` -> `EX_<target>` -> a unique case-insensitive
    substring match among EXCHANGE reactions. Raises with candidates rather than guessing ambiguously.
    """
    for candidate in (target, f"EX_{target}_e", f"EX_{target}"):
        if candidate in model.reactions:
            return model.reactions.get_by_id(candidate)
    t = target.lower()
    hits = [r for r in model.reactions if r.id.startswith("EX_") and t in r.id.lower()]
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise ValueError(f"no exchange reaction matches target '{target}'")
    raise ValueError(
        f"target '{target}' is ambiguous across {len(hits)} exchanges: {[r.id for r in hits[:8]]}"
    )


def product_range_at_growth(model, target_rxn, growth_frac: float, wt_growth: float) -> tuple:
    """Two-sided product flux with growth pinned at >= growth_frac x wild-type. Returns (min,max,growth).

    Caller is responsible for any knockouts already applied. Restores all bounds on exit.
    """
    bio = biomass_reaction(model)
    floor = growth_frac * wt_growth
    with model:
        bio.lower_bound = floor
        model.objective = target_rxn
        model.objective_direction = "min"
        mn = model.slim_optimize()
        model.objective_direction = "max"
        mx = model.slim_optimize()
        # the growth actually achieved while producing at the ceiling
        model.objective = bio
        model.objective_direction = "max"
        gr = model.slim_optimize()
    return mn, mx, gr


def apply_knockout(model, ident: str, level: str = "gene") -> None:
    """Knock out a gene or a reaction, in place (caller supplies the `with model:` context).

    **GENE level is strictly weaker than REACTION level, and not by a little.** iML1515's GPRs carry
    isozymes — `PFL` is `(b0902 and b0903) or (b0902 and b3114) or (b3951 and b3952) or ...` — so
    deleting one gene leaves the reaction fully open. Verified 2026-08-07: knocking out `b0903` (pflB)
    left PFL at bounds (0, 1000) and anaerobic growth unchanged at 0.15754 /h. Any search bounded to
    1-2 GENE deletions therefore cannot reach designs that require disabling an isozyme-backed reaction,
    which is most of the classic fermentation-route designs.
    """
    if level == "reaction":
        model.reactions.get_by_id(ident).knock_out()
    else:
        model.genes.get_by_id(ident).knock_out()


def evaluate_knockouts(
    model, target_rxn, knockouts, growth_frac: float, wt_growth: float, level: str = "gene"
) -> Design:
    """Apply `knockouts` and evaluate coupling of `target_rxn` at NEAR-OPTIMAL growth. Non-mutating.

    The growth floor is `growth_frac` x **the mutant's OWN maximum growth**, not the wild type's. This
    is the standard growth-coupling definition and it is load-bearing (verified 2026-08-07): a floor
    taken relative to the wild type is unreachable for any knockout that slows the cell, and a LOW floor
    leaves so much metabolic slack that the minimum product is ~0 for every design. On anaerobic
    succinate, the literature design `PFL+LDH_D+ALCD2x` guarantees 0.0027 at 10% of the mutant optimum
    and **12.28 at 99%** — same strain, same model; only the floor differs.
    """
    with model:
        for k in knockouts:
            apply_knockout(model, k, level)
        mutant_max = model.slim_optimize()
        if mutant_max is None or mutant_max != mutant_max or mutant_max <= 0:
            return Design(tuple(knockouts), 0.0, 0.0, 0.0, INFEASIBLE, ["knockout is lethal"])
        mn, mx, gr = product_range_at_growth(model, target_rxn, growth_frac, mutant_max)
    return Design(
        knockouts=tuple(knockouts),
        min_product=0.0 if (mn is None or mn != mn) else mn,
        max_product=0.0 if (mx is None or mx != mx) else mx,
        growth=0.0 if (gr is None or gr != gr) else gr,
        coupling=classify_coupling(mn, mx),
        notes=[f"floor = {growth_frac:.0%} of this strain's own max growth ({mutant_max:.4f} /h)"],
    )


def set_anaerobic(model) -> bool:
    """Close O2 uptake in place. Returns False if the model has no O2 exchange.

    Condition matters enormously for coupling: with respiration available a cell can oxidise substrate
    fully to CO2 and never has to secrete a fermentation product, so aerobic coupling is rare. Anaerobic
    growth forces redox balancing through secreted products, which is why the classic growth-coupled
    designs in the literature are anaerobic or micro-aerobic.
    """
    if "EX_o2_e" not in model.reactions:
        return False
    model.reactions.get_by_id("EX_o2_e").lower_bound = 0.0
    return True


def find_coupled_designs(
    model,
    target: str,
    *,
    growth_frac: float = 0.9,
    max_knockouts: int = 3,
    wt_growth: float | None = None,
    pair_pool: int = 40,
    triple_pool: int = 18,
    gene_ids: list[str] | None = None,
    anaerobic: bool = False,
    level: str = "reaction",
    min_growth_frac_of_wt: float = 0.05,
) -> dict:
    """Search knockout sets that make `target` production obligatory for growth.

    Strategy (explicit, and reported so it is never a silent truncation):
      1. Baseline: the wild-type coupling. If already OBLIGATORY there is nothing to design.
      2. Exhaustive SINGLE knockouts over every gene that leaves the cell viable at the growth floor.
      3. PAIRS over the top `pair_pool` singles, ranked by guaranteed product. This is a bounded
         heuristic, NOT an exhaustive double search — the returned dict says so in `search`.
    """
    target_rxn = resolve_target(model, target)
    # The whole search runs inside ONE context so a condition change (anaerobic) applies to every
    # knockout evaluation and is unwound exactly once, on exit.
    with model:
        condition = "model default (aerobic if O2 uptake is open)"
        if anaerobic:
            if not set_anaerobic(model):
                raise ValueError("anaerobic requested but the model has no EX_o2_e exchange")
            condition = "ANAEROBIC (O2 uptake closed)"
            wt_growth = None  # a passed-in aerobic WT would be the wrong reference under anaerobiosis
        if wt_growth is None:
            wt_growth = model.slim_optimize()
        if wt_growth is None or wt_growth != wt_growth or wt_growth <= 0:
            raise ValueError(f"model cannot grow under this condition ({condition}); no design possible")
        floor = growth_frac * wt_growth

        base = evaluate_knockouts(model, target_rxn, [], growth_frac, wt_growth, level)
        base.notes.append("wild type (no knockouts)")

        if gene_ids is not None:
            candidates = list(gene_ids)
        elif level == "reaction":
            bio_id = biomass_reaction(model).id
            # A strain design must be a GENETIC edit, so every candidate must have an associated gene.
            # Requiring a non-empty GPR excludes, in one rule, all the things that are not constructible:
            # exchanges/demands/sinks (boundary pseudo-reactions), spontaneous reactions, and -- the one
            # that actually bit on 2026-08-07 -- ATPM, the non-growth-associated ATP MAINTENANCE
            # pseudo-reaction, which the unguarded search ranked as its TOP succinate "design". You
            # cannot delete a maintenance requirement with a knockout; there is no gene for it.
            candidates = [
                r.id for r in model.reactions
                if not r.id.startswith(("EX_", "DM_", "SK_"))
                and r.id not in (bio_id, target_rxn.id)
                and (r.gene_reaction_rule or "").strip()
            ]
        else:
            candidates = [g.id for g in model.genes]

        viable_floor = min_growth_frac_of_wt * wt_growth   # a design that barely grows is not a design
        singles: list[Design] = []
        n_nonviable = 0
        for cid in candidates:
            with model:
                apply_knockout(model, cid, level)
                g = model.slim_optimize()
            if g is None or g != g or g < viable_floor:
                n_nonviable += 1
                continue
            singles.append(evaluate_knockouts(model, target_rxn, [cid], growth_frac, wt_growth, level))

        ranked_singles = rank_designs(singles)
        pairs: list[Design] = []
        triples: list[Design] = []
        n_pairs_tested = n_triples_tested = 0
        pool = [d.knockouts[0] for d in ranked_singles[:pair_pool]]
        if max_knockouts >= 2:
            for i, a in enumerate(pool):
                for b in pool[i + 1:]:
                    n_pairs_tested += 1
                    pairs.append(
                        evaluate_knockouts(model, target_rxn, [a, b], growth_frac, wt_growth, level)
                    )
        if max_knockouts >= 3:
            # Triples are where the classic fermentation-route designs live (anaerobic succinate needs
            # PFL + LDH_D + ALCD2x -- unreachable at 2). Pool is trimmed because C(n,3) grows fast.
            tri_pool = pool[:triple_pool]
            for i, a in enumerate(tri_pool):
                for j in range(i + 1, len(tri_pool)):
                    for c in tri_pool[j + 1:]:
                        n_triples_tested += 1
                        triples.append(evaluate_knockouts(
                            model, target_rxn, [a, tri_pool[j], c], growth_frac, wt_growth, level))

    all_designs = rank_designs(singles + pairs + triples)
    # A design must BEAT the wild-type guarantee, not merely inherit it (see `improves_on_baseline`).
    coupled = [
        d for d in all_designs
        if d.coupling == OBLIGATORY and improves_on_baseline(d.min_product, base.min_product)
    ]
    n_inherit = sum(
        1 for d in all_designs
        if d.coupling == OBLIGATORY and not improves_on_baseline(d.min_product, base.min_product)
    )

    def _with_gain(d: Design) -> dict:
        rec = d.as_dict()
        rec["improvement_over_wildtype"] = round(d.min_product - base.min_product, 6)
        return rec

    return {
        "target_reaction": target_rxn.id,
        "target_name": getattr(target_rxn, "name", "") or target_rxn.id,
        "condition": condition,
        "knockout_level": level,
        "growth_floor_frac": growth_frac,
        "growth_floor_basis": "fraction of EACH strain's OWN maximum growth (standard coupling definition)",
        "min_growth_frac_of_wt": min_growth_frac_of_wt,
        "growth_floor_per_h": round(floor, 6),
        "wildtype_growth_per_h": round(wt_growth, 6),
        "baseline": base.as_dict(),
        "n_candidates_scanned": len(candidates),
        "n_nonviable_at_floor": n_nonviable,
        "n_singles_evaluated": len(singles),
        "n_pairs_evaluated": n_pairs_tested,
        "n_triples_evaluated": n_triples_tested,
        "search": (
            f"exhaustive single {level} knockouts over {len(candidates)} candidates; "
            + (
                f"pairs are a BOUNDED heuristic over the top {min(pair_pool, len(ranked_singles))} "
                f"singles: {n_pairs_tested} pairs, {n_triples_tested} triples. NOT exhaustive at depth>1"
                if max_knockouts >= 2 else "only single knockouts searched (max_knockouts=1)"
            )
        ),
        "wildtype_already_coupled": base.coupling == OBLIGATORY,
        "n_coupled_designs": len(coupled),
        "n_inheriting_baseline_only": n_inherit,
        "designs": [_with_gain(d) for d in coupled[:20]],
        "best_uncoupled": [_with_gain(d) for d in all_designs[:5]] if not coupled else [],
        "scope": (
            "STOICHIOMETRIC prediction. A coupled design is a HYPOTHESIS FOR THE BENCH, not a validated "
            "strain: FBA does not model regulation, enzyme kinetics, toxicity, metabolic burden, or "
            "whether the knockout strain is constructible."
        ),
    }
