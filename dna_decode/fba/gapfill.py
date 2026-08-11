"""Find and repair the MISSING parts of a metabolic model — the honest form of "predict what's absent".

The design cell asks "given an edit, what trait?" and "given a product, which edits?". This module asks the
third question, the one that bounds both: **what is the model missing?** A genome-scale model is a curated
guess about an organism's biochemistry, and where it is wrong the whole cell is wrong with it — today's
worked example is a measured FALSE NEGATIVE, not a hypothetical.

Two capabilities, deliberately separated because they carry very different evidential weight:

1. **`dead_end_metabolites` — a STRUCTURAL diagnostic that needs no donor, no labels and no network.**
   A metabolite that is produced but never consumed (or consumed but never produced) cannot carry steady-state
   flux, so every reaction that makes it is dead weight. This is a fact about the model, not a claim about
   biology, and it is where the missing parts actually show up.

2. **`propose_repair` — a HYPOTHESIS generator.** Given a donor model (a related organism's reconstruction),
   ask which reactions would restore a trait the model wrongly predicts as absent. cobrapy's `gapfill` does
   the search; this wraps it with the verification that makes the answer trustworthy.

**The honesty rail that matters most here:** a gap is NOT automatically a defect. A model that says "this
organism cannot eat X" may be exactly right — organisms genuinely lack capabilities, and "repairing" a
correct model would be fabricating biology. So a proposal is only ever a candidate, and it is worth acting on
only when there is INDEPENDENT evidence the organism really has the capability. `verify_repair` reports the
restored trait AND a specificity check, because a repair that makes the model grow on everything has
explained nothing.

Worked example (measured end to end, 2026-08-07 — see `wiki/fba_gapfill_2026-08-07.md`):
iML1515 predicts NO growth on sucrose, but BW25113 has a sucrose carbon-source experiment in the
Wetmore/Keio RB-TnSeq set (that assay only runs sources the organism grows on). The structural diagnostic
localises it exactly — `suc6p_c` is a DEAD END, produced by `SUCptspp` and consumed by nothing, so the model
carries a sucrose transporter that leads nowhere. Gap-filling against *Salmonella* iYS1720 proposes a single
reaction, `FFSD` (`h2o_c + suc6p_c --> fru_c + g6p_c`), and adding it takes sucrose growth from **0.000 to
1.7798 /h** — about twice the glucose rate (0.877), which is what a disaccharide should give.
"""
from __future__ import annotations

from dataclasses import dataclass

# A stoichiometric coefficient below this is numerically zero.
COEF_TOL = 1e-9


@dataclass(frozen=True)
class DeadEnd:
    """A metabolite that cannot carry steady-state flux, and which side of it is missing."""
    metabolite: str
    kind: str          # "no_consumer" (produced, never used) | "no_producer" (used, never made)
    reactions: tuple[str, ...]

    def as_dict(self) -> dict:
        return {"metabolite": self.metabolite, "kind": self.kind, "reactions": list(self.reactions)}


def find_dead_ends(stoichiometry: dict[str, dict[str, float]], reversible: set[str] | None = None) -> list[DeadEnd]:
    """PURE: {reaction_id: {metabolite_id: coefficient}} -> the dead-end metabolites.

    A metabolite is a dead end when every reaction touching it sits on the same side: all produce it and
    none consume it, or the reverse. **A REVERSIBLE reaction counts as both** a producer and a consumer —
    it can run either way, so treating it as one-directional would report dead ends that are not dead.

    Kept pure over a plain dict so the logic is testable without cobra and without a model download.
    """
    reversible = reversible or set()
    producers: dict[str, set[str]] = {}
    consumers: dict[str, set[str]] = {}
    for rid, mets in stoichiometry.items():
        for met, coef in mets.items():
            if abs(coef) <= COEF_TOL:
                continue
            if coef > 0 or rid in reversible:
                producers.setdefault(met, set()).add(rid)
            if coef < 0 or rid in reversible:
                consumers.setdefault(met, set()).add(rid)

    out: list[DeadEnd] = []
    for met in sorted(set(producers) | set(consumers)):
        p, c = producers.get(met, set()), consumers.get(met, set())
        if p and not c:
            out.append(DeadEnd(met, "no_consumer", tuple(sorted(p))))
        elif c and not p:
            out.append(DeadEnd(met, "no_producer", tuple(sorted(c))))
    return out


def model_dead_ends(model, exclude_boundary: bool = True) -> list[DeadEnd]:
    """Dead-end metabolites of a live cobra model (thin adapter over the pure `find_dead_ends`).

    Boundary reactions (exchanges / demands / sinks) are excluded by default: they exist precisely to be
    a one-sided source or sink, so counting them would hide the real dead ends behind hundreds of
    uninteresting ones.
    """
    stoich, rev = {}, set()
    for r in model.reactions:
        if exclude_boundary and r.id.startswith(("EX_", "DM_", "SK_")):
            continue
        stoich[r.id] = {m.id: c for m, c in r.metabolites.items()}
        if r.lower_bound < 0 < r.upper_bound:
            rev.add(r.id)
    return find_dead_ends(stoich, rev)


def orphan_uptake_targets(model) -> list[DeadEnd]:
    """Dead ends that a TRANSPORT reaction produces — the structurally suspicious ones.

    A model carrying machinery to import a nutrient, whose product nothing then consumes, is making a claim
    it cannot cash: either the transporter should not be there, or the downstream enzyme is missing. That
    pattern is exactly how the sucrose gap presents (`SUCptspp` makes `suc6p_c`; nothing uses it), so it is
    a good first place to look before reaching for a donor model.
    """
    out = []
    for de in model_dead_ends(model):
        if de.kind != "no_consumer":
            continue
        if any("t" in rid.lower()[1:] or "pts" in rid.lower() for rid in de.reactions):
            out.append(de)
    return out


def growth_on_sole_carbon(model, exchange: str, base_medium: dict, default_carbon: str = "EX_glc__D_e",
                          uptake: float = 10.0) -> float:
    """Growth (/h) with `exchange` swapped in as the sole carbon source. 0.0 for no-growth/infeasible."""
    if exchange not in model.reactions:
        return 0.0
    med = {k: v for k, v in base_medium.items() if k != default_carbon}
    med[exchange] = uptake
    with model:
        model.medium = med
        v = model.slim_optimize()
    return 0.0 if (v is None or v != v) else float(v)


def propose_repair(model, donor, target_exchange: str, *, lower_bound: float = 0.05,
                   base_medium: dict | None = None, default_carbon: str = "EX_glc__D_e") -> dict:
    """Which donor reactions would let `model` grow on `target_exchange`? Returns candidates + evidence.

    The search is cobrapy's `gapfill` over the donor reactions the model lacks. `demand_reactions` and
    `exchange_reactions` are BOTH off: letting the solver invent a demand or an exchange lets it satisfy the
    objective by inventing a sink rather than by finding the missing biochemistry, which is not a repair.

    Every proposal is a HYPOTHESIS. `verify_repair` is what turns it into a measured claim.
    """
    from cobra import Model  # noqa: PLC0415 (lazy heavy import)
    from cobra.flux_analysis import gapfill  # noqa: PLC0415

    base = dict(base_medium if base_medium is not None else model.medium)
    before = growth_on_sole_carbon(model, target_exchange, base, default_carbon)

    have = {r.id for r in model.reactions}
    universal = Model("universal_donor")
    universal.add_reactions([r.copy() for r in donor.reactions if r.id not in have])

    med = {k: v for k, v in base.items() if k != default_carbon}
    med[target_exchange] = 10.0
    with model:
        model.medium = med
        try:
            solutions = gapfill(model, universal, lower_bound=lower_bound,
                                demand_reactions=False, exchange_reactions=False, iterations=1)
        except Exception as e:  # infeasible / no solution -- report, never fabricate one
            return {"target_exchange": target_exchange, "growth_before": round(before, 6),
                    "status": f"no_solution: {type(e).__name__}", "candidates": [],
                    "n_donor_reactions": len(universal.reactions)}

    cands = [[{"id": r.id, "reaction": r.reaction, "name": getattr(r, "name", "")} for r in sol]
             for sol in solutions]
    return {
        "target_exchange": target_exchange,
        "growth_before": round(before, 6),
        "status": "ok" if cands else "no_solution",
        "n_donor_reactions": len(universal.reactions),
        "candidates": cands,
        "claim_status": "HYPOTHESIS — a gap may be correct biology; verify the organism truly has the trait",
    }


def verify_repair(model, donor, reaction_ids: list[str], target_exchange: str, *,
                  specificity_exchanges: tuple[str, ...] = (),
                  base_medium: dict | None = None, default_carbon: str = "EX_glc__D_e") -> dict:
    """Add `reaction_ids` from the donor and MEASURE what changed — the trait, and the collateral.

    Reports growth on the target before and after, plus growth on `specificity_exchanges` before and after.
    The specificity half is the point: a repair that also lifts growth on unrelated carbon sources has made
    the model permissive rather than correct, and that is indistinguishable from success if you only look at
    the target.
    """
    base = dict(base_medium if base_medium is not None else model.medium)
    targets = (target_exchange,) + tuple(specificity_exchanges)
    before = {ex: round(growth_on_sole_carbon(model, ex, base, default_carbon), 6) for ex in targets}
    with model:
        model.add_reactions([donor.reactions.get_by_id(r).copy() for r in reaction_ids])
        after = {ex: round(growth_on_sole_carbon(model, ex, base, default_carbon), 6) for ex in targets}

    changed = [ex for ex in specificity_exchanges if abs(after[ex] - before[ex]) > 1e-6]
    return {
        "added_reactions": list(reaction_ids),
        "target_exchange": target_exchange,
        "growth_before": before[target_exchange],
        "growth_after": after[target_exchange],
        "repaired": after[target_exchange] > 1e-4 >= before[target_exchange],
        "specificity_before": {k: before[k] for k in specificity_exchanges},
        "specificity_after": {k: after[k] for k in specificity_exchanges},
        "specificity_unchanged": not changed,
        "specificity_changed_exchanges": changed,
        "scope": (
            "A repair is a HYPOTHESIS about the organism's biochemistry, not a validated fact. It is worth "
            "acting on only with INDEPENDENT evidence the organism has the trait; otherwise the model may "
            "have been right and the 'repair' fabricates biology."
        ),
    }
