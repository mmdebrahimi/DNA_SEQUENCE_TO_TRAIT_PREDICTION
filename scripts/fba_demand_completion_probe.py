"""For each gene the model never reacts to at all: WHICH missing biomass demand would make it essential?

WHERE THIS SITS
`wiki/fba_missed_gene_partition_2026-08-22.md` partitioned the 131 conditionally-essential genes and found
**34 with no growth effect anywhere** -- deleting them changes growth by nothing, in all 25 labelled
conditions. They are not unreachable (that is a separate, proved set of 31); the model is structurally
capable of caring about them and simply never does. The named cause was "objective incompleteness".

That was a label, not a measurement. This script measures it.

THE TEST (the constructive inverse of a knockout)
For gene *g* and each candidate metabolite *m* its reactions touch:

    maximise  DM_m          subject to   biomass >= 10% of wildtype

computed twice -- once with *g* present, once with *g* deleted. If the model can make *m* normally but
**cannot make it at all without g**, then *g* is the sole route to *m*, and a biomass equation that
demanded *m* would call *g* essential. That names the missing demand, per gene, from the reconstruction
alone.

The growth floor is what makes it meaningful: without it the answer would include routes that are only
available to a dead cell.

DIAGNOSTIC ONLY -- AND THAT SCOPE IS LOAD-BEARING
This identifies *which* demand would flip *which* gene. It does **NOT** add demands to the objective and
re-score, and it does **NOT** claim a recall improvement. That would be a NEW ENDPOINT on the same data,
and this project retracted a result once for exactly that move. If the diagnostic looks promising, the
scoring run gets its own pre-registration first.

CURRENCY METABOLITES ARE EXCLUDED. "Can the cell still make ATP without gene g" is a question about
central metabolism, not about a missing biomass component, and including them would swamp the signal
with trivially-shared cofactors.

Model-only for the probe; the gene list comes from a committed artifact.

Usage:
    uv run python scripts/fba_demand_completion_probe.py
    uv run python scripts/fba_demand_completion_probe.py --growth-floor 0.10
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.model import load_model  # noqa: E402

PARTITION_ARTIFACT = Path("wiki/fba_missed_partition_2026-08-22.json")
RATIOS_ARTIFACT = Path("wiki/fba_eflux_bridge_2026-08-17_ratios.json")
PANEL_ARTIFACT = Path("wiki/fba_conditional_carbon_2026-08-13.json")

EPS = 1e-6

#: Shared cofactors / currency. A demand on these asks "can the cell still run central metabolism",
#: which is not the question -- and every gene would answer it the same way.
CURRENCY_BASES = {
    "h", "h2o", "atp", "adp", "amp", "pi", "ppi", "pppi", "nad", "nadh", "nadp", "nadph",
    "co2", "o2", "nh4", "coa", "accoa", "gtp", "gdp", "gmp", "ctp", "cdp", "cmp", "utp",
    "udp", "ump", "itp", "idp", "imp", "datp", "dadp", "damp", "dgtp", "dctp", "dttp",
    "so4", "so3", "h2s", "fad", "fadh2", "fmn", "fmnh2", "q8", "q8h2", "mql8", "mqn8",
    "2dmmql8", "2dmmq8", "thf", "mlthf", "methf", "10fthf", "5mthf", "nadph", "trdrd",
    "trdox", "grxrd", "grxox", "gthrd", "gthox", "na1", "k", "cl", "ca2", "mg2", "fe2",
    "fe3", "zn2", "mn2", "cu2", "cobalt2", "ni2", "mobd", "cbl1", "adocbl", "amet",
    "ahcys", "met__L", "hcys__L", "glu__L", "gln__L", "akg", "pyr", "acald", "actp",
    "ac", "for", "etoh", "lac__D", "succ", "fum", "mal__L", "oaa", "cit", "icit", "g3p",
    "pep", "3pg", "2pg", "13dpg", "f6p", "g6p", "fdp", "dhap", "r5p", "ru5p__D", "xu5p__D",
    "s7p", "e4p", "6pgc", "6pgl", "ade", "adn", "gua", "gsn", "ura", "urate", "hxan",
}


def base_name(met_id: str) -> str:
    """`adphep_D_c` -> `adphep_D`. Strips only the trailing compartment tag."""
    for suf in ("_c", "_p", "_e"):
        if met_id.endswith(suf):
            return met_id[: -len(suf)]
    return met_id


def is_currency(met_id: str) -> bool:
    return base_name(met_id) in CURRENCY_BASES


def never_fires_genes() -> list[str]:
    """NEVER_FIRES genes inside the prior 131-gene denominator (the arc's consistent base)."""
    part = json.loads(PARTITION_ARTIFACT.read_text(encoding="utf-8"))
    prior: set[str] = set()
    for cells in json.loads(RATIOS_ARTIFACT.read_text(encoding="utf-8"))["arms"]["baseline"].values():
        prior |= set(cells)
    return sorted(r["gene"] for r in part["genes"]
                  if r["class"] == "NEVER_FIRES" and r["gene"] in prior)


def split_no_effect(model, genes: list[str], conds: dict[str, str],
                    single_gene_deletion) -> tuple[list[str], dict[str, float]]:
    """Split NEVER_FIRES into NO-EFFECT vs SUB-THRESHOLD by worst deletion ratio over the panel.

    The 34/29 split was measured inline when the partition memo was corrected but never persisted, so it
    is recomputed here rather than hardcoded -- a hardcoded list would silently rot if the panel changed.
    """
    from dna_decode.fba.model import wildtype_growth

    all_ex = tuple(conds.values())
    worst = {g: 1.0 for g in genes}
    for cond, ex in sorted(conds.items()):
        with model:
            medium = dict(model.medium)
            for e in set(all_ex) | {"EX_glc__D_e"}:
                medium.pop(e, None)
            medium[ex] = 10.0
            model.medium = medium
            wt = wildtype_growth(model)
            if not (wt > 1e-9):
                continue
            res = single_gene_deletion(
                model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
            for _, row in res.iterrows():
                gid = next(iter(row["ids"]))
                v = row["growth"]
                worst[gid] = min(worst[gid], 0.0 if v != v else float(v) / wt)
    no_effect = sorted(g for g, r in worst.items() if r > 0.999)
    return no_effect, worst


def max_demand(model, met_id: str, biomass_id: str, floor: float, knock: str | None) -> float:
    """Max producible flux of `met_id` while the cell still grows at >= `floor`. PURE (context-managed)."""
    import cobra

    with model:
        if knock is not None:
            model.genes.get_by_id(knock).knock_out()
        model.reactions.get_by_id(biomass_id).lower_bound = floor
        dm = cobra.Reaction("ZZ_demand_probe")
        dm.lower_bound, dm.upper_bound = 0.0, 1000.0
        model.add_reactions([dm])
        dm.add_metabolites({model.metabolites.get_by_id(met_id): -1.0})
        model.objective = dm
        v = model.slim_optimize()
    return 0.0 if (v is None or v != v) else float(v)


def candidate_metabolites(model, gene_id: str) -> list[str]:
    """Non-currency metabolites the gene's reactions touch."""
    out: set[str] = set()
    for r in model.genes.get_by_id(gene_id).reactions:
        for m in r.metabolites:
            if not is_currency(m.id):
                out.add(m.id)
    return sorted(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--growth-floor", type=float, default=0.10,
                    help="biomass floor as a fraction of wildtype while maximising the demand")
    ap.add_argument("--condition", default="D-Glucose")
    ap.add_argument("--max-mets", type=int, default=40, help="cap candidates per gene")
    ap.add_argument("--out", default=f"wiki/fba_demand_completion_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from cobra.flux_analysis import single_gene_deletion

    model = load_model()
    biomass_id = [r.id for r in model.reactions if r.objective_coefficient != 0][0]
    conds = json.loads(PANEL_ARTIFACT.read_text(encoding="utf-8"))["conditions"]
    all_ex = tuple(conds.values())

    nf = never_fires_genes()
    print(f"model {model.id} | biomass {biomass_id} | NEVER_FIRES in the 131: {len(nf)}")
    print("splitting NO-EFFECT vs SUB-THRESHOLD over the panel ...", flush=True)
    genes, worst = split_no_effect(model, nf, conds, single_gene_deletion)
    print(f"  NO-EFFECT (worst ratio > 0.999): {len(genes)}   "
          f"SUB-THRESHOLD: {len(nf) - len(genes)}")

    medium = dict(model.medium)
    for ex in set(all_ex) | {"EX_glc__D_e"}:
        medium.pop(ex, None)
    medium[conds[a.condition]] = 10.0
    model.medium = medium

    wt = float(model.slim_optimize())
    floor = a.growth_floor * wt
    print(f"{a.condition} wt={wt:.6f} | growth floor {floor:.6f} ({a.growth_floor:.0%})")
    print(f"NO-EFFECT genes to probe: {len(genes)}\n")

    rows, n_named = [], 0
    for n, g in enumerate(genes, 1):
        name = model.genes.get_by_id(g).name
        mets = candidate_metabolites(model, g)[: a.max_mets]
        required = []
        for m in mets:
            with_g = max_demand(model, m, biomass_id, floor, None)
            if with_g <= EPS:
                continue                      # the model cannot make it anyway -- says nothing about g
            without_g = max_demand(model, m, biomass_id, floor, g)
            if without_g <= EPS:
                required.append({"metabolite": m, "max_with_gene": round(with_g, 6),
                                 "max_without_gene": round(without_g, 6)})
        if required:
            n_named += 1
        rows.append({"gene": g, "name": name, "n_candidates": len(mets),
                     "n_sole_route_metabolites": len(required),
                     "sole_route_metabolites": required})
        tag = ", ".join(r["metabolite"] for r in required[:4]) or "-"
        print(f"  [{n:2}/{len(genes)}] {g} {name:8} candidates={len(mets):3} "
              f"sole-route={len(required):3}  {tag}", flush=True)

    # ---- VERIFY. "sole route to m" only IMPLIES "essential once biomass demands m". Force the demand
    # and confirm the flip, rather than asserting the implication.
    import cobra

    print("\nverifying: force each named demand, confirm the gene flips to essential ...", flush=True)
    n_flip = 0
    for r in rows:
        if not r["sole_route_metabolites"]:
            continue
        met = r["sole_route_metabolites"][0]["metabolite"]
        with model:
            dm = cobra.Reaction("ZZ_force_demand")
            dm.lower_bound, dm.upper_bound = 0.01 * wt, 1000.0
            model.add_reactions([dm])
            dm.add_metabolites({model.metabolites.get_by_id(met): -1.0})
            w2 = model.slim_optimize()
            w2 = 0.0 if (w2 is None or w2 != w2) else float(w2)
            with model:
                model.genes.get_by_id(r["gene"]).knock_out()
                v = model.slim_optimize()
                v = 0.0 if (v is None or v != v) else float(v)
        ratio = (v / w2) if w2 > 1e-9 else None
        flipped = bool(w2 > 1e-9 and ratio is not None and ratio <= 0.01)
        n_flip += int(flipped)
        r["verification"] = {"forced_metabolite": met,
                             "wildtype_with_demand": round(w2, 6),
                             "growth_cost_of_demand": round(wt - w2, 6),
                             "ko_ratio_under_demand": None if ratio is None else round(ratio, 6),
                             "flips_to_essential": flipped}
    print(f"  flipped to essential under their own demand: {n_flip}/{n_named}")

    print(f"\n--- {n_named}/{len(genes)} no-effect genes have at least one metabolite they are the "
          f"SOLE route to ---")
    none_named = [r for r in rows if not r["sole_route_metabolites"]]
    print(f"    {len(none_named)} have NONE: {[r['name'] for r in none_named]}")

    out = {
        "record": "fba-demand-completion-probe-v1",
        "date": date.today().isoformat(),
        "model": model.id, "condition": a.condition,
        "biomass_reaction": biomass_id,
        "wildtype_growth": round(wt, 6),
        "growth_floor_fraction": a.growth_floor,
        "n_genes_probed": len(genes),
        "n_genes_with_a_sole_route_metabolite": n_named,
        "n_verified_flip_to_essential_under_demand": n_flip,
        "genes": rows,
        "scope": ("DIAGNOSTIC ONLY. Names which demand would flip which gene. Does NOT add demands to "
                  "the objective and re-score, and makes NO recall claim -- that is a new endpoint and "
                  "would need its own pre-registration."),
        "method": ("maximise a demand reaction on each non-currency metabolite the gene's reactions "
                   "touch, subject to biomass >= the growth floor, with and without the gene"),
        "caveats": [
            "Currency metabolites are excluded by an explicit list; a mis-classified entry there would "
            "hide or invent a result, so the list is pinned by test.",
            "Sole-route is computed in ONE condition (glucose by default). A route available only on "
            "another carbon source would read as sole-route here.",
            "The growth floor matters: without it, routes available only to a non-growing cell would "
            "count. 10% of wildtype is a choice, reported so it can be varied.",
            "'Sole route to m' does not prove biomass SHOULD demand m -- it identifies the candidate, "
            "which is a modelling question a human answers, not a fact the reconstruction settles.",
        ],
        "needs_fitness_browser_db": False,
    }
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
