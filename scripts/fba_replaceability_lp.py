"""Can the rest of the network do a reaction's job WITHOUT it -- in any medium at all?

WHY THIS EXISTS
`scripts/fba_orphan_protection_screen.py` proves a gene structurally uncallable when every reaction its
knockout disables is universally blocked or replaced by an EXACT stoichiometric duplicate. Exact
duplicates barely exist in iML1515 (12 transformation groups with >1 reaction), so that screen's 572-gene
answer is an explicit LOWER bound: it cannot see FUNCTIONAL redundancy, which is precisely what the iron
case turned out to be (`FE3abcpp` moves Fe(III) with ATP where `FE2tpp` symports Fe(II) -- same job,
different transformation). This script closes that declared gap.

THE TEST
For reaction r, delete r (and every other reaction the same knockout disables), CLOSE EVERY EXCHANGE, and
add a shadow reaction Z whose stoichiometry is the NEGATIVE of r's. At steady state, Z carrying flux f
forces the remaining internal network to perform r's exact transformation at rate f. Maximise f.

  f >= |r's capacity|   r is FULLY replaceable -- the internal network reproduces r's conversion at any
                        rate r itself could reach.
  0 < f < capacity      PARTIALLY replaceable: a route exists but is narrower than r.
  f == 0                r is IRREPLACEABLE by internal metabolism.

CLOSING THE EXCHANGES IS THE POINT. A bypass that draws on an exchange is medium-dependent, so it cannot
support a claim quantified over all media. Excluding them makes the test conservative in the safe
direction: routes are missed, never invented.

THE LOOP CAVEAT, MEASURED NOT ASSUMED
With every exchange shut, the only fluxes left are internal cycles -- exactly where thermodynamically
infeasible loops live. A spurious loop could make a reaction look replaceable when no real route exists.
So every FULLY-replaceable verdict is re-tested with `loopless_solution`, and the artifact reports the
loop-free count separately. An unaudited number here would be worth nothing.

Model-only -- no Fitness Browser DB.

Usage:
    uv run python scripts/fba_replaceability_lp.py                  # genes the exact-duplicate screen missed
    uv run python scripts/fba_replaceability_lp.py --max-genes 40   # quick pass
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.model import load_model  # noqa: E402

SCREEN_ARTIFACT = Path("wiki/fba_orphan_protection_2026-08-21.json")
RATIOS_ARTIFACT = Path("wiki/fba_eflux_bridge_2026-08-17_ratios.json")
EPS = 1e-7


def replaceable_rate(model, rxn_id: str, also_disabled: set[str]) -> float:
    """Max rate at which the network MINUS `also_disabled` can perform `rxn_id`'s transformation,
    using internal reactions only (every exchange closed). PURE w.r.t. the passed model (context-managed).
    """
    import cobra

    r = model.reactions.get_by_id(rxn_id)
    stoich = {m: -c for m, c in r.metabolites.items()}     # the SHADOW: undoes r's net conversion
    if not stoich:
        return 0.0
    with model:
        for ex in model.exchanges:
            ex.lower_bound = ex.upper_bound = 0.0
        for rid in also_disabled | {rxn_id}:
            x = model.reactions.get_by_id(rid)
            x.lower_bound = x.upper_bound = 0.0
        z = cobra.Reaction("ZZ_shadow_probe")
        z.lower_bound, z.upper_bound = 0.0, 1000.0
        model.add_reactions([z])
        z.add_metabolites(stoich)
        model.objective = z
        v = model.slim_optimize()
    return 0.0 if (v is None or v != v) else float(v)


def loop_free_rate(model, rxn_id: str, also_disabled: set[str]) -> float:
    """Same probe, but the solution must be loop-free -- kills replaceability via a futile cycle."""
    import cobra
    from cobra.flux_analysis import loopless_solution

    r = model.reactions.get_by_id(rxn_id)
    stoich = {m: -c for m, c in r.metabolites.items()}
    if not stoich:
        return 0.0
    with model:
        for ex in model.exchanges:
            ex.lower_bound = ex.upper_bound = 0.0
        for rid in also_disabled | {rxn_id}:
            x = model.reactions.get_by_id(rid)
            x.lower_bound = x.upper_bound = 0.0
        z = cobra.Reaction("ZZ_shadow_probe")
        z.lower_bound, z.upper_bound = 0.0, 1000.0
        model.add_reactions([z])
        z.add_metabolites(stoich)
        model.objective = z
        try:
            sol = loopless_solution(model)
        except Exception:
            return float("nan")
        v = sol.fluxes.get("ZZ_shadow_probe", 0.0)
    return 0.0 if (v != v) else float(v)


def capacity(model, rxn_id: str) -> float:
    r = model.reactions.get_by_id(rxn_id)
    return max(abs(r.lower_bound), abs(r.upper_bound))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-genes", type=int, default=0, help="0 = all CALLABLE genes")
    ap.add_argument("--no-loopless", action="store_true")
    ap.add_argument("--out", default=f"wiki/fba_replaceability_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from scripts.fba_orphan_protection_screen import gpr_disabled_reactions, load_blocked

    model = load_model()
    scr = json.loads(SCREEN_ARTIFACT.read_text(encoding="utf-8"))
    uncallable = set(scr["structurally_uncallable"]["genes"])
    blocked, _ = load_blocked(model)

    # Only CALLABLE genes can move: the uncallable ones are already proved.
    targets = [g.id for g in model.genes if g.id not in uncallable]
    if a.max_genes:
        targets = targets[: a.max_genes]
    print(f"model {model.id}: probing {len(targets)} CALLABLE genes "
          f"(the {len(uncallable)} uncallable ones are already proved)")

    cache: dict[tuple, tuple[float, float]] = {}
    newly: list[dict] = []
    partial = 0
    for n, gid in enumerate(targets, 1):
        disabled = [r for r in gpr_disabled_reactions(model, gid) if r not in blocked]
        if not disabled:
            continue
        key = tuple(sorted(disabled))
        if key not in cache:
            others = set(disabled)
            rates = [(replaceable_rate(model, r, others - {r}), capacity(model, r)) for r in disabled]
            worst = min((f / c if c else 1.0) for f, c in rates)
            anyrate = min(f for f, _ in rates)
            cache[key] = (worst, anyrate)
        worst, anyrate = cache[key]
        if worst >= 1.0 - EPS:
            newly.append({"gene": gid, "name": model.genes.get_by_id(gid).name,
                          "n_disabled": len(disabled), "reactions": sorted(disabled),
                          "coverage": round(worst, 6)})
        elif anyrate > EPS:
            partial += 1
        if n % 150 == 0:
            print(f"  {n}/{len(targets)}  fully-replaceable so far: {len(newly)}", flush=True)

    print(f"\nfully replaceable by internal metabolism (LP): {len(newly)} genes")
    print(f"partially replaceable (route exists but narrower): {partial} genes")

    # ---- loop audit: a futile cycle can fake replaceability with every exchange shut.
    audited = []
    if newly and not a.no_loopless:
        print(f"\nloop audit on all {len(newly)} candidates ...", flush=True)
        for i, rec in enumerate(newly, 1):
            others = set(rec["reactions"])
            ok = True
            for r in rec["reactions"]:
                lf = loop_free_rate(model, r, others - {r})
                if lf != lf or lf < capacity(model, r) - EPS:
                    ok = False
                    break
            rec["loop_free"] = ok
            if ok:
                audited.append(rec["gene"])
            if i % 20 == 0:
                print(f"  {i}/{len(newly)}  surviving: {len(audited)}", flush=True)
        print(f"  survive the loop audit: {len(audited)}/{len(newly)}")

    out = {
        "record": "fba-replaceability-lp-v1",
        "date": date.today().isoformat(),
        "model": model.id,
        "n_callable_probed": len(targets),
        "n_uncallable_already_proved": len(uncallable),
        "fully_replaceable": {"n": len(newly), "genes": [r["gene"] for r in newly], "detail": newly},
        "n_partially_replaceable": partial,
        "loop_audit": {"run": bool(newly) and not a.no_loopless,
                       "n_surviving": len(audited), "genes": sorted(audited)},
        "method": ("shadow-reaction LP with every exchange closed: max rate at which internal reactions "
                   "alone reproduce the target reaction's net conversion"),
        "caveats": [
            "Exchanges are closed so the claim is medium-independent; a bypass that needs an exchange is "
            "deliberately not counted. Conservative -- routes are missed, never invented.",
            "Full replaceability is required at the reaction's FULL capacity bound, which for a 1000-bound "
            "reaction is far more than any medium would ever demand. Conservative again.",
            "Only the loop-audited subset is trustworthy: with exchanges shut, futile cycles can fake a "
            "route. Unaudited counts are reported separately and should not be quoted.",
            "Replaceability of each disabled reaction is tested with the others also removed, so a gene's "
            "reactions cannot rescue each other.",
        ],
        "needs_fitness_browser_db": False,
    }

    if RATIOS_ARTIFACT.exists() and audited:
        d = json.loads(RATIOS_ARTIFACT.read_text(encoding="utf-8"))
        ess: set[str] = set()
        for cells in d["arms"]["baseline"].values():
            ess |= set(cells)
        hit = sorted(ess & set(audited))
        out["impact_on_experimental_deficit"] = {
            "n_essential_genes": len(ess),
            "n_newly_uncallable_and_essential": len(hit),
            "genes": hit,
        }
        print(f"\n[impact] of {len(ess)} experimentally-essential genes, {len(hit)} are newly proved "
              f"uncallable by the LP test")

    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
