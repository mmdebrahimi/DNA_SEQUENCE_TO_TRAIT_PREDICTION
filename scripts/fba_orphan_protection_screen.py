"""Which genes can FBA NEVER call essential -- proved from the reconstruction, before any label.

MOTIVATION
`wiki/fba_orphan_redundancy_2026-08-21.md` measured that iron acquisition in iML1515 is protected by
route redundancy, and that one of the parallel routes (`FE3abcpp`) has an EMPTY GPR -- no gene deletion
of any kind can close it. That memo left an honest gap: 113 orphan reactions is a COUNT, not an impact
estimate. This script closes it by turning the observation into a screen over every gene.

THE THEOREM (this is what makes the output a declaration rather than a guess)
Knock out gene g and let D(g) be the reactions that ACTUALLY become non-functional (full-GPR aware, so
an isozyme `or` keeps its reaction alive). If EVERY r in D(g) satisfies one of

  (B) UNIVERSALLY BLOCKED  r carries no flux even with every exchange open. Any medium is a SUBSET of
      all-open, so tighter bounds keep it blocked: r is blocked in every medium. Deleting a reaction
      that can only ever carry zero removes nothing from the feasible set.

  (C) DUPLICATE-RESCUED    some reaction r' survives the knockout, is not itself blocked, performs the
      SAME transformation (identical stoichiometry up to a scalar, reversal included), and has bounds at
      least as permissive as r in canonical units. Every flux through r can be rerouted through r'.

then the knockout leaves the feasible flux set over the retained reactions UNCHANGED. Hence the optimum
is unchanged -- for ANY objective over the retained reactions, in ANY medium -- so the growth ratio is
exactly 1 and the model can never call g essential. D(g) empty (case A, fully isozyme-buffered) is the
degenerate instance of the same statement.

WHAT THE CLAIM IS NOT
- Not "the gene is dispensable in the cell." It is a statement about the MODEL's reach, so a gene in
  this set that is experimentally essential is a declared blind spot, not a prediction.
- Not closed under changing INTERNAL bounds. Medium here means exchange bounds, which is the axis every
  condition in this project varies.
- Objectives are assumed linear over the retained reactions (biomass is). An objective with a coefficient
  ON a deleted reaction is outside the claim.

RESCUER PERMANENCE
A rescuer with an EMPTY GPR is PERMANENT: no deletion anywhere, of any gene, in any combination, can
remove it. A gene-carried rescuer only protects against SINGLE deletions. The screen reports both, so
"orphan-protected" is separated from "isozyme-protected" rather than blurred together.

VERIFICATION IS PART OF THE RUN
The theorem is checked against the thing it predicts: a full single-gene deletion is run over every gene
in each verification medium, and ANY structurally-uncallable gene with growth ratio < 1 falsifies the
screen. That check is reported, not assumed.

Model-only -- no Fitness Browser DB.

Usage:
    uv run python scripts/fba_orphan_protection_screen.py
    uv run python scripts/fba_orphan_protection_screen.py --no-verify        # skip the deletion check
    uv run python scripts/fba_orphan_protection_screen.py --verify-conditions D-Glucose,Acetate
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.model import load_model  # noqa: E402

PANEL_ARTIFACT = Path("wiki/fba_conditional_carbon_2026-08-13.json")
BLINDSPOT_ARTIFACT = Path("wiki/fba_structural_blindspot_2026-08-21.json")
#: committed experimental essentiality calls, so the impact cross-tab needs no Fitness Browser DB
RATIOS_ARTIFACT = Path("wiki/fba_eflux_bridge_2026-08-17_ratios.json")

COEF_ROUND = 9
BOUND_EPS = 1e-9
#: a deletion ratio below this falsifies the "can never be called essential" claim
RATIO_TOL = 1e-6


# --------------------------------------------------------------------------------------- pure helpers

def canon_key(stoich: dict[str, float]) -> tuple[tuple[str, float], ...] | None:
    """Canonical identity of a transformation: stoichiometry divided by the coefficient of its
    lexicographically-smallest metabolite.

    Dividing by a SIGNED reference coefficient normalises direction as well as scale, so `A -> B` and
    `B -> A` collapse to the same key (their difference is carried by the canonical interval instead).
    Returns None for an empty stoichiometry (nothing to compare).
    """
    if not stoich:
        return None
    ref = min(stoich)
    c = stoich[ref]
    if abs(c) <= BOUND_EPS:
        return None
    return tuple(sorted((m, round(v / c, COEF_ROUND)) for m, v in stoich.items()))


def canon_interval(stoich: dict[str, float], lo: float, hi: float) -> tuple[float, float] | None:
    """The interval of CANONICAL rates this reaction can realise.

    With canonical rate u = v * c (v the reaction's own flux, c its reference coefficient), the image of
    [lo, hi] is [lo*c, hi*c] -- whose endpoints swap when c < 0. Returning it sorted keeps the caller's
    containment test direction-agnostic.
    """
    if not stoich:
        return None
    c = stoich[min(stoich)]
    if abs(c) <= BOUND_EPS:
        return None
    a, b = lo * c, hi * c
    return (min(a, b), max(a, b))


def covers(outer: tuple[float, float], inner: tuple[float, float]) -> bool:
    """Does `outer` contain `inner` (canonical units)? Sufficient for rerouting, deliberately strict."""
    return outer[0] <= inner[0] + BOUND_EPS and outer[1] >= inner[1] - BOUND_EPS


def classify_gene(disabled: list[str],
                  blocked: set[str],
                  rescued: dict[str, list[str]]) -> tuple[str, bool]:
    """(class, is_uncallable) for one gene. PURE -- all model access is done by the caller.

    `rescued` maps a disabled reaction to its surviving duplicate rescuers (may be empty).
    Classes are ordered most-specific-first so the label names the REASON, not just the verdict.
    """
    if not disabled:
        return "NO_KO_EFFECT", True
    unexplained = [r for r in disabled if r not in blocked and not rescued.get(r)]
    if unexplained:
        return "CALLABLE", False
    if all(r in blocked for r in disabled):
        return "ALL_DISABLED_BLOCKED", True
    if all(rescued.get(r) for r in disabled):
        return "DUPLICATE_RESCUED", True
    return "MIXED_BLOCKED_AND_RESCUED", True


# ----------------------------------------------------------------------------------- model-facing bits

def gpr_disabled_reactions(model, gene_id: str) -> list[str]:
    """Reactions that ACTUALLY become non-functional when this gene is knocked out.

    Same contract as `scripts/fba_fva_requirement_class.py` -- cobrapy zeroes a reaction only when the
    FULL GPR evaluates false, so a reaction behind an isozyme `or` survives a single-gene deletion.
    """
    with model:
        model.genes.get_by_id(gene_id).knock_out()
        return [r.id for r in model.genes.get_by_id(gene_id).reactions if not r.functional]


def build_duplicate_index(model) -> dict:
    """canonical transformation key -> [reaction ids] performing it."""
    idx: dict[tuple, list[str]] = {}
    for r in model.reactions:
        k = canon_key({m.id: c for m, c in r.metabolites.items()})
        if k is not None:
            idx.setdefault(k, []).append(r.id)
    return idx


def find_rescuers(model, rxn_id: str, disabled: set[str], blocked: set[str], dup_index: dict) -> list[str]:
    """Surviving reactions that can carry everything `rxn_id` could.

    A rescuer must (a) perform the same transformation, (b) survive THIS knockout, (c) not be blocked --
    a dead reaction rescues nothing -- and (d) have a canonical interval covering the target's.
    """
    r = model.reactions.get_by_id(rxn_id)
    st = {m.id: c for m, c in r.metabolites.items()}
    key, want = canon_key(st), canon_interval(st, r.lower_bound, r.upper_bound)
    if key is None or want is None:
        return []
    out = []
    for cand in dup_index.get(key, ()):
        if cand == rxn_id or cand in disabled or cand in blocked:
            continue
        c = model.reactions.get_by_id(cand)
        have = canon_interval({m.id: v for m, v in c.metabolites.items()},
                              c.lower_bound, c.upper_bound)
        if have is not None and covers(have, want):
            out.append(cand)
    return out


def masking_genes(model, gene_id: str) -> set[str]:
    """Every gene sharing a reaction with this one -- the set whose presence can mask its deletion.

    Deleting this SET and re-optimising is the honest unmask test. Zeroing the gene's reactions instead
    is NOT equivalent: for a gene mapped to many reactions (ompC covers 285 diffusion reactions) that
    removes far more capability than the redundancy under test, and it silently mis-scores.
    """
    return {x.id for r in model.genes.get_by_id(gene_id).reactions for x in r.genes}


def unmask_ratio(model, gene_id: str, wt: float) -> float:
    """Growth ratio after deleting the gene AND everything that shares its reactions."""
    with model:
        for gid in masking_genes(model, gene_id):
            model.genes.get_by_id(gid).knock_out()
        v = model.slim_optimize()
    if v is None or v != v:
        return 0.0
    return float(v) / wt


def apply_carbon(model, exchange: str, all_carbon: tuple[str, ...]) -> None:
    """Sole-carbon medium (mirrors `fba_structural_blindspot.apply_carbon`)."""
    medium = dict(model.medium)
    for ex in set(all_carbon) | {"EX_glc__D_e"}:
        medium.pop(ex, None)
    medium[exchange] = 10.0
    model.medium = medium


def load_blocked(model) -> tuple[set[str], str]:
    """Universally-blocked reactions: reuse the committed artifact when present, else recompute."""
    if BLINDSPOT_ARTIFACT.exists():
        d = json.loads(BLINDSPOT_ARTIFACT.read_text(encoding="utf-8"))
        rxns = d.get("blocked", {}).get("blocked_reactions")
        if rxns:
            return set(rxns), f"committed artifact {BLINDSPOT_ARTIFACT.name}"
    from cobra.flux_analysis import find_blocked_reactions
    with model:
        for r in model.exchanges:
            r.lower_bound = min(r.lower_bound, -1000.0)
            r.upper_bound = max(r.upper_bound, 1000.0)
        return set(find_blocked_reactions(model)), "recomputed (all exchanges open)"


# --------------------------------------------------------------------------------------------- driver

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-verify", action="store_true",
                    help="skip the single-gene-deletion falsification check")
    ap.add_argument("--verify-conditions", default="D-Glucose,Potassium acetate",
                    help="comma-separated conditions to run the falsification check in")
    ap.add_argument("--no-unmask", action="store_true",
                    help="skip the isozyme-unmask test on the buffered genes")
    ap.add_argument("--unmask-condition", default="D-Glucose")
    ap.add_argument("--out", default=f"wiki/fba_orphan_protection_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    model = load_model()
    n_gene, n_rxn = len(model.genes), len(model.reactions)
    print(f"model {model.id}: {n_rxn} reactions, {n_gene} genes")

    orphans = {r.id for r in model.reactions if not r.gene_reaction_rule.strip()}
    non_ex_orphans = orphans - {r.id for r in model.exchanges}
    print(f"orphan (empty-GPR) reactions: {len(orphans)} total, {len(non_ex_orphans)} non-exchange")

    blocked, blocked_src = load_blocked(model)
    print(f"universally-blocked reactions: {len(blocked)}  [{blocked_src}]")

    dup_index = build_duplicate_index(model)
    n_dup_groups = sum(1 for v in dup_index.values() if len(v) > 1)
    print(f"transformation groups with >1 reaction: {n_dup_groups}")

    print("\nscreening every gene ...", flush=True)
    classes: dict[str, str] = {}
    detail: dict[str, dict] = {}
    for n, g in enumerate(model.genes, 1):
        disabled = gpr_disabled_reactions(model, g.id)
        dset = set(disabled)
        rescued = {r: find_rescuers(model, r, dset, blocked, dup_index)
                   for r in disabled if r not in blocked}
        cls, uncallable = classify_gene(disabled, blocked, rescued)
        classes[g.id] = cls
        if uncallable and cls != "NO_KO_EFFECT":
            all_resc = sorted({x for v in rescued.values() for x in v})
            detail[g.id] = {
                "name": g.name, "class": cls,
                "n_disabled": len(disabled),
                "n_blocked": sum(1 for r in disabled if r in blocked),
                "rescuers": all_resc,
                "orphan_rescuers": [x for x in all_resc if x in orphans],
                "permanently_protected": bool(all_resc) and all(x in orphans for x in all_resc),
            }
        if n % 300 == 0:
            print(f"  {n}/{n_gene} ...", flush=True)

    counts: dict[str, int] = {}
    for c in classes.values():
        counts[c] = counts.get(c, 0) + 1
    uncallable = sorted(g for g, c in classes.items() if c != "CALLABLE")

    print("\n--- screen ---")
    for c in sorted(counts, key=lambda k: -counts[k]):
        print(f"  {c:28} {counts[c]:5}  ({counts[c]/n_gene:6.1%})")
    print(f"  {'STRUCTURALLY UNCALLABLE':28} {len(uncallable):5}  ({len(uncallable)/n_gene:6.1%})")

    perm = sorted(g for g, d in detail.items() if d["permanently_protected"])
    orphan_touched = sorted(g for g, d in detail.items() if d["orphan_rescuers"])
    print(f"\n  rescued genes whose rescuers are ALL orphan (permanent): {len(perm)}")
    print(f"  rescued genes with at least one orphan rescuer          : {len(orphan_touched)}")

    out = {
        "record": "fba-orphan-protection-screen-v1",
        "date": date.today().isoformat(),
        "model": model.id,
        "n_genes": n_gene,
        "n_reactions": n_rxn,
        "orphan_reactions": {"n_total": len(orphans), "n_non_exchange": len(non_ex_orphans)},
        "blocked_source": blocked_src,
        "n_blocked_reactions": len(blocked),
        "class_counts": counts,
        "structurally_uncallable": {
            "n": len(uncallable), "fraction": round(len(uncallable) / n_gene, 4),
            "genes": uncallable,
        },
        "orphan_protection": {
            "n_genes_all_rescuers_orphan": len(perm),
            "n_genes_any_orphan_rescuer": len(orphan_touched),
            "genes_permanently_protected": perm,
        },
        "detail": detail,
        "theorem": (
            "If every reaction a knockout disables is either universally blocked (no flux with all "
            "exchanges open, hence none in any subset medium) or duplicate-rescued by a surviving "
            "reaction with a covering canonical interval, the knockout leaves the feasible set over the "
            "retained reactions unchanged -- so growth is unchanged for any objective in any medium and "
            "the model can never call the gene essential."
        ),
        "needs_fitness_browser_db": False,
    }

    # ---- falsification: the theorem predicts ratio == 1 for every uncallable gene, everywhere.
    if not a.no_verify:
        from cobra.flux_analysis import single_gene_deletion
        conds = json.loads(PANEL_ARTIFACT.read_text(encoding="utf-8"))["conditions"]
        all_ex = tuple(conds.values())
        want = [c.strip() for c in a.verify_conditions.split(",") if c.strip()]
        checks = {}
        for cond in want:
            if cond not in conds:
                print(f"\n  [verify] unknown condition {cond!r}, skipping")
                continue
            print(f"\n[verify] single-gene deletion, all {n_gene} genes, {cond} ...", flush=True)
            with model:
                apply_carbon(model, conds[cond], all_ex)
                wt = float(model.slim_optimize())
                if not (wt > 1e-9) or wt != wt:
                    print(f"  wildtype infeasible on {cond}; skipping")
                    continue
                res = single_gene_deletion(model, gene_list=list(model.genes), processes=1)
            ratios = {}
            for _, row in res.iterrows():
                gid = next(iter(row["ids"]))
                gr = row["growth"]
                ratios[gid] = 0.0 if gr != gr else float(gr) / wt
            viol = sorted((g, round(ratios[g], 6)) for g in uncallable
                          if ratios.get(g, 1.0) < 1.0 - RATIO_TOL)
            # NON-VACUITY: a screen that admits everything would also show zero violations. The CALLABLE
            # complement must actually contain genes the model DOES call essential, or the split is idle.
            callable_genes = [g for g, c in classes.items() if c == "CALLABLE"]
            hit = [g for g in callable_genes if ratios.get(g, 1.0) <= 0.01]
            red = [g for g in callable_genes if ratios.get(g, 1.0) < 1.0 - RATIO_TOL]
            checks[cond] = {
                "wildtype_growth": round(wt, 6),
                "n_violations": len(viol),
                "violations": viol[:25],
                "verdict": "THEOREM_HELD" if not viol else "THEOREM_FALSIFIED",
                "non_vacuity": {
                    "n_callable": len(callable_genes),
                    "n_callable_called_essential": len(hit),
                    "n_callable_growth_reduced": len(red),
                    "verdict": "SCREEN_IS_SHARP" if hit else "SCREEN_IS_VACUOUS_HERE",
                },
            }
            print(f"  wildtype {wt:.6f} | uncallable genes with ratio < 1: {len(viol)} "
                  f"-> {checks[cond]['verdict']}")
            print(f"  non-vacuity: of {len(callable_genes)} CALLABLE genes, {len(hit)} are called "
                  f"essential and {len(red)} growth-reducing -> {checks[cond]['non_vacuity']['verdict']}")
            for g, r in viol[:10]:
                print(f"    VIOLATION {g} ({classes[g]}) ratio={r}")
        out["falsification"] = checks
        out["falsification_verdict"] = (
            "THEOREM_HELD" if checks and all(v["verdict"] == "THEOREM_HELD" for v in checks.values())
            else ("NOT_RUN" if not checks else "THEOREM_FALSIFIED"))
        print(f"\noverall falsification verdict: {out['falsification_verdict']}")

    # ---- unmask: is a NO_KO_EFFECT gene's FUNCTION essential once the redundancy is removed?
    #      This separates "the model carries an isozyme the strain does not" (curable by GPR curation)
    #      from "the function is dispensable in the model regardless" (an objective/network problem).
    if not a.no_unmask:
        cond = a.unmask_condition
        conds = json.loads(PANEL_ARTIFACT.read_text(encoding="utf-8"))["conditions"]
        if cond not in conds:
            print(f"\n[unmask] unknown condition {cond!r}, skipping")
        else:
            buffered = [g for g, c in classes.items() if c == "NO_KO_EFFECT"]
            print(f"\n[unmask] testing {len(buffered)} isozyme-buffered genes in {cond} ...", flush=True)
            with model:
                apply_carbon(model, conds[cond], tuple(conds.values()))
                wt = float(model.slim_optimize())
                masked = [g for g in buffered if unmask_ratio(model, g, wt) <= 0.01]
            out["unmask"] = {
                "condition": cond, "n_tested": len(buffered),
                "n_isozyme_masked": len(masked), "genes_isozyme_masked": sorted(masked),
                "meaning": ("masked = the FUNCTION is essential and only GPR redundancy hides it, so "
                            "strain-aware GPR curation could recover it; not-masked = the function is "
                            "dispensable in the model regardless, which GPR work will not touch"),
            }
            print(f"  isozyme-masked: {len(masked)}/{len(buffered)}")

    # ---- impact: how much of the EXPERIMENTAL deficit lands inside the declared blind spot
    if RATIOS_ARTIFACT.exists():
        d = json.loads(RATIOS_ARTIFACT.read_text(encoding="utf-8"))
        ess: set[str] = set()
        for cells in d["arms"]["baseline"].values():
            ess |= set(cells)
        inter = sorted(ess & set(uncallable))
        by_cls: dict[str, int] = {}
        for g in inter:
            by_cls[classes[g]] = by_cls.get(classes[g], 0) + 1
        masked_ess = sorted(set(inter) & set(out.get("unmask", {}).get("genes_isozyme_masked", [])))
        out["impact_on_experimental_deficit"] = {
            "source": RATIOS_ARTIFACT.name,
            "n_essential_genes": len(ess),
            "n_essential_structurally_uncallable": len(inter),
            "fraction": round(len(inter) / len(ess), 4) if ess else None,
            "by_class": by_cls,
            "genes": inter,
            "n_uncallable_essential_isozyme_masked": len(masked_ess),
            "genes_isozyme_masked": masked_ess,
            "note": ("a hard floor on the recall of ANY constraint-based method on this gene set: no "
                     "medium, objective or constraint layer can move these cells"),
        }
        print(f"\n[impact] {len(inter)}/{len(ess)} experimentally-essential genes are structurally "
              f"uncallable ({len(inter)/len(ess):.1%})  {by_cls}")
        print(f"         of those, isozyme-masked (GPR-curable): {len(masked_ess)}")

    out["caveats"] = [
        "Model reach, not biology: an experimentally-essential gene in this set is a DECLARED blind "
        "spot, not a prediction that it is dispensable.",
        "Medium means exchange bounds. Changing INTERNAL bounds is outside the claim.",
        "Objectives are assumed linear over the retained reactions; an objective weighting a deleted "
        "reaction is outside the claim.",
        "Duplicate rescue requires the rescuer's canonical interval to COVER the target's, which is "
        "sufficient but conservative -- genes are missed, never falsely admitted.",
        "A gene-carried rescuer protects only against SINGLE deletions; only an orphan rescuer is "
        "permanent. The two are reported separately.",
    ]
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
