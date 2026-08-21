"""Bypass closure: does a gene go from dispensable to essential when its free shortcut is closed?

THE DISCRIMINATION THIS BUYS
Two live hypotheses both predict an idle gene, so neither is separable by looking at idleness:
  MEDIUM REALISM        the model is handed a nutrient for free, so the acquisition machinery is idle.
  OBJECTIVE INCOMPLETE  biomass never demands the product, so the pathway is idle.

The separating observable is a DELETION-RATIO FLIP. Close the free shortcut, keep the wildtype feasible,
then re-run the gene deletion:
  ratio ~1.0 -> <1.0  => MEDIUM REALISM implicated (the gene was only idle because of the free tap)
  ratio stays ~1.0    => OBJECTIVE INCOMPLETENESS implicated (nothing demands the product either way)

FIRST CASE: iron/zinc. `EX_fe2_e`, `EX_fe3_e`, `EX_zn2_e` all sit in the default medium at
(-1000, 1000), and all 12 named Fe/Zn acquisition genes (fepA fes fepB fepC fepD fepG tonB exbB exbD
znuA znuB znuC) carry zero flux at the optimum -- the whole siderophore system idles because iron is a
free tap, while in a real cell it is essential.

WHAT COUNTS AS A FAILED EXPERIMENT (stated before running)
If closing the shortcut makes the WILDTYPE infeasible, the test is INCONCLUSIVE for that nutrient, not a
result: it means the model has no reconstructed physiological route, so "close the tap" is starvation
rather than realism. That is reported as INCONCLUSIVE_WT_INFEASIBLE, never as a flip.

Deletions are GPR-aware by construction (cobrapy's single_gene_deletion evaluates the full GPR), which
is the same correctness point that moved 16 genes out of REQUIRED in the FVA classifier.

Model-only -- no Fitness Browser DB.

Usage:
    uv run python scripts/fba_bypass_closure.py
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

#: bypass name -> (exchanges to close, genes expected to depend on the closed route)
BYPASSES: dict[str, dict] = {
    "iron_zinc": {
        "close": ["EX_fe2_e", "EX_fe3_e", "EX_zn2_e"],
        "genes": ["fepA", "fes", "fepB", "fepC", "fepD", "fepG",
                  "tonB", "exbB", "exbD", "znuA", "znuB", "znuC"],
        "rationale": "direct metal uptake is a free tap; real acquisition needs siderophore/ABC machinery",
    },
    "iron_only": {
        "close": ["EX_fe2_e", "EX_fe3_e"],
        "genes": ["fepA", "fes", "fepB", "fepC", "fepD", "fepG", "tonB", "exbB", "exbD"],
        "rationale": "isolates iron from zinc so a joint infeasibility can be attributed",
    },
}

FLIP_EPS = 1e-6


def genes_by_name(model, names: list[str]) -> dict[str, str]:
    """{gene name -> gene id} for names present in the model."""
    lookup = {(g.name or "").lower(): g.id for g in model.genes}
    return {n: lookup[n.lower()] for n in names if n.lower() in lookup}


def deletion_ratios(model, gene_ids: list[str], single_gene_deletion) -> tuple[float, dict[str, float]]:
    """(wildtype growth, {gene_id: growth ratio}) in the model's CURRENT medium."""
    wt = float(model.slim_optimize())
    if not (wt > 1e-9) or wt != wt:
        return wt, {}
    res = single_gene_deletion(
        model, gene_list=[model.genes.get_by_id(g) for g in gene_ids], processes=1)
    out = {}
    for _, row in res.iterrows():
        gid = next(iter(row["ids"]))
        g = row["growth"]
        out[gid] = 0.0 if g != g else float(g) / wt
    return wt, out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--condition", default="D-Glucose")
    ap.add_argument("--out", default=f"wiki/fba_bypass_closure_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from cobra.flux_analysis import single_gene_deletion

    model = load_model()
    conds = json.loads(PANEL_ARTIFACT.read_text(encoding="utf-8"))["conditions"]
    carbon = conds[a.condition]
    all_ex = tuple(conds.values())

    results = {}
    for name, spec in BYPASSES.items():
        gmap = genes_by_name(model, spec["genes"])
        gids = sorted(set(gmap.values()))
        print(f"\n=== bypass '{name}' — closing {spec['close']} ===")
        print(f"    genes resolved: {len(gmap)}/{len(spec['genes'])}")
        if not gids:
            results[name] = {"verdict": "NO_GENES_RESOLVED"}
            continue

        # --- baseline: the tap is open
        with model:
            medium = dict(model.medium)
            for ex in set(all_ex) | {"EX_glc__D_e"}:
                medium.pop(ex, None)
            medium[carbon] = 10.0
            model.medium = medium
            wt_open, ratios_open = deletion_ratios(model, gids, single_gene_deletion)

        # --- closed: the tap is shut, everything else identical
        with model:
            medium = dict(model.medium)
            for ex in set(all_ex) | {"EX_glc__D_e"}:
                medium.pop(ex, None)
            medium[carbon] = 10.0
            for ex in spec["close"]:
                medium.pop(ex, None)
            model.medium = medium
            wt_closed, ratios_closed = deletion_ratios(model, gids, single_gene_deletion)

        print(f"    wildtype growth  open={wt_open:.6f}  closed={wt_closed:.6f}")

        if not (wt_closed > 1e-9) or wt_closed != wt_closed:
            print("    VERDICT: INCONCLUSIVE_WT_INFEASIBLE — closing the tap starves the model, so the "
                  "reconstruction has no physiological route. Not a flip.")
            results[name] = {
                "verdict": "INCONCLUSIVE_WT_INFEASIBLE",
                "wt_open": round(wt_open, 6), "wt_closed": None,
                "closed": spec["close"], "rationale": spec["rationale"],
                "note": "the free tap is load-bearing: without it there is no route to the nutrient at all",
            }
            continue

        flips = {g: (ratios_open.get(g, 1.0), ratios_closed.get(g, 1.0))
                 for g in gids
                 if ratios_open.get(g, 1.0) > 1 - FLIP_EPS
                 and ratios_closed.get(g, 1.0) < 1 - FLIP_EPS}
        verdict = ("MEDIUM_REALISM_IMPLICATED" if flips else "OBJECTIVE_INCOMPLETENESS_IMPLICATED")
        print(f"    genes flipping dispensable -> growth-reducing: {len(flips)}/{len(gids)}")
        for g, (o, c) in list(flips.items())[:6]:
            nm = model.genes.get_by_id(g).name
            print(f"      {g} {nm:8} {o:.4f} -> {c:.4f}")
        print(f"    VERDICT: {verdict}")
        results[name] = {
            "verdict": verdict,
            "closed": spec["close"], "rationale": spec["rationale"],
            "wt_open": round(wt_open, 6), "wt_closed": round(wt_closed, 6),
            "wt_growth_cost_of_closing": round(wt_open - wt_closed, 6),
            "n_genes": len(gids), "n_flipped": len(flips),
            "flips": {g: [round(o, 4), round(c, 4)] for g, (o, c) in flips.items()},
            "ratios_open": {g: round(v, 4) for g, v in ratios_open.items()},
            "ratios_closed": {g: round(v, 4) for g, v in ratios_closed.items()},
        }

    out = {
        "record": "fba-bypass-closure-v1",
        "date": date.today().isoformat(),
        "model": model.id,
        "condition": a.condition,
        "discrimination": {
            "flip": "MEDIUM_REALISM_IMPLICATED — the gene was idle only because of a free shortcut",
            "no_flip": "OBJECTIVE_INCOMPLETENESS_IMPLICATED — nothing demands the product either way",
            "wt_infeasible": "INCONCLUSIVE — closing the tap is starvation, not realism",
        },
        "results": results,
        "caveats": [
            "One condition. A bypass that matters on glucose may not on another carbon source.",
            "Gene sets are curated by name, not exhaustive: other genes may also depend on the route.",
            "A flip shows the gene BECOMES growth-reducing, not that it reaches the 1% essentiality "
            "cutoff. Magnitude is reported so a graded readout can use it.",
            "Model-only: this does not consult experimental labels, so it identifies a MECHANISM, not a "
            "recovered true positive.",
        ],
    }
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
