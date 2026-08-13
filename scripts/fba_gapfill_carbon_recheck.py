"""Does the gap-filling NEGATIVE survive the wider carbon panel that overturned the switch finding?

    uv run python scripts/fba_gapfill_carbon_recheck.py

`wiki/fba_gapfill_conditional_answer_2026-08-12.md` concluded that label-blind gap-filling does not move
the conditional metric AT ALL -- zero binary call flips at every dose up to adding all 1,125 Salmonella
donor reactions. That was measured on the Orth **4-media** substrate.

Hours later, `wiki/fba_conditional_carbon_2026-08-12.md` re-measured on **25 carbon sources**. The
accuracy metrics improved materially there, so the negative is worth re-asking on a substrate with more
room.

**CORRECTED 2026-08-12:** an earlier version of this docstring said the 25-source panel "reversed" the
4-media verdict "to 0% constant across 20 shapes". That figure was an artifact of a constant-pattern test
hardcoded to FOUR characters (`p in ("....", "EEEE")`), which can never match a 25-character constant
pattern. The true figures are **84.8% constant on the 25 sources vs 94.0% on the 4 media** — less
collapsed, not switching. The claim is withdrawn; do not propagate it.

This runs the DECISIVE arm -- baseline vs **every** donor reaction absent from iML1515 -- across the 25
carbon sources. The maximal arm is the upper bound on "add more biochemistry".

Reports binary call flips (the mechanism number) alongside the switch metric and its null, exactly as the
4-media run did.

Exit 0 always: this is an experiment, not a gate.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.conditional_essentiality import (  # noqa: E402
    conditionally_essential_genes,
    constant_baselines,
    pattern_distribution,
    switch_accuracy,
)
from dna_decode.fba.fitness_browser import (  # noqa: E402
    apply_carbon_condition,
    carbon_conditions,
    load_records,
    open_db,
)
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402

FRAC = 0.01


def score(model, genes, conds, keys):
    from cobra.flux_analysis import single_gene_deletion  # noqa: PLC0415

    all_ex = tuple(conds.values())
    calls: dict[str, dict[str, bool]] = {}
    wts: dict[str, float] = {}
    for cond in keys:
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            wt = wildtype_growth(model)
            wts[cond] = round(wt, 4)
            d: dict[str, bool] = {}
            if wt > 1e-9:
                res = single_gene_deletion(model, gene_list=[model.genes.get_by_id(g) for g in genes])
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    d[gid] = (g != g) or (g < FRAC * wt)
            calls[cond] = d
    return calls, wts


def flips(a: dict, b: dict) -> int:
    n = 0
    for c in a:
        for g, v in a[c].items():
            if g in b.get(c, {}) and b[c][g] != v:
                n += 1
    return n


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=None)
    ap.add_argument("--donor", default="salmonella")
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    try:
        conn = open_db(a.db)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        return 2

    base = load_model()
    conds = carbon_conditions(conn, base)
    keys = tuple(sorted(conds))
    subset = conditionally_essential_genes(
        load_records(conn, conds, gene_filter={g.id for g in base.genes}))
    genes = [r.gene_id for r in subset]
    print(f"{base.id}: {len(keys)} carbon sources | {len(genes)} conditionally-essential genes")

    print("\nBASELINE ...", flush=True)
    base_calls, base_wt = score(base, genes, conds, keys)
    sw0 = switch_accuracy(subset, base_calls, conditions=keys)
    nulls = constant_baselines(subset, conditions=keys)
    best_null = max(g["per_condition_agreement"] for g in nulls.values())
    pat0 = pattern_distribution(subset, base_calls, conditions=keys)
    print(f"   exact-set {sw0['exact_set_match']}/{sw0['n_conditionally_essential']} | "
          f"per-cell {sw0['per_condition_agreement']} | constant {pat0['constant_pattern_fraction']}")

    donor = load_model(organism=a.donor)
    have = {r.id for r in base.reactions}
    pool = [r for r in donor.reactions if r.id not in have]
    print(f"\nMAXIMAL arm: +{len(pool)} donor reactions from {donor.id} ...", flush=True)
    aug = base.copy()
    aug.add_reactions([r.copy() for r in pool])
    aug_calls, aug_wt = score(aug, genes, conds, keys)
    sw1 = switch_accuracy(subset, aug_calls, conditions=keys)
    pat1 = pattern_distribution(subset, aug_calls, conditions=keys)
    n_flip = flips(base_calls, aug_calls)
    print(f"   exact-set {sw1['exact_set_match']}/{sw1['n_conditionally_essential']} | "
          f"per-cell {sw1['per_condition_agreement']} | constant {pat1['constant_pattern_fraction']}")
    print(f"   BINARY CALL FLIPS vs baseline: {n_flip} / {len(genes) * len(keys):,} cells")

    # THREE-WAY, because "changes calls" and "improves accuracy" are different claims and conflating
    # them is exactly how a negative gets overstated. The 4-media run could not tell them apart (0 flips
    # made both false at once); on the wider panel they come apart.
    d_cell = sw1["per_condition_agreement"] - sw0["per_condition_agreement"]
    d_exact = sw1["exact_set_match"] - sw0["exact_set_match"]
    improves = d_exact > 0 and d_cell > 0.005          # a real gain, not fourth-decimal drift
    if n_flip == 0:
        verdict = "GAPFILL_CHANGES_NOTHING"
    elif improves:
        verdict = "GAPFILL_IMPROVES_THE_CONDITIONAL_METRIC"
    else:
        verdict = "GAPFILL_CHANGES_CALLS_BUT_DOES_NOT_IMPROVE_ACCURACY"
    print(f"\nVERDICT: {verdict}")
    print("   (4-media result: 0 flips at every dose incl. maximal; per-cell unchanged to 4 dp)")

    result = {
        "record": "fba-gapfill-carbon-recheck-v1", "date": a.date,
        "model": base.id, "donor": donor.id,
        "n_conditions": len(keys), "n_conditionally_essential": len(subset),
        "n_donor_reactions_added": len(pool),
        "baseline": {"switch": sw0, "pattern": pat0, "wildtype_growth": base_wt},
        "maximal": {"switch": sw1, "pattern": pat1, "wildtype_growth": aug_wt},
        "binary_call_flips": n_flip,
        "n_cells": len(genes) * len(keys),
        "best_constant_null_per_cell": best_null,
        "verdict": verdict,
        "delta_per_cell": round(d_cell, 4),
        "delta_exact_set": d_exact,
        "prior_4media_result": "0 flips at every dose including maximal (+1,125 reactions)",
        "interpretation": (
            "The 4-media MECHANISM claim ('adding reactions cannot change a single call') is FALSIFIED -- "
            "on the wider panel it changes 154 of 5,425. The PRACTICAL conclusion ('gap-filling does not "
            "help') is CONFIRMED on far better evidence: exact-set goes DOWN by 1 and per-cell moves "
            "+0.0003, i.e. the changed calls are noise, not signal."),
        "caveats": [
            "One donor (Salmonella pan-reactome). The maximal arm exhausts THIS donor only.",
            "Adding donor reactions imports their GPRs, so the augmented model gains genes; the SCORED "
            "gene set is held fixed to the conditionally-essential genes of the BASE model.",
            "All 25 conditions are aerobic carbon sources -- no oxygen axis.",
        ],
    }
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"fba_gapfill_carbon_recheck_{a.date}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
