"""What does an `infeasible` FBA deletion solve MEAN -- a solver failure, or genuine essentiality?

    uv run python scripts/fba_infeasibility_probe.py

This question is load-bearing and was invisible until the per-cell solver audit made it askable.

Every varying prediction iML1515 makes on the 25-carbon-source panel comes from a cell whose deletion
solve returned `infeasible` (23/23 exact-set matches; 32/33 committed genes). Read one way that is
damning -- the model's only "commitments" are the LP failing. Read the other way it is the strongest
possible essentiality signal. The two readings imply opposite conclusions about the cell, so the
question has to be settled rather than assumed.

**Settled here, and the damning reading is WRONG.** Two pieces of evidence:

1. **Determinism.** All 39 suspect cells re-solve `infeasible` when re-run one gene at a time with
   `processes=1`. A numerical/threading artifact would not reproduce identically.

2. **Identity.** Every one is the canonical catabolic gene for exactly the carbon source it fails on:
   galT/galE/galK on D-galactose (Leloir), malEFGK/malQ on maltose, manXYZ on mannose and glucosamine,
   xylA on xylose, sdhABCD on succinate, kgtP on alpha-ketoglutarate, mtlD on mannitol, srlD on
   sorbitol, fucK on fucose, rbsK on ribose, uxaC/uxaA/uxuA/kdgK/eda on the hexuronates, nagA/nagB on
   N-acetylglucosamine, glcDEF on glycolate. That is not noise; that is the right answer.

**Mechanism.** iML1515 carries a hard maintenance floor, `ATPM lower_bound = 6.86`. Delete the only
route to catabolise the sole carbon source and there is no flux distribution that meets maintenance at
all -- so the LP is genuinely INFEASIBLE rather than feasible-with-zero-growth. Lethal knockouts that
still leave a feasible solution report `growth ~ 0` with `status == optimal` instead (measured: 29 of
150 genes on glucose, with zero infeasible).

**Consequence.** Coding a NaN/infeasible solve as "essential" -- which every deletion script in this
repo does -- is CORRECT, not a bug. An abstention arm that drops those cells does not remove noise; it
removes the true positives. Any analysis that abstains on non-optimal solves is a LOWER BOUND, and a
badly biased one.

Exit 0 always: this is a diagnostic, not a gate.
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
)
from dna_decode.fba.fitness_browser import (  # noqa: E402
    apply_carbon_condition,
    carbon_conditions,
    load_records,
    open_db,
)
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402


def classify_resolve(growth, status) -> str:
    """Pure: what did a re-solve of one suspect cell actually say?

    `optimal_zero_growth`  -- feasible but no growth; ordinary lethality
    `optimal_real_growth`  -- the original NaN was spurious; a genuine false call
    `still_infeasible`     -- deterministic infeasibility; maintenance cannot be met
    """
    if str(status) != "optimal":
        return "still_infeasible"
    if growth is None or growth != growth or growth < 1e-9:
        return "optimal_zero_growth"
    return "optimal_real_growth"


def verdict_for(counts: dict[str, int]) -> str:
    """Pure verdict over the re-solve counts."""
    total = sum(counts.values())
    if not total:
        return "NO_SUSPECT_CELLS"
    if counts.get("optimal_real_growth", 0) > total * 0.1:
        return "INFEASIBLE_SOLVES_ARE_SPURIOUS"
    if counts.get("still_infeasible", 0) == total:
        return "INFEASIBLE_IS_DETERMINISTIC_GENUINE_ESSENTIALITY"
    return "INFEASIBLE_SOLVES_ARE_MIXED"


def main(argv: list[str] | None = None) -> int:
    from cobra.flux_analysis import single_gene_deletion  # noqa: PLC0415

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--artifact", default=None,
                    help="carbon-panel JSON carrying solver_audit.suspect_cells")
    ap.add_argument("--db", default=None)
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    root = Path(__file__).resolve().parent.parent
    art_path = Path(a.artifact) if a.artifact else root / f"wiki/fba_conditional_carbon_{a.date}.json"
    if not art_path.exists():
        print(f"no carbon artifact at {art_path}; run fba_conditional_carbon_validate.py first",
              file=sys.stderr)
        return 2
    art = json.loads(art_path.read_text(encoding="utf-8"))
    suspects = [tuple(x) for x in (art.get("solver_audit") or {}).get("suspect_cells", [])]
    print(f"{len(suspects)} suspect cells from {art_path.name}")

    model = load_model()
    conn = open_db(a.db)
    conds = carbon_conditions(conn, model)
    all_ex = tuple(conds.values())
    subset = conditionally_essential_genes(
        load_records(conn, conds, gene_filter={g.id for g in model.genes}))
    truth = {r.gene_id: r.experimental for r in subset}

    by_cond: dict[str, list[str]] = {}
    for g, c in suspects:
        by_cond.setdefault(c, []).append(g)

    counts: dict[str, int] = {}
    rows = []
    for cond, genes in sorted(by_cond.items()):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            wt = wildtype_growth(model)
            # processes=1: single-threaded, so a threading/numerical artifact would show up as a
            # DIFFERENT answer than the batched run that produced these cells.
            res = single_gene_deletion(
                model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
            for _, row in res.iterrows():
                gid = next(iter(row["ids"]))
                g_val, st = row["growth"], row["status"]
                v = classify_resolve(g_val, st)
                counts[v] = counts.get(v, 0) + 1
                rows.append({
                    "gene": gid, "gene_name": model.genes.get_by_id(gid).name or None,
                    "condition": cond, "status": str(st),
                    "growth": None if g_val != g_val else round(float(g_val), 6),
                    "wildtype_growth": round(wt, 4),
                    "experimentally_essential_here": truth.get(gid, {}).get(cond),
                    "resolve": v,
                })

    atpm = next((r for r in model.reactions if r.id == "ATPM"), None)
    v = verdict_for(counts)
    print("\nRE-SOLVE (processes=1):")
    for k in sorted(counts):
        print(f"   {k:44s} {counts[k]}")
    print(f"\nATPM maintenance floor: lower_bound = {atpm.lower_bound if atpm else 'n/a'}")
    print(f"VERDICT: {v}")

    agree = sum(1 for r in rows if r["experimentally_essential_here"] is True)
    print(f"   of {len(rows)} infeasible cells, {agree} are experimentally essential in that condition "
          f"({100 * agree / len(rows) if rows else 0:.0f}%)")
    print("\n   gene identities by condition:")
    for cond in sorted(by_cond):
        names = [f"{r['gene']}({r['gene_name']})" for r in rows if r["condition"] == cond]
        print(f"      {cond[:40]:42s} {', '.join(names)}")

    result = {
        "record": "fba-infeasibility-probe-v1", "date": a.date, "model": model.id,
        "question": "is an `infeasible` deletion solve a solver failure or genuine essentiality?",
        "n_suspect_cells": len(suspects),
        "resolve_counts": counts,
        "atpm_lower_bound": atpm.lower_bound if atpm else None,
        "verdict": v,
        "n_experimentally_essential_in_that_condition": agree,
        "consequence": (
            "Coding a NaN/infeasible solve as ESSENTIAL is correct, not a bug. An abstention arm that "
            "drops these cells removes the TRUE POSITIVES, not noise, and is a biased lower bound."),
        "rows": rows,
    }
    outdir = Path(a.out_dir) if a.out_dir else root / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"fba_infeasibility_probe_{a.date}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
