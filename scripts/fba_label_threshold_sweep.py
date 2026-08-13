"""Does the carbon-panel conditional result survive a different LABEL cutoff? And what carries it?

    uv run python scripts/fba_label_threshold_sweep.py

Two questions in one pass, because they share an expensive FBA step.

**Q1 — label sensitivity (the named open gap).** Every number in
`wiki/fba_conditional_carbon_2026-08-12.md` rests on ONE inherited cutoff: a gene is essential in a
condition iff mean RB-TnSeq fitness `< -2` (Bernstein 2023, reused for comparability). Nothing has ever
tested whether the headline moves under a different bar. This sweeps the fitness threshold AND the
per-measurement confidence bar `min_abs_t` (`GeneFitness.t`, reachable since 2026-08-13) and reports the
headline at every setting.

**Q2 — what mechanism actually carries the conditional signal?** The 2026-08-13 infeasibility finding
showed 32 of 33 commitments on the shipped panel came from `infeasible` solves (the ATPM floor unmet),
not from a growth ratio dipping below the 1% cutoff. If that holds across label settings, then the model
is effectively a BOOLEAN can-grow/cannot-grow predictor for conditional purposes, and the entire
continuous-threshold line (`continuous_readout`'s oracle, `deployable_threshold`'s retune) cannot help
this metric no matter how it is tuned. That is worth knowing before anyone tunes it again.

**Efficiency note that makes the sweep cheap:** the FBA deletion calls do NOT depend on the label
threshold — only which genes are SCORED does. So the deletions run ONCE over the union of every gene that
is conditionally essential at any setting, and each setting is then re-scored from that cache in
milliseconds.

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
from dna_decode.fba.solver_audit import audit_deletion_frame, merge_audits  # noqa: E402

FRAC = 0.01
FIT_THRESHOLDS = (-1.0, -1.5, -2.0, -2.5, -3.0)      # -2.0 is the shipped inherited cutoff
T_BARS = (None, 2.0, 3.0, 4.0)
T_MODES = ("per_cell", "all_conditions")


def call_mechanism(ratio: float | None, wt: float, frac: float = FRAC) -> str:
    """Why is this cell called essential? PURE.

    `infeasible`      -- no growth value at all (NaN): the LP could not be solved, which after the
                         2026-08-13 probe means the ATPM maintenance floor cannot be met without this
                         gene. Genuine essentiality, and a BOOLEAN signal -- no threshold involved.
    `sub_threshold`   -- a real, finite growth value that happens to fall below `frac * wt`. This is the
                         ONLY kind of call a threshold retune could ever move.
    `not_essential`   -- above the bar.
    """
    if ratio is None or ratio != ratio:
        return "infeasible"
    return "sub_threshold" if ratio < frac * wt else "not_essential"


def score_setting(subset, calls, keys, growth, wts) -> dict:
    """Headline metrics + the commitment-mechanism breakdown for one label setting. PURE."""
    sw = switch_accuracy(subset, calls, conditions=keys)
    nulls = constant_baselines(subset, conditions=keys)
    best_null = max(g["per_condition_agreement"] for g in nulls.values())
    pat = pattern_distribution(subset, calls, conditions=keys)

    mech = {"infeasible": 0, "sub_threshold": 0}
    committed = 0
    committed_all_infeasible = 0
    for r in subset:
        ess = [c for c in keys if calls.get(c, {}).get(r.gene_id, False)]
        if not ess or len(ess) == len(keys):
            continue                                   # constant prediction: not a commitment
        committed += 1
        kinds = [call_mechanism(growth.get(c, {}).get(r.gene_id), wts[c]) for c in ess]
        for k in kinds:
            if k in mech:
                mech[k] += 1
        if kinds and all(k == "infeasible" for k in kinds):
            committed_all_infeasible += 1

    total_mech = mech["infeasible"] + mech["sub_threshold"]
    return {
        "n_conditionally_essential": len(subset),
        "exact_set_match": sw["exact_set_match"],
        "exact_set_match_rate": sw["exact_set_match_rate"],
        "per_condition_agreement": sw["per_condition_agreement"],
        "best_constant_null": best_null,
        "lift_over_null": (round(sw["per_condition_agreement"] - best_null, 4)
                           if sw["per_condition_agreement"] is not None else None),
        "constant_pattern_fraction": pat["constant_pattern_fraction"],
        "n_committed": committed,
        "n_committed_all_infeasible": committed_all_infeasible,
        "essential_calls_by_mechanism": mech,
        "infeasible_share_of_commitment_calls": (round(mech["infeasible"] / total_mech, 4)
                                                 if total_mech else None),
    }


def verdict_for_sweep(cells: list[dict]) -> dict:
    """Axis-aware verdict. PURE.

    The axes must NOT be pooled. A blind 5/10 count over the whole grid reads as "the headline is
    cutoff-dependent" when in fact every failure sits on ONE axis whose INSTRUMENT is broken
    (`all_conditions` t-mode is anti-selective for switchers). Pooling a broken instrument with a sound
    one manufactures a false negative -- the mirror of the retracted-headline mistake.
    """

    def axis(cs):
        powered = [c for c in cs if c["n_conditionally_essential"] >= 30]
        if not powered:
            return {"status": "UNDERPOWERED", "n_powered": 0, "n_beating_null": 0,
                    "lifts": [], "n_committing": 0}
        beats = [c for c in powered if (c["lift_over_null"] or 0) > 0]
        committing = [c for c in powered if c["n_committed"] > 0]
        if not committing:
            status = "INSTRUMENT_DEGENERATE_NO_COMMITMENTS"
        elif len(beats) == len(powered):
            status = "SURVIVES_EVERY_POWERED_SETTING"
        elif len(beats) >= 0.7 * len(powered):
            status = "SURVIVES_MOST_SETTINGS"
        else:
            status = "CUTOFF_DEPENDENT"
        return {"status": status, "n_powered": len(powered), "n_beating_null": len(beats),
                "lifts": [c["lift_over_null"] for c in powered], "n_committing": len(committing)}

    fit_axis = axis([c for c in cells if c["min_abs_t"] is None])
    per_cell = axis([c for c in cells if c["min_abs_t"] is not None
                     and c.get("min_abs_t_mode") == "per_cell"])
    all_cond = axis([c for c in cells if c["min_abs_t"] is not None
                     and c.get("min_abs_t_mode") == "all_conditions"])

    shares = [c["infeasible_share_of_commitment_calls"] for c in cells
              if c["n_conditionally_essential"] >= 30
              and c["infeasible_share_of_commitment_calls"] is not None]
    if shares and min(shares) >= 0.9:
        mech = "CONDITIONAL_SIGNAL_IS_BINARY_FEASIBILITY"
    elif shares and min(shares) >= 0.5:
        mech = "CONDITIONAL_SIGNAL_IS_MOSTLY_FEASIBILITY"
    elif shares:
        mech = "THRESHOLD_CROSSING_CARRIES_REAL_SIGNAL"
    else:
        mech = "UNDETERMINED"

    return {
        "fit_threshold_axis": fit_axis,
        "t_bar_axis_per_cell": per_cell,
        "t_bar_axis_all_conditions": all_cond,
        "headline_label_sensitivity": fit_axis["status"],
        "mechanism": mech,
        "min_infeasible_share": min(shares) if shares else None,
        "max_infeasible_share": max(shares) if shares else None,
        "note_on_all_conditions_mode": (
            "An `all_conditions` t-bar is ANTI-SELECTIVE for conditional essentiality: a switcher is "
            "confidently essential in ONE condition and confidently NEUTRAL (t~0) in the rest, so it "
            "fails a require-|t|-everywhere bar by construction. Its settings are reported but must "
            "NOT be pooled with the others -- the instrument, not the claim, is what failed there."),
    }


def main(argv: list[str] | None = None) -> int:
    from cobra.flux_analysis import single_gene_deletion  # noqa: PLC0415

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=None)
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    try:
        conn = open_db(a.db)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        return 2

    model = load_model()
    conds = carbon_conditions(conn, model)
    keys = tuple(sorted(conds))
    model_genes = {g.id for g in model.genes}
    print(f"{model.id}: {len(keys)} carbon sources")

    # ---- resolve every label setting FIRST, so the deletion pass covers their union exactly ----
    settings = {}
    union: set[str] = set()
    for fit in FIT_THRESHOLDS:
        for tbar in T_BARS:
            for mode in (T_MODES if tbar is not None else ("per_cell",)):
                recs = load_records(conn, conds, gene_filter=model_genes, threshold=fit,
                                    min_abs_t=tbar, min_abs_t_mode=mode)
                sub = conditionally_essential_genes(recs)
                settings[(fit, tbar, mode)] = sub
                union |= {r.gene_id for r in sub}
                print(f"   fit<{fit:<5} t>={str(tbar):<5} {mode:<15} -> {len(recs):5d} rows, "
                      f"{len(sub):4d} conditionally essential", flush=True)
    genes = sorted(union)
    print(f"\nunion to delete: {len(genes)} genes x {len(keys)} conditions = "
          f"{len(genes) * len(keys):,} knockouts (ONE pass, reused by all "
          f"{len(settings)} settings)", flush=True)

    # ---- the one expensive pass ----
    calls: dict[str, dict[str, bool]] = {}
    growth: dict[str, dict[str, float]] = {}
    wts: dict[str, float] = {}
    audits = {}
    all_ex = tuple(conds.values())
    for n, cond in enumerate(keys, 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            wt = wildtype_growth(model)
            wts[cond] = wt
            d, g_by_gene = {}, {}
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in genes])
                audits[cond] = audit_deletion_frame(res, cond)
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    gv = row["growth"]
                    g_by_gene[gid] = None if gv != gv else float(gv)
                    d[gid] = (gv != gv) or (gv < FRAC * wt)
            calls[cond] = d
            growth[cond] = g_by_gene
        print(f"   [{n:2d}/{len(keys)}] {cond[:38]:40s} wt {wt:.4f} | "
              f"{sum(d.values()):4d} essential", flush=True)

    # ---- re-score every setting from the cache ----
    cells = []
    print("\n=== LABEL-THRESHOLD SWEEP ===")
    print(f"{'fit':<7}{'t_bar':<6}{'t_mode':<10}{'n_ce':>6}{'exact':>8}{'per_cell':>10}{'null':>8}"
          f"{'lift':>9}{'const':>8}{'commit':>8}{'infeas%':>9}")
    for (fit, tbar, mode), sub in settings.items():
        c = score_setting(sub, calls, keys, growth, wts)
        c["fit_threshold"] = fit
        c["min_abs_t"] = tbar
        c["min_abs_t_mode"] = None if tbar is None else mode
        c["is_shipped_setting"] = (fit == -2.0 and tbar is None)
        cells.append(c)
        share = c["infeasible_share_of_commitment_calls"]
        mark = "  <- SHIPPED" if c["is_shipped_setting"] else ""
        print(f"{fit:<7}{str(tbar):<6}{(mode[:9] if tbar else '-'):<10}"
              f"{c['n_conditionally_essential']:>6}"
              f"{c['exact_set_match']:>8}{c['per_condition_agreement']!s:>10}"
              f"{c['best_constant_null']!s:>8}{c['lift_over_null']!s:>9}"
              f"{c['constant_pattern_fraction']!s:>8}{c['n_committed']:>8}"
              f"{'-' if share is None else f'{100 * share:.0f}%':>9}{mark}")

    v = verdict_for_sweep(cells)
    for name, key in (("fit-threshold axis    ", "fit_threshold_axis"),
                      ("t-bar axis (per_cell) ", "t_bar_axis_per_cell"),
                      ("t-bar axis (all_conds)", "t_bar_axis_all_conditions")):
        ax = v[key]
        print(f"\n{name}: {ax['status']} "
              f"({ax['n_beating_null']}/{ax['n_powered']} powered beat null; "
              f"{ax['n_committing']}/{ax['n_powered']} commit at all)")
    print(f"\nHEADLINE label sensitivity: {v['headline_label_sensitivity']}")
    print(f"MECHANISM:                  {v['mechanism']} "
          f"(infeasible share {v['min_infeasible_share']}-{v['max_infeasible_share']})")

    result = {
        "record": "fba-label-threshold-sweep-v1", "date": a.date, "model": model.id,
        "n_conditions": len(keys),
        "fit_thresholds": list(FIT_THRESHOLDS),
        "t_bars": list(T_BARS),
        "t_modes": list(T_MODES),
        "n_genes_deleted_union": len(genes),
        "essentiality_frac": FRAC,
        "cells": cells,
        "verdict": v,
        "solver_audit": merge_audits(audits) if audits else None,
        "caveats": [
            "The FBA calls are IDENTICAL across settings by construction -- only which genes are SCORED "
            "changes. So this measures label-cutoff sensitivity of the METRIC, not of the model.",
            "A lower |fit| bar admits noisier labels and a higher one shrinks the set; both move the "
            "constant null, so `lift_over_null` is the comparable column, not raw per-cell agreement.",
            "min_abs_t drops a gene failing the confidence bar in ANY condition (the complete-row rule), "
            "so it shrinks the set faster than a per-cell mask would.",
            "All 25 conditions are aerobic carbon sources; no oxygen axis.",
        ],
    }
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"fba_label_threshold_sweep_{a.date}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
