"""Why do 147 of 2,015 nitrogen cells differ between two identical passes?

WHY THIS EXISTS
---------------
`scripts/fba_conditional_nitrogen.py` ran its pre-registered determinism gate and FAILED it:

    "deterministic": false,
    "n_cells_differing_between_passes": 147

The prereg (`wiki/fba_nitrogen_prereg_2026-08-17.md`) made that gate mandatory precisely because of
today's retraction (`wiki/fba_eflux_graded_RETRACTION_2026-08-17.md`), where a pre-registered bar was
cleared by a run whose run-to-run variance exceeded the effect. So a determinism failure here must be
DIAGNOSED, not tolerated and not silently relaxed.

There are two very different explanations, and they are NOT equally bad:

  (H1) FLOAT NOISE MISREAD AS NON-DETERMINISM. The gate tests `abs(delta_ratio) > 1e-12`. A ratio is an
       LP objective divided by another LP objective; drifts of 1e-13..1e-9 are ordinary float64
       behaviour, not a different answer. If every drift is tiny and NO cell crosses the FRAC=0.01 call
       line, then the RESULT is reproducible even though the bits are not.

  (H2) GENUINE RUN-TO-RUN VARIANCE. Solver state carried over from pass 1 changes pass 2's answers by
       amounts that flip essentiality calls and move the reported metrics. That is the SAME failure
       class as the retraction and would invalidate the nitrogen numbers outright.

DISCRIMINATION
--------------
H1 and H2 predict different magnitude distributions and different call behaviour, so ONE measurement
separates them. This probe reports, for each pass-pair:

  * the magnitude distribution of the disagreements (max / p99 / median), not just a count
  * how many cells CROSS the FRAC=0.01 threshold (the only drift that can change a claim)
  * whether the per-cell agreement metric itself changes

It also runs a third pass against a FRESHLY LOADED model, which separates "solver state carried within
one process" from "the panel is just not reproducible".

This probe deliberately does NOT decide the tolerance. It measures first; the tolerance question is
answered by what the numbers turn out to be.

Run: .venv/Scripts/python.exe scripts/fba_nitrogen_determinism_probe.py
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.conditional_essentiality import conditionally_essential_genes  # noqa: E402
from dna_decode.fba.fitness_browser import ESSENTIAL_FITNESS, open_db  # noqa: E402
from dna_decode.fba.model import load_model  # noqa: E402
from dna_decode.fba.nitrogen import (  # noqa: E402
    load_nitrogen_records,
    nitrogen_conditions,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fba_conditional_nitrogen import FRAC, run_panel  # noqa: E402

GATE_TOL = 1e-12          # the tolerance the failing gate used


def compare(label, ra, rb, keys, genes):
    """Magnitude + call-level comparison of two passes. Counts are over the SAME cells."""
    deltas, crossers = [], []
    for c in keys:
        for g in genes:
            a = ra[c].get(g, 0.0)
            b = rb[c].get(g, 0.0)
            d = abs(a - b)
            if d > GATE_TOL:
                deltas.append(d)
                if (a <= FRAC) != (b <= FRAC):        # the ONLY drift that can change a claim
                    crossers.append((c, g, a, b))
    n_cells = len(keys) * len(genes)
    out = {
        "pair": label,
        "n_cells": n_cells,
        "n_differ_above_gate_tol": len(deltas),
        "max_abs_delta": max(deltas) if deltas else 0.0,
        "median_abs_delta": statistics.median(deltas) if deltas else 0.0,
        "n_call_flips": len(crossers),
        "call_flip_examples": [(c, g, round(a, 6), round(b, 6)) for c, g, a, b in crossers[:5]],
    }
    print(f"\n--- {label} ---")
    print(f"  cells differing by > {GATE_TOL:g} : {len(deltas)} of {n_cells}")
    if deltas:
        ds = sorted(deltas)
        print(f"  |delta|  max {ds[-1]:.3e}   p99 {ds[int(0.99*(len(ds)-1))]:.3e}   "
              f"median {statistics.median(ds):.3e}   min {ds[0]:.3e}")
    print(f"  cells CROSSING the FRAC={FRAC} call line : {len(crossers)}   <-- the number that matters")
    for c, g, a, b in crossers[:5]:
        print(f"      {c:28} {g:8} {a:.6g} vs {b:.6g}")
    return out


def agreement(ratios, subset, keys, genes):
    """The reported headline metric, recomputed from a given pass -- does it move between passes?"""
    right = total = 0
    for c in keys:
        for r in subset:
            pred = ratios[c].get(r.gene_id, 0.0) <= FRAC
            right += int(pred == r.experimental[c])
            total += 1
    return round(right / total, 6) if total else None


def main():
    from cobra.flux_analysis import single_gene_deletion

    model = load_model()
    conn = open_db()
    conds = nitrogen_conditions(conn, model)
    keys = sorted(conds)
    records = load_nitrogen_records(conn, conds, gene_filter={g.id for g in model.genes},
                                    threshold=ESSENTIAL_FITNESS)
    subset = conditionally_essential_genes(records)
    genes = [r.gene_id for r in subset]

    print("NITROGEN DETERMINISM PROBE - is the 147-cell disagreement float noise or real variance?")
    print("=" * 94)
    print(f"{len(genes)} two-sided genes x {len(keys)} conditions = {len(genes)*len(keys)} cells per pass")

    print("\n=== PASS A ===", flush=True)
    ra, _ = run_panel(model, conds, keys, genes, single_gene_deletion)
    print("\n=== PASS B (same process, same model object -- reproduces the reported gate) ===", flush=True)
    rb, _ = run_panel(model, conds, keys, genes, single_gene_deletion)
    print("\n=== PASS C (FRESHLY loaded model -- separates carried solver state) ===", flush=True)
    model2 = load_model()
    rc, _ = run_panel(model2, conds, keys, genes, single_gene_deletion)

    cmps = [compare("A vs B (in-process)", ra, rb, keys, genes),
            compare("A vs C (fresh model)", ra, rc, keys, genes)]

    agrees = {"A": agreement(ra, subset, keys, genes),
              "B": agreement(rb, subset, keys, genes),
              "C": agreement(rc, subset, keys, genes)}
    print(f"\nHEADLINE METRIC per pass: {agrees}")

    max_delta = max(c["max_abs_delta"] for c in cmps)
    flips = sum(c["n_call_flips"] for c in cmps)
    stable = len(set(agrees.values())) == 1

    print("\nVERDICT")
    if flips == 0 and stable:
        print(f"  H1 SUPPORTED: drifts are numerical only (max {max_delta:.3e}), ZERO cells cross the")
        print(f"  {FRAC} call line, and the headline metric is identical across all three passes.")
        print("  The RESULT is reproducible; the 1e-12 bit-equality gate was measuring the wrong thing.")
    else:
        print(f"  H2 SUPPORTED: {flips} call flip(s), metric stable={stable} (max delta {max_delta:.3e}).")
        print("  This is real run-to-run variance -- the nitrogen numbers must NOT be reported.")
    return {"comparisons": cmps, "agreement_per_pass": agrees, "n_genes": len(genes),
            "n_conditions": len(keys), "max_abs_delta": max_delta, "n_call_flips": flips,
            "metric_stable_across_passes": stable}


def _verify(res):
    """The probe must be non-vacuous: it has to actually observe the failure it was built to explain."""
    ab = res["comparisons"][0]
    assert ab["n_cells"] > 1000, f"expected a full panel, got {ab['n_cells']} cells"
    assert ab["n_differ_above_gate_tol"] > 0, (
        "the probe observed ZERO disagreements, so it did NOT reproduce the reported gate failure -- "
        "its verdict would be vacuous. Do not trust this run.")
    assert res["max_abs_delta"] > 0.0, "no drift measured at all"
    # A call flip and a stable metric are contradictory only if the flip is in a scored cell; report both
    # rather than asserting, so the probe cannot 'pass' by hiding a contradiction.
    print(f"\n[verify-in-batch] PASS: reproduced the failure ({ab['n_differ_above_gate_tol']} cells "
          f"differ above the {GATE_TOL:g} gate tolerance, so the verdict is not vacuous); measured the "
          f"largest drift anywhere at {res['max_abs_delta']:.3e}; counted {res['n_call_flips']} cell(s) "
          f"crossing the FRAC={FRAC} line (the only drift that can change a claim); and recomputed the "
          f"headline metric independently in all three passes "
          f"(stable={res['metric_stable_across_passes']}).")


if __name__ == "__main__":
    r = main()
    _verify(r)
    Path("wiki/fba_nitrogen_determinism_probe_2026-08-17.json").write_text(
        json.dumps(r, indent=2), encoding="utf-8")
    print("wrote wiki/fba_nitrogen_determinism_probe_2026-08-17.json")
