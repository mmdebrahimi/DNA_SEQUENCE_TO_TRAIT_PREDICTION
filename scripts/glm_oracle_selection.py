"""Measure what the expression oracle is worth to the GENERATIVE LOOP, not to a correlation table.

    uv run python scripts/glm_oracle_selection.py [--trials 2000] [--axes Minus10,Minus35]

The oracle reports Spearman rho ~0.59 on held-out promoter elements. The loop does not consume rho — it
consumes a CHOICE: propose k candidates, score them, keep the best. This script simulates exactly that
choice against real measured expression, on promoters built from element variants the oracle never saw, and
reports how much of the AVAILABLE selection value it captures.

**Pre-registered reading, frozen in `PREREGISTERED` below before any number was looked at.** The decision it
informs is whether the oracle is good enough to be the loop's scoring head, or whether the scoring head is
itself the bottleneck. Both thresholds are ASSERTED, not derived — said plainly in the artifact rather than
dressed up; what IS derived is the normaliser (the per-trial ceiling arm), so "captures half the headroom"
is measured against the achievable maximum rather than against a number pulled from the air.

**The honest scope limit, and it is severe.** The candidate pool here is MEASURED GRID sequences. The real
loop's candidates come from a GENERATOR and will sit off this distribution, where the oracle is less
reliable and where no measured label exists to check it. So every number here is an **UPPER BOUND** on the
loop's realisable selection value. It bounds the scoring head on its home turf; it does not simulate the
loop.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

from dna_decode.glm.expression import (
    DEFAULT_CACHE,
    ELEMENT_COLUMNS,
    fit_predict_ridge,
    leave_element_out,
    load_pairs,
)
from dna_decode.glm.selection import (
    assert_k1_is_degenerate,
    permutation_null,
    simulate_selection,
)

#: Frozen BEFORE looking at any result. Both thresholds asserted; the normaliser is derived.
PREREGISTERED = {
    "decision": "is the expression oracle good enough to be the generative loop's scoring head?",
    "at_k": 10,
    "learned_scorer": "onehot",
    "baseline_scorer": "gc",
    "bar_a_headroom_captured": 0.50,
    "bar_a_meaning": "the learned scorer realises at least half the AVAILABLE selection value",
    "bar_b_margin_over_baseline": 0.10,
    "bar_b_margin_is_noise_aware": "the effective margin is max(0.10, 2 * permutation-null stdev at that "
                                   "axis and k) -- a flat 0.10 is the SAME ORDER as the measured "
                                   "scorer-realisation spread, so on its own it could be cleared by noise",
    "bar_b_fixed_at": "instrument-build time, BEFORE any real number was computed, after measuring the null "
                      "spread on synthetic data (60 useless scorers: mean headroom 0.0006, stdev 0.096). "
                      "This is not a post-hoc re-size.",
    "bar_b_meaning": "the learned features earn their keep IN THE LOOP, not only on correlation",
    "min_axes_passing": 3,
    "n_axes": 5,
    "verdicts": {
        "ORACLE_READY_AS_SCORING_HEAD": "bar A and bar B both pass on a majority of axes",
        "CRUDE_BASELINE_SUFFICES": "bar A passes, bar B fails -- GC does the work; the learned oracle is "
                                   "not what the loop is waiting for",
        "SELECTION_IS_THE_BOTTLENECK": "bar A fails -- the scoring head cannot convert its correlation "
                                       "into a choice",
    },
    "thresholds_are": "ASSERTED, not derived; the per-trial ceiling arm that normalises them IS derived",
}

SCORERS = ("onehot", "kmer6", "gc")
K_VALUES = (1, 2, 5, 10, 20)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--trials", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--axes", default=",".join(ELEMENT_COLUMNS))
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    axes = [a.strip() for a in args.axes.split(",") if a.strip()]
    pairs, rep = load_pairs(args.cache_dir, allow_download=not args.no_download)
    print(f"loaded {rep.n_joined} measured promoters (join rate {rep.join_rate:.4f})")

    per_axis: list[dict] = []
    for axis in axes:
        train, test, held = leave_element_out(pairs, axis, seed=args.seed)
        if len(test) < max(K_VALUES):
            print(f"  {axis}: SKIP -- held-out pool {len(test)} < max k {max(K_VALUES)}")
            continue
        meas = [p.expression for p in test]
        spread = max(meas) / min(meas) if min(meas) > 0 else float("nan")
        print(f"  {axis}: held {held} -> train {len(train)} / test {len(test)} "
              f"(measured spread {spread:.1f}x)")

        # ONE null per (axis, k): a shuffled scorer is a uniform permutation whatever it came from.
        nulls = {k: permutation_null(meas, k=k, trials=max(200, args.trials // 5), seed=args.seed + k)
                 for k in K_VALUES}
        for k in K_VALUES:
            nb = nulls[k]
            if nb.get("n_perms"):
                print(f"    null k={k:<3d} mean {nb['mean']:+.4f}  stdev {nb['stdev']:.4f}  "
                      f"p95 {nb['p95']:+.4f}")

        cells: list[dict] = []
        spearmans: dict[str, float] = {}
        for scorer in SCORERS:
            fit = fit_predict_ridge(train, test, feature=scorer, return_pred=True)
            spearmans[scorer] = fit["spearman"]
            for k in K_VALUES:
                r = simulate_selection(fit["pred"], fit["measured"], k=k, scorer=scorer,
                                       trials=args.trials, seed=args.seed + k)
                if k == 1:
                    assert_k1_is_degenerate(r)      # the harness's own correctness self-test
                cells.append(r.as_dict())
            bits = []
            for c in cells:
                if c["scorer"] != scorer:
                    continue
                h = c["headroom_captured"]
                bits.append(f"k={c['k']}:" + ("--" if h is None else f"{h:.2f}"))
            print(f"    {scorer:8s} rho={fit['spearman']:+.4f}  " + "  ".join(bits))

        per_axis.append({
            "axis": axis, "held_out_values": held,
            "n_train": len(train), "n_test": len(test),
            "measured_spread_fold": round(spread, 2),
            "spearman": spearmans,
            "permutation_null_by_k": {str(k): nulls[k] for k in K_VALUES},
            "cells": cells,
        })

    verdict, detail = _verdict(per_axis)
    artifact = {
        "schema": "glm-oracle-selection-v1",
        "date": str(date.today()),
        "substrate": "GEO GSE108535 sigma70 promoter MPRA (E. coli), leave-element-out",
        "load_report": rep.as_dict(),
        "trials_per_cell": args.trials,
        "k_values": list(K_VALUES),
        "scorers": list(SCORERS),
        "preregistered": PREREGISTERED,
        "per_axis": per_axis,
        "verdict": verdict,
        "verdict_detail": detail,
        "honest_limits": [
            "UPPER BOUND on the loop's selection value: the candidate pool is MEASURED GRID sequences, "
            "while the real loop's candidates come from a generator and sit off this distribution.",
            "Candidates within one held-out element set share that element, so they vary only on the other "
            "axes -- a narrower pool than the full grid (per-axis measured spread is reported).",
            "Both pre-registered thresholds are asserted, not derived; the ceiling normaliser is derived.",
            "A designed fully-crossed grid is not a natural promoter distribution.",
            "Ridge on fixed features is the trivial-head half of the design -- a FLOOR, not a ceiling.",
        ],
    }
    out = Path(args.out) if args.out else Path(f"wiki/glm_oracle_selection_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nVERDICT: {verdict}\n  {detail}\nwrote {out}")
    return 0


def _cell(axis: dict, scorer: str, k: int) -> dict | None:
    for c in axis["cells"]:
        if c["scorer"] == scorer and c["k"] == k:
            return c
    return None


def _verdict(per_axis: list[dict]) -> tuple[str, str]:
    """Mechanical application of the frozen bar. No post-hoc re-sizing."""
    if not per_axis:
        return "INDETERMINATE_NO_AXES_SCORED", "no element axis produced a large enough held-out pool"
    k = PREREGISTERED["at_k"]
    learned, base = PREREGISTERED["learned_scorer"], PREREGISTERED["baseline_scorer"]
    a_pass = b_pass = 0
    rows = []
    for ax in per_axis:
        cl, cb = _cell(ax, learned, k), _cell(ax, base, k)
        if cl is None or cb is None or cl["headroom_captured"] is None:
            continue
        hl = cl["headroom_captured"]
        hb = cb["headroom_captured"] if cb["headroom_captured"] is not None else 0.0
        nb = ax.get("permutation_null_by_k", {}).get(str(k), {})
        sd = nb.get("stdev") or 0.0
        margin = max(PREREGISTERED["bar_b_margin_over_baseline"], 2.0 * sd)
        a = hl >= PREREGISTERED["bar_a_headroom_captured"]
        b = (hl - hb) >= margin
        a_pass += int(a); b_pass += int(b)
        rows.append(f"{ax['axis']} learned {hl:.3f} vs gc {hb:.3f} (margin needed {margin:.3f}; "
                    f"A {'y' if a else 'n'} B {'y' if b else 'n'})")
    need = PREREGISTERED["min_axes_passing"]
    detail = f"at k={k}, {a_pass}/{len(rows)} axes pass bar A, {b_pass}/{len(rows)} pass bar B | " + "; ".join(rows)
    # A bar requiring `need` passing axes is UNREACHABLE on fewer than `need` scored axes, so a substantive
    # verdict there would be an artifact of the run's scope, not a finding. Caught on the first smoke run:
    # a single-axis invocation printed SELECTION_IS_THE_BOTTLENECK while reporting that every scored axis
    # passed BOTH bars. Refuse instead -- a partial run must never emit a headline.
    if len(rows) < need:
        return ("INDETERMINATE_INSUFFICIENT_AXES",
                f"only {len(rows)} axes scored but the frozen bar needs {need} passing; "
                f"the bar is unreachable by construction. {detail}")
    if a_pass < need:
        return "SELECTION_IS_THE_BOTTLENECK", detail
    if b_pass < need:
        return "CRUDE_BASELINE_SUFFICES", detail
    return "ORACLE_READY_AS_SCORING_HEAD", detail


if __name__ == "__main__":
    raise SystemExit(main())
