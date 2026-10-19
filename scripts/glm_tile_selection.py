"""What is the tile rung's best scorer worth to a CHOICE, not to a correlation table?

    uv run python scripts/glm_tile_selection.py [--trials 2000] [--seeds 10]

**The question, and why it is the one that matters.** The GLM portfolio's live `ceiling:` criterion asks
for the best achievable PREDICTION of promoter activity. Nine arms are declared against it. But the
project already measured, on 2026-10-07, that prediction is not the binding constraint: on the designed
sigma70 grid the learned scorer at rho ~0.59 realised a cross-axis mean **0.7937** of the AVAILABLE
selection value at k=10 (GC: **0.1328**), verdict `ORACLE_READY_AS_SCORING_HEAD`, and the recorded
consequence was *"the scoring half is adequate and the GENERATOR is the bottleneck."*

That measurement was taken on a DIFFERENT substrate. Every arm since has run on the GSE144621 peak-tile
rung, where the best scorer is the 346-parameter Promoter Calculator at rho **0.3968** against a 0.941
replicate ceiling -- 42% of it. **Nobody has measured what 0.3968 is worth to a choice.** Until that
number exists, "keep pushing prediction" and "the criterion is the wrong target" are both unevidenced.

This script produces that number. It is a DECISION INSTRUMENT, not an attempt to raise anything --
which is exactly why it is deliberately NOT declared as an arm in the ceiling criterion's arm set. An
arm is a candidate approach to improving the metric; this asks whether the metric is worth improving.

**It reuses `dna_decode/glm/selection.py` unchanged.** No new statistics are written here: `picked` vs
a PAIRED `random` from the same candidate set vs the measured `ceiling` of that set, with
`headroom_captured = (picked - random)/(ceiling - random)`. The k=1 degeneracy self-test and the
uniform tie-breaking are properties of that module, already tested.

**Every headroom is reported beside its own permutation null, and that is load-bearing.** A single
cell's headroom carries scorer-realisation spread that `trials` cannot reduce -- measured on synthetic
data, 60 useless scorers over a 300-element pool at k=10 gave mean 0.0006 with **stdev 0.096**. A
headroom quoted without its band is uninterpretable.

**The honest scope limit, and it is the same one the designed-grid run carried.** The candidate pool is
MEASURED peak tiles. The real generative loop's candidates come from a GENERATOR and sit off this
distribution, where no measured label exists. So every number here is an UPPER BOUND on the loop's
realisable selection value. It is a DIFFERENT upper bound from the grid's (natural, peak-enriched
sequence rather than a fully-crossed design), not a stricter one.
"""
from __future__ import annotations

import argparse
import json
import statistics
from datetime import date
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from glm_tile_predictions import DEFAULT_OUT, load_table, require_real  # noqa: E402

from dna_decode.glm.selection import (  # noqa: E402
    assert_k1_is_degenerate,
    permutation_null,
    simulate_selection,
)

#: FROZEN before any tile-rung number was read, and deliberately MIRRORING
#: `scripts/glm_oracle_selection.py::PREREGISTERED` so the two substrates are directly comparable.
#: The primary scorer is `max_tx` because it is the best MEASURED scorer on this rung (0.3968) -- the
#: question is whether the scorer we would actually deploy converts its correlation into a choice.
PREREGISTERED = {
    "decision": "is the tile rung's best scorer good enough to be the generative loop's scoring head?",
    "at_k": 10,
    "primary_scorer": "max_tx",
    "baseline_scorer": "gc",
    "bar_a_headroom_captured": 0.50,
    "bar_a_meaning": "the primary scorer realises at least half the AVAILABLE selection value",
    "bar_b_margin_over_baseline": 0.10,
    "bar_b_margin_is_noise_aware": "the effective margin is max(0.10, 2 * permutation-null stdev at "
                                   "that seed and k) -- a flat 0.10 is the SAME ORDER as the measured "
                                   "scorer-realisation spread (0.096), so alone it could be cleared "
                                   "by noise",
    "bar_b_meaning": "the scorer earns its keep IN THE LOOP, not only on correlation",
    "min_seeds_passing": 6,
    "n_seeds_expected": 10,
    "thresholds_are": "ASSERTED, not derived; the per-trial ceiling arm that normalises them IS derived",
    "mirrors": "scripts/glm_oracle_selection.py (designed grid): bar A 0.50, bar B max(0.10, 2*sd)",
    "comparison_anchor": {
        "artifact": "wiki/glm_oracle_selection_2026-10-07.json",
        "substrate": "GSE108535 designed sigma70 grid, leave-element-out",
        "verdict": "ORACLE_READY_AS_SCORING_HEAD",
        "onehot_mean_headroom_at_k10": 0.7937,
        "gc_mean_headroom_at_k10": 0.1328,
        "note": "cross-AXIS means over 5 element axes; this run reports cross-SEED means over 10 "
                "leave-peak-out splits. Different pooling units -- compare with that in mind.",
    },
    "verdicts": {
        "ORACLE_READY_AS_SCORING_HEAD": "bar A and bar B both pass on a majority of seeds -- further "
                                        "prediction work has measured diminishing returns FOR THE "
                                        "GENERATIVE GOAL",
        "CRUDE_BASELINE_SUFFICES": "bar A passes, bar B fails -- GC does the work; the biophysical "
                                   "scorer is not what the loop is waiting for",
        "SELECTION_IS_THE_BOTTLENECK": "bar A fails -- the scoring head cannot convert its "
                                       "correlation into a choice, so prediction work IS justified",
    },
}

K_VALUES = (1, 2, 5, 10, 20)


def cell_for(cells: list[dict], scorer: str, k: int) -> dict | None:
    for c in cells:
        if c["scorer"] == scorer and c["k"] == k:
            return c
    return None


def verdict(per_seed: list[dict], pre: dict = PREREGISTERED) -> tuple[str, str]:
    """Mechanical application of the frozen bar. No post-hoc re-sizing.

    Refuses with `INDETERMINATE_INSUFFICIENT_SEEDS` when fewer seeds scored than the bar needs -- the
    bar is then UNREACHABLE by construction and any substantive verdict would be an artifact of the
    run's scope. The sibling script earned this branch: a single-axis invocation once printed
    `SELECTION_IS_THE_BOTTLENECK` while reporting that every scored axis passed BOTH bars.
    """
    if not per_seed:
        return "INDETERMINATE_NO_SEEDS_SCORED", "no seed produced a scoreable pool"
    k = pre["at_k"]
    prim, base = pre["primary_scorer"], pre["baseline_scorer"]
    rows, a_pass, b_pass = [], 0, 0
    for s in per_seed:
        cp, cb = cell_for(s["cells"], prim, k), cell_for(s["cells"], base, k)
        if cp is None or cb is None or cp["headroom_captured"] is None:
            continue
        hp = cp["headroom_captured"]
        hb = cb["headroom_captured"] if cb["headroom_captured"] is not None else 0.0
        sd = (s.get("permutation_null_by_k", {}).get(str(k), {}) or {}).get("stdev") or 0.0
        margin = max(pre["bar_b_margin_over_baseline"], 2.0 * sd)
        a, b = hp >= pre["bar_a_headroom_captured"], (hp - hb) >= margin
        a_pass += int(a)
        b_pass += int(b)
        rows.append(f"seed {s['seed']} {prim} {hp:.3f} vs {base} {hb:.3f} "
                    f"(margin needed {margin:.3f}; A {'y' if a else 'n'} B {'y' if b else 'n'})")
    need = pre["min_seeds_passing"]
    detail = (f"at k={k}, {a_pass}/{len(rows)} seeds pass bar A, {b_pass}/{len(rows)} pass bar B | "
              + "; ".join(rows))
    if len(rows) < need:
        return ("INDETERMINATE_INSUFFICIENT_SEEDS",
                f"only {len(rows)} seeds scored but the frozen bar needs {need} passing; the bar is "
                f"unreachable by construction. {detail}")
    if a_pass < need:
        return "SELECTION_IS_THE_BOTTLENECK", detail
    if b_pass < need:
        return "CRUDE_BASELINE_SUFFICES", detail
    return "ORACLE_READY_AS_SCORING_HEAD", detail


def render_md(art: dict) -> str:
    pre = art["preregistered"]
    k = pre["at_k"]
    L = [f"# Tile-rung selection value — {art['date']}", "",
         f"**VERDICT: `{art['verdict']}`**", "", art["verdict_detail"], "",
         "## What this measures", "",
         "Not a correlation — a CHOICE: propose k candidates, score them, keep the best. "
         f"`headroom_captured` is the fraction of the ACHIEVABLE selection value the scorer realises "
         f"at k={k}, normalised by the measured best of the same candidate set.", "",
         f"Substrate: {art['substrate']}", "",
         "## Cross-seed mean headroom captured", "", "| scorer | " +
         " | ".join(f"k={kk}" for kk in art["k_values"]) + " |",
         "|---|" + "---|" * len(art["k_values"])]
    for sc in art["scorers"]:
        cells = []
        for kk in art["k_values"]:
            vals = [h for s in art["per_seed"]
                    if (c := cell_for(s["cells"], sc, kk)) and (h := c["headroom_captured"]) is not None]
            cells.append("—" if not vals else f"{statistics.fmean(vals):.4f}")
        L.append(f"| `{sc}` | " + " | ".join(cells) + " |")
    anchor = pre["comparison_anchor"]
    L += ["", "## Against the designed-grid anchor", "",
          f"- {anchor['substrate']}: learned **{anchor['onehot_mean_headroom_at_k10']}**, "
          f"gc **{anchor['gc_mean_headroom_at_k10']}** (verdict `{anchor['verdict']}`)",
          f"- {anchor['note']}", "",
          "## Null band (scorer-realisation spread, per seed and k)", "",
          "A headroom without this band is uninterpretable: `trials` reduces Monte-Carlo error in "
          "estimating ONE scorer's headroom and says nothing about how far it sits from zero.", ""]
    for s in art["per_seed"][:3]:
        nb = s["permutation_null_by_k"].get(str(k), {})
        if nb.get("n_perms"):
            L.append(f"- seed {s['seed']} k={k}: mean {nb['mean']:+.4f}, stdev {nb['stdev']:.4f}, "
                     f"p95 {nb['p95']:+.4f}")
    L += ["", "## Honest limits", ""] + [f"- {x}" for x in art["honest_limits"]]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", default=str(DEFAULT_OUT))
    ap.add_argument("--trials", type=int, default=2000)
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    in_dir = Path(a.in_dir)
    paths = sorted(in_dir.glob("seed-*.json.gz"))[: a.seeds]
    if not paths:
        raise SystemExit(f"no seed tables in {in_dir}; run scripts/glm_tile_predictions.py first")
    tables = [load_table(p) for p in paths]
    require_real(tables)                 # a smoke table skipped the reconcile gate and cannot score

    man_path = in_dir / "manifest.json"
    if not man_path.is_file():
        raise SystemExit(f"{man_path} is absent -- without the manifest there is no record of WHICH "
                         f"table this is (tile count, split, protocol, anchors). Refusing.")
    man = json.loads(man_path.read_text(encoding="utf-8"))

    scorers = sorted({s for t in tables for s in t["predictions"]})
    print(f"{len(tables)} seed tables | scorers {scorers} | "
          f"pool ~{statistics.fmean([t['n_test'] for t in tables]):.0f} tiles", flush=True)

    per_seed: list[dict] = []
    for t in tables:
        meas = t["measured"]
        nulls = {k: permutation_null(meas, k=k, trials=max(200, a.trials // 5), seed=t["seed"] + k)
                 for k in K_VALUES}
        cells: list[dict] = []
        for sc in scorers:
            if sc not in t["predictions"]:
                continue
            for k in K_VALUES:
                r = simulate_selection(t["predictions"][sc], meas, k=k, scorer=sc,
                                       trials=a.trials, seed=t["seed"] + k)
                if k == 1:
                    assert_k1_is_degenerate(r)     # the harness's own correctness self-test
                cells.append(r.as_dict())
            bits = []
            for c in cells:
                if c["scorer"] != sc:
                    continue
                h = c["headroom_captured"]
                bits.append(f"k={c['k']}:" + ("--" if h is None else f"{h:.2f}"))
            print(f"  seed {t['seed']} {sc:14s} rho={t['spearman_recomputed'][sc]:+.4f}  "
                  + "  ".join(bits), flush=True)
        per_seed.append({
            "seed": t["seed"], "n_test": t["n_test"],
            "spearman": t["spearman_recomputed"],
            "permutation_null_by_k": {str(k): nulls[k] for k in K_VALUES},
            "cells": cells,
        })

    v, detail = verdict(per_seed)
    art = {
        "schema": "glm-tile-selection-v1",
        "date": str(date.today()),
        "substrate": f"GEO GSE144621 peak tiles, {man['n_tiles']} tiles, {man['split']}, "
                     f"ceiling 0.941",
        "prediction_table_manifest": man,
        "trials_per_cell": a.trials,
        "k_values": list(K_VALUES),
        "scorers": scorers,
        "preregistered": PREREGISTERED,
        "per_seed": per_seed,
        "verdict": v,
        "verdict_detail": detail,
        "honest_limits": [
            "UPPER BOUND on the loop's selection value: the candidate pool is MEASURED peak tiles, "
            "while a generator's candidates sit off this distribution where no label exists.",
            "A DIFFERENT upper bound from the designed-grid run's, not a stricter one -- natural "
            "peak-enriched sequence versus a fully-crossed design.",
            "Both bars are ASSERTED, not derived; the per-trial ceiling arm that normalises them IS "
            "derived.",
            "Cross-seed means pool over leave-peak-out splits of ONE experiment; the designed-grid "
            "anchor pooled over element axes. Different units.",
            "Peak tiles are selected from ChIP peaks, so this pool is enriched for real promoters "
            "relative to the genome -- selection is easier here than on random sequence.",
        ],
    }
    out = Path(a.out) if a.out else Path(f"wiki/glm_tile_selection_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(art, indent=1), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(art), encoding="utf-8")
    print(f"\nVERDICT: {v}\n  {detail[:400]}\nwrote {out} + {out.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
