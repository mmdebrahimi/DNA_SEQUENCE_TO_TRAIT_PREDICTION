"""Is being WEAK and DIFFERENTLY-SOURCED enough to lift, with no learned component at all?

    uv run python scripts/glm_tile_maxtx_gc_arm.py [--seeds 10]

**Its own arm, deliberately.** `scripts/glm_tile_hybrid_arm.py` froze the question "does rank-averaging
the learned conv encoder with the biophysical max_tx scalar beat the scalar alone?" -- conv-specific by
construction. This pairing has NO conv in it, so filing it under that arm would retroactively widen a
frozen question. Declared separately, before measuring, per R6 (a declared arm cannot evaporate).

**The question this isolates.** The hybrid arm measured that `conv` alone LOSES to `max_tx` (-0.0767,
0/10) and still LIFTS it (+0.0339, 10/10), with a shuffled-conv control collapsing to -0.1232. Two
accounts survive that result and the arm cannot separate them:

* **source-diversity** -- any weak, differently-sourced scorer lifts, because rank-averaging cancels
  error that is not shared; the conv's own quality is beside the point.
* **conv content** -- the conv specifically carries something `max_tx` lacks, and a different weak
  scorer would not do.

`max_tx + gc` discriminates. `gc` is weak standalone (0.3112, the weakest measured arm on this rung)
and differently sourced (base composition -- neither alignment-derived like `max_tx` nor learned like
`conv`). Under source-diversity it should lift; under conv-content it should not.

**There is already a PRIOR STRIKE against source-diversity, and this bar is written knowing it.** The
same artifact carries a three-way arm: `hybrid_conv_maxtx_gc` at **+0.0302** against the two-way's
**+0.0339** -- adding `gc` to the pair made it slightly WORSE. So a null result here is the OUTCOME
THE EXISTING EVIDENCE POINTS AT, and the bar below is NOT written in the expectation of a lift. What
the three-way cannot do is isolate: with conv present, `gc` is a third channel competing against an
already-orthogonal pair, not a lone second modality against the comparator.

**The control is what makes a positive interpretable.** `max_tx` rank-averaged with a SHUFFLED `gc`
vector: if that lifts too, the gain came from rank-averaging two vectors rather than from `gc`'s
content, and the headline is an artifact. Same shape as the hybrid arm's shuffled-conv control.

**This does NOT discharge `orthogonal-modality-stack-promcalc-plus-N-channels`** (status `planned` in
`soraya_runs/arms-ledger.json`). That arm stacks N channels onto promcalc; this isolates ONE weak
differently-sourced channel against one comparator. A planned arm blocks its criterion and this run
leaves it blocked.

**No fitting, no calibrator, no second-level model.** Both prediction vectors are already out-of-sample
on each seed's held-out peaks, and `rank_average_hybrid` has no parameters, so the combine introduces
no leakage. It RANKS; it does not dose.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from glm_tile_predictions import (  # noqa: E402
    DEFAULT_OUT,
    load_table,
    reconcile_problems,
    require_real,
)
from promcalc_decomposition_arm import paired  # noqa: E402

from dna_decode.forward.variant_effect import rank_average_hybrid  # noqa: E402
from dna_decode.glm.genomewide import spearman  # noqa: E402

#: FROZEN before any number was read -- and, unusually, before the data was even REACHABLE: the
#: prediction tables live on a volume that is detached from this host, so this bar could not have been
#: tuned to a result that nobody can compute yet.
PREREGISTERED = {
    "question": "does rank-averaging the WEAK, differently-sourced gc scalar with the biophysical "
                "max_tx scalar beat the scalar alone on the tile rung, with NO learned component?",
    "isolates": "source-diversity (any weak orthogonal scorer lifts) vs conv-content (the conv "
                "specifically carries something max_tx lacks)",
    "primary_metric": "spearman",
    "comparator": "max_tx",
    "margin": 0.02,
    "min_seeds_winning": 8,
    "n_seeds_expected": 10,
    "mirrors": "scripts/glm_tile_hybrid_arm.py PREREGISTERED (margin 0.02, 8/10 seeds) so the two "
               "combine results are directly comparable",
    "prior_strike": {
        "fact": "the three-way hybrid_conv_maxtx_gc scores +0.0302 against the two-way "
                "hybrid_conv_maxtx's +0.0339 -- adding gc to the pair made it slightly worse; and "
                "gc-single-feature is the weakest measured arm on this rung at 0.3112",
        "consequence": "a NO_GAIN verdict is the outcome the existing evidence points at. This bar is "
                       "NOT written in the expectation of a lift, and NO_GAIN is a real result here, "
                       "not a failed run",
        "why_still_worth_running": "the three-way cannot ISOLATE -- with conv present, gc is a third "
                                   "channel against an already-orthogonal pair, never a lone second "
                                   "modality against the comparator",
    },
    "precondition": {
        "name": "orthogonality",
        "rule": "a naive rank-average lifts only when the modalities are ORTHOGONAL; a near-chance "
                "second modality REDUCES the signal (measured: ESM2+ProSST 0.204 vs ESM2 0.352 on "
                "CTX-M-14, because ProSST alone was at chance)",
        "measured_first": True,
        "not_a_gate": "the orthogonality number is REPORTED, not used to skip the run -- a prediction "
                      "that a combine will fail is not a substitute for measuring it",
        "note": "gc sits at 0.3112 standalone against max_tx's 0.3968. That is WEAK, not at chance, "
                "which is the regime where the published rule says a combine can still pay",
    },
    "control": "max_tx rank-averaged with a SHUFFLED gc vector; if this lifts, the gain is an "
               "artifact of rank-averaging rather than of gc's content",
    "does_not_discharge": {
        "arm": "orthogonal-modality-stack-promcalc-plus-N-channels",
        "status_when_frozen": "planned",
        "why": "that arm stacks N channels onto promcalc; this isolates ONE weak differently-sourced "
               "channel against one comparator. This run leaves that arm's criterion BLOCKED",
    },
    "verdicts": {
        "SOURCE_DIVERSITY_SUFFICES": "paired median >= +0.02 over max_tx AND >= 8/10 seeds win AND "
                                     "control clean -- a weak differently-sourced scorer lifts with "
                                     "no learned component, so component standalone score is the "
                                     "wrong selection criterion",
        "WEAK": "gain > 0 but below the margin or the seed count",
        "NO_GAIN": "no positive paired gain -- being weak and differently-sourced is NOT sufficient, "
                   "and the hybrid arm's lift is better explained by the conv's CONTENT",
        "INDETERMINATE_CONTROL_ALSO_LIFTS": "the shuffled-gc control lifts too; the combine, not gc, "
                                            "is producing the number",
        "INDETERMINATE_RECONCILE_FAILED": "the prediction table's published anchors drift",
        "INDETERMINATE_MISSING_SCORER": "gc or max_tx absent from the table",
        "INSUFFICIENT_DATA": "the prediction tables are not reachable on this host",
    },
}

COMBO = ("max_tx", "gc")


def as_table(pred: list[float]) -> dict[str, float]:
    """A prediction vector as the `{key: score}` table `rank_average_hybrid` consumes.

    Positional keys, exactly as the hybrid arm does it: every vector in one seed's table is aligned to
    the same `tile_indices`, and a positional key makes a length mismatch impossible to paper over.
    """
    return {str(i): float(v) for i, v in enumerate(pred)}


def combine(preds: dict[str, list[float]], names: tuple[str, ...]) -> list[float]:
    """Rank-average `names` into one score vector, in the original index order."""
    tables = [as_table(preds[n]) for n in names]
    combined = rank_average_hybrid(tables)
    n = len(preds[names[0]])
    return [combined[str(i)] for i in range(n)]


def shuffled(vec: list[float], seed: int) -> list[float]:
    """A permutation of `vec` -- the control's second modality. Same distribution, no content."""
    out = list(vec)
    random.Random(2000 + seed).shuffle(out)
    return out


def orthogonality(a: list[float], b: list[float]) -> float:
    """Spearman between two scorers' predictions. Near 1 = redundant; near 0 = orthogonal."""
    return spearman(a, b)


def verdict(arms: dict, pre: dict = PREREGISTERED) -> tuple[str, str]:
    """Mechanical application of the frozen bar, ordered so a refusal wins over a headline."""
    if arms.get("_missing_scorer"):
        return ("INDETERMINATE_MISSING_SCORER",
                f"the table lacks {arms['_missing_scorer']}; nothing to combine")
    if arms.get("_reconcile_problems"):
        return ("INDETERMINATE_RECONCILE_FAILED", "; ".join(arms["_reconcile_problems"]))
    p = arms["_paired_hybrid_vs_comparator"]
    ctrl = arms.get("_paired_control_vs_comparator") or {}
    ctrl_lifts = (ctrl.get("median") or 0.0) >= pre["margin"] and \
                 (ctrl.get("wins") or 0) >= pre["min_seeds_winning"]
    if ctrl_lifts:
        return ("INDETERMINATE_CONTROL_ALSO_LIFTS",
                f"the shuffled-gc control gains {ctrl['median']:+.4f} ({ctrl['wins']}/{ctrl['n']}) "
                f"-- indistinguishable from the real combine's {p['median']:+.4f}; the number is an "
                f"artifact of rank-averaging, not of gc")
    if p["median"] >= pre["margin"] and p["wins"] >= pre["min_seeds_winning"]:
        return ("SOURCE_DIVERSITY_SUFFICES",
                f"max_tx+gc {p['median']:+.4f} over {pre['comparator']}, {p['wins']}/{p['n']} seeds, "
                f"CI95 {p['ci95']} -- a weak differently-sourced scorer lifts with no learned "
                f"component")
    if p["median"] > 0:
        return ("WEAK",
                f"max_tx+gc {p['median']:+.4f} over {pre['comparator']}, {p['wins']}/{p['n']} seeds, "
                f"CI95 {p['ci95']} -- positive but short of the {pre['margin']:+.2f} margin or the "
                f"{pre['min_seeds_winning']}-seed bar")
    return ("NO_GAIN",
            f"max_tx+gc {p['median']:+.4f} over {pre['comparator']}, {p['wins']}/{p['n']} seeds, "
            f"CI95 {p['ci95']} -- being weak and differently-sourced is NOT sufficient; the hybrid "
            f"arm's lift is better explained by the conv's CONTENT")


def render_md(art: dict) -> str:
    a = art["arms"]
    L = [f"# Tile-rung source-diversity arm (max_tx ⊕ gc, no conv) — {art['date']}", "",
         f"**VERDICT: `{art['verdict']}`**", "", art["verdict_detail"], "",
         "## What this isolates", "",
         "`conv` alone LOSES to `max_tx` and still LIFTS it. Two accounts survive that and the hybrid "
         "arm cannot separate them: **source-diversity** (any weak orthogonal scorer lifts) and "
         "**conv content** (the conv specifically carries something). `gc` is weak (0.3112) and "
         "differently sourced (composition), so this pairing discriminates.", "",
         "**Prior strike, declared before measuring:** the three-way `hybrid_conv_maxtx_gc` scores "
         "+0.0302 against the two-way's +0.0339 — adding `gc` made it slightly worse. `NO_GAIN` is "
         "the outcome the existing evidence points at, and is a real result here.", "",
         "## The precondition, measured first", "",
         f"Spearman between the `gc` and `max_tx` prediction vectors, per seed: "
         f"median **{art['orthogonality']['median']:.4f}** "
         f"(min {art['orthogonality']['min']:.4f}, max {art['orthogonality']['max']:.4f}).", "",
         "## Arms", "", "| arm | median spearman | per-seed |", "|---|---|---|"]
    for name in sorted(k for k in a if not k.startswith("_")):
        arm = a[name]
        L.append(f"| `{name}` | {arm['median']:+.4f} | "
                 + ", ".join(f"{v:+.3f}" for v in arm["per_seed"]) + " |")
    L += ["", "## Paired margins (same splits, per seed)", ""]
    for key, label in (("_paired_hybrid_vs_comparator", "max_tx+gc vs max_tx"),
                       ("_paired_control_vs_comparator", "CONTROL (shuffled gc) vs max_tx"),
                       ("_paired_gc_vs_comparator", "gc alone vs max_tx")):
        p = a.get(key)
        if p:
            L.append(f"- **{label}**: median {p['median']:+.4f}, {p['wins']}/{p['n']} seeds, "
                     f"CI95 {p['ci95']}")
    L += ["", "## Honest limits", ""] + [f"- {x}" for x in art["honest_limits"]]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", default=str(DEFAULT_OUT))
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    in_dir = Path(a.in_dir)
    if not in_dir.is_dir():
        print(f"VERDICT: INSUFFICIENT_DATA\n  prediction tables not reachable at {in_dir}; "
              f"the substrate volume is not attached to this host")
        return 2
    paths = sorted(in_dir.glob("seed-*.json.gz"))[: a.seeds]
    if not paths:
        print(f"VERDICT: INSUFFICIENT_DATA\n  no seed tables in {in_dir}; "
              f"run scripts/glm_tile_predictions.py first")
        return 2
    tables = [load_table(p) for p in paths]
    require_real(tables)

    arms: dict = {}
    recon: list[str] = []
    for t in tables:
        recon += reconcile_problems(t["seed"], t["spearman_recomputed"])
    if recon:
        arms["_reconcile_problems"] = recon

    missing = [s for s in COMBO if any(s not in t["predictions"] for t in tables)]
    if missing:
        arms["_missing_scorer"] = missing
        v, detail = verdict(arms)
        print(f"VERDICT: {v}\n  {detail}")
        return 1

    runs: dict[str, list[float]] = {k: [] for k in
                                    ("hybrid_maxtx_gc", "max_tx", "gc", "control_shuffled_gc")}
    orth: list[float] = []
    for t in tables:
        preds, meas = t["predictions"], t["measured"]
        orth.append(orthogonality(preds["gc"], preds["max_tx"]))
        runs["max_tx"].append(spearman(preds["max_tx"], meas))
        runs["gc"].append(spearman(preds["gc"], meas))
        runs["hybrid_maxtx_gc"].append(spearman(combine(preds, COMBO), meas))
        ctrl = dict(preds)
        ctrl["gc"] = shuffled(preds["gc"], t["seed"])
        runs["control_shuffled_gc"].append(spearman(combine(ctrl, COMBO), meas))
        print(f"  seed {t['seed']}: orth {orth[-1]:+.4f} | max_tx {runs['max_tx'][-1]:+.4f} | "
              f"gc {runs['gc'][-1]:+.4f} | max_tx+gc {runs['hybrid_maxtx_gc'][-1]:+.4f} | "
              f"ctrl {runs['control_shuffled_gc'][-1]:+.4f}", flush=True)

    for k, vals in runs.items():
        if vals:
            arms[k] = {"median": round(statistics.median(vals), 4), "n": len(vals),
                       "per_seed": [round(v, 4) for v in vals]}
    arms["_paired_hybrid_vs_comparator"] = paired(runs["hybrid_maxtx_gc"], runs["max_tx"])
    arms["_paired_control_vs_comparator"] = paired(runs["control_shuffled_gc"], runs["max_tx"])
    arms["_paired_gc_vs_comparator"] = paired(runs["gc"], runs["max_tx"])

    v, detail = verdict(arms)
    art = {
        "schema": "glm-tile-maxtx-gc-arm-v1",
        "date": str(date.today()),
        "substrate": "GEO GSE144621 peak tiles, leave_peak_out, ceiling 0.941",
        "arm_name": "promcalc-plus-gc-rank-average-source-diversity-150bp",
        "n_seeds": len(tables),
        "preregistered": PREREGISTERED,
        "orthogonality": {"per_seed": [round(o, 4) for o in orth],
                          "median": round(statistics.median(orth), 4),
                          "min": round(min(orth), 4), "max": round(max(orth), 4)},
        "arms": arms,
        "verdict": v,
        "verdict_detail": detail,
        "honest_limits": [
            "This does NOT discharge `orthogonal-modality-stack-promcalc-plus-N-channels` (status "
            "`planned`). That arm stacks N channels onto promcalc; this isolates ONE weak "
            "differently-sourced channel against one comparator.",
            "It does NOT overturn or revisit the G-A encoder-gate negative, which asked whether a "
            "learned encoder beats one GC number as a STANDALONE representation. Different question.",
            "One rung, one pairing, one combine rule. A negative does not show that no weak "
            "differently-sourced scorer can lift -- only that THIS one, rank-averaged with THIS "
            "scalar, does not.",
            "`max_tx` and `gc` are both computed on the same split, so their errors may be correlated "
            "through the split rather than through the biology.",
            "A positive would show source-diversity is SUFFICIENT on this rung; it would still not "
            "show it is sufficient in general, nor that the conv's content is irrelevant.",
        ],
    }
    out = Path(a.out) if a.out else Path(f"wiki/glm_tile_maxtx_gc_arm_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(art, indent=1), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(art), encoding="utf-8")
    print(f"\nVERDICT: {v}\n  {detail}\nwrote {out} + {out.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
