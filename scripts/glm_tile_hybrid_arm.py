"""Does combining the LEARNED encoder with the BIOPHYSICAL scalar beat the scalar alone?

    uv run python scripts/glm_tile_hybrid_arm.py [--seeds 10]

**The one untried combination on this rung.** The Promoter Calculator was measured AGAINST the conv
encoder (+0.0401 over it, 4/4 seeds) and the 55-feature decomposition was measured against the scalar
(`NO_GAIN`, paired -0.0040, CI95 includes zero). Nobody combined the learned and the biophysical
scorer. That gap is worth one run because this repo has a MEASURED rule about when such a combine pays:

> a naive rank-average of ORTHOGONAL modalities beats the strongest single modality on 84-90% of
> ProteinGym proteins -- `ESM2+GEMME` lifts **even though GEMME alone LOSES to ESM2**.

So "the conv loses to max_tx" is not an argument against trying. **Orthogonality is.**

**And the precondition is measured FIRST, not assumed, because the failure mode is also measured.** On
CTX-M-14 the same rank-average LOST -- ESM2+ProSST 0.204 against ESM2's 0.352 -- because ProSST alone
sat at chance (-0.040), so averaging it in halved the signal. The conv sits at 0.3234 against GC's
0.3112 with a paired CI that includes zero. A negative here is therefore a legitimate, expected-possible
outcome and is recorded as a measured negative (R6/R7), never as a deleted arm.

**No fitting, no calibrator, no second-level model.** Both prediction vectors are already out-of-sample
on each seed's held-out peaks, and `rank_average_hybrid` has no parameters, so the combine introduces
no leakage. It RANKS; it does not dose.

**The control is what makes a positive interpretable.** `max_tx` rank-averaged with a SHUFFLED conv
vector: if that lifts too, the gain came from rank-averaging two vectors rather than from the conv's
content, and the headline is an artifact.
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

#: FROZEN before any hybrid number was read, MIRRORING the decomposition arm's bar so the two
#: "does more information help?" results are directly comparable.
PREREGISTERED = {
    "question": "does rank-averaging the learned conv encoder with the biophysical max_tx scalar beat "
                "the scalar alone on the tile rung?",
    "primary_metric": "spearman",
    "comparator": "max_tx",
    "margin": 0.02,
    "min_seeds_winning": 8,
    "n_seeds_expected": 10,
    "mirrors": "scripts/promcalc_decomposition_arm.py PREREGISTERED (margin 0.02, 8/10 seeds)",
    "precondition": {
        "name": "orthogonality",
        "rule": "a naive rank-average lifts only when the modalities are ORTHOGONAL; a near-chance "
                "second modality REDUCES the signal (measured: ESM2+ProSST 0.204 vs ESM2 0.352 on "
                "CTX-M-14, because ProSST alone was at chance)",
        "measured_first": True,
        "not_a_gate": "the orthogonality number is REPORTED, not used to skip the run -- a prediction "
                      "that a combine will fail is not a substitute for measuring it",
    },
    "control": "max_tx rank-averaged with a SHUFFLED conv vector; if this lifts, the gain is an "
               "artifact of rank-averaging rather than of the conv's content",
    "verdicts": {
        "HYBRID_ADDS": "paired median >= +0.02 over max_tx AND >= 8/10 seeds win AND control clean",
        "WEAK": "gain > 0 but below the margin or the seed count",
        "NO_GAIN": "no positive paired gain -- the learned encoder carries nothing orthogonal here",
        "INDETERMINATE_CONTROL_ALSO_LIFTS": "the shuffled-conv control lifts too; the combine, not "
                                            "the conv, is producing the number",
        "INDETERMINATE_RECONCILE_FAILED": "the prediction table's published anchors drift",
        "INDETERMINATE_MISSING_SCORER": "conv or max_tx absent from the table",
    },
}

COMBOS = {
    "hybrid_conv_maxtx": ("conv", "max_tx"),
    "hybrid_conv_maxtx_gc": ("conv", "max_tx", "gc"),
}


def as_table(pred: list[float]) -> dict[str, float]:
    """A prediction vector as the `{key: score}` table `rank_average_hybrid` consumes.

    Keys are the POSITION in this seed's test vector as a string. Position, not tile id: every vector
    in one seed's table is already aligned to the same `tile_indices`, and using a positional key makes
    a length mismatch impossible to paper over -- `rank_average_hybrid` raises on an empty intersection.
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
    random.Random(1000 + seed).shuffle(out)
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
        return ("INDETERMINATE_RECONCILE_FAILED",
                "; ".join(arms["_reconcile_problems"]))
    p = arms["_paired_hybrid_vs_comparator"]
    ctrl = arms.get("_paired_control_vs_comparator") or {}
    ctrl_lifts = (ctrl.get("median") or 0.0) >= pre["margin"] and \
                 (ctrl.get("wins") or 0) >= pre["min_seeds_winning"]
    if ctrl_lifts:
        return ("INDETERMINATE_CONTROL_ALSO_LIFTS",
                f"the shuffled-conv control gains {ctrl['median']:+.4f} ({ctrl['wins']}/{ctrl['n']}) "
                f"-- indistinguishable from the real combine's {p['median']:+.4f}; the number is an "
                f"artifact of rank-averaging, not of the conv")
    if p["median"] >= pre["margin"] and p["wins"] >= pre["min_seeds_winning"]:
        return ("HYBRID_ADDS",
                f"hybrid {p['median']:+.4f} over {pre['comparator']}, {p['wins']}/{p['n']} seeds, "
                f"CI95 {p['ci95']}")
    if p["median"] > 0:
        return ("WEAK",
                f"hybrid {p['median']:+.4f} over {pre['comparator']}, {p['wins']}/{p['n']} seeds, "
                f"CI95 {p['ci95']} -- positive but short of the {pre['margin']:+.2f} margin or the "
                f"{pre['min_seeds_winning']}-seed bar")
    return ("NO_GAIN",
            f"hybrid {p['median']:+.4f} over {pre['comparator']}, {p['wins']}/{p['n']} seeds, "
            f"CI95 {p['ci95']} -- the learned encoder carries nothing orthogonal the scalar lacks")


def render_md(art: dict) -> str:
    a = art["arms"]
    L = [f"# Tile-rung hybrid arm (conv ⊕ max_tx) — {art['date']}", "",
         f"**VERDICT: `{art['verdict']}`**", "", art["verdict_detail"], "",
         "## The precondition, measured first", "",
         f"Spearman between the conv and `max_tx` prediction vectors, per seed: "
         f"median **{art['orthogonality']['median']:.4f}** "
         f"(min {art['orthogonality']['min']:.4f}, max {art['orthogonality']['max']:.4f}).", "",
         "Near 1 would mean the two scorers are redundant and a rank-average cannot add; near 0 means "
         "they are orthogonal, which is the condition under which the combine is measured to pay.", "",
         "## Arms", "", "| arm | median spearman | per-seed |", "|---|---|---|"]
    for name in sorted(k for k in a if not k.startswith("_")):
        arm = a[name]
        L.append(f"| `{name}` | {arm['median']:+.4f} | "
                 + ", ".join(f"{v:+.3f}" for v in arm["per_seed"]) + " |")
    L += ["", "## Paired margins (same splits, per seed)", ""]
    for key, label in (("_paired_hybrid_vs_comparator", "hybrid vs max_tx"),
                       ("_paired_hybrid3_vs_comparator", "hybrid+gc vs max_tx"),
                       ("_paired_control_vs_comparator", "CONTROL (shuffled conv) vs max_tx")):
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
    paths = sorted(in_dir.glob("seed-*.json.gz"))[: a.seeds]
    if not paths:
        raise SystemExit(f"no seed tables in {in_dir}; run scripts/glm_tile_predictions.py first")
    tables = [load_table(p) for p in paths]
    require_real(tables)

    arms: dict = {}
    recon: list[str] = []
    for t in tables:
        recon += reconcile_problems(t["seed"], t["spearman_recomputed"])
    if recon:
        arms["_reconcile_problems"] = recon

    missing = [s for s in ("conv", "max_tx") if any(s not in t["predictions"] for t in tables)]
    if missing:
        arms["_missing_scorer"] = missing
        v, detail = verdict(arms)
        print(f"VERDICT: {v}\n  {detail}")
        return 1

    runs: dict[str, list[float]] = {k: [] for k in
                                    (*COMBOS, "max_tx", "conv", "control_shuffled_conv")}
    orth: list[float] = []
    for t in tables:
        preds, meas = t["predictions"], t["measured"]
        orth.append(orthogonality(preds["conv"], preds["max_tx"]))
        runs["max_tx"].append(spearman(preds["max_tx"], meas))
        runs["conv"].append(spearman(preds["conv"], meas))
        for name, names in COMBOS.items():
            if all(n in preds for n in names):
                runs[name].append(spearman(combine(preds, names), meas))
        ctrl = dict(preds)
        ctrl["conv"] = shuffled(preds["conv"], t["seed"])
        runs["control_shuffled_conv"].append(
            spearman(combine(ctrl, ("conv", "max_tx")), meas))
        print(f"  seed {t['seed']}: orth {orth[-1]:+.4f} | max_tx {runs['max_tx'][-1]:+.4f} | "
              f"conv {runs['conv'][-1]:+.4f} | hybrid {runs['hybrid_conv_maxtx'][-1]:+.4f} | "
              f"ctrl {runs['control_shuffled_conv'][-1]:+.4f}", flush=True)

    for k, vals in runs.items():
        if vals:
            arms[k] = {"median": round(statistics.median(vals), 4), "n": len(vals),
                       "per_seed": [round(v, 4) for v in vals]}
    arms["_paired_hybrid_vs_comparator"] = paired(runs["hybrid_conv_maxtx"], runs["max_tx"])
    if runs["hybrid_conv_maxtx_gc"]:
        arms["_paired_hybrid3_vs_comparator"] = paired(runs["hybrid_conv_maxtx_gc"], runs["max_tx"])
    arms["_paired_control_vs_comparator"] = paired(runs["control_shuffled_conv"], runs["max_tx"])
    arms["_paired_conv_vs_comparator"] = paired(runs["conv"], runs["max_tx"])

    v, detail = verdict(arms)
    art = {
        "schema": "glm-tile-hybrid-arm-v1",
        "date": str(date.today()),
        "substrate": "GEO GSE144621 peak tiles, leave_peak_out, ceiling 0.941",
        "arm_name": "promcalc-plus-conv-rank-average-hybrid-150bp",
        "n_seeds": len(tables),
        "preregistered": PREREGISTERED,
        "orthogonality": {"per_seed": [round(o, 4) for o in orth],
                          "median": round(statistics.median(orth), 4),
                          "min": round(min(orth), 4), "max": round(max(orth), 4)},
        "arms": arms,
        "verdict": v,
        "verdict_detail": detail,
        "honest_limits": [
            "The combine is label-free and parameter-free (a rank-average), so it RANKS and does not "
            "dose -- the same deployability class as the inverse cell.",
            "Both inputs are out-of-sample on each seed's held-out peaks, so no second-level leakage "
            "is introduced; but the pool is MEASURED peak tiles, not generator output.",
            "One rung, one experiment, one combine rule. A negative does not show that no learned "
            "representation can add here -- only that THIS encoder, rank-averaged with THIS scalar, "
            "does not.",
            "`max_tx` and `conv` are both fitted by ridge on the same split, so their errors may be "
            "correlated through the split rather than through the biology.",
        ],
    }
    out = Path(a.out) if a.out else Path(f"wiki/glm_tile_hybrid_arm_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(art, indent=1), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(art), encoding="utf-8")
    print(f"\nVERDICT: {v}\n  {detail}\nwrote {out} + {out.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
