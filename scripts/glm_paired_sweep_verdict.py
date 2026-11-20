"""Recompute the paired-sweep verdict from the committed grid, with a corrected null gate.

WHY THIS SCRIPT EXISTS. The in-kernel reading wrote INSTRUMENT_INVALID, and that was a defect in
the gate rather than a statement about the data. The gate required EVERY cell's Markov null to sit
within 1e-6 of BAR_BOTH = 0.6679 -- but 0.6679 is the null's value AT SEED 7. The falsifier takes a
single `--seed` (scripts/glm_generator_falsifier.py:124, default 7) and threads it into
`build_natural(..., a.seed)` at line 185, which at line 57 does `random.Random(seed).shuffle(wins)`
and therefore SELECTS the natural windows and the fit/eval split, and into
`MarkovGenerator(order=..., seed=a.seed).fit(fit_set)` at line 284. So the null is refit per seed
and legitimately varies; in the committed grid it takes 10 distinct values from 0.5674 to 0.7024.

THE CORRECTED GATE anchors reproduction at seed 7 ONLY, where 0.6679 is the right constant, and
records every other seed's null as observed. Each cell carries its OWN seed-matched null, which is
what makes the model-versus-null comparison paired as well.

WHAT THE SPREAD MEASURES -- a correction to the pre-registration's wording. The bar
(wiki/glm_generator_paired_noise_bar_2026-11-20.json) calls the across-seed spread
"generation-sampling variation". Because the seed also moves the corpus split and the null fit,
the spread is really MEASUREMENT variation: generation plus corpus-split plus null-fit. The PAIRED
delta is unaffected -- within one seed both models see the same windows, the same split and the
same null, so the confound cancels inside each delta -- but the band's name was too narrow. The
rules themselves are applied exactly as frozen.

Read-only apart from its own output artifact. Run:
    uv run python scripts/glm_paired_sweep_verdict.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID = ROOT / "wiki" / "glm_generator_paired_sweep_2026-11-20.json"
BAR = ROOT / "wiki" / "glm_generator_paired_noise_bar_2026-11-20.json"
OUT = ROOT / "wiki" / "glm_generator_paired_sweep_verdict_2026-11-20.json"

#: The committed arms, mode-both, lower is better. Both ran at the falsifier's default seed 7.
COMMITTED_ZEROSHOT = 0.7787
COMMITTED_FINETUNED = 0.8058
COMMITTED_DELTA = COMMITTED_FINETUNED - COMMITTED_ZEROSHOT  # +0.0271
#: The null's value at seed 7, and ONLY at seed 7.
NULL_AT_SEED7 = 0.6679
ANCHOR_SEED = 7
TOL = 1e-6


def load_grid(path: Path = GRID) -> list[dict]:
    return json.loads(path.read_text(encoding="utf-8"))["grid"]


def cell(grid: list[dict], model: str, seed: int) -> float | None:
    for c in grid:
        if c["model"] == model and c["seed"] == seed:
            return c["both"]
    return None


def null_for(grid: list[dict], seed: int) -> float | None:
    """Each cell carries the null from its OWN falsifier run, so this is seed-matched."""
    for c in grid:
        if c["seed"] == seed and c.get("null"):
            v = c["null"].get("both")
            if v is not None:
                return v
    return None


def instrument_ok(grid: list[dict], seeds: list[int]) -> tuple[bool, str]:
    """Corrected gate: reproduction is checked at the anchor seed only."""
    n7 = null_for(grid, ANCHOR_SEED)
    if n7 is None:
        return False, "no null recorded at the anchor seed %d" % ANCHOR_SEED
    if abs(n7 - NULL_AT_SEED7) > TOL:
        return False, ("the null did not reproduce AT SEED 7: %.4f vs %.4f"
                       % (n7, NULL_AT_SEED7))
    for m in ("zeroshot", "epoch1", "epoch2", "epoch3"):
        missing = [s for s in seeds if cell(grid, m, s) is None]
        if missing:
            return False, "missing cells for %s at seeds %s" % (m, missing)
    return True, "null reproduces at seed 7 (%.4f) and the grid is complete" % n7


def delta_verdict(grid: list[dict], seeds: list[int]) -> dict:
    """The FROZEN delta rule, applied verbatim: same sign as the committed delta on every
    seed AND mean greater than the deltas' own spread."""
    d = {s: cell(grid, "epoch3", s) - cell(grid, "zeroshot", s) for s in seeds}
    v = [d[s] for s in seeds]
    mean = sum(v) / len(v)
    spread = max(v) - min(v)
    same_sign = all((x > 0) == (COMMITTED_DELTA > 0) for x in v)
    verdict = "DELTA_ROBUST" if (same_sign and mean > spread) else "DELTA_NOT_ESTABLISHED"
    return {"verdict": verdict,
            "per_seed_delta": {str(s): round(d[s], 6) for s in seeds},
            "mean_delta": round(mean, 6),
            "spread_delta": round(spread, 6),
            "n_against_committed_sign": sum(1 for x in v if (x > 0) != (COMMITTED_DELTA > 0)),
            "committed_delta": round(COMMITTED_DELTA, 6),
            "delta_at_anchor_seed": round(d[ANCHOR_SEED], 6),
            "sign_flipped_at_anchor_seed": (d[ANCHOR_SEED] > 0) != (COMMITTED_DELTA > 0),
            "spread_over_committed_effect": round(spread / abs(COMMITTED_DELTA), 2)}


def trajectory_verdict(grid: list[dict], seeds: list[int], spread: float) -> dict:
    """The FROZEN trajectory rule on paired consecutive differences.

    PRECEDENCE NOTE: the bar defines FLAT_NOT_EPOCHS ("no consecutive difference is
    sign-consistent") and TRAJECTORY_UNRESOLVED (the same, AND the delta spread exceeds the total
    epoch-1-to-epoch-3 range) without stating which wins when both hold. TRAJECTORY_UNRESOLVED is
    strictly the more specific of the two, so it is checked first; the bar's silence on precedence
    is recorded as a defect in the bar rather than resolved silently.
    """
    steps = [("zeroshot", "epoch1"), ("epoch1", "epoch2"), ("epoch2", "epoch3")]
    per_step = {}
    for a, b in steps:
        vals = [cell(grid, b, s) - cell(grid, a, s) for s in seeds]
        pos = sum(1 for x in vals if x > 0)
        per_step["%s->%s" % (a, b)] = {
            "mean": round(sum(vals) / len(vals), 6),
            "n_positive": pos, "n_negative": sum(1 for x in vals if x < 0),
            "sign_consistent": pos == len(vals) or pos == 0}
    any_consistent = any(v["sign_consistent"] for v in per_step.values())
    totals = [cell(grid, "epoch3", s) - cell(grid, "epoch1", s) for s in seeds]
    rng = max(totals) - min(totals)

    if not any_consistent and spread > rng:
        verdict = "TRAJECTORY_UNRESOLVED"
    elif not any_consistent:
        verdict = "FLAT_NOT_EPOCHS"
    elif all(per_step[k]["mean"] > 0 and per_step[k]["sign_consistent"] for k in per_step):
        verdict = "OVERFIT_MONOTONE"
    else:
        verdict = "FLAT_NOT_EPOCHS"
    return {"verdict": verdict, "per_step": per_step,
            "ep1_to_ep3_range": round(rng, 6), "delta_spread": round(spread, 6),
            "precedence_note": ("FLAT_NOT_EPOCHS and TRAJECTORY_UNRESOLVED both hold; the bar did "
                                "not state precedence, so the more specific one is reported."),
            "suggestive_but_not_a_finding": (
                "Step means rise monotonically (%s) and positive counts rise (%s), which is what "
                "progressive over-training would look like. The frozen rule says UNRESOLVED and it "
                "is NOT upgraded: its own text prescribes more draws, not a conclusion."
                % ([per_step[k]["mean"] for k in per_step],
                   [per_step[k]["n_positive"] for k in per_step]))}


def null_comparison(grid: list[dict], seeds: list[int]) -> dict:
    """Model versus its OWN seed-matched null. This is the claim the sweep does NOT overturn."""
    gaps = {}
    for s in seeds:
        z, n = cell(grid, "zeroshot", s), null_for(grid, s)
        gaps[s] = z - n
    v = list(gaps.values())
    return {"per_seed_zeroshot_minus_null": {str(s): round(gaps[s], 6) for s in seeds},
            "mean_gap": round(sum(v) / len(v), 6),
            "n_seeds_model_beats_null": sum(1 for x in v if x < 0),
            "n_seeds": len(v),
            "distinct_null_values": sorted({round(null_for(grid, s), 4) for s in seeds}),
            "reading": ("The generator is worse than its own seed-matched 3-mer Markov null at "
                        "%d of %d seeds. This claim is NOT affected by the delta result."
                        % (len(v) - sum(1 for x in v if x < 0), len(v)))}


def compute() -> dict:
    grid = load_grid()
    seeds = sorted({c["seed"] for c in grid})
    ok, why = instrument_ok(grid, seeds)
    rec = {"record": "glm-generator-paired-sweep-verdict-v1",
           "date": "2026-11-20",
           "bar": "wiki/glm_generator_paired_noise_bar_2026-11-20.json",
           "grid": "wiki/glm_generator_paired_sweep_2026-11-20.json",
           "claim_status": "DIAGNOSTIC_METHOD_VARIANT",
           "seeds": seeds,
           "instrument_ok": ok, "instrument_note": why,
           "supersedes_in_kernel_verdict": (
               "The grid artifact's own `verdict` field reads INSTRUMENT_INVALID. That was a gate "
               "defect: it required the null within 1e-6 of 0.6679 at EVERY seed, but 0.6679 is "
               "the null's value at seed 7 only, because --seed also drives the corpus split and "
               "the null refit. The grid data itself is sound and is used unmodified.")}
    if not ok:
        rec["delta"] = {"verdict": "INSTRUMENT_INVALID"}
        rec["trajectory"] = {"verdict": "INSTRUMENT_INVALID"}
        return rec
    d = delta_verdict(grid, seeds)
    rec["delta"] = d
    rec["trajectory"] = trajectory_verdict(grid, seeds, d["spread_delta"])
    rec["null_comparison"] = null_comparison(grid, seeds)
    rec["anchor"] = {
        "zeroshot_at_seed7": cell(grid, "zeroshot", ANCHOR_SEED),
        "committed_zeroshot": COMMITTED_ZEROSHOT,
        "zeroshot_reproduces_exactly":
            abs(cell(grid, "zeroshot", ANCHOR_SEED) - COMMITTED_ZEROSHOT) < TOL,
        "epoch3_at_seed7": cell(grid, "epoch3", ANCHOR_SEED),
        "committed_finetuned": COMMITTED_FINETUNED,
        "epoch3_drift": round(abs(cell(grid, "epoch3", ANCHOR_SEED) - COMMITTED_FINETUNED), 6),
        "reading": ("Zero-shot reproduces the committed number EXACTLY at the same seed, so the "
                    "scoring path is deterministic and the instrument is validated. Epoch 3 does "
                    "NOT, and cannot be expected to: it is a different training run. The loss "
                    "curves differ (7.43638 / 7.16777 / 6.91858 here vs 7.4361 / 7.1655 / 6.9148 "
                    "committed), so training is nondeterministic -- which is exactly the variance "
                    "the bar recorded this design does NOT bound.")}
    rec["what_this_does_not_bound"] = (
        "TRAINING variance, and it turned out to dominate. Every scored model comes from ONE "
        "training trajectory; the seeds vary generation and the corpus split. The epoch-3 drift "
        "from the committed arm at the SAME seed is larger than the effect under test.")
    rec["honest_limits"] = [
        "One organism (E. coli K-12), one corpus, ONE training trajectory.",
        "An upstream window is a PROXY for a promoter, not a mapped TSS.",
        "A k-mer discriminator sees LOCAL composition; the published failures are LONG-RANGE, so "
        "a low score is necessary, NOT sufficient -- and this run produced no low scores.",
        "n=120 generated sequences per cell, the committed geometry, kept for comparability.",
        "The across-seed spread is MEASUREMENT variation (generation + corpus split + null refit), "
        "not generation alone; the bar's wording was narrower than what the seed actually moves.",
    ]
    return rec


def main() -> int:
    rec = compute()
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    d, t = rec["delta"], rec["trajectory"]
    print("instrument : %s -- %s" % (rec["instrument_ok"], rec["instrument_note"]))
    print("DELTA      : %s" % d["verdict"])
    if "mean_delta" in d:
        print("  mean %+.4f | spread %.4f | %d of %d seeds against the committed sign"
              % (d["mean_delta"], d["spread_delta"], d["n_against_committed_sign"], len(rec["seeds"])))
        print("  spread is %.1fx the committed effect of %+.4f"
              % (d["spread_over_committed_effect"], d["committed_delta"]))
        print("  at seed 7: %+.4f (sign flipped vs committed: %s)"
              % (d["delta_at_anchor_seed"], d["sign_flipped_at_anchor_seed"]))
    print("TRAJECTORY : %s" % t["verdict"])
    if "null_comparison" in rec:
        nc = rec["null_comparison"]
        print("NULL       : model worse than its own seed-matched null at %d of %d seeds, mean gap %+.4f"
              % (nc["n_seeds"] - nc["n_seeds_model_beats_null"], nc["n_seeds"], nc["mean_gap"]))
    print("wrote %s" % OUT.relative_to(ROOT).as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
