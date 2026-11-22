"""Is the kmer mode's UNRESOLVED winnable region a sample-size problem, and does more data fix it?

    uv run python scripts/glm_margin_sample_size_sweep.py      # read-only; writes one artifact

**The question this answers, and why it is the one that matters.**
`scripts/glm_grid_reconcile.py` measures the winnable region against a floor in the null's own geometry
and lands on two different verdicts: `both` mode resolves (+0.1376, CI95 [0.0928, 0.1824]) while `kmer`
mode does NOT (+0.0329, CI95 [-0.0131, 0.0788], sd nearly twice the mean). An unresolved margin is not an
absent one, and the distinction is the whole of the open decision about kmer mode: if the margin is real
but merely unmeasured at the committed geometry, then the fix is more data and the mode is salvageable; if
it stays unresolved as the data grows, the mode genuinely has nothing to offer against this null.

The committed geometry uses **1200 of 4340 available promoter windows and generates 120 sequences**. There
is therefore a large amount of free headroom, and nobody has looked at what it buys.

**PRE-REGISTERED PREDICTION, written before the sweep was run.** Both terms of the margin should move in
the direction of resolvability, for different reasons:

  * The **null** should RISE. A held-out discriminator with 352 features trained on ~360 rows at the
    committed geometry is close to p=n; more rows let it actually learn the natural-versus-Markov
    difference rather than fit split noise, so the measured separation should increase toward its true
    value.
  * The **floor** should stay near 0.50 but with SHRINKING spread. It is real sequence against real
    sequence, so its true AUROC is 0.5 at every sample size; only its sampling variation falls.

So the margin's mean should rise and its SD should fall, and `kmer` should become resolvable at some
geometry if the margin is real at all. The informative outcome is the one where it does NOT: a margin that
stays unresolved while its error bar visibly shrinks is evidence the margin is genuinely near zero, and
that would settle the kmer question in the opposite direction. **Both outcomes are named here before any
number exists so neither can be retro-fitted.**

**This runs OFF the committed geometry by design, so its numbers must never be joined to the grid's.**
Changing `limit_natural` moves the Markov fit size, the discriminator's train and test sizes and the class
balance together; changing `n_generated` moves the balance again. That is the point of the sweep and it is
also why the output is a separate artifact that reports the committed cell explicitly as one row among
several, rather than an amendment to any generator figure.

Free: pure Python plus sklearn, no torch, no GPU, no network, no D: volume.
"""

from __future__ import annotations

import json
import math
import statistics
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFSEQ = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
OUT = ROOT / "wiki" / "glm_margin_sample_size_sweep_2026-11-21.json"

LENGTH = 150
KIND = "promoter"
KMER_K = 4
POSITIONAL_BINS = 6
MARKOV_ORDER = 3
MODES = ("kmer", "both")
SEEDS = tuple(range(10))

# (limit_natural, n_generated). The first row IS the committed kernel geometry; the rest spend the
# headroom. Each needs limit + n_generated spare windows for the floor, so the largest is bounded by the
# 4340 available promoters.
CONFIGS = (
    (1200, 120),     # the committed cell -- must reproduce the grid's nulls
    (1200, 400),
    (2400, 400),
    (3600, 400),
    (3600, 600),
)

PREDICTION = {
    "registered_before_any_number": True,
    "null": "RISES with more natural windows, because the committed geometry trains ~352 features on "
            "~360 rows and is close to p=n, so added rows convert fitted split noise into real signal",
    "floor": "stays near 0.50 with SHRINKING spread, because real-versus-real has a true AUROC of 0.5 at "
             "every sample size and only its sampling variation depends on n",
    "margin": "mean rises and sd falls, so kmer becomes resolvable at some geometry IF the margin is real",
    "the_informative_failure": "a kmer margin that stays unresolved while its error bar visibly shrinks "
                               "is evidence the margin is genuinely near zero, which would settle the "
                               "kmer question in the opposite direction",
}


def require_sklearn() -> None:
    import numpy  # noqa: F401
    import sklearn  # noqa: F401


def full_draw(seed: int):
    """The whole shuffled promoter set at one seed, built ONCE and sliced.

    `build_natural_set` shuffles all 4340 windows with `Random(seed)` and truncates afterwards, so
    `build(limit=L)` equals `build(limit=None)[:L]` for every L. Building once per seed instead of once
    per (seed, limit) is a 5x saving here, and the equality is ASSERTED at runtime rather than trusted --
    if it ever broke, every row of this sweep would silently describe a different corpus than its label.
    """
    from dna_decode.glm.corpus import build_natural_set

    wins, _ = build_natural_set(REFSEQ, length=LENGTH, limit=None, seed=seed, kind=KIND)
    return wins


def assert_slice_equivalence(seed: int, limit: int) -> None:
    from dna_decode.glm.corpus import build_natural_set

    direct, _ = build_natural_set(REFSEQ, length=LENGTH, limit=limit, seed=seed, kind=KIND)
    sliced = full_draw(seed)[:limit]
    if [w.seq for w in direct] != [w.seq for w in sliced]:
        raise SystemExit(f"REFUSE: slicing the full draw does not equal build(limit={limit}) at "
                         f"seed {seed}; every row of this sweep would be mislabelled.")


def cell(wins, *, limit: int, n_gen: int, seed: int) -> dict:
    """One (limit, n_generated, seed) cell: the null, the perfect-generator floor, and the margin.

    The floor uses REAL windows held out of both halves of the draw, filtered by SEQUENCE because the
    corpus contains 7 sequences duplicated across 24 windows -- an index-disjoint spare can be
    sequence-identical to an eval window and the discriminator sees only sequence.
    """
    from dna_decode.glm.baselines import MarkovGenerator
    from dna_decode.glm.corpus import fit_eval_halves
    from dna_decode.glm.falsifier import heldout_auroc

    committed = wins[:limit]
    fit, ev = fit_eval_halves(committed)
    used = set(fit) | set(ev)
    spare, skipped = [], 0
    for w in wins[limit:]:
        if len(spare) >= n_gen:
            break
        if w.seq in used:
            skipped += 1
            continue
        spare.append(w.seq)
        used.add(w.seq)
    if len(spare) < n_gen:
        raise SystemExit(f"REFUSE: only {len(spare)} sequence-disjoint spares for limit={limit} "
                         f"n_gen={n_gen} at seed {seed}")
    gen = MarkovGenerator(order=MARKOV_ORDER, seed=seed).fit(fit).generate(n_gen, LENGTH)

    out = {"seed": seed, "limit_natural": limit, "n_generated": n_gen,
           "n_fit": len(fit), "n_eval": len(ev), "n_spare_skipped_as_duplicates": skipped}
    for mode in MODES:
        kw = {"k": KMER_K, "seed": seed, "mode": mode}
        if mode == "both":
            kw["bins"] = POSITIONAL_BINS
        null = heldout_auroc(ev, gen, **kw)
        floor = heldout_auroc(ev, spare, **kw)
        out[mode] = {
            "null": round(null, 4) if null is not None else None,
            "floor": round(floor, 4) if floor is not None else None,
            "margin": round(null - floor, 4) if (null is not None and floor is not None) else None,
        }
    return out


def resolve(margins: list[float]) -> dict:
    """Mean margin with a 95% interval, and whether it clears zero. Paired by construction."""
    n = len(margins)
    if n < 3:
        return {"n": n, "resolvable": None}
    mean = statistics.fmean(margins)
    sd = statistics.stdev(margins)
    se = sd / math.sqrt(n)
    crit = 2.262 if n == 10 else 1.96          # t(0.975, df=9) at n=10
    lo, hi = mean - crit * se, mean + crit * se
    return {"n": n, "mean": round(mean, 4), "sd": round(sd, 4), "stderr": round(se, 4),
            "ci95": [round(lo, 4), round(hi, 4)],
            "t": round(mean / se, 3) if se else None,
            "resolvable_above_zero": bool(lo > 0),
            "n_positive": sum(1 for v in margins if v > 0)}


def feature_count(mode: str) -> int:
    """How many features the discriminator fits, so p-against-n is reported rather than asserted."""
    kmer = 4 ** KMER_K
    pos = POSITIONAL_BINS * (4 ** min(KMER_K, 2))
    return kmer if mode == "kmer" else kmer + pos


def main() -> int:
    if not (REFSEQ / "genome.fna").exists():
        print(f"REFUSE: no local reference at {REFSEQ}")
        return 2
    require_sklearn()
    t0 = time.time()

    assert_slice_equivalence(SEEDS[0], CONFIGS[0][0])
    print("slice equivalence: build(limit=L) == build(limit=None)[:L]")
    print(f"configs: {list(CONFIGS)}   seeds: {list(SEEDS)}")
    print(f"features: kmer {feature_count('kmer')}, both {feature_count('both')}")
    print()

    draws = {s: full_draw(s) for s in SEEDS}
    rows, by_config = [], {}
    for limit, n_gen in CONFIGS:
        cells = [cell(draws[s], limit=limit, n_gen=n_gen, seed=s) for s in SEEDS]
        rows.extend(cells)
        key = f"limit{limit}_n{n_gen}"
        entry = {"limit_natural": limit, "n_generated": n_gen,
                 "n_eval": cells[0]["n_eval"], "n_fit": cells[0]["n_fit"],
                 "train_rows_approx": (cells[0]["n_eval"] + n_gen) // 2,
                 "is_committed_geometry": (limit, n_gen) == CONFIGS[0]}
        for mode in MODES:
            nulls = [c[mode]["null"] for c in cells if c[mode]["null"] is not None]
            floors = [c[mode]["floor"] for c in cells if c[mode]["floor"] is not None]
            margins = [c[mode]["margin"] for c in cells if c[mode]["margin"] is not None]
            entry[mode] = {
                "null_mean": round(statistics.fmean(nulls), 4),
                "null_sd": round(statistics.stdev(nulls), 4),
                "floor_mean": round(statistics.fmean(floors), 4),
                "floor_sd": round(statistics.stdev(floors), 4),
                "margin": resolve(margins),
                "features_over_train_rows": round(
                    feature_count(mode) / max(1, (cells[0]["n_eval"] + n_gen) // 2), 3),
            }
        by_config[key] = entry
        for mode in MODES:
            e = entry[mode]
            print(f"  limit {limit:5d} n {n_gen:4d} [{mode:5s}] null {e['null_mean']:.4f} "
                  f"floor {e['floor_mean']:.4f} (sd {e['floor_sd']:.4f}) -> margin "
                  f"{e['margin']['mean']:+.4f} CI {e['margin']['ci95']} "
                  f"resolvable={e['margin']['resolvable_above_zero']}  p/n {e['features_over_train_rows']}")

    # ---- did the prediction hold? Checked mechanically, not narrated. ----
    first, last = list(by_config.values())[0], list(by_config.values())[-1]
    checks = {}
    for mode in MODES:
        checks[mode] = {
            "null_rose": last[mode]["null_mean"] > first[mode]["null_mean"],
            "null_change": round(last[mode]["null_mean"] - first[mode]["null_mean"], 4),
            "floor_near_half_everywhere": all(
                abs(v[mode]["floor_mean"] - 0.5) < 0.05 for v in by_config.values()),
            "floor_spread_shrank": last[mode]["floor_sd"] < first[mode]["floor_sd"],
            "floor_sd_change": round(last[mode]["floor_sd"] - first[mode]["floor_sd"], 4),
            "margin_stderr_shrank": last[mode]["margin"]["stderr"] < first[mode]["margin"]["stderr"],
            "resolvable_at_committed": first[mode]["margin"]["resolvable_above_zero"],
            "resolvable_at_largest": last[mode]["margin"]["resolvable_above_zero"],
            "became_resolvable": (not first[mode]["margin"]["resolvable_above_zero"]
                                  and last[mode]["margin"]["resolvable_above_zero"]),
            "n_configs_resolvable": sum(
                1 for v in by_config.values() if v[mode]["margin"]["resolvable_above_zero"]),
        }

    def verdict(mode: str) -> str:
        c = checks[mode]
        if c["resolvable_at_committed"] and c["resolvable_at_largest"]:
            return "RESOLVED_AT_EVERY_GEOMETRY_TESTED"
        if c["became_resolvable"]:
            return "SAMPLE_SIZE_FIXES_IT__MARGIN_IS_REAL_AND_WAS_MERELY_UNMEASURED"
        if c["n_configs_resolvable"] == 0 and c["margin_stderr_shrank"]:
            return "STAYS_UNRESOLVED_AS_PRECISION_IMPROVES__MARGIN_IS_NEAR_ZERO"
        if c["n_configs_resolvable"] == 0:
            return "UNRESOLVED_AND_PRECISION_DID_NOT_IMPROVE__SWEEP_UNINFORMATIVE"
        return "RESOLVABILITY_IS_GEOMETRY_DEPENDENT"

    rec = {
        "record": "glm-margin-sample-size-sweep-v1",
        "date": date.today().isoformat(),
        "question": ("Is the kmer mode's unresolved winnable region a sample-size problem that the "
                     "unused 3140 promoter windows would fix?"),
        "prediction": PREDICTION,
        "geometry_note": ("DELIBERATELY OFF the committed geometry. limit_natural moves the Markov fit "
                          "size, the discriminator's train and test sizes and the class balance "
                          "together; n_generated moves the balance again. These numbers must NEVER be "
                          "joined to the committed grid's generator figures. The first config IS the "
                          "committed cell and is flagged as such."),
        "fixed": {"window_length_bp": LENGTH, "natural_set_arg": KIND, "kmer_k": KMER_K,
                  "positional_bins": POSITIONAL_BINS, "markov_order": MARKOV_ORDER,
                  "feature_modes": list(MODES), "seeds": list(SEEDS),
                  "n_promoter_windows_available": len(draws[SEEDS[0]])},
        "configs": [{"limit_natural": a, "n_generated": b} for a, b in CONFIGS],
        "rows": rows,
        "by_config": by_config,
        "prediction_checks": checks,
        "verdict_by_mode": {m: verdict(m) for m in MODES},
        "honest_limits": [
            "Five geometries on a coarse grid, not a fitted curve. It shows the DIRECTION the margin "
            "moves with data; it does not identify the sample size at which a mode becomes usable.",
            "limit_natural and n_generated are not separated: each config moves fit size, train size, "
            "test size and class balance at once. A clean separation would need a factorial over fit-n "
            "and eval-n, which this does not run.",
            "The floor is EMPIRICAL at one draw per seed. Its true value is 0.5 by construction for "
            "real-versus-real, so its mean is a check on the instrument and its SD is the quantity of "
            "interest here.",
            "Ten seeds per config, so each interval is itself imprecise; the resolvability flag is a "
            "reading of an interval, not a power calculation.",
            "The null is the Markov chain, which is refitted at each limit on a larger fit half, so a "
            "rising null mixes 'the discriminator sees more' with 'the null is fitted better'. Both push "
            "in the same direction here, and the sweep does not separate them.",
            "A winnable region is necessary, NOT sufficient: these features see LOCAL composition only.",
            "One organism, one null, one window length.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print()
    for mode in MODES:
        print(f"[{mode}] {verdict(mode)}")
        c = checks[mode]
        print(f"      null {c['null_change']:+.4f} across the sweep | floor sd "
              f"{c['floor_sd_change']:+.4f} | floor ~0.5 everywhere: "
              f"{c['floor_near_half_everywhere']} | resolvable at "
              f"{c['n_configs_resolvable']}/{len(CONFIGS)} geometries")
    print(f"\n{rec['wall_seconds']}s -> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
