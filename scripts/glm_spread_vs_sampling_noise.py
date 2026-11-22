"""Is the across-seed spread a seed-role effect at all, or just the sampling noise of a held-out AUROC?

    uv run python scripts/glm_spread_vs_sampling_noise.py      # read-only; writes one artifact

**Why this has to run before any seed-role conclusion stands.** A held-out AUROC between `n1` and `n2`
items has a known sampling standard deviation under the null that the two classes are exchangeable --
the Mann-Whitney variance, `sqrt((n1 + n2 + 1) / (12 * n1 * n2))`. At the committed geometry the
discriminator's test half holds roughly 300 natural and 60 generated sequences, for which that formula
gives **0.0409**. The committed grid's ten nulls have an observed across-seed SD of **0.0410**. At
`n_generated=400` the formula gives 0.0264 and the off-geometry v1 probe measured 0.0266.

Two numbers agreeing to three decimal places twice, across a geometry change that moved them by a third,
is not something a seed-role mechanism predicts. It is what a score dominated by test-set sampling noise
looks like. If that holds up, then **no seed role carries the spread**, the only lever on the error bar is
test-set size, and a `(B + D + E) / A` ratio near 3 is simply what a noise-dominated score produces rather
than evidence of roles cancelling.

**The floor is the clean test and that is the point of using it.** The perfect-generator floor scores real
sequence against real sequence, so its true AUROC is exactly 0.5 at every seed and every sample size --
which is precisely the null the formula assumes. Its across-seed SD therefore has no real signal mixed in
and can be compared to the benchmark directly. The null's SD is a weaker test, because a true AUROC above
0.5 shifts the variance away from the H0 value.

**The one resolvable seed-role result is also the one that rests on a single arbitrary choice.** In
`glm_seed_role_variance_v2` the stream-decoupled arm A2 came in at SD ratio 1.515, CI [1.191, 1.912] in
`both` mode -- but that is one interval out of ten, in one of two modes (kmer gave 1.049), from ONE pair
of stream offsets picked for no reason beyond being round numbers. Measured against the sampling
benchmark it is A2, not the all-coupled arm, that looks anomalous. So A2 is re-run here at several offset
pairs, and the all-coupled arm is extended to many more replicates so its own SD stops being the
imprecise term in every ratio.

**Both outcomes are named before the run.** If the benchmark fits and A2 does not replicate, the
suppression was an artifact of one offset choice and the honest conclusion is that the spread is sampling
noise -- which argues against a seed-role split on the grounds that there is nothing for it to remove. If
A2 DOES replicate across offsets, the runner's reported error bars are biased low by stream sharing and a
role split becomes a correctness fix rather than a precision one, which is the opposite recommendation.

Free: pure Python plus sklearn. No torch, no GPU, no network, no D: volume.
"""

from __future__ import annotations

import json
import math
import random
import statistics
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFSEQ = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
GRID = ROOT / "wiki" / "glm_generator_paired_sweep_2026-11-20.json"
SWEEP = ROOT / "wiki" / "glm_margin_sample_size_sweep_2026-11-21.json"
OUT = ROOT / "wiki" / "glm_spread_vs_sampling_noise_2026-11-21.json"

LENGTH = 150
LIMIT = 1200
N_GENERATED = 120
KIND = "promoter"
KMER_K = 4
POSITIONAL_BINS = 6
MARKOV_ORDER = 3
MODES = ("kmer", "both")

N_MAIN = 240              # replicates for the all-coupled arm, so ITS sd stops being the loose term
N_OFFSET = 60             # replicates per stream-offset variant
OFFSET_PAIRS = ((10_000, 20_000),      # the pair v2 used -- must be included to test replication
                (1_000, 2_000),
                (123_457, 987_653),
                (7, 13),
                (500_000, 900_000))
N_BOOTSTRAP = 2000

PREDICTION = {
    "registered_before_any_number": True,
    "benchmark": ("a held-out AUROC's null sampling sd is sqrt((n1+n2+1)/(12*n1*n2)); at the committed "
                  "geometry that is about 0.0409 and the grid's observed across-seed sd is 0.0410"),
    "if_the_benchmark_fits_the_floor": ("the across-seed spread is test-set sampling noise, no seed role "
                                       "carries it, and the only lever on the error bar is test-set "
                                       "size. A (B+D+E)/A ratio near 3 is then what a noise-dominated "
                                       "score gives, NOT role cancellation."),
    "if_A2_replicates_across_offsets": ("stream sharing really does suppress the reported spread, the "
                                        "runner's error bars are biased LOW, and a seed-role split "
                                        "becomes a correctness fix -- the opposite recommendation"),
    "if_A2_does_not_replicate": ("the single resolvable seed-role result was an artifact of one "
                                 "arbitrary offset choice and must be withdrawn"),
}


def require_sklearn() -> None:
    import numpy  # noqa: F401
    import sklearn  # noqa: F401


def mann_whitney_null_sd(n1: int, n2: int) -> float:
    """Sampling sd of an AUROC between `n1` and `n2` items when the two classes are exchangeable.

    The standard Mann-Whitney U null variance. It is the right benchmark for the perfect-generator floor,
    whose true AUROC is 0.5 by construction, and only an approximation for the Markov null, whose true
    AUROC is above 0.5.
    """
    if n1 < 1 or n2 < 1:
        return float("nan")
    return math.sqrt((n1 + n2 + 1) / (12.0 * n1 * n2))


def expected_test_split(n_natural: int, n_generated: int) -> tuple[int, int]:
    """How many of each class land in the held-out half.

    `heldout_auroc` shuffles the combined list and trains on the first half, so the test half is about
    half of everything and keeps the overall class proportion in expectation.
    """
    total = n_natural + n_generated
    test = total - total // 2
    frac = n_natural / total
    n1 = int(round(test * frac))
    return n1, test - n1


def sd_ci(values: list[float], *, n: int = N_BOOTSTRAP, seed: int = 999) -> dict:
    """Bootstrap interval for one sample's sd, so an observed sd can be compared to a benchmark fairly.

    v1 and v2 both compared sd values without an interval on the sd ITSELF; at ten replicates the sd of
    an sd estimate is roughly 24% of its value, which is wide enough to swallow most of these effects.
    """
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        b = [rng.choice(values) for _ in values]
        try:
            out.append(statistics.stdev(b))
        except statistics.StatisticsError:
            continue
    out.sort()
    lo = out[int(0.025 * len(out))]
    hi = out[min(len(out) - 1, int(0.975 * len(out)))]
    return {"sd": round(statistics.stdev(values), 4), "n": len(values),
            "ci95": [round(lo, 4), round(hi, 4)]}


def benchmark_against_sweep() -> dict:
    """Free check: does the formula predict the FLOOR's observed spread at all five swept geometries?

    The floor is the clean case -- real against real, true AUROC 0.5 -- so this is the strongest available
    test of whether the across-seed spread is just sampling noise, and it costs nothing because the sweep
    already recorded every number.
    """
    if not SWEEP.exists():
        return {"status": "sweep artifact absent"}
    rec = json.loads(SWEEP.read_text(encoding="utf-8"))
    rows = {}
    for key, cfg in rec["by_config"].items():
        n1, n2 = expected_test_split(cfg["n_eval"], cfg["n_generated"])
        pred = mann_whitney_null_sd(n1, n2)
        entry = {"n_eval": cfg["n_eval"], "n_generated": cfg["n_generated"],
                 "test_n_natural": n1, "test_n_generated": n2,
                 "predicted_sd": round(pred, 4)}
        for mode in MODES:
            obs_floor = cfg[mode]["floor_sd"]
            obs_null = cfg[mode]["null_sd"]
            entry[mode] = {
                "observed_floor_sd": obs_floor,
                "floor_over_predicted": round(obs_floor / pred, 2),
                "observed_null_sd": obs_null,
                "null_over_predicted": round(obs_null / pred, 2),
            }
        rows[key] = entry
    ratios = [rows[k][m]["floor_over_predicted"] for k in rows for m in MODES]
    return {"by_config": rows,
            "floor_over_predicted_min": min(ratios),
            "floor_over_predicted_max": max(ratios),
            "floor_over_predicted_mean": round(statistics.fmean(ratios), 2),
            "note": ("the floor is the clean test: real against real has a true AUROC of 0.5, which is "
                     "exactly the null the benchmark assumes. The null's ratio is reported too but is a "
                     "weaker test because a true AUROC above 0.5 moves the variance."),
            "caveat": ("each observed sd comes from 10 seeds, where the sd of an sd estimate is roughly "
                       "24% of its value, so a ratio inside about 0.76-1.24 is not a mismatch")}


def all_grid_pairs(grid: dict) -> dict:
    """EVERY model-vs-null and model-vs-model pair, not a chosen three.

    Stating a pattern from three hand-picked pairs is how a sign becomes a mechanism. The paired-delta
    benefit is reported per pair with a Fisher-z interval on the correlation, and the pairs are grouped by
    KIND (against the null, or against another checkpoint) only after all of them are computed.
    """
    series: dict[str, dict[int, dict[str, float]]] = {}
    for c in grid["grid"]:
        series.setdefault(c["model"], {})[c["seed"]] = {m: c[m] for m in MODES}
        series.setdefault("null", {})[c["seed"]] = {m: c["null"][m] for m in MODES}
    names = sorted(series)
    out = {}
    for mode in MODES:
        res = {}
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                seeds = sorted(set(series[a]) & set(series[b]))
                if len(seeds) < 4:
                    continue
                xa = [series[a][s][mode] for s in seeds]
                xb = [series[b][s][mode] for s in seeds]
                d = [p - q for p, q in zip(xa, xb)]
                sa, sb = statistics.stdev(xa), statistics.stdev(xb)
                unp = math.sqrt(sa ** 2 + sb ** 2)
                r = statistics.correlation(xa, xb)
                z = 0.5 * math.log((1 + r) / (1 - r)) if abs(r) < 1 else None
                se = 1.0 / math.sqrt(len(seeds) - 3)
                ci = ([round(math.tanh(z - 1.96 * se), 3), round(math.tanh(z + 1.96 * se), 3)]
                      if z is not None else None)
                res[f"{a}__vs__{b}"] = {
                    "kind": "vs_null" if "null" in (a, b) else "vs_checkpoint",
                    "n_seeds": len(seeds),
                    "mean_delta": round(statistics.fmean(d), 4),
                    "sd_paired_delta": round(statistics.stdev(d), 4),
                    "sd_unpaired_sum": round(unp, 4),
                    "paired_over_unpaired": round(statistics.stdev(d) / unp, 3) if unp else None,
                    "correlation": round(r, 3),
                    "correlation_ci95": ci,
                    "pairing_benefit_resolvable": bool(statistics.stdev(d) < unp and ci and ci[0] > 0),
                }
        grouped = {}
        for kind in ("vs_null", "vs_checkpoint"):
            sel = [v for v in res.values() if v["kind"] == kind]
            grouped[kind] = {
                "n_pairs": len(sel),
                "n_with_resolvable_benefit": sum(1 for v in sel if v["pairing_benefit_resolvable"]),
                "median_ratio": round(statistics.median([v["paired_over_unpaired"] for v in sel]), 3),
                "ratios": sorted(v["paired_over_unpaired"] for v in sel),
            }
        out[mode] = {"pairs": res, "grouped": grouped}
    return out


def main() -> int:
    if not (REFSEQ / "genome.fna").exists():
        print(f"REFUSE: no local reference at {REFSEQ}")
        return 2
    if not GRID.exists():
        print(f"REFUSE: grid absent at {GRID}")
        return 2
    require_sklearn()
    t0 = time.time()

    from dna_decode.glm.baselines import MarkovGenerator
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves
    from dna_decode.glm.falsifier import heldout_auroc

    n1, n2 = expected_test_split(LIMIT // 2, N_GENERATED)
    pred = mann_whitney_null_sd(n1, n2)
    print(f"committed geometry: test half ~{n1} natural vs {n2} generated "
          f"-> benchmark sd {pred:.4f}")

    bench = benchmark_against_sweep()
    if "by_config" in bench:
        print("\nbenchmark against the swept FLOOR spreads (true AUROC is 0.5 there):")
        for key, row in bench["by_config"].items():
            print(f"  {key:16s} predicted {row['predicted_sd']:.4f} | "
                  f"kmer floor {row['kmer']['observed_floor_sd']:.4f} "
                  f"({row['kmer']['floor_over_predicted']}x) | "
                  f"both floor {row['both']['observed_floor_sd']:.4f} "
                  f"({row['both']['floor_over_predicted']}x)")
        print(f"  floor/predicted across all ten cells: "
              f"{bench['floor_over_predicted_min']}-{bench['floor_over_predicted_max']} "
              f"(mean {bench['floor_over_predicted_mean']})")

    # ---- arm A at many replicates, and A2 at several offset pairs ----
    cache: dict[int, tuple[list[str], list[str]]] = {}

    def corpus(seed: int):
        if seed not in cache:
            wins, _ = build_natural_set(REFSEQ, length=LENGTH, limit=LIMIT, seed=seed, kind=KIND)
            cache[seed] = fit_eval_halves(wins)
        return cache[seed]

    def score(corpus_seed: int, gen_seed: int, split_seed: int) -> dict[str, float | None]:
        fit, ev = corpus(corpus_seed)
        gen = MarkovGenerator(order=MARKOV_ORDER, seed=gen_seed).fit(fit).generate(N_GENERATED, LENGTH)
        o = {}
        for mode in MODES:
            kw = {"k": KMER_K, "seed": split_seed, "mode": mode}
            if mode == "both":
                kw["bins"] = POSITIONAL_BINS
            o[mode] = heldout_auroc(ev, gen, **kw)
        return o

    print(f"\narm A at {N_MAIN} replicates...")
    a_vals = {m: [] for m in MODES}
    for s in range(N_MAIN):
        r = score(s, s, s)
        for m in MODES:
            if r[m] is not None:
                a_vals[m].append(r[m])
    a_stat = {m: sd_ci(a_vals[m]) for m in MODES}
    for m in MODES:
        print(f"  A [{m:5s}] sd {a_stat[m]['sd']:.4f} CI {a_stat[m]['ci95']} "
              f"(benchmark {pred:.4f})")

    print(f"\narm A2 at {len(OFFSET_PAIRS)} offset pairs, {N_OFFSET} replicates each...")
    a2 = {}
    for og, osp in OFFSET_PAIRS:
        vals = {m: [] for m in MODES}
        for s in range(N_OFFSET):
            r = score(s, s + og, s + osp)
            for m in MODES:
                if r[m] is not None:
                    vals[m].append(r[m])
        key = f"gen+{og}_split+{osp}"
        a2[key] = {m: sd_ci(vals[m]) for m in MODES}
        a2[key]["sorted_values_both"] = sorted(round(v, 4) for v in vals["both"])
        for m in MODES:
            print(f"  A2 {key:28s} [{m:5s}] sd {a2[key][m]['sd']:.4f} "
                  f"CI {a2[key][m]['ci95']}")

    # ---- did A2's inflation replicate? Compared against arm A at N_MAIN replicates. ----
    replication = {}
    for m in MODES:
        sd_a = a_stat[m]["sd"]
        hits = []
        for key, v in a2.items():
            lo, hi = v[m]["ci95"]
            hits.append({"offsets": key, "sd": v[m]["sd"],
                         "ci95": v[m]["ci95"],
                         "above_A_resolvably": bool(lo > sd_a),
                         "ratio_to_A": round(v[m]["sd"] / sd_a, 3)})
        n_above = sum(1 for h in hits if h["above_A_resolvably"])
        replication[m] = {
            "sd_A": sd_a, "sd_A_ci95": a_stat[m]["ci95"], "n_A": a_stat[m]["n"],
            "benchmark_sd": round(pred, 4),
            "A_consistent_with_benchmark": bool(a_stat[m]["ci95"][0] <= pred <= a_stat[m]["ci95"][1]),
            "by_offset": hits,
            "n_offsets_resolvably_above_A": n_above,
            "n_offsets": len(hits),
            "replicates": n_above >= len(hits) - 1,
        }

    def verdict(m: str) -> str:
        r = replication[m]
        if r["replicates"]:
            return "STREAM_SHARING_SUPPRESSES_THE_SPREAD__ERROR_BARS_BIASED_LOW"
        if r["n_offsets_resolvably_above_A"] == 0:
            return "NO_STREAM_EFFECT__V2_RESULT_WAS_ONE_OFFSET_CHOICE"
        return "STREAM_EFFECT_IS_OFFSET_DEPENDENT__NOT_A_MECHANISM"

    rec = {
        "record": "glm-spread-vs-sampling-noise-v1",
        "date": date.today().isoformat(),
        "question": ("Is the across-seed spread a seed-role effect, or the sampling noise of a held-out "
                     "AUROC at this test-set size and class balance?"),
        "prediction": PREDICTION,
        "geometry": {"window_length_bp": LENGTH, "limit_natural": LIMIT,
                     "n_generated_per_candidate": N_GENERATED, "markov_order": MARKOV_ORDER,
                     "kmer_k": KMER_K, "positional_bins": POSITIONAL_BINS,
                     "test_half_n_natural": n1, "test_half_n_generated": n2},
        "benchmark": {
            "formula": "sqrt((n1 + n2 + 1) / (12 * n1 * n2))  -- Mann-Whitney U null variance",
            "committed_geometry_sd": round(pred, 4),
            "grid_observed_sd_both": 0.0410,
            "grid_observed_sd_kmer": 0.0493,
            "v1_offgeometry_n400_observed_kmer": 0.0266,
            "v1_offgeometry_n400_predicted": round(
                mann_whitney_null_sd(*expected_test_split(LIMIT // 2, 400)), 4),
        },
        "benchmark_vs_sweep_floors": bench,
        "arm_A_many_replicates": a_stat,
        "arm_A2_by_offset_pair": a2,
        "replication": replication,
        "verdict_by_mode": {m: verdict(m) for m in MODES},
        "all_grid_pairs": all_grid_pairs(json.loads(GRID.read_text(encoding="utf-8"))),
        "honest_limits": [
            "The Mann-Whitney benchmark is the null variance for exchangeable classes. It applies exactly "
            "to the perfect-generator floor, whose true AUROC is 0.5 by construction, and only "
            "approximately to the Markov null, whose true AUROC is above 0.5.",
            "It also assumes the discriminator's score is a fixed function. Here the classifier is "
            "refitted per split, which adds variance the formula does not model, so the benchmark is a "
            "LOWER bound on what to expect rather than an exact prediction.",
            "The expected test-half class counts are used; the actual split is random, so n1 and n2 vary "
            "slightly between replicates.",
            "The swept floor spreads come from 10 seeds each, where the sd of an sd estimate is roughly "
            "24% of its value. Ratios inside about 0.76-1.24 are not mismatches, and that band is wide.",
            "Arm A2 is tested at five offset pairs; offsets are an unbounded space and five is not a "
            "proof that none behaves differently.",
            "Measured on the MARKOV candidate at one geometry, one organism, one null.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print()
    for m in MODES:
        r = replication[m]
        print(f"[{m}] {verdict(m)}")
        print(f"      A sd {r['sd_A']} CI {r['sd_A_ci95']} (n={r['n_A']}) | benchmark "
              f"{r['benchmark_sd']} | A consistent with benchmark: "
              f"{r['A_consistent_with_benchmark']}")
        print(f"      A2 resolvably above A at {r['n_offsets_resolvably_above_A']}/{r['n_offsets']} "
              f"offset pairs")
    print()
    for m in MODES:
        g = rec["all_grid_pairs"][m]["grouped"]
        print(f"[{m}] pairing benefit resolvable: vs_null "
              f"{g['vs_null']['n_with_resolvable_benefit']}/{g['vs_null']['n_pairs']} "
              f"(median ratio {g['vs_null']['median_ratio']}) | vs_checkpoint "
              f"{g['vs_checkpoint']['n_with_resolvable_benefit']}/{g['vs_checkpoint']['n_pairs']} "
              f"(median ratio {g['vs_checkpoint']['median_ratio']})")
    print(f"\n{rec['wall_seconds']}s -> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
