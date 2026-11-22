"""Decompose the across-seed variance AT THE COMMITTED GEOMETRY, with the mechanism arms v1 lacked.

    uv run python scripts/glm_seed_role_variance_v2.py      # read-only; writes one artifact

**Why v2 exists.** `scripts/glm_seed_role_variance.py` answered the right question with the wrong
geometry. It hardcoded `N_GENERATED = 400`; the committed sweep kernel froze `n = 120`
(`scripts/kaggle/glm_paired_sweep_kernel.py:45`). Its all-coupled arm was documented as "reproduces the
status quo" and was never reconciled against the status quo, so nothing caught it: at n=400 that arm
reports `both` sd 0.0140, while the grid's own ten nulls have sd **0.0410**. Every arm in v1 is therefore
off-geometry, and because the comparison between arms was the entire finding, the finding does not survive
the correction. `scripts/glm_grid_reconcile.py` establishes that at n=120 the hand-written null path is
byte-identical to the runner at 20 of 20 seed-mode cells, and this script GATES on that same reproduction
before reporting anything.

**Why v1's mechanism claim was unsupported even on its own numbers.** v1 explained a coupled variance
below the sum of the decoupled ones by saying the Markov null is fitted on the fit half of the very corpus
draw it is scored against. That IS true of the runner -- and it is true in EVERY v1 arm, because every arm
took `fit, ev` from one `get_corpus(corpus_seed)` call. A property held constant across all arms cannot
explain a difference between them. Testing it needs an arm where the fit corpus and the eval corpus are
different draws, which is arm F here.

**Seven arms. The plan is a pure function of (arm, seed) returning four roles, so the whole design is one
table and is testable without running a discriminator.**

    arm                      fit-corpus   eval-corpus   generation   split
    A_all_coupled            s            s             s            s      <- the runner. GATED.
    A2_streams_decoupled     s            s             s+10000      s+20000
    B_generation_only        anchor       anchor        s            anchor
    C_corpus_and_split       s            s             anchor       s
    D_corpus_only            s            s             anchor       anchor
    E_split_only             anchor       anchor        anchor       s
    F_fit_decoupled          s+500        s             s            s

A2 is the rival v1 never tested. In arm A one integer seeds three separate `random.Random` instances, so
the corpus shuffle, the generator and the discriminator split all consume the SAME pseudo-random stream.
Any cancellation that produces could be an artifact of re-using one integer rather than a statistical
relationship between the roles. A2 keeps one logical seed per replicate while giving the three roles
independent streams; if the suppression is real it should survive, and if it is a stream artifact it
should vanish.

D and E split what v1's arm C fused. v1 moved corpus and split together and then reported a decomposition,
so it could not say which of the two carried the variance. With D and E the additivity of corpus against
split is checkable in its own right (`Var(C)` against `Var(D) + Var(E)`).

**The claim is about a paired DELTA, and v1 measured no delta.** The precision that matters is the
precision of `score(candidate) - score(comparator)`, because that is what an arm set is ranked on.
`Var(A - B) = Var(A) + Var(B) - 2Cov(A, B)`, so a design can leave each arm's own spread untouched and
still tighten or loosen the delta. The delta experiment here scores TWO candidates per replicate (Markov
order 3 against order 1, a materially worse null) under three designs -- shared-varying, shared-pinned,
and independent -- which is the direct test of whether sharing a seed buys precision ON THE DELTA.

Free, pure Python, seconds per replicate: no torch, no GPU, no model download, no D: volume. sklearn is
asserted on entry because `heldout_auroc` catches ImportError and returns None, which is also its
"insufficient data" return -- a missing dependency must never be able to pose as a measurement.
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
OUT = ROOT / "wiki" / "glm_seed_role_variance_v2_2026-11-21.json"

# THE KERNEL'S FROZEN GEOMETRY (scripts/kaggle/glm_paired_sweep_kernel.py:45).
LENGTH = 150
LIMIT = 1200
N_GENERATED = 120        # v1 used 400. That single value is why v1 does not reproduce the runner.
KIND = "promoter"
KMER_K = 4
POSITIONAL_BINS = 6
MARKOV_ORDER = 3
DELTA_RIVAL_ORDER = 1    # a materially worse null, so the delta is not pure noise
MODES = ("kmer", "both")

ANCHOR_SEED = 7
N_REPLICATES = 60        # v1 used 10; a ratio of two SDs at n=10 cannot separate these hypotheses
STREAM_OFFSET_GEN = 10_000
STREAM_OFFSET_SPLIT = 20_000
FIT_DECOUPLE_OFFSET = 500
N_BOOTSTRAP = 2000
GATE_TOL = 5e-4

ARMS = ("A_all_coupled", "A2_streams_decoupled", "B_generation_only", "C_corpus_and_split",
        "D_corpus_only", "E_split_only", "F_fit_decoupled")


def seed_plan(arm: str, seed: int, anchor: int = ANCHOR_SEED) -> tuple[int, int, int, int]:
    """`(fit_corpus_seed, eval_corpus_seed, gen_seed, split_seed)` for one replicate of one arm.

    Pure, so the experiment's entire design is checkable without a discriminator. v1 kept this isolated
    for the same reason and it is worth restating: two arms that silently shared a plan would produce two
    identical SDs and a decomposition that looked perfectly reasonable.
    """
    a = anchor
    if arm == "A_all_coupled":
        return seed, seed, seed, seed
    if arm == "A2_streams_decoupled":
        return seed, seed, seed + STREAM_OFFSET_GEN, seed + STREAM_OFFSET_SPLIT
    if arm == "B_generation_only":
        return a, a, seed, a
    if arm == "C_corpus_and_split":
        return seed, seed, a, seed
    if arm == "D_corpus_only":
        return seed, seed, a, a
    if arm == "E_split_only":
        return a, a, a, seed
    if arm == "F_fit_decoupled":
        return seed + FIT_DECOUPLE_OFFSET, seed, seed, seed
    raise ValueError(f"unknown arm {arm!r}")


def arms_are_distinct(seeds: list[int], anchor: int = ANCHOR_SEED) -> bool:
    """Checks the PLANS, not the outputs: identical plans give identical numbers and still look fine."""
    plans = {arm: tuple(seed_plan(arm, s, anchor) for s in seeds) for arm in ARMS}
    return len(set(plans.values())) == len(ARMS)


def arm_varies_only(arm: str, seeds: list[int], anchor: int = ANCHOR_SEED) -> set[str]:
    """Which roles actually MOVE in this arm. A role named in the arm's title must appear here.

    This is the non-vacuity check for the design: an arm called `split_only` whose split seed happened to
    be pinned would report a spread caused by something else entirely, under a name asserting otherwise.
    """
    names = ("fit_corpus", "eval_corpus", "generation", "split")
    plans = [seed_plan(arm, s, anchor) for s in seeds]
    return {n for i, n in enumerate(names) if len({p[i] for p in plans}) > 1}


def require_sklearn() -> None:
    import numpy  # noqa: F401
    import sklearn  # noqa: F401


def _corpus_cache():
    """Build each distinct corpus draw once; many arms reuse the anchor draw."""
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves

    cache: dict[int, tuple[list[str], list[str]]] = {}

    def get(corpus_seed: int) -> tuple[list[str], list[str]]:
        if corpus_seed not in cache:
            wins, _ = build_natural_set(REFSEQ, length=LENGTH, limit=LIMIT,
                                        seed=corpus_seed, kind=KIND)
            cache[corpus_seed] = fit_eval_halves(wins)
        return cache[corpus_seed]

    return get


def run_replicate(get_corpus, fit_seed: int, eval_seed: int, gen_seed: int, split_seed: int,
                  *, order: int = MARKOV_ORDER) -> dict[str, float | None]:
    """Fit the null on ONE draw's fit half, score it against ANOTHER draw's eval half.

    When `fit_seed == eval_seed` this is exactly the runner's path. When they differ, the fit-set matching
    is broken while every other role is held as the arm's plan specifies -- which is the only way to test
    whether that matching is what produces the cancellation.
    """
    from dna_decode.glm.baselines import MarkovGenerator
    from dna_decode.glm.falsifier import heldout_auroc

    fit, _ = get_corpus(fit_seed)
    _, ev = get_corpus(eval_seed)
    gen = MarkovGenerator(order=order, seed=gen_seed).fit(fit).generate(N_GENERATED, LENGTH)
    out: dict[str, float | None] = {}
    for mode in MODES:
        kw = {"k": KMER_K, "seed": split_seed, "mode": mode}
        if mode == "both":
            kw["bins"] = POSITIONAL_BINS
        out[mode] = heldout_auroc(ev, gen, **kw)
    return out


def summarise(values: list[float | None]) -> dict:
    clean = [v for v in values if v is not None]
    return {
        "n": len(clean),
        "n_distinct": len({round(v, 6) for v in clean}),
        "mean": round(statistics.fmean(clean), 4) if clean else None,
        "sd": round(statistics.stdev(clean), 4) if len(clean) > 1 else None,
        "min": round(min(clean), 4) if clean else None,
        "max": round(max(clean), 4) if clean else None,
    }


def bootstrap_sd_ratio(vals_x: list[float], vals_y: list[float], *, n: int = N_BOOTSTRAP,
                       seed: int = 12345) -> dict:
    """CI on `sd(x) / sd(y)` by resampling replicates.

    A point ratio of two SDs at modest n is a coarse reading -- v1 said so and then leaned on the point
    value anyway. Reporting whether the interval excludes 1.0 is the difference between "these differ" and
    "these might differ".
    """
    rng = random.Random(seed)
    ratios = []
    for _ in range(n):
        bx = [rng.choice(vals_x) for _ in vals_x]
        by = [rng.choice(vals_y) for _ in vals_y]
        try:
            sx, sy = statistics.stdev(bx), statistics.stdev(by)
        except statistics.StatisticsError:
            continue
        if sy > 0:
            ratios.append(sx / sy)
    if not ratios:
        return {"point": None, "ci95": None, "excludes_one": None}
    ratios.sort()
    lo = ratios[int(0.025 * len(ratios))]
    hi = ratios[min(len(ratios) - 1, int(0.975 * len(ratios)))]
    point = statistics.stdev(vals_x) / statistics.stdev(vals_y)
    return {"point": round(point, 3), "ci95": [round(lo, 3), round(hi, 3)],
            "excludes_one": bool(lo > 1.0 or hi < 1.0), "n_bootstrap": len(ratios)}


def delta_designs(get_corpus, seeds: list[int], anchor: int = ANCHOR_SEED) -> dict:
    """Does sharing a seed tighten the DELTA? Three designs, same two candidates.

    shared_varying  -- both candidates see the SAME corpus/split, which varies per replicate. The runner's
                       paired design.
    shared_pinned   -- both candidates see the anchor corpus/split; only their generation RNGs vary.
    independent     -- each candidate gets its OWN corpus/split seed, so nothing is shared.

    `Var(A - B) = Var(A) + Var(B) - 2Cov(A, B)`: only the independent design has Cov = 0, so comparing
    the three isolates what the pairing actually buys. Reported with the per-candidate SDs so a reader can
    see whether a tighter delta came from a tighter candidate or from the covariance.
    """
    out = {}
    for mode in MODES:
        res = {}
        for design in ("shared_varying", "shared_pinned", "independent"):
            da, db, dd = [], [], []
            for i, s in enumerate(seeds):
                if design == "shared_varying":
                    ca = cb = s
                    pa = pb = s
                elif design == "shared_pinned":
                    ca = cb = anchor
                    pa = pb = anchor
                else:
                    ca, pa = s, s
                    cb, pb = seeds[(i + 1) % len(seeds)], seeds[(i + 1) % len(seeds)]
                a = run_replicate(get_corpus, ca, ca, s, pa, order=MARKOV_ORDER)[mode]
                b = run_replicate(get_corpus, cb, cb, s, pb, order=DELTA_RIVAL_ORDER)[mode]
                if a is None or b is None:
                    continue
                da.append(a)
                db.append(b)
                dd.append(b - a)
            if len(dd) < 3:
                continue
            sd_a, sd_b = statistics.stdev(da), statistics.stdev(db)
            unpaired = math.sqrt(sd_a ** 2 + sd_b ** 2)
            res[design] = {
                "n": len(dd),
                "mean_delta": round(statistics.fmean(dd), 4),
                "sd_delta": round(statistics.stdev(dd), 4),
                "sd_candidate_order3": round(sd_a, 4),
                "sd_candidate_order1": round(sd_b, 4),
                "sd_unpaired_sum": round(unpaired, 4),
                "delta_over_unpaired": round(statistics.stdev(dd) / unpaired, 3) if unpaired else None,
                "correlation": (round(statistics.correlation(da, db), 3) if len(dd) > 2 else None),
                "sign_consistent": all(x > 0 for x in dd) or all(x < 0 for x in dd),
            }
        out[mode] = res
    return out


def main() -> int:
    if not (REFSEQ / "genome.fna").exists():
        print(f"REFUSE: no local reference at {REFSEQ}")
        return 2
    if not GRID.exists():
        print(f"REFUSE: committed grid absent at {GRID}")
        return 2
    require_sklearn()

    grid = json.loads(GRID.read_text(encoding="utf-8"))
    grid_nulls: dict[int, dict[str, float]] = {}
    for c in grid["grid"]:
        grid_nulls.setdefault(c["seed"], {"kmer": c["null"]["kmer"], "both": c["null"]["both"]})
    grid_seeds = sorted(grid_nulls)

    seeds = list(range(N_REPLICATES))
    print(f"replicates: {len(seeds)} (seeds {seeds[0]}..{seeds[-1]}), anchor {ANCHOR_SEED}")
    print(f"geometry: n={N_GENERATED} limit={LIMIT} order={MARKOV_ORDER} k={KMER_K} "
          f"bins={POSITIONAL_BINS}")

    if not arms_are_distinct(seeds):
        print("REFUSE: two arms share a seed plan, so the decomposition would measure nothing.")
        return 3
    # Every arm's name must describe a role that actually moves in it.
    expect = {
        "A_all_coupled": {"fit_corpus", "eval_corpus", "generation", "split"},
        "A2_streams_decoupled": {"fit_corpus", "eval_corpus", "generation", "split"},
        "B_generation_only": {"generation"},
        "C_corpus_and_split": {"fit_corpus", "eval_corpus", "split"},
        "D_corpus_only": {"fit_corpus", "eval_corpus"},
        "E_split_only": {"split"},
        "F_fit_decoupled": {"fit_corpus", "eval_corpus", "generation", "split"},
    }
    for arm, want in expect.items():
        got = arm_varies_only(arm, seeds)
        if got != want:
            print(f"REFUSE: arm {arm} varies {sorted(got)}, its design says {sorted(want)}")
            return 3
    print("design: 7 arms distinct, every arm varies exactly the roles its name claims")

    t0 = time.time()
    get_corpus = _corpus_cache()

    # ---- THE GATE: arm A must reproduce the committed grid at the grid's own seeds. ----
    gate_rows, gate_hits = [], {m: 0 for m in MODES}
    for s in grid_seeds:
        fs, es, gs, ps = seed_plan("A_all_coupled", s)
        got = run_replicate(get_corpus, fs, es, gs, ps)
        for m in MODES:
            ok = got[m] is not None and abs(got[m] - grid_nulls[s][m]) <= GATE_TOL
            gate_hits[m] += int(ok)
            gate_rows.append({"seed": s, "mode": m, "grid": grid_nulls[s][m],
                              "arm_A": round(got[m], 4) if got[m] is not None else None,
                              "reconciles": ok})
    gate_ok = all(gate_hits[m] == len(grid_seeds) for m in MODES)
    print(f"GATE: arm A reconciles {gate_hits} of {len(grid_seeds)} grid seeds per mode -> "
          f"{'PASS' if gate_ok else 'FAIL'}")
    if not gate_ok:
        print("REFUSE: the all-coupled arm does not reproduce the runner, so no arm built on this "
              "harness describes the committed geometry. This is exactly the v1 defect.")
        OUT.write_text(json.dumps({"record": "glm-seed-role-variance-v2",
                                   "date": date.today().isoformat(),
                                   "status": "GATE_FAILED_NO_NUMBERS_REPORTED",
                                   "gate_rows": gate_rows}, indent=2), encoding="utf-8")
        return 3

    per_arm: dict[str, dict[str, list]] = {a: {m: [] for m in MODES} for a in ARMS}
    for arm in ARMS:
        for s in seeds:
            fs, es, gs, ps = seed_plan(arm, s)
            r = run_replicate(get_corpus, fs, es, gs, ps)
            for m in MODES:
                per_arm[arm][m].append(r[m])
        got = {m: summarise(per_arm[arm][m]) for m in MODES}
        print(f"  {arm:22s} kmer sd {got['kmer']['sd']}  both sd {got['both']['sd']}"
              f"   (distinct {got['both']['n_distinct']}/{got['both']['n']})")

    summary = {a: {m: summarise(per_arm[a][m]) for m in MODES} for a in ARMS}
    for a in ARMS:
        for m in MODES:
            if summary[a][m]["n_distinct"] < 2:
                print(f"REFUSE: arm {a} mode {m} gave {summary[a][m]['n_distinct']} distinct value(s); "
                      "an SD over a repeated number is not a measurement.")
                return 3

    # ---- the comparisons that decide the question ----
    comparisons = {}
    for m in MODES:
        va = [v for v in per_arm["A_all_coupled"][m] if v is not None]
        get = lambda arm: [v for v in per_arm[arm][m] if v is not None]  # noqa: E731
        comparisons[m] = {
            "sd_by_arm": {a: summary[a][m]["sd"] for a in ARMS},
            "vs_A": {a: bootstrap_sd_ratio(get(a), va) for a in ARMS if a != "A_all_coupled"},
            "additivity_corpus_vs_split": {
                "var_C_corpus_and_split": round(statistics.variance(get("C_corpus_and_split")), 6),
                "var_D_corpus_only": round(statistics.variance(get("D_corpus_only")), 6),
                "var_E_split_only": round(statistics.variance(get("E_split_only")), 6),
                "sum_D_plus_E": round(statistics.variance(get("D_corpus_only"))
                                      + statistics.variance(get("E_split_only")), 6),
                "sum_over_C_ratio": round((statistics.variance(get("D_corpus_only"))
                                           + statistics.variance(get("E_split_only")))
                                          / statistics.variance(get("C_corpus_and_split")), 2),
            },
            "additivity_full": {
                "var_A_all_coupled": round(statistics.variance(va), 6),
                "sum_B_plus_D_plus_E": round(statistics.variance(get("B_generation_only"))
                                             + statistics.variance(get("D_corpus_only"))
                                             + statistics.variance(get("E_split_only")), 6),
                "sum_over_A_ratio": round((statistics.variance(get("B_generation_only"))
                                           + statistics.variance(get("D_corpus_only"))
                                           + statistics.variance(get("E_split_only")))
                                          / statistics.variance(va), 2),
                "reading": ("independent roles would make this about 1.0; well above 1 means the roles "
                            "cancel inside one seed, well below means they reinforce"),
            },
        }

    deltas = delta_designs(get_corpus, seeds[:20])

    def mech_verdict(m: str) -> str:
        """Ordered, first match wins. Names which rival the numbers actually support."""
        c = comparisons[m]
        a2 = c["vs_A"]["A2_streams_decoupled"]
        f = c["vs_A"]["F_fit_decoupled"]
        sd_a = c["sd_by_arm"]["A_all_coupled"]
        sd_b = c["sd_by_arm"]["B_generation_only"]
        sd_d = c["sd_by_arm"]["D_corpus_only"]
        if sd_b >= sd_a and sd_d >= sd_a:
            base = "COUPLED_BELOW_BOTH_DECOUPLED_ARMS"
        elif sd_b <= sd_a and sd_d <= sd_a:
            base = "COUPLED_ABOVE_BOTH_DECOUPLED_ARMS"
        else:
            base = "MIXED_NO_CONSISTENT_DIRECTION"
        if a2["excludes_one"]:
            base += "__STREAM_SHARING_MATTERS"
        if f["excludes_one"]:
            base += "__FIT_SET_MATCHING_MATTERS"
        return base

    rec = {
        "record": "glm-seed-role-variance-v2",
        "date": date.today().isoformat(),
        "supersedes": "wiki/glm_seed_role_variance_2026-11-21.json",
        "why_it_supersedes": (
            "v1 ran at n_generated=400; the committed kernel froze n=120. v1's all-coupled arm was "
            "documented as reproducing the status quo and was never reconciled against it -- at n=400 it "
            "reports `both` sd 0.0140 against the grid's own 0.0410. Since the comparison BETWEEN arms "
            "was the whole finding, the finding does not survive the geometry correction."),
        "question": ("At the committed geometry, which seed roles carry the across-seed spread, does the "
                     "coupling suppress or inflate it, and is any suppression a statistical relationship "
                     "or an artifact of three RNGs sharing one integer?"),
        "probe_candidate": f"MarkovGenerator(order={MARKOV_ORDER}) -- free to replicate; NOT the HF "
                           "generator",
        "geometry": {"window_length_bp": LENGTH, "limit_natural": LIMIT,
                     "n_generated_per_candidate": N_GENERATED, "natural_set_arg": KIND,
                     "kmer_k": KMER_K, "positional_bins": POSITIONAL_BINS,
                     "markov_order": MARKOV_ORDER, "feature_modes": list(MODES),
                     "source": "scripts/kaggle/glm_paired_sweep_kernel.py:45 RECORDED"},
        "n_replicates": len(seeds),
        "anchor_seed": ANCHOR_SEED,
        "arms": {
            "A_all_coupled": "one seed drives fit corpus, eval corpus, generation and split (the runner)",
            "A2_streams_decoupled": (f"same logical coupling, but generation seed +{STREAM_OFFSET_GEN} "
                                     f"and split seed +{STREAM_OFFSET_SPLIT} so the three RNG streams "
                                     "are independent"),
            "B_generation_only": "corpus and split pinned at the anchor; only generation varies",
            "C_corpus_and_split": "generation pinned; corpus and split vary together (v1's arm C)",
            "D_corpus_only": "generation and split pinned; only the corpus draw varies",
            "E_split_only": "corpus and generation pinned; only the discriminator split varies",
            "F_fit_decoupled": (f"fit corpus taken from seed+{FIT_DECOUPLE_OFFSET} while the eval corpus "
                                "stays at seed, breaking the fit-set matching v1 named as the mechanism"),
        },
        "seed_plan_table": {a: [dict(zip(("fit_corpus", "eval_corpus", "generation", "split"),
                                         seed_plan(a, s))) for s in seeds[:5]] for a in ARMS},
        "gate": {"rows": gate_rows, "hits": gate_hits, "n_grid_seeds": len(grid_seeds),
                 "tol": GATE_TOL, "passed": gate_ok,
                 "meaning": "arm A is byte-equivalent to the runner at every committed seed and mode"},
        "summary": summary,
        "comparisons": comparisons,
        "mechanism_verdict_by_mode": {m: mech_verdict(m) for m in MODES},
        "paired_delta_designs": deltas,
        "honest_limits": [
            "Measured on the MARKOV candidate. It bounds the INSTRUMENT's nuisance variance and says "
            "nothing directly about the HF generator's own sampling variance. The generator's paired "
            "delta IS measurable from the committed grid and is reported in "
            "wiki/glm_grid_reconcile_2026-11-21.json instead of being estimated here.",
            f"{len(seeds)} replicates per arm. Every arm-to-arm claim is reported as a bootstrap CI on "
            "the SD ratio rather than a point value, and a CI containing 1.0 is reported as such.",
            "The anchor is one seed. Arms B, D and E pin the corpus or split at it, so an anchor-specific "
            "property would show up as an arm effect; v1 swept the anchor for arm B and found the effect "
            "held at 0, 3, 7 and 9, but that sweep was at the wrong geometry and is NOT inherited here.",
            "Arm F changes the fit corpus by a fixed offset, so its fit and eval draws are different but "
            "not independent in any stronger sense; it tests whether the matching matters, not how.",
            "The delta experiment uses Markov order 1 as the rival candidate. It establishes what the "
            "pairing does to a delta between two CHEAP candidates, which is the instrument question; it "
            "is not a statement about a delta between two fine-tuned checkpoints.",
            "One geometry, one organism, one null.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print()
    for m in MODES:
        c = comparisons[m]
        print(f"[{m}] sd by arm: " + "  ".join(f"{a.split('_')[0]}={c['sd_by_arm'][a]}" for a in ARMS))
        for a in ("A2_streams_decoupled", "B_generation_only", "D_corpus_only", "E_split_only",
                  "F_fit_decoupled"):
            r = c["vs_A"][a]
            print(f"      {a:22s} sd ratio vs A {r['point']} CI {r['ci95']} "
                  f"{'(excludes 1)' if r['excludes_one'] else '(includes 1)'}")
        print(f"      additivity B+D+E over A: {c['additivity_full']['sum_over_A_ratio']}x")
        print(f"      -> {mech_verdict(m)}")
        d = deltas.get(m, {})
        for design, v in d.items():
            print(f"      delta[{design:15s}] sd {v['sd_delta']}  unpaired {v['sd_unpaired_sum']}  "
                  f"ratio {v['delta_over_unpaired']}  corr {v['correlation']}")
    print(f"\n{rec['wall_seconds']}s -> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
