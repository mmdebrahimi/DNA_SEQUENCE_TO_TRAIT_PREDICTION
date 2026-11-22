"""Reconcile every seed-role claim against the grid that produced it, and build a LIKE-FOR-LIKE floor.

    uv run python scripts/glm_grid_reconcile.py        # read-only; writes one artifact

**Why this exists: two probes were scored off-geometry and nothing caught it.** The committed sweep
kernel froze `{"n": 120, "limit_natural": 1200, "bins": 6, "length": 150}`
(`scripts/kaggle/glm_paired_sweep_kernel.py:45`). Both seed-role probes hardcoded `N_GENERATED = 400`.
Everything else matched, so both produced entirely plausible numbers that described a geometry no committed
artifact uses. The all-coupled arm's docstring claimed it "reproduces the status quo" and it was never
checked against the status quo; at n=400 it reports `both` sd 0.0140 where the grid's own ten nulls have sd
0.0410. **A harness that re-implements a runner's path by hand must be gated against that runner's output
before any arm built on it is read.** This script is that gate.

**The second, deeper defect it closes: the committed margin subtracted two different measurements.**
`margin = null - ceiling` paired

  * a NULL measured as 600 natural EVAL windows against 120 GENERATED sequences (5:1, test half ~360), with
  * a CEILING measured as `self_control_ceiling`, which splits the 600 eval windows 300-vs-300 (1:1).

Different class balance and different test-set size, so the difference of the two AUROCs is not a margin in
any geometry. A held-out AUROC's sampling variation depends on both, which is exactly why they cannot be
subtracted. The fix is a floor measured in the NULL's own geometry: score the same 600 eval windows against
**120 REAL natural windows held out of both halves**. That is the score a PERFECT generator would earn --
real sequence against real sequence at the null's own balance -- so

    like_for_like_margin = null - perfect_generator_floor

is a subtraction of two commensurable numbers. The 120 extra windows come free: `build_natural_set`
shuffles all 4340 promoters with `Random(seed)` and THEN truncates, so building at `limit=1320` leaves the
first 1200 byte-identical to the committed `limit=1200` draw and hands back 120 disjoint spares. That
prefix-stability is asserted at runtime, not assumed -- if it ever broke, the "held-out" windows could
overlap the eval set and the floor would be measured against itself.

Also reported, all free from artifacts already on disk:
  * the MEAN-based ceiling, which the committed artifact listed as an unquantified limit while shipping
    every one of its five per-seed values;
  * per-seed `model - null` in kmer mode, which is what Decision One actually turns on: a mode can carry a
    registerable distance from the generator to the null even where the null itself is near its floor.

No network, no GPU, no torch, no D: volume. sklearn is asserted on entry because `heldout_auroc` catches
ImportError and returns None, which a caller reads as "insufficient data" -- a missing dependency must
never be able to pose as a measurement.
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
GRID = ROOT / "wiki" / "glm_generator_paired_sweep_2026-11-20.json"
CEILINGS = ROOT / "wiki" / "glm_seed_role_ceilings_2026-11-21.json"
OUT = ROOT / "wiki" / "glm_grid_reconcile_2026-11-21.json"

# THE KERNEL'S FROZEN GEOMETRY, transcribed from scripts/kaggle/glm_paired_sweep_kernel.py:45.
LENGTH = 150
LIMIT = 1200
N_GENERATED = 120          # the value both earlier probes got wrong (they used 400)
KIND = "promoter"
KMER_K = 4
POSITIONAL_BINS = 6
MARKOV_ORDER = 3
MODES = ("kmer", "both")

RECONCILE_TOL = 5e-4       # the grid rounds to 4 dp


def require_sklearn() -> None:
    """Fail loudly on a missing discriminator dependency.

    `falsifier.heldout_auroc` wraps its sklearn import in `except ImportError: return None`, and None is
    the same value it returns when a class is too small to split. So on an interpreter without sklearn
    every measurement in this file would come back None and the run would read as a data problem. Measured
    once for real: under bare `python` on this host every AUROC was None; under `uv run` all of them
    resolve.
    """
    import sklearn  # noqa: F401
    import numpy  # noqa: F401


def corpus_with_spares(seed: int, *, n_spare: int):
    """The committed draw, plus `n_spare` REAL windows SEQUENCE-disjoint from both halves of it.

    Returns `(fit, eval, spare, n_skipped)`. The committed draw is `wins[:LIMIT]` halved -- re-truncated
    BY HAND rather than letting `fit_eval_halves` halve the longer list, because halving 1320 would give
    660/660 and silently move the eval set off the committed geometry, the same class of error this script
    exists to catch.

    **Spares are filtered by SEQUENCE, not by index, and that is not pedantry: the first version of this
    function refused at seed 1 because a spare window sequence genuinely appears in the eval half.** The
    corpus is 4340 windows over 4323 distinct sequences -- 7 sequences duplicated across 24 windows, one
    of them 6 times -- because distinct genes can share an identical upstream region (operon structure,
    duplicated genes, repeated elements). An index-disjoint spare can therefore be sequence-identical to
    an eval window, and the discriminator sees only sequence. Leaving one in would put the same feature
    vector in both classes and depress the floor, so the floor would read as cleaner than it is.
    """
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves

    wins, _stats = build_natural_set(REFSEQ, length=LENGTH, limit=None, seed=seed, kind=KIND)
    if len(wins) < LIMIT + n_spare:
        raise SystemExit(f"REFUSE: only {len(wins)} windows available; need {LIMIT + n_spare}")
    committed = wins[:LIMIT]
    fit, ev = fit_eval_halves(committed)
    used = set(fit) | set(ev)
    spare, skipped = [], 0
    for w in wins[LIMIT:]:
        if len(spare) >= n_spare:
            break
        if w.seq in used:
            skipped += 1
            continue
        spare.append(w.seq)
        used.add(w.seq)          # also keeps the spare set internally duplicate-free
    if len(spare) < n_spare:
        raise SystemExit(f"REFUSE: only {len(spare)} sequence-disjoint spares at seed {seed}")
    return fit, ev, spare, skipped


def corpus_redundancy() -> dict:
    """Measure how far the fit/eval split actually delivers the disjointness the runner promises.

    `scripts/glm_generator_falsifier.py` states that the Markov null is fitted on FIT only and every
    discrimination is against EVAL only, so that no generator is ever compared against sequences it was
    trained on. Duplicate upstream windows mean that is not exactly true. Measured rather than assumed,
    and reported in proportion: the leak is a handful of windows in 600, which is a real caveat on the
    promise and far too small to move any score.
    """
    from collections import Counter

    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves

    wins, _ = build_natural_set(REFSEQ, length=LENGTH, limit=None, seed=0, kind=KIND)
    seqs = [w.seq for w in wins]
    c = Counter(seqs)
    dups = {k: v for k, v in c.items() if v > 1}
    per_seed = {}
    for s in range(10):
        ww, _ = build_natural_set(REFSEQ, length=LENGTH, limit=LIMIT, seed=s, kind=KIND)
        fit, ev = fit_eval_halves(ww)
        fs = set(fit)
        per_seed[str(s)] = sum(1 for e in ev if e in fs)
    return {
        "n_windows": len(seqs),
        "n_distinct_sequences": len(c),
        "n_duplicated_sequences": len(dups),
        "n_windows_involved_in_a_duplicate": sum(dups.values()),
        "max_copies_of_one_sequence": max(c.values()),
        "eval_windows_also_in_fit_by_seed": per_seed,
        "worst_case_leak_fraction_of_eval": round(max(per_seed.values()) / (LIMIT // 2), 5),
        "reading": ("The runner claim that no generator is scored against sequences it was trained on "
                    "holds up to duplicate upstream windows: at most 3 of 600 eval windows (0.5%) also "
                    "appear in the fit half. Recorded as a real caveat on the promise, and far too small "
                    "to move a score."),
    }


def assert_prefix_stable(seed: int) -> None:
    """The spares are only 'held out' if over-building does not reshuffle the committed prefix.

    `build_natural_set` shuffles the full 4340-window list with `Random(seed)` and truncates afterwards, so
    the first LIMIT entries must be identical at any larger limit. If that ever stopped holding, the
    'disjoint' spares could land inside the eval set and the perfect-generator floor would be a set scored
    against itself -- which would read as a suspiciously clean floor, not as an error.
    """
    from dna_decode.glm.corpus import build_natural_set

    a, _ = build_natural_set(REFSEQ, length=LENGTH, limit=LIMIT, seed=seed, kind=KIND)
    b, _ = build_natural_set(REFSEQ, length=LENGTH, limit=LIMIT + N_GENERATED, seed=seed, kind=KIND)
    if [w.seq for w in a] != [w.seq for w in b[:LIMIT]]:
        raise SystemExit(
            f"REFUSE: corpus prefix is NOT stable under a larger limit at seed {seed}. The spare windows "
            "cannot be assumed disjoint from the eval half, so no like-for-like floor can be built.")


def nulls_from_grid(grid: dict) -> dict[int, dict[str, float]]:
    """Per-seed null pair, read FROM the grid. Refuses if a seed carries two distinct pairs."""
    per: dict[int, set] = {}
    for c in grid["grid"]:
        per.setdefault(c["seed"], set()).add((c["null"]["kmer"], c["null"]["both"]))
    out: dict[int, dict[str, float]] = {}
    for seed, pairs in sorted(per.items()):
        if len(pairs) != 1:
            raise SystemExit(f"REFUSE: seed {seed} carries {len(pairs)} distinct null pairs")
        kmer, both = next(iter(pairs))
        out[seed] = {"kmer": kmer, "both": both}
    return out


def model_minus_null_from_grid(grid: dict) -> dict:
    """Decision One's actual premise, free from the grid: how far is the generator FROM the null?

    A mode whose null sits at its own noise floor can still be perfectly usable, because the quantity an
    arm set is scored on is the generator's DISTANCE from the null, not the null's distance from the floor.
    If the generator sits far above the null in kmer mode then movement toward the null is registerable
    there, and only the last stretch -- actually beating the null -- is unresolvable. Retiring the mode on
    the margin alone would discard a working measurement.
    """
    out = {}
    for mode in MODES:
        per_model: dict[str, list[float]] = {}
        for c in grid["grid"]:
            gap = c[mode] - c["null"][mode]
            per_model.setdefault(c["model"], []).append(gap)
        allv = [g for v in per_model.values() for g in v]
        out[mode] = {
            "by_model": {m: {"mean": round(statistics.fmean(v), 4),
                             "sd": round(statistics.stdev(v), 4) if len(v) > 1 else None,
                             "min": round(min(v), 4), "max": round(max(v), 4),
                             "n_positive": sum(1 for g in v if g > 0), "n": len(v)}
                         for m, v in sorted(per_model.items())},
            "all_models_mean": round(statistics.fmean(allv), 4),
            "all_models_min": round(min(allv), 4),
            "n_cells_generator_above_null": sum(1 for g in allv if g > 0),
            "n_cells": len(allv),
        }
    return out


def paired_delta_from_grid(grid: dict) -> dict:
    """The one GPU-free view of the HF generator's own variance: is its delta tighter than the unpaired sum?

    `Var(A - B) = Var(A) + Var(B) - 2Cov(A, B)`. Measuring the per-seed delta's SD against
    `sqrt(Var(A) + Var(B))` says directly whether the shared seed buys precision ON THE DELTA, which is the
    quantity the arm is actually scored on. This is the real generator, not the Markov probe, so it does not
    inherit the probe's scope limit.
    """
    out = {}
    for mode in MODES:
        by_model: dict[str, dict[int, float]] = {}
        nulls: dict[int, float] = {}
        for c in grid["grid"]:
            by_model.setdefault(c["model"], {})[c["seed"]] = c[mode]
            nulls[c["seed"]] = c["null"][mode]
        res = {}
        pairs = [("zeroshot", "null"), ("epoch3", "zeroshot"), ("epoch3", "null")]
        for a, b in pairs:
            va = by_model.get(a)
            vb = nulls if b == "null" else by_model.get(b)
            if not va or not vb:
                continue
            seeds = sorted(set(va) & set(vb))
            if len(seeds) < 3:
                continue
            xa = [va[s] for s in seeds]
            xb = [vb[s] for s in seeds]
            d = [p - q for p, q in zip(xa, xb)]
            sd_a, sd_b = statistics.stdev(xa), statistics.stdev(xb)
            unpaired = math.sqrt(sd_a ** 2 + sd_b ** 2)
            res[f"{a}_minus_{b}"] = {
                "n_seeds": len(seeds),
                "mean_delta": round(statistics.fmean(d), 4),
                "sd_paired_delta": round(statistics.stdev(d), 4),
                "sd_unpaired_sum": round(unpaired, 4),
                "paired_over_unpaired": round(statistics.stdev(d) / unpaired, 3) if unpaired else None,
                "correlation": round(statistics.correlation(xa, xb), 3) if len(seeds) > 2 else None,
                "pairing_helps": statistics.stdev(d) < unpaired,
            }
        out[mode] = res
    return out


def mean_based_ceiling(ceil_rec: dict) -> dict:
    """The limit the committed ceilings artifact NAMED and never computed.

    It states that the ceiling is the MAX over five self-control seeds, "so a mean-based floor would be
    lower and every margin correspondingly more positive", and ships all five values per seed -- so the
    quantity was one line of arithmetic away the whole time. Naming a limit is not the same as bounding it.
    """
    out = {}
    for mode in MODES:
        rows = [r for r in ceil_rec["rows"] if r["mode"] == mode]
        per = []
        for r in rows:
            vals = r["ceiling_values"]
            per.append({"seed": r["seed"], "max": r["ceiling"],
                        "mean": round(statistics.fmean(vals), 4),
                        "min": round(min(vals), 4),
                        "null": r["null"],
                        "margin_vs_max": r["margin"],
                        "margin_vs_mean": round(r["null"] - statistics.fmean(vals), 4)})
        mm = [p["margin_vs_mean"] for p in per]
        mx = [p["margin_vs_max"] for p in per]
        out[mode] = {
            "per_seed": per,
            "ceiling_max_mean": round(statistics.fmean([p["max"] for p in per]), 4),
            "ceiling_mean_mean": round(statistics.fmean([p["mean"] for p in per]), 4),
            "margin_vs_max_mean": round(statistics.fmean(mx), 4),
            "margin_vs_mean_mean": round(statistics.fmean(mm), 4),
            "positive_vs_max": sum(1 for v in mx if v > 0),
            "positive_vs_mean": sum(1 for v in mm if v > 0),
            "n_seeds": len(per),
            "verdict_flips_on_floor_rule": (statistics.fmean(mx) <= 0) != (statistics.fmean(mm) <= 0),
        }
    return out


def margin_significance(margins: list[float]) -> dict:
    """Is the mean like-for-like margin resolvably above zero at this replicate count?

    The margin is already a PAIRED difference -- null and floor are measured at the same seed, on the same
    eval half, through the same split -- so a one-sample test on the ten margins is the right shape. This
    matters because a count like "positive at 7 of 10 seeds" reads as a finding while being entirely
    consistent with a true margin of zero: both terms carry real sampling variation at 120 generated
    sequences, and in kmer mode the margin's own SD comes out roughly twice its mean.

    Reported as a t statistic with a normal-approximation two-sided p and a 95% interval. No scipy
    dependency is taken for this; the approximation is stated rather than hidden, and with n=10 the
    t-vs-normal difference does not change any reading here.
    """
    n = len(margins)
    if n < 3:
        return {"n": n, "resolvable": None, "note": "too few replicates to test"}
    mean = statistics.fmean(margins)
    sd = statistics.stdev(margins)
    se = sd / math.sqrt(n)
    tstat = mean / se if se else None
    # Two-sided normal approximation via the error function; stated, not hidden.
    p = (2 * (1 - 0.5 * (1 + math.erf(abs(tstat) / math.sqrt(2))))) if tstat is not None else None
    crit = 2.262 if n == 10 else 1.96          # t(0.975, df=9) when n is exactly 10
    return {
        "n": n,
        "mean": round(mean, 4),
        "sd": round(sd, 4),
        "stderr": round(se, 4),
        "t": round(tstat, 3) if tstat is not None else None,
        "p_two_sided_normal_approx": round(p, 4) if p is not None else None,
        "ci95": [round(mean - crit * se, 4), round(mean + crit * se, 4)],
        "resolvable_above_zero": bool(tstat is not None and (mean - crit * se) > 0),
        "sd_over_mean": round(sd / mean, 2) if mean else None,
    }


def power_control(seeds: list[int]) -> dict:
    """Can this discriminator detect a BADLY wrong generator at all, in each mode?

    A winnable region is only meaningful if the instrument has power, and a floor that sits at chance
    could in principle mean the discriminator is simply blind. The uniform i.i.d. generator is the
    standard control: it is obviously not DNA-like, so a working discriminator must separate it easily.
    This also re-derives the committed CPU reference (kmer 0.9408 / both 0.9641 at seed 7 from
    wiki/glm_generator_falsifier_REAL_2026-10-07.json via the kernel's RECORDED block), so the power claim
    is reproduced here rather than cited.

    It is deliberately NOT a floor: uniform sequence is far more different from natural than any real
    generator, so this bounds the instrument's CEILING of detectability, not the target a generator must
    beat.
    """
    from dna_decode.glm.baselines import UniformGenerator
    from dna_decode.glm.falsifier import heldout_auroc

    out = {}
    for mode in MODES:
        vals = []
        for s in seeds:
            _fit, ev, _spare, _sk = corpus_with_spares(s, n_spare=1)
            uni = UniformGenerator(seed=s).generate(N_GENERATED, LENGTH)
            kw = {"k": KMER_K, "seed": s, "mode": mode}
            if mode == "both":
                kw["bins"] = POSITIONAL_BINS
            v = heldout_auroc(ev, uni, **kw)
            if v is not None:
                vals.append(v)
        out[mode] = {
            "n": len(vals),
            "mean": round(statistics.fmean(vals), 4),
            "min": round(min(vals), 4),
            "max": round(max(vals), 4),
            "seed7": round(vals[seeds.index(7)], 4) if 7 in seeds else None,
        }
    return out


def main() -> int:
    if not (REFSEQ / "genome.fna").exists():
        print(f"REFUSE: no local reference at {REFSEQ}")
        return 2
    for p in (GRID, CEILINGS):
        if not p.exists():
            print(f"REFUSE: required artifact absent at {p}")
            return 2
    require_sklearn()

    from dna_decode.glm.baselines import MarkovGenerator
    from dna_decode.glm.falsifier import heldout_auroc

    grid = json.loads(GRID.read_text(encoding="utf-8"))
    ceil_rec = json.loads(CEILINGS.read_text(encoding="utf-8"))
    nulls = nulls_from_grid(grid)
    seeds = sorted(nulls)
    print(f"grid: {len(grid['grid'])} cells, seeds {seeds}")
    print(f"geometry under test: n={N_GENERATED} limit={LIMIT} order={MARKOV_ORDER} "
          f"k={KMER_K} bins={POSITIONAL_BINS}")

    t0 = time.time()
    assert_prefix_stable(seeds[0])
    print("prefix stability: committed draw is byte-identical under the larger limit")

    rows = []
    recon_hits = {m: 0 for m in MODES}
    for s in seeds:
        fit, ev, spare, n_skipped = corpus_with_spares(s, n_spare=N_GENERATED)
        gen = MarkovGenerator(order=MARKOV_ORDER, seed=s).fit(fit).generate(N_GENERATED, LENGTH)
        # Disjointness is CHECKED after filtering, not argued: a leftover overlap would put the same
        # feature vector in both classes and make the floor a partial self-comparison.
        overlap = len(set(spare) & (set(ev) | set(fit)))
        if overlap:
            print(f"REFUSE: {overlap} spare windows still overlap the committed draw at seed {s}")
            return 3
        for mode in MODES:
            kw = {"k": KMER_K, "seed": s, "mode": mode}
            if mode == "both":
                kw["bins"] = POSITIONAL_BINS
            repro = heldout_auroc(ev, gen, **kw)
            floor = heldout_auroc(ev, spare, **kw)
            committed = nulls[s][mode]
            ok = repro is not None and abs(repro - committed) <= RECONCILE_TOL
            recon_hits[mode] += int(ok)
            rows.append({
                "seed": s, "mode": mode,
                "null_committed": committed,
                "null_reproduced": round(repro, 4) if repro is not None else None,
                "reconciles": ok,
                "perfect_generator_floor": round(floor, 4) if floor is not None else None,
                "like_for_like_margin": (round(committed - floor, 4)
                                         if floor is not None else None),
                "n_eval": len(ev), "n_generated": len(gen), "n_spare": len(spare),
                "n_spare_skipped_as_duplicates": n_skipped,
            })
        print(f"  seed {s}: null {nulls[s]['kmer']:.4f}/{nulls[s]['both']:.4f}  "
              f"floor {rows[-2]['perfect_generator_floor']:.4f}/{rows[-1]['perfect_generator_floor']:.4f}"
              f"  reconciles {rows[-2]['reconciles']}/{rows[-1]['reconciles']}")

    # ---- THE GATE. Nothing downstream is readable if the harness does not reproduce the grid. ----
    n = len(seeds)
    gate_ok = all(recon_hits[m] == n for m in MODES)
    print(f"\nGATE: reconciled {recon_hits} of {n} seeds per mode -> "
          f"{'PASS' if gate_ok else 'FAIL'}")

    summary = {}
    for mode in MODES:
        sel = [r for r in rows if r["mode"] == mode and r["like_for_like_margin"] is not None]
        mg = [r["like_for_like_margin"] for r in sel]
        fl = [r["perfect_generator_floor"] for r in sel]
        nu = [r["null_committed"] for r in sel]
        summary[mode] = {
            "n_seeds": len(sel),
            "null_mean": round(statistics.fmean(nu), 4),
            "null_sd": round(statistics.stdev(nu), 4),
            "perfect_floor_mean": round(statistics.fmean(fl), 4),
            "perfect_floor_sd": round(statistics.stdev(fl), 4),
            "perfect_floor_min": min(fl), "perfect_floor_max": max(fl),
            "like_for_like_margin_mean": round(statistics.fmean(mg), 4),
            "like_for_like_margin_sd": round(statistics.stdev(mg), 4),
            "seeds_positive": sum(1 for v in mg if v > 0),
            "seeds_negative": sum(1 for v in mg if v < 0),
        }

    sig = {m: margin_significance([r["like_for_like_margin"] for r in rows
                                   if r["mode"] == m and r["like_for_like_margin"] is not None])
           for m in MODES}
    pw = power_control(seeds)
    mb = mean_based_ceiling(ceil_rec)
    mmn = model_minus_null_from_grid(grid)
    pd = paired_delta_from_grid(grid)

    def mode_verdict(mode: str) -> str:
        """Ordered, first match wins.

        `UNRESOLVED_AT_THIS_N` sits ABOVE the seed-dependent label deliberately. A positive mean whose
        confidence interval straddles zero is not a small winnable region, it is an unmeasured one, and
        reporting "positive at 7 of 10 seeds" without that distinction is how a count becomes a finding.
        """
        s = summary[mode]
        if s["like_for_like_margin_mean"] <= 0:
            return "NO_WINNABLE_REGION_ON_AVERAGE"
        if not sig[mode]["resolvable_above_zero"]:
            return "WINNABLE_REGION_UNRESOLVED_AT_THIS_N"
        if s["seeds_negative"] > 0:
            return "WINNABLE_REGION_SEED_DEPENDENT"
        return "WINNABLE_REGION_AT_EVERY_SEED"

    rec = {
        "record": "glm-grid-reconcile-v1",
        "date": date.today().isoformat(),
        "question": ("Do the seed-role probes reproduce the grid they were scored against, and what is "
                     "the winnable region when the floor is measured in the NULL's own geometry?"),
        "geometry": {"window_length_bp": LENGTH, "limit_natural": LIMIT,
                     "n_generated_per_candidate": N_GENERATED, "natural_set_arg": KIND,
                     "kmer_k": KMER_K, "positional_bins": POSITIONAL_BINS,
                     "markov_order": MARKOV_ORDER, "feature_modes": list(MODES),
                     "source": "scripts/kaggle/glm_paired_sweep_kernel.py:45 RECORDED"},
        "the_defect_this_found": {
            "what": "Both seed-role probes used n_generated=400; the committed kernel froze n=120.",
            "consequence": ("The all-coupled arm claimed to reproduce the status quo and did not. At "
                            "n=400 it reports `both` sd 0.0140; the grid's own ten nulls have sd 0.0410, "
                            "so the arm it was compared against was never the status quo."),
            "detected_by": "reconciling the arm's per-seed values against the grid, which was never done",
        },
        "the_second_defect": {
            "what": ("`margin = null - ceiling` subtracted a 600-vs-120 AUROC from a 300-vs-300 one."),
            "why_it_matters": ("A held-out AUROC's sampling variation depends on test-set size and class "
                               "balance, so the two numbers are not commensurable and their difference "
                               "is not a margin in any geometry."),
            "fix": ("score the same eval half against 120 REAL held-out natural windows -- the score a "
                    "PERFECT generator earns at the null's own balance."),
        },
        "gate": {"reconciled_per_mode": recon_hits, "n_seeds": n, "tol": RECONCILE_TOL,
                 "passed": gate_ok,
                 "meaning": ("the hand-written null path is byte-equivalent to the runner's at every "
                             "seed, so arms built on it describe the committed geometry")},
        "rows": rows,
        "summary": summary,
        "verdict_by_mode": {m: mode_verdict(m) for m in MODES},
        "margin_significance": sig,
        "power_control_uniform": pw,
        "corpus_redundancy": corpus_redundancy(),
        "floor_rule_sensitivity": mb,
        "generator_distance_to_null": mmn,
        "paired_delta_on_the_real_generator": pd,
        "honest_limits": [
            "The perfect-generator floor uses 120 REAL windows, so it is the floor for a generator that "
            "emits genuine promoter sequence. It is an EMPIRICAL floor at one draw per seed, not an "
            "analytic one, and carries its own sampling variation.",
            "Spares are filtered by SEQUENCE against both halves of the committed draw. The corpus "
            "holds 7 sequences duplicated across 24 of 4340 windows, so an index-disjoint spare can be "
            "sequence-identical to an eval window; the count filtered out is reported per seed.",
            "The spare windows come from the same organism and the same extraction, so they share any "
            "systematic property of the corpus. That is the point -- a floor must not be more different "
            "from the eval set than a perfect generator would be -- but it does mean this floor cannot "
            "detect a corpus-level artifact.",
            "The nulls are READ from the grid and separately REPRODUCED here; the reproduction is what "
            "licenses the comparison, and it is reported per seed rather than asserted.",
            "`model - null` and the paired delta are read from the committed grid, whose own verdict is "
            "INSTRUMENT_INVALID for an unrelated reason (a per-cell null gate that could only ever pass "
            "at seed 7; corrected in scripts/glm_paired_sweep_verdict.py). The per-cell scores are not "
            "in dispute; the verdict label is.",
            "A winnable region is necessary, NOT sufficient: a k-mer-plus-positional discriminator sees "
            "LOCAL composition, so clearing it would not establish functional realism.",
            "The margin significance test uses a normal approximation to the t p-value and the "
            "t(0.975, df=9) critical value for the interval. Stated rather than hidden; at n=10 it "
            "changes no reading here, and the interval is what the verdict depends on.",
            "The uniform power control bounds the instrument's CEILING of detectability, not the target "
            "a generator must beat: i.i.d. sequence is far more different from natural than any real "
            "generator, so high power against it does not imply power against a near-miss.",
            "One organism, one geometry, one null.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print()
    for mode in MODES:
        s = summary[mode]
        print(f"[{mode}] null {s['null_mean']:.4f} vs PERFECT-GENERATOR floor "
              f"{s['perfect_floor_mean']:.4f} -> like-for-like margin "
              f"{s['like_for_like_margin_mean']:+.4f} (sd {s['like_for_like_margin_sd']:.4f}), "
              f"positive {s['seeds_positive']}/{s['n_seeds']} -> {mode_verdict(mode)}")
        print(f"         committed margin vs MAX self-control ceiling was "
              f"{mb[mode]['margin_vs_max_mean']:+.4f}; vs MEAN ceiling "
              f"{mb[mode]['margin_vs_mean_mean']:+.4f} "
              f"(flips: {mb[mode]['verdict_flips_on_floor_rule']})")
        q = sig[mode]
        print(f"         margin t={q['t']} p~{q['p_two_sided_normal_approx']} "
              f"CI95 {q['ci95']} resolvable={q['resolvable_above_zero']} "
              f"(sd/mean {q['sd_over_mean']})")
        print(f"         POWER control vs uniform: mean {pw[mode]['mean']:.4f} "
              f"(range {pw[mode]['min']:.4f}-{pw[mode]['max']:.4f})")
        g = mmn[mode]
        print(f"         generator ABOVE null in {g['n_cells_generator_above_null']}/{g['n_cells']} "
              f"cells, mean distance {g['all_models_mean']:+.4f} (min {g['all_models_min']:+.4f})")
    print(f"\n{rec['wall_seconds']}s -> {OUT.relative_to(ROOT).as_posix()}")
    return 0 if gate_ok else 3


if __name__ == "__main__":
    raise SystemExit(main())
