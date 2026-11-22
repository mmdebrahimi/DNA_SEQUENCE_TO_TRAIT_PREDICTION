"""Sweep the two instrument levers nobody has swept: the FEATURE MODE and the NULL's ORDER.

    uv run python scripts/glm_instrument_levers.py        # read-only; writes one artifact

**Two questions, both decision-relevant and neither previously measured.**

**1. Does `positional` alone match `both`?** The committed grid scores only `kmer` and `both`, and `both`
is literally the concatenation of the two feature blocks (256 k-mer features plus 96 positional ones). So
whenever `both` beats `kmer` it is natural to credit the positional block -- but nobody has run the
positional block on its own. If positional alone matches `both`, the honest recommendation is the smaller
feature set rather than the union, and the k-mer block is carrying nothing. If `both` beats positional
alone, the two blocks are genuinely complementary and the union earns its size. This matters directly to
the open question about which mode to score on, and it cannot be read off the existing table.

**2. Is "the Markov null" even well specified?** Every generator claim in this arm is "the generator is
worse than a 3-mer Markov chain". The ORDER is a free parameter fixed at 3 by a CLI default, and nothing
has measured how much the null's strength depends on it. An order-1 chain knows only base composition; an
order-5 chain memorises 6-mers and, fitted on 600 windows of 150 bp, approaches memorising the training
set. If the null's distance above the floor moves a lot with order, then "beats the Markov null" is a
claim about a tuning choice as much as about the generator, and that needs saying out loud.

**The floor is computed once per (mode, seed) and reused across orders, which is correct and is the point.**
The perfect-generator floor scores real sequence against real sequence; it does not involve the null at
all, so it cannot depend on the null's order. Recomputing it per order would add noise and invite a reader
to think order affects the floor. Any apparent order effect is therefore entirely in the null term.

Geometry is otherwise the committed kernel's. Free: pure Python plus sklearn, no torch, no GPU, no
network, no D: volume.
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
OUT = ROOT / "wiki" / "glm_instrument_levers_2026-11-21.json"

# A SECOND substrate, used only to test whether the order-dependence pattern is STRUCTURAL.
# The mechanism claim is that a Markov chain's sufficient statistic IS the k-mer conditional
# distribution, so raising its order must saturate a k-mer discriminator in ANY substrate while leaving
# a positional one alone. That is a strong prediction and it is cheap to falsify: Campylobacter jejuni
# is GC 0.266 against E. coli's 0.453, a compact genome with different promoter composition. If the
# pattern is an E. coli quirk it should break here; if it is structural it must reproduce.
SECOND = ROOT / "data" / "refseq_cache" / "GCA_060772745.1"
SECOND_LABEL = "Campylobacter jejuni GCA_060772745.1 (GC 0.266)"
REF_LABEL = "Escherichia coli K-12 MG1655 (GC 0.453)"

STRUCTURAL_PREDICTION = {
    "registered_before_the_second_substrate_was_run": True,
    "claim": ("the kmer margin collapses as the Markov null's order rises because the chain's "
              "sufficient statistic IS the k-mer conditional distribution, while the positional margin "
              "is order-invariant because no order of Markov chain models WHERE composition sits"),
    "therefore": ("the pattern must be STRUCTURAL and must reproduce in any substrate: kmer range "
                  "across orders much larger than positional range, and no positional order resolvably "
                  "different from the committed one"),
    "falsified_if": ("the second substrate shows a flat kmer profile, or a positional profile that "
                     "varies with order -- either would make the E. coli result a substrate quirk and "
                     "the mechanism claim unsupported"),
}

LENGTH = 150
LIMIT = 1200
N_GENERATED = 120
KIND = "promoter"
KMER_K = 4
POSITIONAL_BINS = 6
SEEDS = tuple(range(10))

MODES = ("kmer", "positional", "both")
ORDERS = (1, 2, 3, 4, 5)
COMMITTED_ORDER = 3


def require_sklearn() -> None:
    import numpy  # noqa: F401
    import sklearn  # noqa: F401


def feature_count(mode: str) -> int:
    """Feature counts as `build_features` actually produces them.

    `positional_features` is called with `k=min(KMER_K, 2)`, so the positional block is
    `bins * 4**2 = 96` and NOT `bins * 4**4`. Getting this wrong would misreport every p-against-n
    figure, and the clamp is easy to miss because it lives inside the builder rather than at the
    call site.
    """
    kmer = 4 ** KMER_K
    pos = POSITIONAL_BINS * (4 ** min(KMER_K, 2))
    return {"kmer": kmer, "positional": pos, "both": kmer + pos}[mode]


def resolve(margins: list[float]) -> dict:
    n = len(margins)
    if n < 3:
        return {"n": n, "resolvable_above_zero": None}
    mean = statistics.fmean(margins)
    sd = statistics.stdev(margins)
    se = sd / math.sqrt(n)
    crit = 2.262 if n == 10 else 1.96
    lo, hi = mean - crit * se, mean + crit * se
    return {"n": n, "mean": round(mean, 4), "sd": round(sd, 4), "stderr": round(se, 4),
            "ci95": [round(lo, 4), round(hi, 4)], "t": round(mean / se, 3) if se else None,
            "resolvable_above_zero": bool(lo > 0),
            "n_positive": sum(1 for v in margins if v > 0)}


def paired_diff(a: list[float], b: list[float]) -> dict:
    """Seed-matched difference between two measurements, which is the only fair way to compare them.

    Both arms see the same corpus draw and the same split at each seed, so differencing within a seed
    removes the draw-to-draw variation that dominates either arm on its own. Comparing two unpaired
    means at ten seeds would hide a real difference behind corpus noise.
    """
    d = [x - y for x, y in zip(a, b)]
    n = len(d)
    mean = statistics.fmean(d)
    sd = statistics.stdev(d) if n > 1 else 0.0
    se = sd / math.sqrt(n) if n else 0.0
    crit = 2.262 if n == 10 else 1.96
    # A zero-variance difference is MAXIMALLY resolvable, not unresolvable. An earlier version guarded
    # with `bool(se and ...)`, which made sd == 0 report False -- the exact inversion of the truth. It
    # survived its own unit test only because that test built the second arm by subtracting a constant,
    # and floating-point error left se tiny-but-nonzero rather than exactly zero. Real data here never
    # hits sd == 0, so no reported number was affected, but the branch is wrong and a future caller
    # would be silently told the clearest possible result is undetectable.
    if n < 2:
        resolvable = False
    elif se == 0:
        resolvable = mean != 0
    else:
        resolvable = (mean - crit * se > 0) or (mean + crit * se < 0)
    return {"n": n, "mean_diff": round(mean, 4), "sd_diff": round(sd, 4),
            "ci95": [round(mean - crit * se, 4), round(mean + crit * se, 4)],
            "resolvably_nonzero": bool(resolvable),
            "n_positive": sum(1 for v in d if v > 0)}


def main() -> int:
    if not (REFSEQ / "genome.fna").exists():
        print(f"REFUSE: no local reference at {REFSEQ}")
        return 2
    require_sklearn()
    t0 = time.time()

    from dna_decode.glm.baselines import MarkovGenerator
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves
    from dna_decode.glm.falsifier import heldout_auroc

    print(f"geometry: length={LENGTH} limit={LIMIT} n={N_GENERATED} k={KMER_K} "
          f"bins={POSITIONAL_BINS} seeds={list(SEEDS)}")
    print(f"features: " + ", ".join(f"{m} {feature_count(m)}" for m in MODES))
    print(f"modes {list(MODES)} x orders {list(ORDERS)}")
    print()

    def measure(ref_dir: Path):
        """All (mode, order, seed) nulls plus the order-independent floors, for one substrate."""
        return _measure(ref_dir, build_natural_set, fit_eval_halves, MarkovGenerator, heldout_auroc)

    # floors: one per (mode, seed), independent of the null's order by construction
    floors, nulls = measure(REFSEQ)
    return _report(floors, nulls, measure, t0)


def _measure(ref_dir, build_natural_set, fit_eval_halves, MarkovGenerator, heldout_auroc):
    floors: dict[tuple[str, int], float] = {}
    nulls: dict[tuple[str, int, int], float] = {}
    for s in SEEDS:
        wins, _ = build_natural_set(ref_dir, length=LENGTH, limit=None, seed=s, kind=KIND)
        if len(wins) < LIMIT + N_GENERATED:
            raise SystemExit(f"REFUSE: only {len(wins)} windows at {ref_dir}")
        committed = wins[:LIMIT]
        fit, ev = fit_eval_halves(committed)
        used = set(fit) | set(ev)
        spare = []
        for w in wins[LIMIT:]:
            if len(spare) >= N_GENERATED:
                break
            if w.seq in used:
                continue
            spare.append(w.seq)
            used.add(w.seq)
        gens = {o: MarkovGenerator(order=o, seed=s).fit(fit).generate(N_GENERATED, LENGTH)
                for o in ORDERS}
        for mode in MODES:
            kw = {"k": KMER_K, "seed": s, "mode": mode}
            if mode in ("positional", "both"):
                kw["bins"] = POSITIONAL_BINS
            f = heldout_auroc(ev, spare, **kw)
            if f is not None:
                floors[(mode, s)] = f
            for o in ORDERS:
                v = heldout_auroc(ev, gens[o], **kw)
                if v is not None:
                    nulls[(mode, o, s)] = v
    return floors, nulls


def _report(floors, nulls, measure, t0: float):
    """Assemble, cross-check and write. Split out of `main` so the second substrate can reuse `measure`.

    `t0` is threaded in explicitly rather than read from an enclosing scope: these are module-level
    functions, not closures, so a wall-clock reference would raise NameError only at the very END of a
    multi-minute run -- after every measurement had been made and before anything was written.
    """
    # ---- assemble ----
    cells = {}
    for mode in MODES:
        for o in ORDERS:
            m = [nulls[(mode, o, s)] - floors[(mode, s)] for s in SEEDS
                 if (mode, o, s) in nulls and (mode, s) in floors]
            nl = [nulls[(mode, o, s)] for s in SEEDS if (mode, o, s) in nulls]
            cells[f"{mode}_order{o}"] = {
                "mode": mode, "order": o,
                "n_features": feature_count(mode),
                "null_mean": round(statistics.fmean(nl), 4),
                "null_sd": round(statistics.stdev(nl), 4) if len(nl) > 1 else None,
                "floor_mean": round(statistics.fmean([floors[(mode, s)] for s in SEEDS
                                                      if (mode, s) in floors]), 4),
                "margin": resolve(m),
                "is_committed_cell": (o == COMMITTED_ORDER and mode in ("kmer", "both")),
            }
            e = cells[f"{mode}_order{o}"]
            print(f"  {mode:10s} order {o}  null {e['null_mean']:.4f}  floor "
                  f"{e['floor_mean']:.4f} -> margin {e['margin']['mean']:+.4f} "
                  f"CI {e['margin']['ci95']}  resolvable={e['margin']['resolvable_above_zero']}")

    # ---- Q1: does positional alone match both? Seed-paired, at the committed order. ----
    o = COMMITTED_ORDER
    q1 = {}
    for a, b in (("both", "positional"), ("both", "kmer"), ("positional", "kmer")):
        va = [nulls[(a, o, s)] - floors[(a, s)] for s in SEEDS]
        vb = [nulls[(b, o, s)] - floors[(b, s)] for s in SEEDS]
        q1[f"{a}_minus_{b}"] = paired_diff(va, vb)
    both_vs_pos = q1["both_minus_positional"]
    q1["verdict"] = (
        "BOTH_BEATS_POSITIONAL_ALONE__UNION_EARNS_ITS_SIZE"
        if both_vs_pos["resolvably_nonzero"] and both_vs_pos["mean_diff"] > 0 else
        "POSITIONAL_ALONE_MATCHES_BOTH__KMER_BLOCK_ADDS_NOTHING_RESOLVABLE"
        if not both_vs_pos["resolvably_nonzero"] else
        "POSITIONAL_ALONE_BEATS_BOTH__THE_KMER_BLOCK_HURTS")

    # ---- Q2: how much does the null's strength depend on its order? ----
    q2 = {}
    for mode in MODES:
        ms = {o: cells[f"{mode}_order{o}"]["margin"]["mean"] for o in ORDERS}
        best = max(ms, key=lambda k: ms[k])
        worst = min(ms, key=lambda k: ms[k])
        vs_committed = {}
        for o in ORDERS:
            if o == COMMITTED_ORDER:
                continue
            va = [nulls[(mode, o, s)] - floors[(mode, s)] for s in SEEDS]
            vb = [nulls[(mode, COMMITTED_ORDER, s)] - floors[(mode, s)] for s in SEEDS]
            vs_committed[f"order{o}_minus_order{COMMITTED_ORDER}"] = paired_diff(va, vb)
        q2[mode] = {
            "margin_by_order": {str(k): round(v, 4) for k, v in ms.items()},
            "range": round(max(ms.values()) - min(ms.values()), 4),
            "strongest_null_order": worst,     # a STRONGER null sits CLOSER to the floor
            "weakest_null_order": best,
            "vs_committed_order": vs_committed,
            "n_orders_resolvably_different_from_committed": sum(
                1 for v in vs_committed.values() if v["resolvably_nonzero"]),
        }

    # ---- the structural cross-check, on a substrate chosen for being as different as available ----
    second = None
    if SECOND.exists() and (SECOND / "genome.fna").exists():
        print(f"\nstructural cross-check on {SECOND_LABEL}...")
        f2, n2 = measure(SECOND)
        prof = {}
        for mode in MODES:
            ms = {}
            for o in ORDERS:
                m = [n2[(mode, o, s)] - f2[(mode, s)] for s in SEEDS
                     if (mode, o, s) in n2 and (mode, s) in f2]
                ms[o] = round(statistics.fmean(m), 4)
            vs_committed = {}
            for o in ORDERS:
                if o == COMMITTED_ORDER:
                    continue
                va = [n2[(mode, o, s)] - f2[(mode, s)] for s in SEEDS]
                vb = [n2[(mode, COMMITTED_ORDER, s)] - f2[(mode, s)] for s in SEEDS]
                vs_committed[f"order{o}"] = paired_diff(va, vb)
            prof[mode] = {
                "margin_by_order": {str(k): v for k, v in ms.items()},
                "range": round(max(ms.values()) - min(ms.values()), 4),
                "n_orders_resolvably_different_from_committed": sum(
                    1 for v in vs_committed.values() if v["resolvably_nonzero"]),
            }
            print(f"   {mode:10s} " + "  ".join(f"o{o}={ms[o]:+.4f}" for o in ORDERS)
                  + f"   range {prof[mode]['range']:.4f}")
        held = (prof["kmer"]["range"] > 5 * prof["positional"]["range"]
                and prof["positional"]["n_orders_resolvably_different_from_committed"] == 0)
        second = {
            "substrate": SECOND_LABEL,
            "profile": prof,
            "kmer_over_positional_range_ratio": (
                round(prof["kmer"]["range"] / prof["positional"]["range"], 2)
                if prof["positional"]["range"] else None),
            "prediction_held": bool(held),
            "verdict": ("ORDER_DEPENDENCE_IS_STRUCTURAL__REPLICATES_ACROSS_SUBSTRATES" if held else
                        "ORDER_DEPENDENCE_DID_NOT_REPLICATE__MECHANISM_CLAIM_UNSUPPORTED"),
        }
        print(f"   -> {second['verdict']}")

    rec = {
        "record": "glm-instrument-levers-v1",
        "date": date.today().isoformat(),
        "questions": {
            "q1": "Does the positional feature block alone match `both`, or do the two blocks add?",
            "q2": ("How much of the Markov null's strength is a property of its ORDER, which is fixed "
                   "at 3 by a CLI default and has never been swept?"),
        },
        "geometry": {"window_length_bp": LENGTH, "limit_natural": LIMIT,
                     "n_generated_per_candidate": N_GENERATED, "natural_set_arg": KIND,
                     "kmer_k": KMER_K, "positional_bins": POSITIONAL_BINS,
                     "seeds": list(SEEDS), "committed_markov_order": COMMITTED_ORDER,
                     "feature_counts": {m: feature_count(m) for m in MODES}},
        "floor_is_order_independent_by_construction": (
            "The perfect-generator floor scores real sequence against real sequence and never involves "
            "the null, so it cannot depend on the null's order. It is computed once per (mode, seed) and "
            "reused, which keeps any apparent order effect entirely in the null term."),
        "cells": cells,
        "q1_feature_blocks": q1,
        "q2_null_order": q2,
        "structural_prediction": STRUCTURAL_PREDICTION,
        "structural_cross_check_second_substrate": second,
        "primary_substrate": REF_LABEL,
        "honest_limits": [
            "The positional block is `bins * 4**min(k,2)` = 96 features, NOT bins * 4**k. The clamp "
            "lives inside `positional_features` rather than at the call site, so `both` is 352 features "
            "and not 1792. Every p-against-n reading here depends on that.",
            "Ten seeds per cell. Each comparison between cells is seed-PAIRED, which is what makes a "
            "difference of this size readable at n=10; the per-cell margins themselves are noisier.",
            "A higher-order Markov chain fitted on 600 windows of 150bp approaches memorising its "
            "training set, so a high-order null is not simply a 'better' null -- part of any change is "
            "overfitting the fit half rather than modelling the substrate. This sweep measures the "
            "effect and does NOT separate those two causes.",
            "The order sweep itself runs on E. coli K-12; the structural cross-check runs the same "
            "sweep on ONE Campylobacter jejuni draft assembly. Two substrates test whether the pattern "
            "is structural and do not characterise how it varies across organisms.",
            "One window length, one geometry, one null family.",
            "A winnable region is necessary, NOT sufficient: all three feature sets see LOCAL "
            "composition, so the positional block adds WHERE composition sits, not long-range "
            "organisation.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    if GRID.exists():
        rec["note_on_the_committed_grid"] = (
            "The committed grid scores modes kmer and both at order 3 only. The positional-alone column "
            "and every non-3 order here are NEW cells, not a re-reading of that grid.")
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print()
    print(f"Q1 (order {COMMITTED_ORDER}, seed-paired margins): {q1['verdict']}")
    for k, v in q1.items():
        if k == "verdict":
            continue
        print(f"   {k:26s} {v['mean_diff']:+.4f} CI {v['ci95']} "
              f"resolvable={v['resolvably_nonzero']} ({v['n_positive']}/{v['n']} seeds positive)")
    print()
    print("Q2 null strength against its order (margin above floor; SMALLER = stronger null):")
    for mode in MODES:
        q = q2[mode]
        print(f"   {mode:10s} " + "  ".join(f"o{o}={q['margin_by_order'][str(o)]:+.4f}"
                                            for o in ORDERS))
        print(f"              range {q['range']:.4f}, strongest null at order "
              f"{q['strongest_null_order']}, "
              f"{q['n_orders_resolvably_different_from_committed']}/{len(ORDERS) - 1} orders "
              f"resolvably different from the committed {COMMITTED_ORDER}")
    print(f"\n{rec['wall_seconds']}s -> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
