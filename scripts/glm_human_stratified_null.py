"""The ACTIVITY-STRATIFIED null -- and a measurement of the oracle's own ceiling.

WHAT IT IS. An order-3 Markov chain fit SEPARATELY per activity quintile of the generator's training
folds, sampled at the quintile matching each requested target. It is the null that actually bites an
activity-conditioned generator: if a 4-mer model conditioned on the same target reproduces the
activity-sequence relationship, the generator has added nothing over 4-mer statistics however well it
hits the target in absolute terms.

WHY IT IS ALSO A DIAGNOSTIC OF THE JUDGE, WHICH IS THE MORE IMPORTANT USE. The oracle's power IS
k-mer composition -- measured on held-out folds, HepG2 GC 0.2402 / 4-mer 0.4245 / 6-mer 0.4598 /
positional 0.2733. So a generator matching per-stratum 6-mer statistics scores as well as this oracle
can detect, and "designs to spec" AS THIS JUDGE MEASURES IT is substantially a composition-matching
task. The k-mer saturation problem did not disappear in the pivot; it moved from the discriminator to
the scorer.

That makes this run load-bearing in a way a baseline usually is not. If the stratified chain reaches
high target monotonicity, it PROVES the oracle's activity signal is k-mer-expressible, and the
remaining room for the generator is only the gap between an order-3 chain (4-mer statistics) and a
6-mer readout. Narrow, real, and worth knowing BEFORE a generator is trained rather than after.

THREE ARMS, so the result is interpretable in both directions:
  * stratified  -- order-3 chain per activity quintile. The null the bar names.
  * pooled      -- one order-3 chain over all training folds, sampled identically. Has NO access to
                   the target, so its monotonicity is the no-information floor for a sequence model
                   that still produces realistic composition.
  * uniform     -- i.i.d. uniform sequence. The chance floor, and the sanity gate: if the oracle
                   assigns uniform sequence an activity trend, the oracle is reading an artifact.

LEAK DISCIPLINE. Chains are fit on folds 1-6 only. The oracle is fit on folds 7-8 only. Nothing is
FIT on the held-out folds -- but the real_ceiling arm READS them (it scores real held-out sequence
stratified by its own measured activity), so they are not untouched. An earlier version of this
docstring claimed nothing here touches folds 9-10, which its own real_ceiling arm refutes; the
reservation artifact counted that read as a spent query and was right. Pass --held-out-folds 9 to
keep the read off the FINAL_ONLY shard. So the
null is scored by a judge that never saw the null's training data, exactly as the generator will be.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dna_decode.glm.baselines import MarkovGenerator, UniformGenerator  # noqa: E402
from dna_decode.glm.lentimpra import load_cell  # noqa: E402
from scripts.glm_human_activity_oracle import (  # noqa: E402
    ADAPTOR, HELD_OUT_FOLDS, LOCAL_DIR, ORACLE_TRAIN_FOLDS, WORKING_HELD_OUT_FOLDS,
    OraclePreconditionFailed, assert_oracle_fixture, build, fold_tag, parse_fold_set,
    spearman,
)

GENERATOR_FOLDS = {1, 2, 3, 4, 5, 6}
#: NO module-level OUT: the path depends on the held-out fold set and is built in main() via
#: fold_tag(), so a fold-9 result cannot land on the folds-9+10 record. A default constant here
#: was printed while a different path was written -- removed rather than left as a trap.
ELEMENT_BP = 200
ORDER = 3
ORACLE_MODE = "kmer6"      # the winning oracle arm on held-out folds, both cells


def fit_oracle(seqs, ys, *, mode: str = ORACLE_MODE, alpha: float = 1.0):
    """Return a predict(seqs) -> list[float]. Standardisation is fit on TRAIN only."""
    import numpy as np
    A = np.asarray(build(seqs, mode, 10), dtype=np.float64)
    mu, sd = A.mean(0), A.std(0)
    sd[sd == 0] = 1.0
    A = np.hstack([(A - mu) / sd, np.ones((A.shape[0], 1))])
    reg = alpha * np.eye(A.shape[1])
    reg[-1, -1] = 0.0
    w = np.linalg.solve(A.T @ A + reg, A.T @ np.asarray(ys, dtype=np.float64))

    def predict(xs):
        B = np.asarray(build(xs, mode, 10), dtype=np.float64)
        B = np.hstack([(B - mu) / sd, np.ones((B.shape[0], 1))])
        return (B @ w).tolist()
    return predict


def quintiles(recs, n: int = 5):
    """Equal-count strata by activity. Equal-COUNT, not equal-width: activity is skewed, and
    equal-width bins would leave the top stratum nearly empty and its chain unfittable."""
    s = sorted(recs, key=lambda r: r.activity)
    out, step = [], len(s) / n
    for i in range(n):
        out.append(s[int(i * step):int((i + 1) * step)])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cells", default="K562,HepG2")
    ap.add_argument("--per-stratum", type=int, default=240)
    ap.add_argument("--strata", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ceiling-seeds", type=int, default=10,
                    help="how many seeds to re-draw the matched-n ceiling with, to give the "
                         "verdict rule's DENOMINATOR a measured noise scale rather than a "
                         "single point")
    ap.add_argument("--held-out-folds", default=None,
                    help="comma-separated held-out fold set the real_ceiling arm READS. Default "
                         "9,10 reproduces the committed measurement; pass 9 to keep the ceiling "
                         "off the FINAL_ONLY shard, which is required because the verdict rule is "
                         "a fraction OF this ceiling.")
    a = ap.parse_args()

    # Same third precondition: this script FITS and CALLS the oracle, so a changed feature
    # pipeline would move the ceiling -- the denominator of every fraction reported below.
    try:
        precondition = assert_oracle_fixture()
    except OraclePreconditionFailed as exc:
        print(f"PRECONDITION FAILED: {exc}", file=sys.stderr)
        return 5
    print("precondition: oracle-feature-fixture-v1 PASSED")

    try:
        held_folds = parse_fold_set(a.held_out_folds, HELD_OUT_FOLDS)
    except ValueError as exc:
        print(f"REFUSE: {exc}", file=sys.stderr)
        return 3
    out_path = ROOT / "wiki" / f"glm_human_stratified_null{fold_tag(held_folds)}_2026-11-22.json"

    res = {
        "record": "glm-human-stratified-null-v1",
        "date": "2026-11-22",
        "bar": "wiki/human_design_bar_2026-11-22.json",
        "purpose": ("the null an activity-conditioned generator must beat, AND a measurement of "
                    "whether the oracle's activity signal is k-mer-expressible at all"),
        "oracle_mode": ORACLE_MODE,
        "precondition": precondition,
        "geometry": f"element_{ELEMENT_BP}nt_adaptors_stripped",
        "leak": {"chains_fit_on_folds": sorted(GENERATOR_FOLDS),
                 "oracle_fit_on_folds": sorted(ORACLE_TRAIN_FOLDS),
                 "folds_read_by_real_ceiling_arm_never_fitted": sorted(held_folds),
                 "folds_not_read_at_all": sorted(set(range(1, 11)) - GENERATOR_FOLDS
                                                 - ORACLE_TRAIN_FOLDS - held_folds),
                 "note": ("the previous key was named folds_untouched, which was FALSE -- the "
                          "real_ceiling arm reads these folds (scored, never fitted). Renamed "
                          "rather than silently kept.")},
        "cells": {},
    }
    t0 = time.time()

    for cell in [c.strip() for c in a.cells.split(",") if c.strip()]:
        recs, _ = load_cell(cell, LOCAL_DIR)
        gen_recs = [r for r in recs if r.fold in GENERATOR_FOLDS]
        ora_recs = [r for r in recs if r.fold in ORACLE_TRAIN_FOLDS]
        strip = lambda s: s[ADAPTOR:len(s) - ADAPTOR]  # noqa: E731
        print(f"{cell}: chain-fit {len(gen_recs)} | oracle-fit {len(ora_recs)}")

        oracle = fit_oracle([strip(r.seq) for r in ora_recs], [r.activity for r in ora_recs])

        # Sanity: the oracle must rank REAL held-out-fold sequence sensibly, or nothing below means
        # anything. Uses folds 7-8's own in-sample rank only as a smoke check, not as a result.
        strata = quintiles(gen_recs, a.strata)
        real_means = [round(statistics.fmean(r.activity for r in s), 4) for s in strata]
        print(f"   real stratum mean activity: {real_means}")

        arms: dict[str, dict] = {}
        pooled = MarkovGenerator(order=ORDER, seed=a.seed).fit(
            [strip(r.seq) for r in gen_recs])
        uni = UniformGenerator(seed=a.seed)

        held = [r for r in recs if r.fold in held_folds]
        held_strata = quintiles(held, a.strata)

        for arm in ("real_ceiling", "stratified", "pooled", "uniform"):
            targets, preds, counts = [], [], []
            for qi, stratum in enumerate(strata):
                if arm == "real_ceiling":
                    # REAL held-out sequence, stratified by its OWN measured activity. The oracle
                    # never saw the held-out folds, so this is the monotonicity a perfect generator
                    # could reach against this judge -- bounded by the oracle's own ~0.5 precision.
                    #
                    # MATCHED-n RANDOM sample, NOT the first n. `quintiles()` returns each stratum
                    # sorted by activity ASCENDING, so [:n] took the n LOWEST-activity members of
                    # every stratum and truncated the top quintile's tail, where most of the
                    # activity signal lives. Measured 2026-11-22: that understated the ceiling by
                    # 0.1095 in K562 (-4.14 sd of the random-draw spread) and 0.0300 in HepG2
                    # (-2.74 sd), while random-n and ALL-records agree to within 0.4 sd -- so the
                    # sample SIZE was never the problem, the ORDER was. Direction matters: the
                    # ceiling is the DENOMINATOR of stratified_frac_of_ceiling, so understating it
                    # OVERSTATES every generator's fraction of achievable signal.
                    #
                    # n is kept equal to the generated arms' n deliberately, so the ratio compares
                    # like with like; the all-records estimate ships beside it as a diagnostic.
                    sel = random.Random(a.seed * 1000 + qi).sample(
                        held_strata[qi], min(a.per_stratum, len(held_strata[qi])))
                    gen = [strip(r.seq) for r in sel]
                elif arm == "stratified":
                    ch = MarkovGenerator(order=ORDER, seed=a.seed + qi).fit(
                        [strip(r.seq) for r in stratum])
                    gen = ch.generate(a.per_stratum, ELEMENT_BP)
                elif arm == "pooled":
                    gen = pooled.generate(a.per_stratum, ELEMENT_BP)
                else:
                    gen = uni.generate(a.per_stratum, ELEMENT_BP)
                p = oracle(gen)
                counts.append(len(p))
                targets.extend([qi] * len(p))
                preds.extend(p)
            rho = spearman(targets, preds)
            by_q, off = [], 0
            for cnt in counts:
                by_q.append(round(statistics.fmean(preds[off:off + cnt]), 4) if cnt else None)
                off += cnt
            arms[arm] = {
                "target_monotonicity_spearman": round(rho, 4),
                "predicted_activity_by_target": by_q,
                "top_minus_bottom_lift": round(by_q[-1] - by_q[0], 4),
            }
            print(f"   {arm:11s} rho {rho:+.4f} | lift {arms[arm]['top_minus_bottom_lift']:+.4f} "
                  f"| by target {by_q}")

        # ---- the denominator deserves a noise scale, not a single point ----
        def _ceiling_rho(pick):
            t, p = [], []
            for qi, st in enumerate(held_strata):
                pr = oracle([strip(r.seq) for r in pick(st, qi)])
                t.extend([qi] * len(pr))
                p.extend(pr)
            return round(spearman(t, p), 4)

        ceil_all = _ceiling_rho(lambda st, qi: st)
        ceil_spread = [
            _ceiling_rho(lambda st, qi, s=s: random.Random(s * 1000 + qi).sample(
                st, min(a.per_stratum, len(st))))
            for s in range(max(0, a.ceiling_seeds))
        ]
        ceil_diag = {
            "matched_n_rho_by_seed": ceil_spread,
            "matched_n_mean": round(statistics.fmean(ceil_spread), 4) if ceil_spread else None,
            "matched_n_sd": (round(statistics.stdev(ceil_spread), 4)
                             if len(ceil_spread) > 1 else None),
            "all_held_out_records_rho": ceil_all,
            "n_held_out_records": len(held),
            "why": ("the ceiling is the DENOMINATOR of stratified_frac_of_ceiling, so its own "
                    "sampling noise propagates into every fraction quoted against it. Reported as "
                    "a measured spread because the bar refuses frozen numeric tolerances."),
            "superseded_estimator_note": (
                "a first-n subsample of an activity-SORTED stratum was used before 2026-11-22; it "
                "understated this quantity (K562 by 0.1095 = -4.14 sd, HepG2 by 0.0300 = -2.74 sd) "
                "by truncating each stratum's top edge"),
        }

        res["cells"][cell] = {
            "n_chain_fit": len(gen_recs), "n_oracle_fit": len(ora_recs),
            "ceiling_diagnostics": ceil_diag,
            "real_stratum_mean_activity": real_means,
            "real_top_minus_bottom": round(real_means[-1] - real_means[0], 4),
            "arms": arms,
        }

    # ---- verdict: what did we learn about the JUDGE, not about a generator ----
    strat = [res["cells"][c]["arms"]["stratified"]["target_monotonicity_spearman"]
             for c in res["cells"]]
    pool = [res["cells"][c]["arms"]["pooled"]["target_monotonicity_spearman"]
            for c in res["cells"]]
    unif = [abs(res["cells"][c]["arms"]["uniform"]["target_monotonicity_spearman"])
            for c in res["cells"]]

    # CEILING-RELATIVE, not absolute. Measured 2026-11-22: the frozen absolute bar of 0.30 sat at
    # 80% of a ceiling of 0.3729, so it labelled a null capturing 54% of achievable signal "weak".
    # This is the FOURTH threshold this session to misdescribe its data for the same reason -- an
    # absolute bar set before the scale was measured. A fraction of the measured ceiling has no such
    # failure mode, and the real-sequence arm is what makes it computable.
    def _pcts(d):
        """'K562 53.0%, HepG2 55.6%' -- named per cell and rounded. Cells are NEVER pooled on this
        substrate (they have already been measured contradicting each other), so a bare list hides
        exactly the thing a reader needs: which cell each number belongs to."""
        return ", ".join(f"{c} {v:.1%}" if v is not None else f"{c} n/a" for c, v in d.items())

    def _nums(d):
        return ", ".join(f"{c} {v:.4f}" for c, v in d.items())

    def _ceil_point(cell):
        """The ALL-HELD-OUT-RECORDS ceiling, not a single matched-n draw.

        Measured 2026-11-22: matched-n draws average 0.4923 (sd 0.0254) against an all-records
        0.4912 in K562 -- agreement to 0.04 sd -- so sample SIZE does not bias this statistic,
        while a single draw carries ~5% relative noise. Every fraction below is divided by this
        number, so a seed-dependent denominator would make the verdict rule seed-dependent.
        """
        return res["cells"][cell]["ceiling_diagnostics"]["all_held_out_records_rho"]

    def _ceil_sd(cell):
        return res["cells"][cell]["ceiling_diagnostics"]["matched_n_sd"]

    def frac(cell, arm):
        c = _ceil_point(cell)
        v = res["cells"][cell]["arms"][arm]["target_monotonicity_spearman"]
        return round(v / c, 4) if c and c > 0 else None

    fr = {c: frac(c, "stratified") for c in res["cells"]}
    ceil = {c: _ceil_point(c) for c in res["cells"]}
    pl = {c: frac(c, "pooled") for c in res["cells"]}

    def _frac_sd(cell):
        f, c, s = fr[cell], ceil[cell], _ceil_sd(cell)
        if f is None or not c or s is None:
            return None
        return round(abs(f) * (s / c), 4)        # first-order, denominator only

    res["ceiling_relative"] = {
        "ceiling_spearman": ceil,
        "ceiling_estimator": ("all held-out records per stratum -- the lowest-variance estimate. "
                              "Per-cell matched-n draws and their spread are in "
                              "cells.*.ceiling_diagnostics."),
        "ceiling_sd_from_matched_n_draws": {c: _ceil_sd(c) for c in res["cells"]},
        "stratified_frac_of_ceiling": fr,
        "pooled_frac_of_ceiling": pl,
        "frac_sd_propagated_from_ceiling_only": {c: _frac_sd(c) for c in res["cells"]},
        "frac_sd_note": ("first-order propagation of the DENOMINATOR's spread only. The null arm's "
                         "own seed spread is NOT included, so this is a floor on each fraction's "
                         "uncertainty and not its total."),
    }
    vals = [v for v in fr.values() if v is not None]

    if max(unif) > 0.10:
        res["verdict"] = "ORACLE_READS_AN_ARTIFACT"
        res["verdict_detail"] = (
            f"Uniform i.i.d. sequence shows a target trend (|rho| up to {max(unif):.4f}) although it "
            "carries no target information. The oracle is reading the sampling rather than activity. "
            "REFUSING to use this null until that is explained.")
    elif min(ceil.values()) <= 0.10:
        res["verdict"] = "NO_USABLE_CEILING"
        res["verdict_detail"] = (
            f"Real sequence stratified by its OWN measured activity reaches only {_nums(ceil)} against "
            "this "
            "oracle, so there is almost no achievable signal to compete for and no arm is "
            "interpretable. The oracle, not the null, is the thing to fix.")
    elif min(vals) >= 0.80:
        res["verdict"] = "ORACLE_SIGNAL_IS_KMER_EXPRESSIBLE"
        res["verdict_detail"] = (
            f"An order-3 chain conditioned only on activity quintile captures {_pcts(fr)} of the "
            f"ACHIEVABLE monotonicity (ceiling {_nums(ceil)}). The oracle's activity signal is "
            "essentially "
            "4-mer-expressible and a generator has almost no room against this judge.")
    elif min(vals) >= 0.25:
        res["verdict"] = "KMER_NULL_IS_SUBSTANTIAL_ROOM_REMAINS"
        res["verdict_detail"] = (
            f"The stratified 4-mer null captures {_pcts(fr)} of the achievable monotonicity "
            f"(ceiling {_nums(ceil)} from real sequence stratified by its own measured activity), "
            f"while an UNCONDITIONED chain captures {_pcts(pl)} -- so the signal is "
            "genuinely coming from the conditioning and not from the generation process. A real "
            "bar with real room left. A generator must be quoted against BOTH the null and the "
            "ceiling; beating the null while far below the ceiling is a partial result.")
    else:
        res["verdict"] = "KMER_NULL_IS_A_LOW_BAR"
        res["verdict_detail"] = (
            f"The stratified null captures only {_pcts(fr)} of achievable monotonicity "
            f"(ceiling {_nums(ceil)}). "
            "Conditioning an order-3 chain barely moves this oracle, so most of the signal is NOT "
            "4-mer-expressible and the null is a low bar.")

    res["seconds"] = round(time.time() - t0, 1)
    print("\n" + "=" * 74)
    print(f"VERDICT: {res['verdict']}\n  {res['verdict_detail']}")
    out_path.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\n-> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
