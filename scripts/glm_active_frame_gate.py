"""The published conv on the ACTIVE frame, under its own frozen bar.

WHY A SEPARATE BAR AND NOT AN EXTENSION
----------------------------------------
`scripts/glm_encoder_gate.py::PREREGISTERED` is a frozen pre-registration carrying
`reconcile_anchor = {"arm": "gc", "seed": 0, "expected": 0.3000, "tolerance": 0.01}`. The active
frame breaks that anchor BY CONSTRUCTION -- GC scores near zero there, not 0.3000 -- so sharing
the dict would both retroactively rewrite a commitment and force `INDETERMINATE_RECONCILE_FAILED`
on every run. That file is imported READ-ONLY here and a test pins it byte-identical to HEAD.

THE ORDERING IS LOAD-BEARING
----------------------------
Phase 1 scores GC and `encoder-multi` on the POOLED frame and asserts the published anchor
before the frame changes. The pooled frame is the one carrying a published number, so it is the
only place where "am I reading the same data under the same split?" can be answered. Validating
on the frame under question would prove nothing. Phase 2 then scores the active frame. Phase 3
reports the margin inside equal-sized rank bins, so the headline does not rest on the single
0.5778 cut -- but the binned table is a SHAPE, never a second bar, and a test pins that it
cannot move the verdict.

The bar below is frozen in source before the run and it says twice that it is not blind: the
2026-11-16 scratch measurement motivated the thresholds. It is also frozen AFTER the mechanism
analysis, which this project's own standing rule requires -- the two bars mis-specified this
quarter were both frozen before the analysis that said what a correct fix could achieve.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics as st
import sys
import time
from datetime import date
from pathlib import Path

# Thread pinning MUST precede the numeric imports: a recorded BLAS-oversubscription instance
# cost a 2.34x throughput error and a wrong projection.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.glm.active_frame import threshold_audit  # noqa: E402
from dna_decode.glm.encoder import DEFAULT_WIDTHS, REGISTERED_PROTOCOL, fit_encoder  # noqa: E402
from dna_decode.glm.genomewide import leave_peak_out, load_peak_tiles, spearman  # noqa: E402
# single-sourced from Step 2 rather than copied: the repo has a recorded duplication trap, and
# relocating it into dna_decode/glm/active_frame.py would change Step 1's declared contract
# after that step was verified.
from scripts.active_frame_rescore import _binned_margin  # noqa: E402
from scripts.glm_encoder_gate import PREREGISTERED as POOLED_PREREGISTERED  # noqa: E402
from scripts.glm_encoder_gate import ridge_gc  # noqa: E402

VERDICTS = (
    "FRAME_EXPLAINS_THE_GAP",
    "FRAME_CONTRIBUTES",
    "FRAME_DOES_NOT_EXPLAIN",
    "INDETERMINATE_NULL_NOT_CLEAN",
    "INDETERMINATE_POOLED_RECONCILE_FAILED",
)

#: FROZEN IN SOURCE before the run. Thresholds are ASSERTED; the ceilings that contextualise
#: them are DERIVED from the tiles' own replicates.
PREREGISTERED_ACTIVE_FRAME = {
    "question": "on the ACTIVE frame -- the 14,913 tiles above a pure threshold on the label at "
                "0.5778 -- does the published conv encoder (widths 6/20/75, 13,217 params) beat "
                "a single GC number by MORE than the +0.0122 it manages on the pooled frame?",
    "primary_metric": "spearman",
    "comparator": "gc",
    "candidate_arms": ["encoder-multi"],
    "null_arm": "encoder-shuffled",
    "registered_seeds": list(range(10)),
    "frame": "active-only; leave-peak-out applied WITHIN the active subset, so train and test "
             "are both active and the comparison matches the 2026-11-16 scratch design whose "
             "+0.1820 margin calibrated the thresholds below",
    "device": "cpu",
    "protocol": "dna_decode.glm.encoder.REGISTERED_PROTOCOL, unchanged -- architecture-matched "
                "to the published +0.0122 rather than retuned for this frame",
    "pooled_reconcile_anchor": {
        "arm": "gc", "seed": 0, "expected": 0.3000, "tolerance": 0.01,
        "source": "wiki/glm_genomewide_oracle_2026-10-07.json",
        "why_pooled": "the anchor exists only on the pooled frame; asserting it there BEFORE "
                      "changing frames is the only way to prove the data and split are right",
    },
    "expected_n_parameters": 13217,
    "explains_margin": 0.15,
    "explains_min_seeds": 8,
    "contributes_margin": 0.05,
    "null_max_abs": 0.05,
    "not_a_blind_pre_registration": True,
    "why_not_blind": "the 2026-11-16 scratch measurement (4-mer ridge, +0.1820 active-frame "
                     "margin on 10 of 10 seeds) motivated these thresholds, so they were set "
                     "with a directional expectation in view. Recorded in the same field name "
                     "wiki/asgiven_arm_result_2026-11-02.json uses.",
    "frozen_after_the_mechanism_analysis": True,
    "why_frozen_after": "this project's standing rule: a bar frozen BEFORE the analysis that "
                        "says what a correct fix can achieve is mis-specified, and the two bars "
                        "mis-specified that way this quarter were both frozen too early. The "
                        "mechanism here -- GC is a detector, a learned representation is not -- "
                        "was analysed first and the thresholds set against it.",
    "does_not_extend": "scripts/glm_encoder_gate.py::PREREGISTERED is NOT edited. Its reconcile "
                       "anchor is broken by the active frame by construction, so sharing it "
                       "would rewrite a commitment and force INDETERMINATE_RECONCILE_FAILED.",
    "verdicts": {
        "FRAME_EXPLAINS_THE_GAP": "active-frame margin >= +0.15 on at least 8 of 10 seeds -- the "
                                  "pooled +0.0122 was a frame artifact to first order",
        "FRAME_CONTRIBUTES": "margin > +0.05 with a paired CI excluding zero -- the frame "
                             "carries part of the gap but not the whole of it",
        "FRAME_DOES_NOT_EXPLAIN": "margin <= +0.05, or the CI includes zero -- the frame "
                                  "hypothesis is tested and rejected",
        "INDETERMINATE_NULL_NOT_CLEAN": "the shuffled-label arm itself scores on this frame, so "
                                        "the pipeline fabricates signal and NOTHING is graded",
        "INDETERMINATE_POOLED_RECONCILE_FAILED": "gc seed 0 did not reproduce 0.3000 on the "
                                                 "POOLED frame, so this script is not reading "
                                                 "the published data and nothing is interpretable",
    },
    #: Where the parent plan's substrate-re-framing question goes. Fixed before the run so the
    #: framing cannot be read off the result. RECOMMENDATION ONLY -- adopting any of these edits
    #: a ledger criterion, which is a user authority call.
    "framing_recommendation": {
        "FRAME_EXPLAINS_THE_GAP": "DECLARE TWO CRITERIA -- detection and strength -- and score "
                                  "future bacterial arms on both. A single pooled number would "
                                  "keep averaging two tasks whose answers differ in sign.",
        "FRAME_CONTRIBUTES": "KEEP the pooled frame as primary and ADD the active-frame margin "
                             "as a required disclosure beside it. Not enough separation to "
                             "justify re-basing every future number.",
        "FRAME_DOES_NOT_EXPLAIN": "KEEP THE POOLED FRAME AND DECLARE NOTHING. This is the "
                                  "load-bearing entry: re-framing onto a subset with no measured "
                                  "headroom would score every future arm on a frame strictly "
                                  "worse than today's.",
        "INDETERMINATE_NULL_NOT_CLEAN": "NO FRAMING CHANGE. An ungraded run cannot inform a "
                                        "scope decision.",
        "INDETERMINATE_POOLED_RECONCILE_FAILED": "NO FRAMING CHANGE. Fix the read path first.",
    },
    "thresholds_are": "ASSERTED; the active-frame ceiling (0.9690 from the tiles' own "
                      "replicates) is DERIVED and contextualises them without setting them",
}

EXIT_OK = 0
EXIT_INDETERMINATE = 3
EXIT_DATA = 5


def compute_verdict(arms: dict, *, bar: dict = PREREGISTERED_ACTIVE_FRAME) -> tuple[str, str]:
    """Pure, and ORDERED: the two indeterminates take precedence over any graded outcome.

    `arms` needs `pooled_reconcile_ok`, `null_max_abs`, `active_margin_per_seed` and
    `active_margin_ci`. Nothing about the binned table is read, which a test pins.
    """
    if not arms.get("pooled_reconcile_ok", False):
        return ("INDETERMINATE_POOLED_RECONCILE_FAILED",
                f"gc seed 0 on the POOLED frame did not reproduce "
                f"{bar['pooled_reconcile_anchor']['expected']} within "
                f"{bar['pooled_reconcile_anchor']['tolerance']}")
    nm = float(arms.get("null_max_abs", 1.0))
    if nm >= bar["null_max_abs"]:
        return ("INDETERMINATE_NULL_NOT_CLEAN",
                f"shuffled-label arm reached |{nm:.4f}| against a ceiling of "
                f"{bar['null_max_abs']}")

    margins = list(arms["active_margin_per_seed"])
    med = st.median(margins)
    wins_at_bar = sum(1 for m in margins if m >= bar["explains_margin"])
    ci_lo, ci_hi = arms["active_margin_ci"]
    if med >= bar["explains_margin"] and wins_at_bar >= bar["explains_min_seeds"]:
        return ("FRAME_EXPLAINS_THE_GAP",
                f"median active-frame margin {med:+.4f} >= {bar['explains_margin']} and "
                f"{wins_at_bar}/{len(margins)} seeds clear it")
    if med > bar["contributes_margin"] and (ci_lo > 0 or ci_hi < 0):
        return ("FRAME_CONTRIBUTES",
                f"median margin {med:+.4f} > {bar['contributes_margin']} with CI "
                f"[{ci_lo:+.4f}, {ci_hi:+.4f}] excluding zero")
    return ("FRAME_DOES_NOT_EXPLAIN",
            f"median margin {med:+.4f} against {bar['contributes_margin']} with CI "
            f"[{ci_lo:+.4f}, {ci_hi:+.4f}]")


def framing_recommendation(verdict: str, *, bar: dict = PREREGISTERED_ACTIVE_FRAME) -> str:
    """The pre-committed framing recommendation. Refuses an unknown verdict rather than guessing."""
    try:
        return bar["framing_recommendation"][verdict]
    except KeyError:
        raise KeyError(
            f"no pre-committed framing recommendation for verdict {verdict!r}; the map must be "
            f"TOTAL over the verdict set, and adding a verdict without one is the drift this "
            f"refusal exists to catch"
        ) from None


def _paired_bootstrap(margins, n_boot=1000, seed=707):
    rng = np.random.default_rng(seed)
    arr = np.asarray(margins, dtype=float)
    draws = arr[rng.integers(0, arr.size, size=(n_boot, arr.size))].mean(axis=1)
    return (round(float(np.percentile(draws, 2.5)), 4),
            round(float(np.percentile(draws, 97.5)), 4))


def _fit_one(train, test, *, seed, shuffle):
    t0 = time.time()
    res = fit_encoder(train, test, widths=DEFAULT_WIDTHS, protocol=REGISTERED_PROTOCOL,
                      seed=seed, device="cpu", shuffle_labels=shuffle)
    return {
        "spearman": round(float(res["spearman"]), 4),
        "n_parameters": int(res["n_parameters"]),
        "epochs_ran": int(res["epochs_ran"]),
        "device": res["device"],
        "wall_seconds": round(time.time() - t0, 1),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache", default="data/mpra")
    ap.add_argument("--seeds", type=int, default=10,
                    help="how many registered seeds to run; fewer is a SMOKE run and is "
                         "stamped as such rather than graded")
    ap.add_argument("--n-bins", type=int, default=10)
    ap.add_argument("--subsample", type=int, default=0,
                    help="smoke only: cap the active tile count")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    t_start = time.time()
    threads = {v: os.environ.get(v) for v in
               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS")}
    print(f"threads pinned: {threads}")

    try:
        loaded, _st = load_peak_tiles(args.cache, allow_download=False)
    except Exception as e:  # noqa: BLE001
        print(f"INDETERMINATE_SUBSTRATE_UNREADABLE: {e}")
        return EXIT_DATA
    tiles = [t for t in loaded if t.category == "tile" and len(t.seq) == 150]
    audit = threshold_audit(tiles)
    print(f"loaded {len(tiles)} tiles; {audit['n_active']} active / {audit['n_inactive']} "
          f"inactive, cut {audit['implied_cut']}, overlap {audit['overlap_either_side']}")

    seeds = list(range(args.seeds))
    smoke = args.seeds < len(PREREGISTERED_ACTIVE_FRAME["registered_seeds"]) or args.subsample > 0

    # ---------------- PHASE 1: the pooled anchor, before any frame change ----------------
    print("\nPHASE 1 -- pooled frame, asserting the published anchor")
    anchor = PREREGISTERED_ACTIVE_FRAME["pooled_reconcile_anchor"]
    ptr, pte, _h = leave_peak_out(tiles, frac=0.25, seed=anchor["seed"])
    gc0, _pred = ridge_gc(ptr, pte)
    reconcile_ok = abs(gc0 - anchor["expected"]) <= anchor["tolerance"]
    print(f"  gc seed {anchor['seed']} pooled = {gc0:+.4f} "
          f"(expected {anchor['expected']} +/- {anchor['tolerance']}) -> "
          f"{'OK' if reconcile_ok else 'FAILED'}")
    pooled_conv = None
    if reconcile_ok:
        pooled_conv = _fit_one(ptr, pte, seed=anchor["seed"], shuffle=False)
        print(f"  encoder-multi seed {anchor['seed']} pooled = {pooled_conv['spearman']:+.4f} "
              f"({pooled_conv['n_parameters']} params, {pooled_conv['epochs_ran']} epochs, "
              f"{pooled_conv['wall_seconds']}s)")
        if pooled_conv["n_parameters"] != PREREGISTERED_ACTIVE_FRAME["expected_n_parameters"]:
            print(f"  PARAMETER COUNT MISMATCH: {pooled_conv['n_parameters']} != "
                  f"{PREREGISTERED_ACTIVE_FRAME['expected_n_parameters']}")
            reconcile_ok = False

    # ---------------- PHASE 2: the active frame ----------------
    active = [t for t in tiles if t.active == "active"]
    if args.subsample:
        active = active[: args.subsample]
    print(f"\nPHASE 2 -- active frame, {len(active)} tiles, {len(seeds)} seed(s)"
          f"{' [SMOKE]' if smoke else ''}")

    per_seed = {}
    if reconcile_ok:
        for s in seeds:
            tr, te, held = leave_peak_out(active, frac=0.25, seed=s)
            g, gpred = ridge_gc(tr, te)
            cand = _fit_one(tr, te, seed=s, shuffle=False)
            null = _fit_one(tr, te, seed=s, shuffle=True)
            per_seed[s] = {
                "n_train": len(tr), "n_test": len(te), "held_out_peaks": int(held),
                "gc": round(float(g), 4),
                "encoder_multi": cand["spearman"],
                "encoder_shuffled": null["spearman"],
                "margin": round(cand["spearman"] - float(g), 4),
                "epochs_ran": cand["epochs_ran"],
                "epochs_ran_null": null["epochs_ran"],
                "wall_seconds": cand["wall_seconds"] + null["wall_seconds"],
                "n_parameters": cand["n_parameters"],
                "conditions_on_the_outcome": True,
            }
            print(f"  seed {s:<2d} gc {g:+.4f}  conv {cand['spearman']:+.4f}  "
                  f"margin {per_seed[s]['margin']:+.4f}  null {null['spearman']:+.4f}  "
                  f"({cand['epochs_ran']}ep, {per_seed[s]['wall_seconds']:.0f}s)")

    margins = [per_seed[s]["margin"] for s in sorted(per_seed)]
    nulls = [abs(per_seed[s]["encoder_shuffled"]) for s in sorted(per_seed)]
    ci = _paired_bootstrap(margins) if margins else (0.0, 0.0)
    arms_summary = {
        "pooled_reconcile_ok": bool(reconcile_ok),
        "null_max_abs": round(max(nulls), 4) if nulls else 1.0,
        "active_margin_per_seed": margins,
        "active_margin_ci": ci,
    }
    verdict, detail = compute_verdict(arms_summary)
    if smoke:
        verdict, detail = ("SMOKE_NOT_GRADED",
                           f"ran {len(seeds)} of "
                           f"{len(PREREGISTERED_ACTIVE_FRAME['registered_seeds'])} registered "
                           f"seeds"
                           + (f" on a {args.subsample}-tile subsample" if args.subsample else "")
                           + f"; the bar would have returned {verdict}")

    # ---------------- PHASE 3: the shape of the margin ----------------
    binned = None
    if reconcile_ok and per_seed:
        s0 = sorted(per_seed)[0]
        tr, te, _h = leave_peak_out(active, frac=0.25, seed=s0)
        g, gpred = ridge_gc(tr, te)
        conv = fit_encoder(tr, te, widths=DEFAULT_WIDTHS, protocol=REGISTERED_PROTOCOL,
                          seed=s0, device="cpu")
        y = [float(t.expression) for t in te]
        binned = dict(
            _binned_margin([float(t.rep1) for t in te], conv["predictions"], gpred, y,
                           args.n_bins),
            seed=s0,
            binning_variable="RNA_exp_sum_1, group 1's published summary",
            conditions_on_the_outcome=False,
            is_not_a_second_bar=True,
            reads_as="the SHAPE of the margin across the range. The verdict is computed from "
                     "the seed medians alone and a test pins that this table cannot move it",
        )

    out = {
        "schema": "glm-active-frame-gate-v1",
        "date": str(date.today()),
        "preregistered": PREREGISTERED_ACTIVE_FRAME,
        "pooled_preregistration_imported_read_only": {
            "reconcile_anchor": POOLED_PREREGISTERED["reconcile_anchor"],
            "note": "imported from scripts/glm_encoder_gate.py for reference only; that file is "
                    "not edited and a test pins it byte-identical to HEAD",
        },
        "substrate": {"cache": args.cache, "n_tiles": len(tiles),
                      "n_active": audit["n_active"], "threshold_audit": audit},
        "threads_pinned": threads,
        "phase1_pooled": {"gc_seed0": round(float(gc0), 4), "reconcile_ok": bool(reconcile_ok),
                          "encoder_multi_seed0": pooled_conv},
        "phase2_active": {"n_active_used": len(active), "per_seed": per_seed,
                          "gc_median": round(st.median(per_seed[s]["gc"] for s in per_seed), 4)
                          if per_seed else None,
                          "encoder_multi_median": round(
                              st.median(per_seed[s]["encoder_multi"] for s in per_seed), 4)
                          if per_seed else None,
                          "margin_median": round(st.median(margins), 4) if margins else None,
                          "n_seeds_clearing_explains_bar": sum(
                              1 for m in margins
                              if m >= PREREGISTERED_ACTIVE_FRAME["explains_margin"]),
                          "margin_ci": ci,
                          "null_max_abs": arms_summary["null_max_abs"],
                          "conditions_on_the_outcome": True,
                          "not_a_ledger_candidate": True},
        "phase3_binned_margin": binned,
        "verdict": verdict,
        "verdict_detail": detail,
        "framing_recommendation": (framing_recommendation(verdict)
                                   if verdict in VERDICTS else
                                   "none -- a smoke run does not inform a scope decision"),
        "smoke": bool(smoke),
        "wall_seconds": round(time.time() - t_start, 1),
    }
    dest = Path(args.out) if args.out else Path(f"wiki/glm_active_frame_gate_{date.today()}.json")
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nVERDICT {verdict}\n  {detail}")
    print(f"framing recommendation: {out['framing_recommendation']}")
    print(f"wrote {dest} in {out['wall_seconds']}s")
    return EXIT_OK if verdict in ("FRAME_EXPLAINS_THE_GAP", "FRAME_CONTRIBUTES",
                                  "FRAME_DOES_NOT_EXPLAIN", "SMOKE_NOT_GRADED") \
        else EXIT_INDETERMINATE


if __name__ == "__main__":
    raise SystemExit(main())
