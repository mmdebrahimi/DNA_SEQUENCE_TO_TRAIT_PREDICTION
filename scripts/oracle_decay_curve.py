"""Does the oracle's accuracy DECAY with distance from its training data? Measure the curve.

    uv run python scripts/oracle_decay_curve.py --self-check
    uv run python scripts/oracle_decay_curve.py --ckpt data/promoteratlas/finetune_checkpoints/seed_0.pt \
        --pa-repo <PromoterAtlas clone> --out wiki/oracle_decay_curve_<date>.json

**WHY THIS EXISTS, and the motivation is WEAKER than the plan originally claimed.** The plan asserted
that published work ESTABLISHES this decay. Reading the primary sources (2026-11-18,
`wiki/generator_hypotheses_source_notes_2026-11-18.md`) refuted that:

  * The edit-distance result cited for it measures training-data INFLUENCE (via Kronfluence), not
    accuracy, and that paper never plots performance against distance at all.
  * The other support -- extrapolation collapsing by the K=5 mutational regime -- is on SIMULATED NK
    landscapes, scored with PEARSON r, counting regimes from a seed, and explicitly NOT distance to
    the nearest training sequence.

So the decay is PLAUSIBLE and NOT established for the quantity measured here. That makes this
measurement MORE novel, not less: two papers cited as establishing it instead confirm it is un-done.
**Do not report this as confirming a known effect.**

**WHAT IS MEASURED.** For each held-out variant, the MINIMUM Hamming distance to any sequence in the
TRAIN split -- i.e. how far the nearest thing the model saw actually is. Variants are binned by that
distance and Spearman against MEASURED activity is reported PER BIN, with the per-bin n beside it.

**THE FLOOR IS FIXED BEFORE ANY RUN AND REFUSES RATHER THAN REPORTS.** `MIN_BIN_N` bins carrying
fewer than that many variants emit `status: refused:below_floor` with NO correlation, because a
Spearman over a handful of points is a figure, not a measurement, and a decay "trend" assembled from
such bins is exactly the artefact this script exists to avoid producing.

**ONE DELIBERATE DEVIATION FROM THE PLAN'S WORDING.** The plan asked for "a Wilson-style interval".
A Wilson interval is for a PROPORTION; this statistic is a correlation, so the correct analogue is
the **Fisher z-transform** interval, which is what `fisher_z_ci` implements. Using a Wilson interval
on a correlation would have produced a plausible-looking number that means nothing. Recorded here
rather than silently substituted.

**SCOPE, stamped into every artifact.** This measures decay on an MPRA PROMOTER substrate
(GSE108535), not a protein DMS substrate. The encoder was fine-tuned on a DIFFERENT dataset
(GSE144621 peak tiles), so "distance from training data" here means distance to the held-in
GSE108535 split, NOT to the fine-tuning corpus. Those are different quantities and the artifact says
so.

Reuses, never reimplements: `glm/expression.py` owns the GSE108535 substrate (its docstring records
that a naive join silently loses 65% of the data) and `scripts/glm_selfdistill_finetuned_teacher.py`
owns the fine-tuned-encoder scorer.
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import date
from pathlib import Path

import numpy as np

from dna_decode.glm.expression import (
    SEQ_LENGTH,
    leave_element_out,
    load_pairs,
    spearman,
)

#: Distance bins, FIXED before any run. Roughly geometric, because a linear grid over 0..150 would
#: put almost every combinatorial variant in one bucket. Upper bound is the sequence length, so the
#: bins are exhaustive by construction and a distance cannot fall outside them.
BINS: tuple[tuple[int, int], ...] = (
    (0, 1), (2, 3), (4, 7), (8, 15), (16, 31), (32, 63), (64, SEQ_LENGTH),
)

#: A bin with fewer than this many variants REFUSES rather than reporting a correlation. Fixed here,
#: before any data is read, so it cannot be tuned to make a curve look monotone.
MIN_BIN_N = 30

#: Fisher z needs n > 3 for a finite interval.
MIN_N_FOR_CI = 4

NUC = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode(seqs: list[str]) -> np.ndarray:
    """Sequences -> (N, L) uint8. Any non-ACGT character becomes a 5th symbol, so it always MISMATCHES."""
    out = np.full((len(seqs), SEQ_LENGTH), 4, dtype=np.uint8)
    for i, s in enumerate(seqs):
        for j, ch in enumerate(s[:SEQ_LENGTH]):
            out[i, j] = NUC.get(ch, 4)
    return out


def hamming_min_to_set(test_enc: np.ndarray, train_enc: np.ndarray, *, chunk: int = 128) -> np.ndarray:
    """For each test row, the MINIMUM Hamming distance to any train row.

    Chunked on purpose: the full pairwise comparison is n_test x n_train x L booleans, which at the
    real cohort size (2,725 x 8,173 x 150) is billions of elements. Chunking bounds peak memory
    without changing the result, and the result is exact -- this is not an approximate nearest
    neighbour.
    """
    if train_enc.size == 0:
        raise ValueError("empty train set: a minimum distance to nothing is undefined, not infinity")
    mins = np.empty(len(test_enc), dtype=np.int32)
    for i in range(0, len(test_enc), chunk):
        block = test_enc[i:i + chunk]
        # (block, 1, L) vs (1, train, L) -> (block, train, L)
        d = (block[:, None, :] != train_enc[None, :, :]).sum(axis=2)
        mins[i:i + chunk] = d.min(axis=1)
    return mins


def bin_index(distance: int) -> int:
    """Which BINS slot a distance falls in. Raises rather than silently dropping an out-of-range value."""
    for k, (lo, hi) in enumerate(BINS):
        if lo <= distance <= hi:
            return k
    raise ValueError("distance %d falls outside BINS, which are meant to be exhaustive" % distance)


def fisher_z_ci(rho: float, n: int, *, z: float = 1.96) -> tuple[float, float] | None:
    """Fisher z-transform interval for a correlation. None when n is too small for one to exist.

    This is the correct analogue of the plan's "Wilson-style interval" -- Wilson is for a proportion.
    """
    if n < MIN_N_FOR_CI or not (-1.0 < rho < 1.0):
        return None
    se = 1.0 / math.sqrt(n - 3)
    zr = 0.5 * math.log((1.0 + rho) / (1.0 - rho))
    lo, hi = zr - z * se, zr + z * se
    return (math.tanh(lo), math.tanh(hi))


def summarise_bins(distances, predicted, measured, *, spearman_fn=spearman) -> list[dict]:
    """Per-bin Spearman with n and CI, REFUSING below the fixed floor.

    `spearman_fn` is injected so the pure logic is testable without the real metric.
    """
    distances = np.asarray(distances)
    rows: list[dict] = []
    for k, (lo, hi) in enumerate(BINS):
        sel = (distances >= lo) & (distances <= hi)
        n = int(sel.sum())
        row: dict = {"bin": "%d-%d" % (lo, hi), "bin_index": k, "n": n}
        if n < MIN_BIN_N:
            row["status"] = "refused:below_floor"
            row["spearman"] = None
            row["ci95"] = None
            row["why"] = ("n=%d is below the pre-registered floor of %d; a correlation here would be "
                          "a figure, not a measurement" % (n, MIN_BIN_N))
        else:
            p = [float(v) for v, s in zip(predicted, sel) if s]
            m = [float(v) for v, s in zip(measured, sel) if s]
            r = float(spearman_fn(p, m))
            row["status"] = "scored"
            row["spearman"] = round(r, 4)
            ci = fisher_z_ci(r, n)
            row["ci95"] = None if ci is None else [round(ci[0], 4), round(ci[1], 4)]
        rows.append(row)
    return rows


def monotonicity_report(rows: list[dict]) -> dict:
    """REPORT whether scored bins decline with distance -- never assume or assert it."""
    scored = [(r["bin_index"], r["spearman"]) for r in rows if r["status"] == "scored"]
    if len(scored) < 2:
        return {"verdict": "indeterminate:fewer_than_two_scored_bins", "n_scored_bins": len(scored)}
    vals = [v for _, v in scored]
    deltas = [round(b - a, 4) for a, b in zip(vals, vals[1:])]
    strictly_down = all(d < 0 for d in deltas)
    return {
        "verdict": "declines_monotonically" if strictly_down else "not_monotonic",
        "n_scored_bins": len(scored),
        "scored_bin_indices": [i for i, _ in scored],
        "spearman_by_scored_bin": vals,
        "successive_deltas": deltas,
        "caveat": ("REPORTED, not assumed. With few scored bins a monotone run is weak evidence, and "
                   "the bins are not independent samples of one underlying curve."),
    }


def _self_check() -> int:
    """Exercise the pure helpers on synthetic data. No substrate, no checkpoint, no torch."""
    ok = True

    a = encode(["A" * SEQ_LENGTH])
    b = encode(["A" * SEQ_LENGTH, "T" + "A" * (SEQ_LENGTH - 1)])
    d = hamming_min_to_set(a, b)
    print("  identical row present -> min distance %d (expect 0)" % d[0]); ok &= d[0] == 0

    c = encode(["T" * 3 + "A" * (SEQ_LENGTH - 3)])
    d2 = hamming_min_to_set(c, encode(["A" * SEQ_LENGTH]))
    print("  three substitutions   -> min distance %d (expect 3)" % d2[0]); ok &= d2[0] == 3

    print("  bin(0)=%d bin(1)=%d bin(2)=%d bin(7)=%d bin(8)=%d bin(150)=%d"
          % tuple(bin_index(x) for x in (0, 1, 2, 7, 8, SEQ_LENGTH)))
    ok &= (bin_index(0), bin_index(1)) == (0, 0)
    ok &= (bin_index(2), bin_index(3)) == (1, 1)
    ok &= bin_index(7) == 2 and bin_index(8) == 3

    try:
        bin_index(SEQ_LENGTH + 1); print("  out-of-range distance did NOT raise -- BAD"); ok = False
    except ValueError:
        print("  out-of-range distance raises, as intended")

    try:
        hamming_min_to_set(a, np.empty((0, SEQ_LENGTH), dtype=np.uint8))
        print("  empty train set did NOT raise -- BAD"); ok = False
    except ValueError:
        print("  empty train set raises rather than returning infinity")

    # below-floor bin must REFUSE, and a populated bin must score.
    few = np.zeros(MIN_BIN_N - 1, dtype=int)
    rows = summarise_bins(few, list(range(len(few))), list(range(len(few))))
    r0 = rows[0]
    print("  bin of n=%d -> status %s, spearman %s" % (r0["n"], r0["status"], r0["spearman"]))
    ok &= r0["status"] == "refused:below_floor" and r0["spearman"] is None

    many = np.zeros(MIN_BIN_N + 10, dtype=int)
    vals = list(range(len(many)))
    rows2 = summarise_bins(many, vals, vals)
    print("  bin of n=%d -> status %s, spearman %s" % (rows2[0]["n"], rows2[0]["status"], rows2[0]["spearman"]))
    ok &= rows2[0]["status"] == "scored" and rows2[0]["spearman"] == 1.0

    ci = fisher_z_ci(0.5, 100)
    print("  fisher CI(rho=0.5, n=100) = %s" % (None if ci is None else [round(v, 4) for v in ci]))
    ok &= ci is not None and ci[0] < 0.5 < ci[1]
    print("  fisher CI(rho=0.5, n=3)   = %s (expect None)" % (fisher_z_ci(0.5, 3),))
    ok &= fisher_z_ci(0.5, 3) is None

    mono = monotonicity_report(rows2)
    print("  monotonicity on one scored bin -> %s" % mono["verdict"])
    ok &= mono["verdict"] == "indeterminate:fewer_than_two_scored_bins"

    print("\nSELF-CHECK: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-check", action="store_true", help="exercise the pure helpers and exit")
    ap.add_argument("--cache-dir", default="data/mpra_gse108535")
    ap.add_argument("--ckpt", default="data/promoteratlas/finetune_checkpoints/seed_0.pt")
    ap.add_argument("--pa-repo", default=None)
    ap.add_argument("--axis", default="Minus10")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.self_check:
        return _self_check()

    pairs, _rep = load_pairs(args.cache_dir, allow_download=False)
    train, test, held = leave_element_out(pairs, args.axis, seed=args.seed)
    print("%s: held %s -> train %d / test %d" % (args.axis, held, len(train), len(test)))

    tr_seqs = [p.seq for p in train]
    te_seqs = [p.seq for p in test]
    measured = [p.expression for p in test]

    print("computing exact min Hamming distance, %d test x %d train ..." % (len(te_seqs), len(tr_seqs)))
    dists = hamming_min_to_set(encode(te_seqs), encode(tr_seqs))
    print("  distance range %d..%d, median %d" % (dists.min(), dists.max(), int(np.median(dists))))

    from scripts.glm_selfdistill_finetuned_teacher import load_finetuned_teacher
    score, ck_meta = load_finetuned_teacher(
        Path(args.ckpt), Path(args.pa_repo) if args.pa_repo else None)
    preds = [float(v) for v in score(te_seqs)]
    overall = float(spearman(preds, measured))
    print("  overall Spearman (all bins pooled) = %+.4f" % overall)

    rows = summarise_bins(dists, preds, measured)
    for r in rows:
        if r["status"] == "scored":
            print("  bin %-8s n=%-5d rho=%+.4f  ci95=%s" % (r["bin"], r["n"], r["spearman"], r["ci95"]))
        else:
            print("  bin %-8s n=%-5d REFUSED (below floor %d)" % (r["bin"], r["n"], MIN_BIN_N))

    artifact = {
        "schema": "oracle-decay-curve-v1",
        "date": str(date.today()),
        "question": ("does the fine-tuned encoder's rank accuracy decay with MINIMUM Hamming distance "
                     "to the nearest training sequence?"),
        "motivation_is_weaker_than_the_plan_claimed": (
            "The plan asserted published work ESTABLISHES this decay. Primary-source reading refuted "
            "that: the edit-distance result cited measures training-data INFLUENCE not accuracy and "
            "never plots performance against distance, and the K-regime collapse is on SIMULATED NK "
            "landscapes with PEARSON r counting regimes from a seed, explicitly not distance to the "
            "nearest training sequence. The decay is PLAUSIBLE, NOT established for this quantity. Do "
            "not report this as confirming a known effect. See "
            "wiki/generator_hypotheses_source_notes_2026-11-18.md"),
        "substrate": {
            "dataset": "GEO GSE108535 (Urtecho/Kosuri sigma70 variants)",
            "n_pairs": len(pairs), "axis": args.axis, "held_out_values": held,
            "n_train": len(train), "n_test": len(test), "seq_length": SEQ_LENGTH,
        },
        "scope_limits": [
            "MPRA PROMOTER substrate, NOT protein DMS. Does not transfer to the DMS regime.",
            "The encoder was fine-tuned on a DIFFERENT dataset (GSE144621 peak tiles), so distance "
            "here is to the held-in GSE108535 split, NOT to the fine-tuning corpus. Different "
            "quantities; do not quote one as the other.",
            "One checkpoint of ten seeds, one element axis, one substrate.",
            "Distance is EXACT minimum Hamming to the nearest train sequence, not an approximation.",
        ],
        "teacher": ck_meta,
        "bins_fixed_before_run": [list(b) for b in BINS],
        "min_bin_n_fixed_before_run": MIN_BIN_N,
        "interval_method": ("Fisher z-transform. The plan said 'Wilson-style'; Wilson is for a "
                            "PROPORTION and this is a correlation, so Wilson would have produced a "
                            "plausible-looking number that means nothing. Deviation recorded, not "
                            "silently substituted."),
        "distance_summary": {
            "min": int(dists.min()), "max": int(dists.max()),
            "median": int(np.median(dists)), "mean": round(float(dists.mean()), 2),
        },
        "overall_spearman_pooled": round(overall, 4),
        "bins": rows,
        "monotonicity": monotonicity_report(rows),
    }
    out = Path(args.out) if args.out else Path("wiki/oracle_decay_curve_%s.json" % date.today())
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print("\nwrote %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
