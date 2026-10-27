#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Does a 2.22-decade pyramid beat the length-invariant k-mer default on VARIABLE-LENGTH fragments?

    uv run python scripts/glm_fragment_pyramid_arm.py [--seeds 10] [--limit N] [--out PATH]

THE ARM: `pyramid-4scale-2.22decades-499bp-fragments`, declared on the ceiling criterion 2026-10-27
as the reachable re-scope of `multiresolution-pyramid-4-decades-UNRUN`, which its own executed probe
(`wiki/pyramid_span_feasibility_2026-10-27.json`) measured as substrate-gated at 60.1x.

THE BAR BELOW IS FROZEN BEFORE THE ARM RUNS. `predicted_outcome` and `negative_scope` are written
pre-run on purpose: on this rung a negative is the honest prediction, and a bar plus a prediction
recorded after a number exists is not a pre-registration.

================================================================================================
THREE DESIGN FORKS THE SUBSTRATE FORCES, SETTLED HERE RATHER THAN SILENTLY
================================================================================================

(1) THE RAGGED-INPUT FORK -- and `encoder.py` already took a position against the easy answer.
    `require_uniform_length` REFUSES a ragged batch, with this reasoning:

        A ragged batch has exactly two silent resolutions -- pad or truncate -- and both invent
        sequence that was never measured, then return a plausible number.

    Fragments are 48-499 bp over 359 distinct lengths, so that refusal is pointed straight at this
    arm. TRUNCATING to 150 bp would discard measured sequence AND cap the span at 1.70 decades,
    which deletes the question. NAIVE PADDING would do what the module forbids.

    SETTLED: pad to 499 with ALL-ZERO columns **and exclude padded positions from the pooling**
    (`masked_mean_pool`). The distinction is not cosmetic -- under masked pooling no padded position
    can influence any output, so nothing is invented: the pad is scaffolding for a rectangular
    tensor, not data. Two things make this checkable rather than asserted:
      * all-zero already MEANS "no evidence of any base here" in this project -- that is `one_hot`'s
        own documented encoding of an N, so the pad carries no new semantics;
      * `test_padding_cannot_change_a_score` asserts a sequence scores IDENTICALLY whatever it is
        padded to. If that test ever fails, the pad is contributing and this arm is invalid.

(2) THE LEAKAGE-UNIT FORK. `fit_encoder` deliberately does not hardcode the blocking variable:
    *"a new substrate must NAME its leakage unit instead of silently inheriting one that does not
    apply"*, and `tile_peak_key` REFUSES a record with no `peak` -- which every Fragment is.
    321,123 fragments tile a 4.64 Mb genome, so they OVERLAP heavily; a random split would train and
    test on near-duplicate sequence. CLAUDE.md prices this: random splits inflate ~3x.

    SETTLED: the leakage unit is a CONTIGUOUS GENOMIC BLOCK, `start // BLOCK_BP` with
    BLOCK_BP = 50,000 -- 100x the longest fragment, so block-boundary overlap is at most one
    fragment length per boundary and is bounded, reported, and not hidden.

(3) THE COMPARATOR FORK, and it is the one most likely to be got wrong by a reader.
    This arm is on the FRAGMENT rung, whose own ceiling is 0.7928. The measured pyramid sibling
    (0.3167) and promcalc (0.3968) are on the TILE rung, ceiling 0.941. The census `ceiling_rule`
    says each rung carries its OWN ceiling and none is inherited, so those numbers are NOT
    comparators here and a delta against them would be meaningless.

    SETTLED: the comparator is `kmer4` -- 4-mer counts + ridge, the length-invariant default
    `load_fragments`' own docstring names -- recomputed IN THIS RUN on IDENTICAL splits and seeds.
    The census records kmer4 at 0.1313 on this rung at 30,000 rows, but that figure is NOT used as
    the comparator: a paired within-run delta is the only honest comparison, since a cited number
    from a different split cannot be paired.

Read-only with respect to every ledger and every frozen file. Writes one artifact.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.glm.encoder import PYRAMID_WIDTHS, EncoderDataError  # noqa: E402
from dna_decode.glm.masked_pyramid import pool_coverage  # noqa: E402
from dna_decode.glm.genomewide import DEFAULT_CACHE, load_fragments, spearman  # noqa: E402

ARM = "pyramid-4scale-2.22decades-499bp-fragments"

#: The leakage unit. 100x the longest fragment (499 bp), so boundary overlap is bounded.
BLOCK_BP = 50_000

#: Pad target. The MEASURED max over both fragment rungs, not a round number.
PAD_TO = 499

PREREGISTERED = {
    "arm": ARM,
    "rung": "fragments_LB",
    "rung_ceiling": 0.7928,
    "comparator": "kmer4",
    "comparator_note": (
        "4-mer counts + ridge, recomputed IN-RUN on identical splits/seeds. The census's 0.1313 is "
        "context, NOT the comparator: a cited number from a different split cannot be paired."
    ),
    "primary_metric": "spearman on the held-out genomic blocks",
    "margin": 0.02,
    "margin_note": (
        "INHERITED from the tile-rung arms so the bars stay commensurable, and the inheritance is "
        "disclosed because it is NOT scale-free: 0.02 against this rung's incumbent (~0.13) is a "
        "~15% relative gain, whereas against the tile rung's ~0.40 it is ~5%. The verdict keys on "
        "the ABSOLUTE 0.02; the relative figure is reported beside it so no reader is misled."
    ),
    "min_seeds_winning": 8,
    "n_seeds": 10,
    "widths": list(PYRAMID_WIDTHS),
    # CORRECTED 2026-10-27 from a single MAX to the measured DISTRIBUTION. The max-only figure
    # (2.2210) is realised on 49 of 321,123 fragments; the typical row carries 1.9103. Reporting the
    # max alone overstated the re-scope by 2.5x, which is the forgotten-substrate error in another
    # costume: a claim about a distribution taken from its extreme.
    "span_decades_at_max": round(math.log10(PAD_TO / PYRAMID_WIDTHS[0]), 4),
    "span_decades_at_median": 1.9103,
    "span_decades_note": (
        "MEASURED over all 321,123 LB fragments: median 244 bp (1.9103 decades), p99 308 bp "
        "(2.011), max 499 bp (2.2210). Only 4.0% of rows reach 300 bp and 0.015% reach 400 bp. "
        "Against the measured tile sibling's 1.70 decades the honest gain is +0.21 at the MEDIAN "
        "(244/150 = 1.63x more length), NOT the +0.52 / 3.33x the maximum suggests. Quote the "
        "median figure; the max describes 49 fragments."
    ),
    "leakage_unit": "contiguous genomic block, start // %d" % BLOCK_BP,
    "length_rule": "zero-pad to %d with MASKED pooling; padded positions excluded" % PAD_TO,
    "verdicts": {
        "PYRAMID_ADDS": "gain >= margin AND >= min_seeds_winning seeds win -- a reproducible lift",
        "WEAK_PYRAMID": "gain > 0 AND >= min_seeds_winning seeds win, but below the margin -- a "
                        "REPRODUCIBLE gain that fell short of the effect-size bar",
        "INDETERMINATE_GAIN_NOT_REPRODUCIBLE":
            "gain > 0 but FEWER than min_seeds_winning seeds win -- a positive gain that did NOT "
            "reproduce, which is neither a measured lift nor a measured negative and must not be "
            "read as either",
        "NO_GAIN": "gain <= 0 -- the pyramid does not beat the length-invariant default",
        "INVALID_NULL_SCORED": "the shuffled-label control scored; nothing may be graded",
        "INVALID_PAD_LEAKED": "padding changed a score; the masked-pooling premise is broken",
    },
    "validity_controls": {
        "shuffled_labels": "train labels permuted; must score ~0 or the run is INVALID",
        "pad_invariance": "a sequence must score identically whatever length it is padded to",
        "length_matched_subset": (
            "AMENDED 2026-10-27 from >= 400 bp to >= 300 bp, BEFORE any encoder number existed. "
            "The 400 bp threshold was frozen without measuring the distribution and selects 49 of "
            "321,123 rows (0.015%) -- a control too small to falsify anything, which is worse than "
            "no control because it LOOKS like one. 300 bp selects 12,896 rows (4.0%) and is "
            "powered. Still the minority of the substrate, so it is a padding check, NOT a second "
            "headline. If the headline and this subset disagree in SIGN, padding is implicated and "
            "the headline is WITHDRAWN rather than explained."
        ),
        "length_matched_threshold_bp": 300,
    },
    "predicted_outcome": "NO_GAIN or WEAK_PYRAMID",
    "prediction_rationale": (
        "Recorded BEFORE the run. Four reasons a negative is expected: (a) the conv encoder already "
        "failed to beat a SINGLE GC NUMBER on the easier tile rung (+0.0122, paired CI includes "
        "zero); (b) `load_fragments`' own docstring argues k-mer features are the honest default "
        "here precisely BECAUSE they are length-invariant, which is the property a fixed-window "
        "conv lacks; (c) this rung's best measured arm reaches only 0.1313 of a 0.7928 ceiling "
        "(16.6%), so the recoverable signal looks thin for every model class, not just this one; "
        "(d) more SPAN is only useful if regulatory signal actually lives above 150 bp, which is "
        "the hypothesis under test and not an assumption in its favour."
    ),
    "negative_scope": (
        "A negative here means: a 4-scale conv pyramid over zero-padded, masked-pooled fragments, "
        "trained under the registered protocol on block-split LB fragments, did not beat 4-mer "
        "counts + ridge by 0.02 Spearman. It would NOT mean sequence context above 150 bp is "
        "uninformative, NOT that multi-resolution modelling fails in general, and NOT anything "
        "about the tile rung, which has a different ceiling and a different incumbent."
    ),
    "revival_condition": (
        "a length-invariant architecture (dilated conv, attention with length masking) or a rung "
        "whose fragments are long enough that padding is negligible across the whole set."
    ),
}


class ArmDataError(ValueError):
    """Refuse rather than return a plausible number."""


# ------------------------------------------------------------------ pure helpers

def block_key(start: int, block_bp: int = BLOCK_BP) -> str:
    """The declared leakage unit. A pure function of the coordinate, so blocks are disjoint."""
    if start < 0:
        raise ArmDataError("negative start coordinate %d" % start)
    return "blk%06d" % (start // block_bp)


def pad_sequences(seqs: list[str], pad_to: int = PAD_TO) -> tuple[list[str], list[int]]:
    """Right-pad with 'N' to a common length; return the padded seqs and the TRUE lengths.

    'N' is used rather than a novel symbol because `one_hot` already maps an unknown base to an
    ALL-ZERO column -- the project's own documented encoding of "no evidence of any base here".
    The true lengths are returned so pooling can exclude the pad; they are not optional.
    """
    out, lengths = [], []
    for s in seqs:
        if len(s) > pad_to:
            raise ArmDataError(
                "sequence of %d bp exceeds pad target %d; the pad target must be the MEASURED max, "
                "and truncating here would discard measured sequence" % (len(s), pad_to))
        lengths.append(len(s))
        out.append(s + "N" * (pad_to - len(s)))
    return out, lengths


def kmer4_features(seqs: list[str], lengths: list[int] | None = None):
    """256-dim 4-mer counts, L1-normalised. Length-invariant BY CONSTRUCTION.

    If `lengths` is given, only the TRUE sequence is counted -- a padded tail must not dilute the
    composition, which would make the comparator silently length-dependent and hand this arm an
    unearned advantage over its own baseline.
    """
    import numpy as np
    bases = "ACGT"
    index = {}
    for i, a in enumerate(bases):
        for j, b in enumerate(bases):
            for k, c in enumerate(bases):
                for m, d in enumerate(bases):
                    index[a + b + c + d] = ((i * 4 + j) * 4 + k) * 4 + m
    X = np.zeros((len(seqs), 256), dtype=np.float32)
    for n, s in enumerate(seqs):
        true = s[: lengths[n]] if lengths is not None else s
        for p in range(len(true) - 3):
            idx = index.get(true[p:p + 4])
            if idx is not None:
                X[n, idx] += 1.0
        tot = X[n].sum()
        if tot > 0:
            X[n] /= tot
    return X


def split_by_block(blocks: list[str], *, test_fraction: float, seed: int):
    """Hold out WHOLE genomic blocks. Disjointness is a pure function of the block label."""
    import random as _r
    uniq = sorted(set(blocks))
    if len(uniq) < 2:
        raise ArmDataError("need >= 2 genomic blocks to split, got %d" % len(uniq))
    ps = list(uniq)
    _r.Random(7_000 + seed).shuffle(ps)
    k = max(1, min(len(ps) - 1, int(round(len(ps) * test_fraction))))
    held = set(ps[:k])
    tr = [i for i, b in enumerate(blocks) if b not in held]
    te = [i for i, b in enumerate(blocks) if b in held]
    if not tr or not te:
        raise ArmDataError("degenerate block split")
    return tr, te


def paired_stats(deltas: list[float], margin: float, min_seeds: int) -> dict:
    """Median, win count and a bootstrap-free CI95 over seeds. Pure."""
    if not deltas:
        raise ArmDataError("no paired deltas")
    n = len(deltas)
    med = statistics.median(deltas)
    wins = sum(1 for d in deltas if d > 0)
    if n > 1:
        sd = statistics.stdev(deltas)
        half = 1.96 * sd / math.sqrt(n)
        mean = statistics.fmean(deltas)
        ci = [round(mean - half, 4), round(mean + half, 4)]
    else:
        ci = [None, None]
    return {"median": round(med, 4), "wins": wins, "n": n, "ci95": ci,
            "clears_margin": med >= margin, "reproducible": wins >= min_seeds}


def verdict(p: dict, pre: dict = PREREGISTERED) -> tuple[str, str]:
    """The frozen vocabulary. The seed condition is read from the bar, never written as a literal.

    Ordering note: the two halves of the bar are SEPARATE branches, so a reproducible-but-small gain
    and unreproduced positive noise can never return one string -- the defect closed repo-wide on
    2026-10-26.
    """
    m, w, n = p["median"], p["wins"], p["n"]
    need = pre["min_seeds_winning"]
    if m >= pre["margin"] and w >= need:
        return ("PYRAMID_ADDS",
                "pyramid %+.4f over %s, %d/%d seeds, CI95 %s -- clears the %+.2f margin and the "
                "%d-seed bar" % (m, pre["comparator"], w, n, p["ci95"], pre["margin"], need))
    if m > 0 and w >= need:
        return ("WEAK_PYRAMID",
                "pyramid %+.4f over %s, %d/%d seeds, CI95 %s -- REPRODUCIBLE (%d/%d clears the "
                "%d-seed bar) but short of the %+.2f margin"
                % (m, pre["comparator"], w, n, p["ci95"], w, n, need, pre["margin"]))
    if m > 0:
        return ("INDETERMINATE_GAIN_NOT_REPRODUCIBLE",
                "pyramid %+.4f over %s, only %d/%d seeds, CI95 %s -- the gain is positive but did "
                "NOT reproduce: %d of %d seeds win, against the %d the bar requires; this is "
                "neither a measured lift nor a measured negative and must not be read as either"
                % (m, pre["comparator"], w, n, p["ci95"], w, n, need))
    return ("NO_GAIN",
            "pyramid %+.4f over %s, %d/%d seeds, CI95 %s -- does not beat the length-invariant "
            "default" % (m, pre["comparator"], w, n, p["ci95"]))


# ------------------------------------------------------------------ the run

def _ridge_fit_predict(Xtr, ytr, Xte, alpha: float = 1.0):
    """Closed-form ridge. No sklearn dependency, so the baseline cannot drift under a version bump."""
    import numpy as np
    Xtr = np.asarray(Xtr, dtype=np.float64)
    Xte = np.asarray(Xte, dtype=np.float64)
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1.0
    A = (Xtr - mu) / sd
    B = (Xte - mu) / sd
    ym = float(np.mean(ytr))
    y = np.asarray(ytr, dtype=np.float64) - ym
    G = A.T @ A + alpha * np.eye(A.shape[1])
    w = np.linalg.solve(G, A.T @ y)
    return (B @ w) + ym


def run(cache_dir=DEFAULT_CACHE, *, which="LB", seeds=10, limit=None,
        allow_download=False, test_fraction=0.2) -> dict:
    frags, st = load_fragments(cache_dir, which=which, limit=limit,
                               allow_download=allow_download)
    seqs = [f.seq for f in frags]
    ys = [f.expression for f in frags]
    blocks = [block_key(f.start) for f in frags]

    padded, lengths = pad_sequences(seqs)
    Xk = kmer4_features(padded, lengths)

    n_blocks = len(set(blocks))
    long_mask = [i for i, L in enumerate(lengths)
                 if L >= PREREGISTERED["validity_controls"]["length_matched_threshold_bp"]]
    srt = sorted(lengths)
    span_dist = {
        "length_median": srt[len(srt) // 2],
        "length_p99": srt[min(len(srt) - 1, int(0.99 * len(srt)))],
        "length_max": srt[-1],
        "span_at_median": round(math.log10(srt[len(srt) // 2] / PYRAMID_WIDTHS[0]), 4),
        "span_at_max": round(math.log10(srt[-1] / PYRAMID_WIDTHS[0]), 4),
        "frac_ge_300bp": round(sum(1 for L in srt if L >= 300) / len(srt), 4),
    }

    rows = []
    for seed in range(seeds):
        tr, te = split_by_block(blocks, test_fraction=test_fraction, seed=seed)
        base = spearman(list(_ridge_fit_predict(Xk[tr], [ys[i] for i in tr], Xk[te])),
                        [ys[i] for i in te])
        rows.append({"seed": seed, "n_train": len(tr), "n_test": len(te),
                     "kmer4_spearman": round(base, 4)})

    return {
        "schema": "glm-fragment-pyramid-arm-v1",
        "arm": ARM,
        "date": str(date.today()),
        "preregistered": PREREGISTERED,
        "substrate": st.as_dict(),
        "n_sequences": len(frags),
        "n_genomic_blocks": n_blocks,
        "n_fragments_in_length_matched_subset": len(long_mask),
        "span_distribution": span_dist,
        "branch_coverage": pool_coverage(lengths, PYRAMID_WIDTHS),
        "baseline_rows": rows,
        "status": "BASELINE_ONLY_ENCODER_NOT_RUN",
        "status_note": (
            "The k-mer comparator and every pre-registered design decision are in place and the "
            "block split is exercised. The pyramid half needs the masked-pooling encoder, which is "
            "NOT the shipped `fit_encoder` (that one requires uniform length and pools unmasked). "
            "No verdict is emitted: a verdict from a one-sided run would be a figure, not a "
            "measurement."
        ),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--which", default="LB", choices=["LB", "M9"])
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--allow-download", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    try:
        res = run(a.cache_dir, which=a.which, seeds=a.seeds, limit=a.limit,
                  allow_download=a.allow_download)
    except (ArmDataError, EncoderDataError) as exc:
        print("REFUSED: %s" % exc, file=sys.stderr)
        return 3

    base = [r["kmer4_spearman"] for r in res["baseline_rows"]]
    print("arm: %s" % ARM)
    print("  substrate: n=%d, %d-%d bp, %d genomic blocks, %d in the length-matched subset"
          % (res["n_sequences"], res["substrate"]["length_min"],
             res["substrate"]["length_max"], res["n_genomic_blocks"],
             res["n_fragments_in_length_matched_subset"]))
    sd = res["span_distribution"]
    print("  span: %.4f decades at the MEDIAN (%d bp), %.4f at the max (%d bp) over widths %s"
          % (sd["span_at_median"], sd["length_median"], sd["span_at_max"],
             sd["length_max"], PREREGISTERED["widths"]))
    print("  quote the median: only %.1f%% of rows reach 300 bp"
          % (100 * sd["frac_ge_300bp"]))
    print("  comparator kmer4 spearman over %d seeds: median %.4f  min %.4f  max %.4f"
          % (len(base), statistics.median(base), min(base), max(base)))
    print("  STATUS: %s" % res["status"])

    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
        print("  wrote %s" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
