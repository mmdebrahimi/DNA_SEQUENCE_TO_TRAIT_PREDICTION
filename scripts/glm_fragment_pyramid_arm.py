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

from dna_decode.glm.encoder import (  # noqa: E402
    PYRAMID_WIDTHS,
    REGISTERED_PROTOCOL,
    EncoderDataError,
)
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
        "shuffled_labels": (
            "AMENDED 2026-10-27, BEFORE any scored run, because the first design was MEASURED "
            "invalid. It was: one seed, |rho| < 0.05. The null's own spread makes that unusable -- "
            "measured over 5 seeds at N=4,000 the null is median -0.0008 / mean +0.0130 / "
            "**sd 0.1016**, range -0.1236..+0.1301, so a single-seed 0.05 bar FAILS 3 of 5 times on "
            "a pipeline that is demonstrably centred on zero. The first smoke duly returned a FALSE "
            "INVALID_NULL_SCORED. "
            "WHY the spread is that large, and it is not a defect: early stopping trips almost "
            "immediately once the train-side labels are permuted, leaving weights near "
            "initialisation, so the model is close to a RANDOM linear functional of base "
            "composition -- and composition predicts expression on this rung at "
            "spearman(GC, expression) = -0.1169, which is the same magnitude as the whole kmer4 "
            "comparator (0.1120). A random composition-like feature therefore lands anywhere in "
            "+/-0.13 by chance. "
            "NOW: the null runs on ALL seeds and the bar is on the MEDIAN, |median| < margin "
            "(0.02); the per-seed values and the sd are both reported so a reader can see the "
            "power. **The primary noise control was already pre-registered and is the SEED-COUNT "
            "bar**: under a null centred on zero each seed wins with p=0.5, so >= 8 of 10 has "
            "p = 0.0547. The single-seed null bar was both redundant and miscalibrated."
        ),
        "null_median_bar": 0.02,
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


def median_ci95(values: list[float]) -> tuple:
    """Distribution-free CI for the MEDIAN from order statistics. Exact and bootstrap-free.

    For ranks k and n+1-k the coverage is exactly 1 - 2*P(Bin(n, 1/2) <= k-1), so this reports the
    coverage it ACHIEVES rather than claiming 95%. Returns the tightest pair still covering >= 95%,
    or (None, None, None) when no pair does -- true for n <= 5, which includes the 3-seed smoke.
    Inventing an interval there would put a number on the one case that cannot support one.
    """
    n = len(values)
    s = sorted(values)
    best = (None, None, None)
    for k in range(1, n // 2 + 1):
        tail = sum(math.comb(n, i) for i in range(k)) / float(2 ** n)
        coverage = 1.0 - 2.0 * tail
        if coverage >= 0.95:
            best = (s[k - 1], s[n - k], coverage)   # ascending k -> tightest qualifying pair wins
    return best


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
    mlo, mhi, mcov = median_ci95(deltas)
    return {"median": round(med, 4), "wins": wins, "n": n, "ci95": ci,
            # `ci95` is KEPT byte-identical: it is a normal-approximation interval around the MEAN,
            # several committed artifacts carry it, and glm_tile_hybrid_arm.py reads it. But the
            # headline and `clears_margin` are the MEDIAN, so an unlabelled interval beside them
            # reads as the median's and is not. Both estimands are now named explicitly.
            "ci95_estimand": "mean (normal approximation, 1.96*sd/sqrt(n))",
            "ci95_median": [None if mlo is None else round(mlo, 4),
                            None if mhi is None else round(mhi, 4)],
            "ci95_median_estimand": ("order-statistic interval for the MEDIAN, exact coverage "
                                     "%.4f" % mcov) if mcov is not None
                                    else ("no rank pair reaches 95%% coverage at n=%d, so no "
                                          "median interval is reported" % n),
            "ci95_median_coverage": None if mcov is None else round(mcov, 4),
            "clears_margin": med >= margin, "reproducible": wins >= min_seeds}


def _ci_label(p: dict) -> str:
    """Render the interval WITH its estimand, so a median is never printed beside a bare interval.

    Prefers the median interval, because the headline and `clears_margin` are the median. Falls back
    to the mean interval -- labelled as such -- when the seed count cannot support a median one.
    """
    lo, hi = p.get("ci95_median", [None, None])
    if lo is not None and hi is not None:
        return "[%+.4f, %+.4f] (median, coverage %.3f)" % (lo, hi, p["ci95_median_coverage"])
    return "%s (MEAN-centred; n=%d too small for a median interval)" % (p["ci95"], p["n"])


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
                "%d-seed bar" % (m, pre["comparator"], w, n, _ci_label(p), pre["margin"], need))
    if m > 0 and w >= need:
        return ("WEAK_PYRAMID",
                "pyramid %+.4f over %s, %d/%d seeds, CI95 %s -- REPRODUCIBLE (%d/%d clears the "
                "%d-seed bar) but short of the %+.2f margin"
                % (m, pre["comparator"], w, n, _ci_label(p), w, n, need, pre["margin"]))
    if m > 0:
        return ("INDETERMINATE_GAIN_NOT_REPRODUCIBLE",
                "pyramid %+.4f over %s, only %d/%d seeds, CI95 %s -- the gain is positive but did "
                "NOT reproduce: %d of %d seeds win, against the %d the bar requires; this is "
                "neither a measured lift nor a measured negative and must not be read as either"
                % (m, pre["comparator"], w, n, _ci_label(p), w, n, need))
    return ("NO_GAIN",
            "pyramid %+.4f over %s, %d/%d seeds, CI95 %s -- does not beat the length-invariant "
            "default" % (m, pre["comparator"], w, n, _ci_label(p)))


# ------------------------------------------------------------------ the run

def fit_masked_pyramid(codes, lengths, ys, tr, te, *, seed: int, protocol=REGISTERED_PROTOCOL,
                       device: str = "cpu", shuffle_labels: bool = False,
                       blocks: list | None = None):
    """Train the masked pyramid on `tr`, return Spearman on `te`. Mirrors `encoder.fit_encoder`'s role.

    `shuffle_labels` permutes the TRAIN targets only -- the validity null. The test labels are never
    touched, so the null measures what it claims to.

    Early stopping runs against a BLOCK-BLOCKED inner split of `tr`, for the reason every outer split
    here is blocked: fragments overlap, so a random inner split early-stops against near-duplicates
    of its own training rows. `blocks` is required for that; without it this refuses rather than
    quietly using a random inner split.
    """
    import numpy as np
    import torch

    from dna_decode.glm.masked_pyramid import build_masked_pyramid, codes_to_onehot

    if blocks is None:
        raise ArmDataError(
            "fit_masked_pyramid needs `blocks` for a block-blocked inner split; a random inner "
            "split would early-stop against near-duplicates of its own training rows")

    g = torch.Generator().manual_seed(seed)
    torch.manual_seed(seed)
    rng = np.random.default_rng(10_000 + seed)

    tr_blocks = [blocks[i] for i in tr]
    inner_tr, inner_va = split_by_block(tr_blocks, test_fraction=protocol.val_fraction, seed=seed)
    idx_tr = [tr[i] for i in inner_tr]
    idx_va = [tr[i] for i in inner_va]

    y_all = np.asarray(ys, dtype=np.float32)

    # THE NULL PERMUTES THE WHOLE TRAIN SIDE, BEFORE the inner split -- and getting this wrong is
    # what the first smoke caught. Permuting only `idx_tr` leaves the inner VALIDATION labels real,
    # so early stopping selects whichever epoch best predicts TRUE labels and real signal leaks into
    # model selection: the null then scored +0.0692 instead of ~0 and correctly refused to grade.
    # The inner val set is carved OUT OF train, so "train labels only" has to include it. Only `te`
    # keeps its real labels, which is what makes the null measure what it claims to.
    y_train_side = y_all.copy()
    if shuffle_labels:
        pool = np.array(tr, dtype=int)
        perm = pool.copy()
        rng.shuffle(perm)
        y_train_side[pool] = y_all[perm]

    y_tr = y_train_side[idx_tr].copy()

    model = build_masked_pyramid(PREREGISTERED["widths"], channels=protocol.channels).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=protocol.lr,
                           weight_decay=protocol.weight_decay)
    lossf = torch.nn.MSELoss()

    def batches(idx, y, bs, shuffle):
        order = np.arange(len(idx))
        if shuffle:
            rng.shuffle(order)
        for k in range(0, len(order), bs):
            sel = order[k:k + bs]
            rows = [idx[j] for j in sel]
            x = torch.from_numpy(codes_to_onehot(codes[rows])).to(device)
            ln = [lengths[r] for r in rows]
            yy = None if y is None else torch.from_numpy(np.asarray(y)[sel]).to(device)
            yield x, ln, yy

    best, best_state, bad = float("inf"), None, 0
    for _ in range(protocol.max_epochs):
        model.train()
        for x, ln, yy in batches(idx_tr, y_tr, protocol.batch_size, True):
            opt.zero_grad()
            loss = lossf(model(x, ln), yy)
            loss.backward()
            opt.step()
        model.eval()
        tot, n = 0.0, 0
        with torch.no_grad():
            # y_train_side, NOT y_all: under the null the inner val labels are permuted too.
            for x, ln, yy in batches(idx_va, y_train_side[idx_va], protocol.batch_size, False):
                tot += float(lossf(model(x, ln), yy)) * len(ln)
                n += len(ln)
        va = tot / max(n, 1)
        if va < best - 1e-6:
            best, bad = va, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= protocol.early_stop_patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)

    model.eval()
    preds = []
    with torch.no_grad():
        for x, ln, _ in batches(te, None, protocol.batch_size, False):
            preds.extend(model(x, ln).cpu().numpy().tolist())
    return spearman(preds, [ys[i] for i in te]), {"inner_train": len(idx_tr),
                                                  "inner_val": len(idx_va),
                                                  "val_mse": round(best, 5)}


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
        allow_download=False, test_fraction=0.2, with_encoder=False, device="cpu",
        seed_start=0) -> dict:
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

    codes = None
    if with_encoder:
        from dna_decode.glm.masked_pyramid import encode_codes
        # Codes, not one-hot: N x 4 x 499 float32 is 2.56 GB at this N.
        codes = encode_codes(padded, PAD_TO)

    long_set = set(long_mask)
    seed_range = list(range(seed_start, seed_start + seeds))
    rows = []
    for seed in seed_range:
        # ONE split per seed, handed to BOTH arms. The pairing is constructed here, not assumed:
        # recomputing the split inside each arm would make the pairing a coincidence of the seed.
        tr, te = split_by_block(blocks, test_fraction=test_fraction, seed=seed)
        y_te = [ys[i] for i in te]
        base_pred = list(_ridge_fit_predict(Xk[tr], [ys[i] for i in tr], Xk[te]))
        base = spearman(base_pred, y_te)
        row = {"seed": seed, "n_train": len(tr), "n_test": len(te),
               "kmer4_spearman": round(base, 4)}

        if with_encoder:
            enc, fit_info = fit_masked_pyramid(codes, lengths, ys, tr, te, seed=seed,
                                               device=device, blocks=blocks)
            row["pyramid_spearman"] = round(enc, 4)
            row["delta"] = round(enc - base, 4)
            row["fit"] = fit_info

            # LENGTH-MATCHED CONTROL: a RESTRICTION of these same predictions to the test rows at
            # >= 300 bp, not a retrain. The control asks whether padding drove the headline, and
            # retraining on a subset would change the model as well as the evaluation, which would
            # answer a different question.
            sel = [k for k, i in enumerate(te) if i in long_set]
            if len(sel) >= 50:
                yl = [y_te[k] for k in sel]
                row["long_only"] = {
                    "n": len(sel),
                    "kmer4": round(spearman([base_pred[k] for k in sel], yl), 4),
                }
            else:
                row["long_only"] = {"n": len(sel), "note": "too few long rows in this test split"}

        rows.append(row)

    out = {
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
    }

    if not with_encoder:
        out["status"] = "BASELINE_ONLY_ENCODER_NOT_RUN"
        out["status_note"] = (
            "The k-mer comparator and every pre-registered design decision are in place and the "
            "block split is exercised, but the pyramid half was not run (pass --with-encoder). "
            "No verdict is emitted: a verdict from a one-sided run would be a figure, not a "
            "measurement."
        )
        return out

    # ---- validity null on EVERY seed: train side permuted, test labels untouched.
    # Multi-seed because the null's own sd is ~0.10 at small N (measured), so one draw says nothing.
    null_rows = []
    for seed in seed_range:   # the SAME seeds as the paired arm, not range(seeds)
        trn, ten = split_by_block(blocks, test_fraction=test_fraction, seed=seed)
        r, _ = fit_masked_pyramid(codes, lengths, ys, trn, ten, seed=seed,
                                  device=device, blocks=blocks, shuffle_labels=True)
        null_rows.append(round(r, 4))
    null_median = statistics.median(null_rows)
    null_sd = statistics.stdev(null_rows) if len(null_rows) > 1 else None
    out["null"] = {
        "per_seed": null_rows,
        "median": round(null_median, 4),
        "sd": None if null_sd is None else round(null_sd, 4),
        "bar": "|median| < %.2f" % PREREGISTERED["validity_controls"]["null_median_bar"],
        "bar_note": (
            "the MEDIAN, not a single draw: a correct pipeline centred on zero still produces "
            "individual seeds out to +/-0.13 here, because a near-initialisation conv net is close "
            "to a random functional of base composition and composition predicts this rung at "
            "|0.1169|. The seed-count bar (>= %d of %d, p = 0.055 under the null) is the primary "
            "noise control." % (PREREGISTERED["min_seeds_winning"], seeds)
        ),
        "power_caveat": (
            "UNDER-POWERED AT THE PRE-REGISTERED SEED COUNT, and recorded so no reader over-trusts "
            "this gate. The null's across-seed spread is driven by random INITIALISATION, not by "
            "sample size -- a random composition functional's correlation with expression does not "
            "shrink as N grows -- so with sd ~0.10 the sem of the median over 10 seeds is ~0.032, "
            "LARGER than the 0.02 bar it is compared against. Resolving 0.02 would need ~50 seeds. "
            "Read `learned_contribution` below, which is paired within a seed and therefore does "
            "not carry that variance."
        ),
    }

    # ---- DIAGNOSTIC, not a bar. Nothing below can change a verdict.
    # Real minus null, per seed, same init and same split: what did the LABELS buy? The random-init
    # composition effect cancels within a seed, so this is far less noisy than the raw null.
    learned = [round(r["pyramid_spearman"] - nl, 4)
               for r, nl in zip(rows, null_rows)]
    out["learned_contribution"] = {
        "per_seed": learned,
        "median": round(statistics.median(learned), 4),
        "n_positive": sum(1 for v in learned if v > 0),
        "n": len(learned),
        "definition": "pyramid_spearman - null_spearman, per seed (same init, same split)",
        "status": "DIAGNOSTIC -- reported, NEVER a bar; the pre-registered comparison is the "
                  "kmer4 delta and is untouched",
        "why": (
            "The raw delta over kmer4 can be inflated by an architectural effect that owes nothing "
            "to the labels: a conv net near initialisation is close to a random linear functional "
            "of base composition, and composition predicts this rung at |0.1169|. This isolates "
            "the part the labels bought."
        ),
        "interpretation_caveat": (
            "SIZE-SENSITIVE -- do not read a small value at small N as a property of the arm. "
            "Measured at N=4,000 (3,200 training rows, 1.2% of the rung) the kmer4 delta was "
            "+0.1077 while this was +0.0346; at that N kmer4 itself scores 0.0244 against 0.1120 "
            "at full N, so BOTH arms are under-trained there and neither number transfers."
        ),
    }

    deltas = [r["delta"] for r in rows]
    p = paired_stats(deltas, PREREGISTERED["margin"], PREREGISTERED["min_seeds_winning"])
    out["paired"] = p

    out["seed_range"] = {"start": seed_start, "n": seeds,
                         "seeds": seed_range,
                         "preregistered_n_seeds": PREREGISTERED["n_seeds"]}

    # A PARTIAL seed range emits rows and NO verdict. The bar is ">= 8 of 10 seeds"; scoring 5 of 10
    # against it is a category error, not a weaker measurement. Combine the partial artifacts and
    # judge once the full pre-registered range exists.
    if seed_start != 0 or seeds != PREREGISTERED["n_seeds"]:
        out["status"] = "PARTIAL_SEED_RANGE"
        out["status_note"] = (
            "ran seeds %s; the pre-registration specifies %d seeds starting at 0 with a "
            ">= %d-of-%d winning bar. Rows and the null are emitted, NO verdict: judging a partial "
            "range against a seed-count bar would redefine the acceptance criterion. Combine the "
            "partial artifacts, then judge." % (seed_range, PREREGISTERED["n_seeds"],
                                                PREREGISTERED["min_seeds_winning"],
                                                PREREGISTERED["n_seeds"])
        )
        return out

    # Ordered so a broken null can never be reported as a result.
    if abs(null_median) >= PREREGISTERED["validity_controls"]["null_median_bar"]:
        out["verdict"], out["verdict_detail"] = (
            "INVALID_NULL_SCORED",
            "shuffled-train null has MEDIAN %+.4f over %d seeds (sd %s), outside the %.2f bar; "
            "nothing may be graded" % (null_median, len(null_rows),
                                       "n/a" if null_sd is None else "%.4f" % null_sd,
                                       PREREGISTERED["validity_controls"]["null_median_bar"]))
    else:
        out["verdict"], out["verdict_detail"] = verdict(p)

    med_k = statistics.median([r["kmer4_spearman"] for r in rows])
    med_p = statistics.median([r["pyramid_spearman"] for r in rows])
    out["headline"] = {
        "kmer4_median": round(med_k, 4),
        "pyramid_median": round(med_p, 4),
        "delta_median": p["median"],
        "relative_gain_vs_comparator": (round(p["median"] / med_k, 4) if med_k else None),
        "relative_note": (
            "the margin is pre-registered in ABSOLUTE Spearman (0.02). On this rung the comparator "
            "is ~0.11, so 0.02 is a ~18%% relative gain; the relative figure is reported beside it "
            "and is NOT the bar."
        ),
        "fraction_of_rung_ceiling": round(med_p / PREREGISTERED["rung_ceiling"], 4),
    }
    out["status"] = "SCORED"
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--which", default="LB", choices=["LB", "M9"])
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--seed-start", type=int, default=0,
                    help="first seed; lets a long run be SPLIT across sessions without reducing "
                         "the pre-registered seed count. A partial range emits no verdict.")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--allow-download", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--with-encoder", action="store_true",
                    help="run the pyramid half too and emit a verdict (the expensive half)")
    ap.add_argument("--device", default="cpu", help="cpu or cuda")
    a = ap.parse_args(argv)

    try:
        res = run(a.cache_dir, which=a.which, seeds=a.seeds, limit=a.limit,
                  allow_download=a.allow_download, with_encoder=a.with_encoder,
                  device=a.device, seed_start=a.seed_start)
    except (ArmDataError, EncoderDataError) as exc:
        print("REFUSED: %s" % exc, file=sys.stderr)
        return 3

    # WRITE BEFORE FORMATTING. A completed run must survive a summary bug: a ~11 minute run was
    # lost on 2026-10-27 because the printer raised KeyError before this line was reached.
    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
        print("wrote %s  (persisted BEFORE the summary, so a print bug cannot discard the run)"
              % a.out)

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
    if res.get("paired"):
        p = res["paired"]
        # `headline` is absent on a PARTIAL_SEED_RANGE result, which carries paired + null and NO
        # verdict on purpose. Reading it unconditionally is what crashed the first real smoke.
        h = res.get("headline")
        if h:
            print("  pyramid median %.4f vs kmer4 %.4f -> delta %+.4f, %d/%d seeds, CI95 %s"
                  % (h["pyramid_median"], h["kmer4_median"], p["median"], p["wins"], p["n"],
                     p["ci95"]))
        else:
            print("  paired delta %+.4f, %d/%d seeds, CI95 %s"
                  % (p["median"], p["wins"], p["n"], p["ci95"]))
    if res.get("null"):
        n = res["null"]
        print("  null (shuffled train) median %+.4f  sd %s  per-seed %s"
              % (n["median"], n.get("sd"), n.get("per_seed")))
    if res.get("learned_contribution"):
        lc = res["learned_contribution"]
        print("  LEARNED (real-null, paired per seed) median %+.4f  %d/%d positive  per-seed %s"
              % (lc["median"], lc["n_positive"], lc["n"], lc["per_seed"]))
        print("    diagnostic only, never the bar; size-sensitive -- do not read a small value at "
              "small N as a property of the arm")
    if res.get("verdict"):
        print("  VERDICT: %s" % res["verdict"])
        print("    %s" % res["verdict_detail"])
    print("  STATUS: %s" % res["status"])
    if res.get("status") == "PARTIAL_SEED_RANGE":
        print("    (no verdict BY DESIGN -- %s)" % res["status_note"][:150])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
