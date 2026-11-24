"""Does corpus-window POSITION explain PromoterAtlas's 0.5076 on the peak-tile rung? Three cells.

The subtractive half of this investigation removed the control-vs-tile detection direction and the
atlas channel barely moved. This is the ADDITIVE half, which is the one a sceptic asks for: hand a
cheap baseline the exact positional rule that generated the pretraining corpus and see whether the
margin is absorbed.

  C1  Does the encoder KNOW corpus membership? Predict `in_window_any` from the 128-dim embedding,
      held out over the committed splits, against the GC-only equivalent. The positional
      restatement of the already-measured 0.9156 control-vs-tile separability.
  C2  Does handing the frame to the baseline absorb the gain? The identical ridge is fitted on gc,
      gc+dinucleotide, gc+kmer4, gc+frame and atlas, and the atlas margin is reported over EACH.
      The series is the control that separates absorption from "the baseline simply got better".
  C3  The stratified cell, and the strongest: within in-window tiles only, and within out-of-window
      tiles only, does atlas still beat GC? This controls for frame membership BY CONSTRUCTION, so
      it does not require the frame to be linearly modellable. Its null is a SIZE-MATCHED RANDOM
      stratification, because stratifying at all compresses the outcome and attenuates a rank
      correlation for BOTH channels.

**Nothing is re-embedded, and nothing touches a GPU or the network.** The cached pass at
`data/promoteratlas/tile_embeddings_native150_n44106.npy` is reused, so the 15.7-minute embedding
run is not repeated.

**The anchor runs FIRST.** Before any new number is read, the cached embeddings must reproduce the
committed `arms.native150` on the committed splits. A harness that drifted since the preceding run
would otherwise produce a confident frame verdict about a representation that is no longer the one
the arm measured.

**Every frame-bearing cell is measured on the LOCATED SUBSET, and its own baseline is re-measured
there too.** 356 distinct sequences are ambiguous or not found, so the subset is not the full 44,106
rows and its atlas number is NOT 0.5076. Comparing a subset margin against the full-cohort figure
would confound the frame with the change in n and in peak composition, which is the bar's
pre-commitment (c).

**The strengthening check uses ORDINARY least squares on purpose.** Adding a column cannot lower an
OLS training R-squared -- that is exact -- whereas a ridge penalty can move it a hair in either
direction, so a strict inequality on the ridge figure could fire falsely and report a plumbing
fault that is really a rounding artifact. The gate therefore reads the OLS training R-squared as a
pure design-matrix diagnostic ("did the baseline actually receive the frame?"), while every SCORED
number in this run comes from the unchanged `ridge_spearman` the incumbent arm used.

Per-tile vectors and per-seed detail land in gitignored `data/promoteratlas/`; the wiki artifact
carries medians only, because a 102,755-decimal dump broke two contract-number guards on
2026-11-01.

Run:  uv run python -m scripts.promoteratlas_frame_conditioning
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
ARM = REPO / "wiki/promoteratlas_tile_arm_2026-11-01.json"
BAR = REPO / "wiki/promoteratlas_frame_acceptance_bar.json"
CACHE = REPO / "data/mpra"
EMB_DIR = REPO / "data/promoteratlas"
FNA = REPO / "data/cache/refseq/GCF_000005845.2/genome.fna"
GFF = REPO / "data/cache/refseq/GCF_000005845.2/annotations.gff3"
LOCATE_CACHE = EMB_DIR / "frame_locate_cache.json"
OUT = REPO / "wiki/promoteratlas_frame_conditioning_2026-11-01.json"
DETAIL = EMB_DIR / "frame_conditioning_detail_2026-11-01.json"

SCHEMA = "promoteratlas-frame-conditioning-v1"
N_SEEDS = 10
FRAC = 0.25
TILE_LEN = 150
ANCHOR_TOL = 0.02


def _median(xs) -> float:
    return float(np.median(np.asarray(xs, dtype=float)))


def _spread(xs) -> dict:
    a = np.asarray(sorted(float(x) for x in xs), dtype=float)
    return {
        "n": int(a.size),
        "median": round(float(np.median(a)), 4),
        "q1": round(float(np.percentile(a, 25)), 4),
        "q3": round(float(np.percentile(a, 75)), 4),
        "iqr": round(float(np.percentile(a, 75) - np.percentile(a, 25)), 4),
        "sd": round(float(a.std(ddof=1)) if a.size > 1 else 0.0, 4),
    }


def _stability(xs, bar: float) -> dict:
    """Seed spread of a statistic against the bar it is graded on.

    The grader reads POINT MEDIANS. When a statistic's seed spread straddles its bar, that median
    is a coin flip and the verdict it produces is not resolved by the run -- so every statistic the
    verdict turns on ships this block beside it. `n_seeds_below_bar` near half the seeds is the
    signature to look for.
    """
    a = np.asarray([float(x) for x in xs], dtype=float)
    a = a[~np.isnan(a)]
    if a.size == 0:
        return {"n_seeds": 0, "note": "no usable seeds"}
    below = int((a < bar).sum())
    return {
        "n_seeds": int(a.size),
        "median": round(float(np.median(a)), 4),
        "min": round(float(a.min()), 4),
        "max": round(float(a.max()), 4),
        "n_seeds_below_bar": below,
        "n_seeds_negative": int((a < 0).sum()),
        "straddles_bar": 0 < below < a.size,
        "straddles_zero": bool(a.min() < 0 < a.max()),
        "bar": bar,
    }


def _ols_train_r2(X, y_raw) -> float:
    """Training R-squared of ORDINARY least squares on log1p(y). Exactly monotone in columns.

    A design-matrix diagnostic only -- it answers "did the baseline receive the feature?" and is
    never a scored number.
    """
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import StandardScaler

    y = np.log1p(np.asarray(y_raw, dtype=float))
    sc = StandardScaler().fit(X)
    Z = sc.transform(X)
    return float(LinearRegression().fit(Z, y).score(Z, y))


def _logistic_auroc(Xtr, ytr, Xte, yte) -> float:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    from scripts.promoteratlas_subspace_removal import auroc

    if len(set(int(v) for v in ytr)) < 2 or len(set(int(v) for v in yte)) < 2:
        return float("nan")
    sc = StandardScaler().fit(Xtr)
    clf = LogisticRegression(max_iter=1000, C=1.0).fit(sc.transform(Xtr), ytr)
    return auroc(clf.decision_function(sc.transform(Xte)), yte)


def _stratum_margin(Xa, Xg, y, tr_idx, te_idx):
    """Atlas and GC re-measured WITHIN one stratum, on that stratum's own split."""
    from dna_decode.glm.atlas_embed import ridge_spearman

    if len(tr_idx) < 50 or len(te_idx) < 50:
        return None
    a, _ = ridge_spearman(Xa[tr_idx], y[tr_idx], Xa[te_idx], y[te_idx])
    g, _ = ridge_spearman(Xg[tr_idx], y[tr_idx], Xg[te_idx], y[te_idx])
    return float(a), float(g)


def run(verbose: bool = True) -> dict:
    from dna_decode.data.annotations import parse_gff3
    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.atlas_embed import ridge_spearman
    from dna_decode.glm.corpus_window import (
        cds_rows_from_annotation_table,
        corpus_windows,
        window_coverage,
    )
    from dna_decode.glm.frame_features import FEATURE_NAMES, frame_features
    from dna_decode.glm.tile_locate import (
        load_locate_cache,
        locate_report,
        locate_sequences,
        read_fasta_single,
    )

    t0 = time.time()
    arm = json.loads(ARM.read_text(encoding="utf-8"))
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    committed_atlas = float(arm["arms"]["native150"]["median"])
    committed_gc = float(arm["gc_anchor"]["median"])
    margin = float(bar["inherited_margin"])
    primary = bar["primary_frame_granularity"]

    # ---- substrate, embeddings, frame -------------------------------------------------------
    rows, _ = gw.load_peak_tiles(CACHE, allow_download=False)
    tiles = [t for t in rows if t.category == "tile" and len(t.seq) == TILE_LEN]
    n = len(tiles)
    Xa = np.load(EMB_DIR / ("tile_embeddings_native150_n%d.npy" % n))
    if Xa.shape != (n, 128):
        raise SystemExit("embedding cache is %r, expected (%d, 128)" % (Xa.shape, n))

    _contig, genome = read_fasta_single(FNA)
    cds = cds_rows_from_annotation_table(parse_gff3(GFF))
    mask = window_coverage(corpus_windows(cds, len(genome)), len(genome))

    if LOCATE_CACHE.is_file():
        located = load_locate_cache(LOCATE_CACHE)
        if {t.seq for t in tiles} - set(located):
            located = locate_sequences([t.seq for t in tiles], genome)
    else:
        located = locate_sequences([t.seq for t in tiles], genome)
    lrep = locate_report(located, tile_length=TILE_LEN, genome_length=len(genome))
    table = frame_features(located, mask, cds, tile_length=TILE_LEN)
    by_seq = {s: i for i, s in enumerate(table.seqs)}

    y = np.asarray([t.expression for t in tiles], dtype=float)
    gc = np.asarray(gw.gc_features([t.seq for t in tiles]), dtype=float)
    if verbose:
        print("building k-mer baselines ...", flush=True)
    k2 = np.asarray(gw.kmer_features([t.seq for t in tiles], k=2), dtype=float)
    k4 = np.asarray(gw.kmer_features([t.seq for t in tiles], k=4), dtype=float)

    # frame rows per TILE; NaN where the tile's sequence is not uniquely located, so that an
    # accidental inclusion raises in the scaler instead of quietly scoring an invented frame
    F = np.full((n, len(FEATURE_NAMES)), np.nan, dtype=float)
    has_frame = np.zeros(n, dtype=bool)
    for i, t in enumerate(tiles):
        j = by_seq.get(t.seq)
        if j is not None:
            F[i] = table.matrix[j]
            has_frame[i] = True
    n_located_tiles = int(has_frame.sum())
    n_excluded_tiles = n - n_located_tiles

    in_window = np.zeros(n, dtype=float)
    in_full = np.zeros(n, dtype=float)
    pj = list(FEATURE_NAMES).index(primary)
    fj = list(FEATURE_NAMES).index("in_window_full")
    in_window[has_frame] = F[has_frame][:, pj]
    in_full[has_frame] = F[has_frame][:, fj]

    gc_frame = np.hstack([gc, F])
    gc_k2 = np.hstack([gc, k2])
    gc_k4 = np.hstack([gc, k4])

    if verbose:
        print("tiles %d, embeddings %s, located tiles %d, excluded %d (%.0fs)"
              % (n, Xa.shape, n_located_tiles, n_excluded_tiles, time.time() - t0), flush=True)

    per_seed: dict[str, list] = {k: [] for k in (
        "anchor_atlas_full", "gc_full",
        "atlas", "gc", "gc_dinuc", "gc_kmer4", "gc_frame", "frame_alone",
        "dinuc_alone", "kmer4_alone",
        "train_r2_gc", "train_r2_gc_frame",
        "c1_auroc_embedding", "c1_auroc_gc",
        "margin_in_window", "margin_out_window",
        "atlas_in_window", "gc_in_window", "atlas_out_window", "gc_out_window",
        "random_stratum_min_margin",
        "margin_in_full", "margin_out_full",
    )}

    for seed in range(N_SEEDS):
        train, test, _ = gw.leave_peak_out(tiles, frac=FRAC, seed=seed)
        pos = {id(t): i for i, t in enumerate(tiles)}
        tr = np.asarray([pos[id(t)] for t in train])
        te = np.asarray([pos[id(t)] for t in test])

        # ---- anchor on the FULL cohort, before any new number is read ----------------------
        a_full, _ = ridge_spearman(Xa[tr], y[tr], Xa[te], y[te])
        g_full, _ = ridge_spearman(gc[tr], y[tr], gc[te], y[te])
        per_seed["anchor_atlas_full"].append(round(a_full, 4))
        per_seed["gc_full"].append(round(g_full, 4))

        # ---- the located subset, every baseline re-measured there --------------------------
        trL = tr[has_frame[tr]]
        teL = te[has_frame[te]]

        for key, X in (("atlas", Xa), ("gc", gc), ("gc_dinuc", gc_k2),
                       ("gc_kmer4", gc_k4), ("gc_frame", gc_frame), ("frame_alone", F),
                       ("dinuc_alone", k2), ("kmer4_alone", k4)):
            v, _ = ridge_spearman(X[trL], y[trL], X[teL], y[teL])
            per_seed[key].append(round(float(v), 4))

        per_seed["train_r2_gc"].append(round(_ols_train_r2(gc[trL], y[trL]), 6))
        per_seed["train_r2_gc_frame"].append(round(_ols_train_r2(gc_frame[trL], y[trL]), 6))

        # ---- C1: does the embedding know corpus membership? --------------------------------
        lab_tr, lab_te = in_window[trL], in_window[teL]
        per_seed["c1_auroc_embedding"].append(
            round(_logistic_auroc(Xa[trL], lab_tr, Xa[teL], lab_te), 4))
        per_seed["c1_auroc_gc"].append(
            round(_logistic_auroc(gc[trL], lab_tr, gc[teL], lab_te), 4))

        # ---- C3: the real strata, each baseline re-measured ON ITS OWN stratum -------------
        for tag, sel in (("in_window", in_window[trL] > 0.5), ("out_window", in_window[trL] <= 0.5)):
            sel_te = (in_window[teL] > 0.5) if tag == "in_window" else (in_window[teL] <= 0.5)
            got = _stratum_margin(Xa, gc, y, trL[sel], teL[sel_te])
            if got is None:
                per_seed["margin_%s" % tag].append(float("nan"))
                per_seed["atlas_%s" % tag].append(float("nan"))
                per_seed["gc_%s" % tag].append(float("nan"))
            else:
                a, g = got
                per_seed["atlas_%s" % tag].append(round(a, 4))
                per_seed["gc_%s" % tag].append(round(g, 4))
                per_seed["margin_%s" % tag].append(round(a - g, 4))

        # secondary granularity, REPORTED not graded
        for tag, thr in (("in_full", True), ("out_full", False)):
            s_tr = (in_full[trL] > 0.5) if thr else (in_full[trL] <= 0.5)
            s_te = (in_full[teL] > 0.5) if thr else (in_full[teL] <= 0.5)
            got = _stratum_margin(Xa, gc, y, trL[s_tr], teL[s_te])
            per_seed["margin_%s" % tag].append(
                float("nan") if got is None else round(got[0] - got[1], 4))

        # ---- the SIZE-MATCHED RANDOM stratification null -----------------------------------
        rng = np.random.default_rng(770_000 + seed)
        k_tr = int((in_window[trL] > 0.5).sum())
        k_te = int((in_window[teL] > 0.5).sum())
        ptr = rng.permutation(len(trL))
        pte = rng.permutation(len(teL))
        rand_margins = []
        for a_tr, a_te in ((ptr[:k_tr], pte[:k_te]), (ptr[k_tr:], pte[k_te:])):
            got = _stratum_margin(Xa, gc, y, trL[a_tr], teL[a_te])
            if got is not None:
                rand_margins.append(got[0] - got[1])
        per_seed["random_stratum_min_margin"].append(
            round(float(min(rand_margins)), 4) if rand_margins else float("nan"))

        if verbose:
            print("  seed %d done (%.0fs)" % (seed, time.time() - t0), flush=True)

    med = {k: _median(v) for k, v in per_seed.items()}

    # ---- grade --------------------------------------------------------------------------------
    from dna_decode.glm.frame_verdict import FrameInput, grade

    fi = FrameInput(
        harness_reconciled=abs(med["anchor_atlas_full"] - committed_atlas) <= ANCHOR_TOL,
        unique_locate_rate=lrep.unique_rate,
        frame_alone_spearman=med["frame_alone"],
        train_fit_gc=med["train_r2_gc"],
        train_fit_gc_frame=med["train_r2_gc_frame"],
        atlas=med["atlas"],
        gc_baseline=med["gc"],
        gc_plus_frame=med["gc_frame"],
        gc_plus_matched=max(med["gc_dinuc"], med["gc_kmer4"]),
        atlas_in_window=med["atlas_in_window"],
        gc_in_window=med["gc_in_window"],
        atlas_out_window=med["atlas_out_window"],
        gc_out_window=med["gc_out_window"],
        random_stratum_min_margin=med["random_stratum_min_margin"],
        best_matched_alone=max(med["dinuc_alone"], med["kmer4_alone"]),
        margin=margin,
    )
    v = grade(fi)

    expr_in = [float(y[i]) for i in range(n) if has_frame[i] and in_window[i] > 0.5]
    expr_out = [float(y[i]) for i in range(n) if has_frame[i] and in_window[i] <= 0.5]

    rec = {
        "schema": SCHEMA,
        "date": "2026-11-05",
        "verdict": v.verdict,
        "verdict_detail": v.detail,
        "is_valid": v.is_valid,
        "verdict_stability": {
            "the_verdict_turns_on": (
                "C3's weaker stratum margin against the inherited %.2f bar, and C2's surviving "
                "margin against the same bar." % margin
            ),
            "c3_weaker_stratum_straddles_the_bar": bool(
                _stability(per_seed["margin_out_window"], margin)["straddles_bar"]
            ),
            "c2_surviving_margin_clears_the_bar_by": round(fi.margin_over_gc_frame - margin, 4),
            "resolved": not _stability(per_seed["margin_out_window"], margin)["straddles_bar"],
            "reading": (
                "If the weaker stratum's margin straddles the bar across seeds, then the point "
                "median that selected FRAME_PARTIAL over FRAME_EXPLAINS is a coin flip and THIS "
                "RUN DOES NOT RESOLVE THE DISTINCTION. The frozen verdict is still reported as "
                "computed -- re-sizing a bar after seeing the result is the move the bar exists to "
                "prevent -- but it may not be quoted as though the stratified cell came out clean. "
                "Whether the bar should require an interval on the stratum margins rather than a "
                "point median is a change to a frozen bar and therefore a user decision."
            ),
        },
        "bar": "wiki/promoteratlas_frame_acceptance_bar.json",
        "inherited_margin": margin,
        "primary_frame_granularity": primary,
        "n_seeds": N_SEEDS,
        "anchor": {
            "committed_atlas_native150_median": committed_atlas,
            "local_atlas_full_cohort_median": round(med["anchor_atlas_full"], 4),
            "abs_delta": round(abs(med["anchor_atlas_full"] - committed_atlas), 4),
            "tolerance": ANCHOR_TOL,
            "reconciled": fi.harness_reconciled,
            "committed_gc_median": committed_gc,
            "local_gc_full_cohort_median": round(med["gc_full"], 4),
        },
        "accounting": {
            "n_tiles_scored_by_the_arm": n,
            "n_tiles_with_a_frame": n_located_tiles,
            "n_tiles_excluded": n_excluded_tiles,
            "reconciles": n_located_tiles + n_excluded_tiles == n,
            "distinct_sequences": lrep.n_distinct,
            "distinct_unique": lrep.n_unique,
            "distinct_ambiguous": lrep.n_ambiguous,
            "distinct_not_found": lrep.n_not_found,
            "unique_locate_rate": round(lrep.unique_rate, 6),
            "feature_table": table.accounting(),
        },
        "c1_does_the_encoder_know_corpus_membership": {
            "auroc_embedding_median": round(med["c1_auroc_embedding"], 4),
            "auroc_gc_only_median": round(med["c1_auroc_gc"], 4),
            "reading": (
                "The positional restatement of the measured 0.9156 control-vs-tile separability. A "
                "high value says the embedding CARRIES corpus membership; it does not by itself "
                "say the expression gain is made of it, which is what C2 and C3 test."
            ),
        },
        "c2_matched_strength_series": {
            "atlas_median": round(med["atlas"], 4),
            "gc_median": round(med["gc"], 4),
            "gc_plus_dinucleotide_median": round(med["gc_dinuc"], 4),
            "gc_plus_kmer4_median": round(med["gc_kmer4"], 4),
            "gc_plus_frame_median": round(med["gc_frame"], 4),
            "frame_alone_median": round(med["frame_alone"], 4),
            "margin_over_gc": round(fi.margin_over_gc, 4),
            "margin_over_gc_plus_frame": round(fi.margin_over_gc_frame, 4),
            "margin_over_best_matched_nonpositional": round(fi.margin_over_gc_matched, 4),
            "frame_absorption": round(fi.frame_absorption, 4),
            "matched_absorption": round(fi.matched_absorption, 4),
            "attributable_absorption": round(fi.attributable_absorption, 4),
            "train_r2_gc_median": round(med["train_r2_gc"], 6),
            "train_r2_gc_frame_median": round(med["train_r2_gc_frame"], 6),
            "baseline_was_strengthened": med["train_r2_gc_frame"] >= med["train_r2_gc"],
            "matched_control_premise": {
                "frame_alone_median": round(med["frame_alone"], 4),
                "dinucleotide_alone_median": round(med["dinuc_alone"], 4),
                "kmer4_alone_median": round(med["kmer4_alone"], 4),
                "best_matched_alone": round(max(med["dinuc_alone"], med["kmer4_alone"]), 4),
                "control_is_at_least_as_strong_as_the_frame": (
                    max(med["dinuc_alone"], med["kmer4_alone"]) >= med["frame_alone"]
                ),
                "why_this_matters": (
                    "The control's whole premise is MATCHED INFORMATIVENESS: 'shrinkage is not "
                    "absorption unless it exceeds what an equally strong NON-positional feature "
                    "produces'. If the non-positional comparators are WEAKER standalone than the "
                    "frame, then they were never a matched control and the attributable-absorption "
                    "figure is only weakly controlled -- it would partly measure the frame being a "
                    "better feature rather than a POSITIONAL one. These standalone numbers are "
                    "what decides that, and the first run of this experiment did not measure them."
                ),
            },
            "reading": (
                "Every number here is measured on the LOCATED SUBSET with its own baseline, so the "
                "atlas median is not the full-cohort 0.5076 and must not be compared to it. "
                "`attributable_absorption` is what the frame takes IN EXCESS of an equally strong "
                "non-positional feature; shrinkage below that is the baseline getting better."
            ),
        },
        "c3_stratified": {
            "atlas_in_window_median": round(med["atlas_in_window"], 4),
            "gc_in_window_median": round(med["gc_in_window"], 4),
            "margin_in_window": round(med["margin_in_window"], 4),
            "atlas_out_window_median": round(med["atlas_out_window"], 4),
            "gc_out_window_median": round(med["gc_out_window"], 4),
            "margin_out_window": round(med["margin_out_window"], 4),
            "min_stratum_margin": round(fi.min_stratum_margin, 4),
            "random_stratum_min_margin_median": round(med["random_stratum_min_margin"], 4),
            "secondary_granularity_in_window_full": {
                "margin_in_full": round(med["margin_in_full"], 4),
                "margin_out_full": round(med["margin_out_full"], 4),
                "status": "REPORTED, not graded; the bar grades the named primary only",
            },
            "stability_margin_out_window": _stability(per_seed["margin_out_window"], margin),
            "stability_margin_in_window": _stability(per_seed["margin_in_window"], margin),
            "stability_random_stratum_min_margin": _stability(
                per_seed["random_stratum_min_margin"], margin),
            "median_of_differences_vs_difference_of_medians": {
                "median_of_per_seed_out_window_margins": round(med["margin_out_window"], 4),
                "difference_of_out_window_medians": round(fi.margin_out_window, 4),
                "note": (
                    "These are not the same statistic and they do not have to agree. The grader "
                    "reads the difference of medians; both are reported because a gap between "
                    "them is itself evidence the underlying quantity is unstable across seeds."
                ),
            },
            "expression_spread_in_window": _spread(expr_in),
            "expression_spread_out_window": _spread(expr_out),
            "reading": (
                "Each stratum's baseline is re-measured ON THAT STRATUM, never against the "
                "full-cohort figure. The spreads ship beside the margins because the strata are "
                "measurably not exchangeable, and without them a reader cannot tell range "
                "restriction from absorption. The random-stratification median is the null for how "
                "much stratifying alone moves a margin."
            ),
        },
        "honest_limits": bar["honest_limits_fixed_in_advance"],
        "forecast_scored": {
            "forecast": {k: v for k, v in bar["forecast"].items() if k.startswith("p_")},
            "outcome": v.verdict,
            "p_assigned_to_outcome": bar["forecast"].get("p_%s" % v.verdict),
        },
        "wall_seconds": round(time.time() - t0, 1),
    }

    DETAIL.parent.mkdir(parents=True, exist_ok=True)
    DETAIL.write_text(json.dumps({"schema": SCHEMA + "-detail", "per_seed": per_seed},
                                 indent=2, sort_keys=True), encoding="utf-8")
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    for p in (ARM, BAR, FNA, GFF):
        if not p.is_file():
            print("REFUSED: %s absent" % p)
            return 3

    rec = run(verbose=not a.quiet)
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    print()
    print("ANCHOR   local %.4f vs committed %.4f (delta %.4f, tol %.2f) -> reconciled=%s" % (
        rec["anchor"]["local_atlas_full_cohort_median"],
        rec["anchor"]["committed_atlas_native150_median"],
        rec["anchor"]["abs_delta"], rec["anchor"]["tolerance"], rec["anchor"]["reconciled"]))
    c2 = rec["c2_matched_strength_series"]
    print("C1       embedding AUROC %.4f vs GC-only %.4f" % (
        rec["c1_does_the_encoder_know_corpus_membership"]["auroc_embedding_median"],
        rec["c1_does_the_encoder_know_corpus_membership"]["auroc_gc_only_median"]))
    print("C2       atlas %.4f | gc %.4f | +dinuc %.4f | +kmer4 %.4f | +frame %.4f" % (
        c2["atlas_median"], c2["gc_median"], c2["gc_plus_dinucleotide_median"],
        c2["gc_plus_kmer4_median"], c2["gc_plus_frame_median"]))
    print("         frame alone %.4f | margin over gc %.4f -> over gc+frame %.4f" % (
        c2["frame_alone_median"], c2["margin_over_gc"], c2["margin_over_gc_plus_frame"]))
    print("         absorption frame %.4f vs matched %.4f -> attributable %+.4f" % (
        c2["frame_absorption"], c2["matched_absorption"], c2["attributable_absorption"]))
    c3 = rec["c3_stratified"]
    print("C3       in-window margin %.4f | out-of-window %.4f | random null min %.4f" % (
        c3["margin_in_window"], c3["margin_out_window"],
        c3["random_stratum_min_margin_median"]))
    st = c3["stability_margin_out_window"]
    print("         out-of-window per-seed: min %.4f max %.4f, %d of %d seeds below the bar, "
          "%d negative" % (st["min"], st["max"], st["n_seeds_below_bar"], st["n_seeds"],
                           st["n_seeds_negative"]))
    mc = c2["matched_control_premise"]
    print("CONTROL  standalone: frame %.4f | dinuc %.4f | kmer4 %.4f -> control at least as "
          "strong = %s" % (mc["frame_alone_median"], mc["dinucleotide_alone_median"],
                           mc["kmer4_alone_median"],
                           mc["control_is_at_least_as_strong_as_the_frame"]))
    print()
    print("VERDICT  %s   (resolved by this run: %s)" % (
        rec["verdict"], rec["verdict_stability"]["resolved"]))
    print(rec["verdict_detail"])
    print()
    print("-> %s  (%.0fs)" % (OUT, rec["wall_seconds"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
