"""Experiment A-extended — remove the detection SUBSPACE, not one direction, and watch both curves.

**Why this exists: Experiment A's decisive cell was underpowered by construction, and its own numbers
say so.** A1 measured control-vs-tile separability at AUROC 0.9156 against GC's 0.7641, so the
promoter-DETECTION signal is unmistakably present in the embedding. A2 then removed the single best
separating direction and the atlas channel did not move (drop -0.0019, against 0.0001 for a random
direction). Those two facts are only consistent if the detection signal is DISTRIBUTED over many of
the 128 dimensions -- in which case deleting one axis was never going to remove it. Worse, two
defensible definitions of "the separation direction" disagreed by 8x (logistic -0.0019 vs
mean-difference 0.0163), which is what a sensitivity-to-construction looks like.

So a 1-D null is weak evidence, and reporting EXPOSURE_EXCLUDED_LINEAR off it would overstate what
was tested. This runs the version that can actually fail:

  * Build an ORDERED orthonormal basis of discriminative directions by DEFLATION -- fit a separator,
    remove it, refit on the residual, repeat.
  * For each k, project out span(w_1..w_k) from the tile embeddings and refit the identical ridge on
    the identical committed splits.
  * At every k report BOTH curves: the surviving atlas Spearman AND the RESIDUAL control-vs-tile
    separability.

**The two curves together are the whole point, and neither alone decides anything.** If separability
falls to chance while the atlas number holds, the detection signal has been genuinely removed and the
expression signal did not depend on it -- that is real evidence against the exposure explanation. If
the two fall together, exposure explains the gain. If the atlas falls while separability stays high,
the projection is just destroying capacity, and nothing is attributable.

**The random-subspace control is load-bearing at every k.** Removing k of 128 dimensions costs
something by itself, and that cost grows with k. Without subtracting it, a decline in the atlas curve
reads as exposure when it is only lost capacity.

Directions are fit on TRAIN tiles + TRAIN controls only, never on the test fold, and never using the
expression target.
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
CACHE = REPO / "data/mpra"
EMB_DIR = REPO / "data/promoteratlas"
OUT = REPO / "wiki/promoteratlas_subspace_removal_2026-11-01.json"

KS = (1, 2, 4, 8, 16, 32)
N_SEEDS = 10
FRAC = 0.25
FIT_SUBSAMPLE = 8000  # tiles used to FIT directions; the ridge always uses every tile


def auroc(scores, labels) -> float:
    s = np.asarray(scores, dtype=float)
    y = np.asarray(labels, dtype=int)
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[order[j + 1]] == s[order[i]]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def unit(v):
    n = float(np.linalg.norm(v))
    return v / n if n > 0 else v


def project_out_basis(X, B):
    """Remove the span of the orthonormal rows of B from every row of X."""
    if B is None or len(B) == 0:
        return X
    B = np.asarray(B, dtype=float)
    return X - (X @ B.T) @ B


def sep_direction(Xa, Xb, seed=0):
    from sklearn.linear_model import LogisticRegression

    X = np.vstack([Xa, Xb])
    y = np.concatenate([np.zeros(len(Xa)), np.ones(len(Xb))])
    sd = X.std(axis=0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    clf = LogisticRegression(max_iter=1000, C=1.0, random_state=seed).fit(X / sd, y)
    return unit(clf.coef_[0] / sd)


def deflation_basis(Xc, Xt, kmax, seed=0):
    """Ordered orthonormal discriminative basis, by fit-remove-refit."""
    B = []
    A, Bm = Xc.copy(), Xt.copy()
    for _ in range(kmax):
        w = sep_direction(A, Bm, seed=seed)
        # orthogonalise against what is already removed (Gram-Schmidt)
        for b in B:
            w = w - float(w @ b) * b
        n = float(np.linalg.norm(w))
        if n < 1e-8:
            break
        w = w / n
        B.append(w)
        A = A - np.outer(A @ w, w)
        Bm = Bm - np.outer(Bm @ w, w)
    return np.asarray(B)


def run(verbose=True):
    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.atlas_embed import ridge_spearman
    from dna_decode.glm.tile_controls import load_control_tiles

    t0 = time.time()
    arm = json.loads(ARM.read_text(encoding="utf-8"))
    committed = float(arm["arms"]["native150"]["median"])
    gc_committed = float(arm["gc_anchor"]["median"])

    tiles_all, _ = gw.load_peak_tiles(CACHE, allow_download=False)
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    controls = load_control_tiles(CACHE)

    Xt = np.load(EMB_DIR / ("tile_embeddings_native150_n%d.npy" % len(tiles)))
    Xc = np.load(EMB_DIR / ("control_embeddings_native150_n%d.npy" % len(controls)))
    if verbose:
        print("tiles %s controls %s  k sweep %s" % (Xt.shape, Xc.shape, list(KS)), flush=True)

    idx = {id(t): i for i, t in enumerate(tiles)}
    cperm = np.random.default_rng(12345).permutation(len(controls))
    kmax = max(KS)

    rows = {k: {"atlas": [], "sep": [], "atlas_rand": [], "sep_rand": []} for k in KS}
    base_per_seed, gc_per_seed = [], []

    for seed in range(N_SEEDS):
        train, test, _ = gw.leave_peak_out(tiles, frac=FRAC, seed=seed)
        tr = np.asarray([idx[id(t)] for t in train])
        te = np.asarray([idx[id(t)] for t in test])
        ytr = [tiles[i].expression for i in tr]
        yte = [tiles[i].expression for i in te]
        kc = int(round(len(controls) * (1 - FRAC)))
        c_tr, c_te = cperm[:kc], cperm[kc:]

        rng = np.random.default_rng(31337 + seed)
        fit_rows = tr if len(tr) <= FIT_SUBSAMPLE else rng.choice(tr, FIT_SUBSAMPLE, replace=False)

        base, _ = ridge_spearman(Xt[tr], ytr, Xt[te], yte)
        base_per_seed.append(round(base, 4))
        gcm = np.asarray(gw.gc_features([t.seq for t in tiles]), dtype=float)
        g, _ = ridge_spearman(gcm[tr], ytr, gcm[te], yte)
        gc_per_seed.append(round(g, 4))

        B = deflation_basis(Xc[c_tr], Xt[fit_rows], kmax, seed=seed)
        # a RANDOM orthonormal basis of the same rank, as the capacity-cost control
        R = np.linalg.qr(rng.normal(size=(Xt.shape[1], kmax)))[0].T

        lab = np.concatenate([np.zeros(len(c_te)), np.ones(len(te))])
        for k in KS:
            if k > len(B):
                continue
            Bk = B[:k]
            a, _ = ridge_spearman(
                project_out_basis(Xt[tr], Bk), ytr, project_out_basis(Xt[te], Bk), yte
            )
            w = sep_direction(project_out_basis(Xc[c_tr], Bk),
                              project_out_basis(Xt[fit_rows], Bk), seed=seed)
            s = auroc(
                np.concatenate([project_out_basis(Xc[c_te], Bk) @ w,
                                project_out_basis(Xt[te], Bk) @ w]), lab)

            Rk = R[:k]
            ar, _ = ridge_spearman(
                project_out_basis(Xt[tr], Rk), ytr, project_out_basis(Xt[te], Rk), yte
            )
            wr = sep_direction(project_out_basis(Xc[c_tr], Rk),
                               project_out_basis(Xt[fit_rows], Rk), seed=seed)
            sr = auroc(
                np.concatenate([project_out_basis(Xc[c_te], Rk) @ wr,
                                project_out_basis(Xt[te], Rk) @ wr]), lab)

            rows[k]["atlas"].append(a)
            rows[k]["sep"].append(s)
            rows[k]["atlas_rand"].append(ar)
            rows[k]["sep_rand"].append(sr)
        if verbose:
            print("  seed %d done (%.0fs)" % (seed, time.time() - t0), flush=True)

    base_med = float(np.median(base_per_seed))
    gc_med = float(np.median(gc_per_seed))
    curve = {}
    for k in KS:
        if not rows[k]["atlas"]:
            continue
        curve[k] = {
            "atlas_median": round(float(np.median(rows[k]["atlas"])), 4),
            "residual_separability_median": round(float(np.median(rows[k]["sep"])), 4),
            "atlas_median_random_subspace": round(float(np.median(rows[k]["atlas_rand"])), 4),
            "residual_separability_random_subspace": round(float(np.median(rows[k]["sep_rand"])), 4),
            "atlas_drop_vs_baseline": round(base_med - float(np.median(rows[k]["atlas"])), 4),
            "atlas_drop_attributable": round(
                float(np.median(rows[k]["atlas_rand"])) - float(np.median(rows[k]["atlas"])), 4
            ),
        }

    if verbose:
        print("\n%4s %10s %12s %12s %12s %12s" % (
            "k", "atlas", "sep(resid)", "atlas(rand)", "sep(rand)", "attrib.drop"))
        for k, c in curve.items():
            print("%4d %10.4f %12.4f %12.4f %12.4f %12.4f" % (
                k, c["atlas_median"], c["residual_separability_median"],
                c["atlas_median_random_subspace"],
                c["residual_separability_random_subspace"], c["atlas_drop_attributable"]))

    rec = {
        "schema": "promoteratlas-subspace-removal-v1",
        "date": "2026-11-01",
        "device": "cpu",
        "why": (
            "Experiment A's 1-D removal was underpowered BY CONSTRUCTION: separability was 0.9156 "
            "while removing the single best direction cost -0.0019, which is only consistent with a "
            "detection signal distributed over many of the 128 dimensions. Two defensible "
            "definitions of the direction also disagreed 8x. This is the version that can fail."
        ),
        "reading_rule_fixed_in_advance": (
            "BOTH curves or nothing. Separability -> chance while atlas holds = evidence AGAINST the "
            "exposure explanation. Both fall together = exposure explains the gain. Atlas falls while "
            "separability holds = the projection is destroying capacity and nothing is attributable, "
            "which is what the random-subspace control measures at each k."
        ),
        "baseline": {"atlas_median": round(base_med, 4), "gc_median": round(gc_med, 4),
                     "committed_atlas": committed, "committed_gc": gc_committed},
        "n_tiles": len(tiles),
        "n_controls": len(controls),
        "n_seeds": N_SEEDS,
        "fit_subsample_tiles": FIT_SUBSAMPLE,
        "curve": {str(k): v for k, v in curve.items()},
        "honest_limits": [
            "Still LINEAR. A nonlinear detection signal survives any linear projection, at any k.",
            "Deflation gives ONE ordered basis, not the optimal rank-k discriminative subspace.",
            "Directions are fit on a %d-tile subsample; the ridge uses every tile." % FIT_SUBSAMPLE,
            "1,406 controls against 44,106 tiles, and pos_control is geometry-incompatible.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    print("\n-> %s" % OUT, flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    return run(verbose=not ap.parse_args(argv).quiet)


if __name__ == "__main__":
    sys.exit(main())
