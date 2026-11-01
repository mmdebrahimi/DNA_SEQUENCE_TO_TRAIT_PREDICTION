"""How much detection signal is NONLINEAR? Bound the one caveat the linear experiments cannot touch.

Every intervention in `wiki/promoteratlas_exposure_result_2026-11-01.md` was LINEAR, so the result is
pre-committed to the bounded form -- it may never be written as "exposure excluded", because a
nonlinear detection signal survives a linear projection. That caveat is correct and it is also, as
stated, unquantified: it names a loophole without saying how big it is.

This measures the loophole instead of gesturing at it. It does NOT try to project out a nonlinear
direction (that is not a well-defined operation); it asks the prior question:

    is the control-vs-tile detection signal LINEARLY encoded in this embedding at all?

Compare, on identical held-out folds, how well controls separate from tiles under
  * a LINEAR separator (logistic) -- what the removal experiments were able to delete, AUROC 0.9156
  * a NONLINEAR separator (a small MLP) -- the ceiling of what any detection probe could find here

**The reading rule, fixed before the run.** If the nonlinear probe is barely better than the linear
one, then almost all the detectable signal IS linear, the linear subspace removal deleted almost all
of it, and the "nonlinear survives" caveat is TIGHT -- small, bounded, and quantified. If the
nonlinear probe is substantially better, then real detection headroom was left behind by construction
and the caveat is LOOSE, which would materially weaken the main result.

**This is a bound on the PROBE, not a removal.** Even a tight result does not license "exposure
excluded": it says the signal a linear projection could not reach is small on this control cohort and
this embedding, which is narrower than "no nonlinear exposure exists". Not pre-registered in the
frozen bar; corroboration only, and it cannot change the graded verdict label.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "data/mpra"
EMB_DIR = REPO / "data/promoteratlas"
OUT = REPO / "wiki/promoteratlas_nonlinear_bound_2026-11-01.json"

N_SEEDS = 10
FRAC = 0.25
FIT_SUBSAMPLE = 8000

#: Fixed before the run. The linear figure this is compared against is the measured 0.9156.
TIGHT_IF_GAIN_BELOW = 0.02


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


def run(verbose=True):
    from sklearn.linear_model import LogisticRegression
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler

    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.tile_controls import load_control_tiles

    t0 = time.time()
    tiles_all, _ = gw.load_peak_tiles(CACHE, allow_download=False)
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    controls = load_control_tiles(CACHE)
    Xt = np.load(EMB_DIR / ("tile_embeddings_native150_n%d.npy" % len(tiles)))
    Xc = np.load(EMB_DIR / ("control_embeddings_native150_n%d.npy" % len(controls)))
    gc_t = np.asarray(gw.gc_features([t.seq for t in tiles]), dtype=float)
    gc_c = np.asarray(gw.gc_features([c.seq for c in controls]), dtype=float)
    if verbose:
        print("tiles %s controls %s" % (Xt.shape, Xc.shape), flush=True)

    idx = {id(t): i for i, t in enumerate(tiles)}
    cperm = np.random.default_rng(12345).permutation(len(controls))
    lin, nonlin, lin_gc, nonlin_gc = [], [], [], []

    for seed in range(N_SEEDS):
        train, test, _ = gw.leave_peak_out(tiles, frac=FRAC, seed=seed)
        tr = np.asarray([idx[id(t)] for t in train])
        te = np.asarray([idx[id(t)] for t in test])
        kc = int(round(len(controls) * (1 - FRAC)))
        c_tr, c_te = cperm[:kc], cperm[kc:]
        rng = np.random.default_rng(31337 + seed)
        fit_rows = tr if len(tr) <= FIT_SUBSAMPLE else rng.choice(tr, FIT_SUBSAMPLE, replace=False)

        lab_te = np.concatenate([np.zeros(len(c_te)), np.ones(len(te))])

        def pair(A_tr, B_tr, A_te, B_te):
            Xtr = np.vstack([A_tr, B_tr])
            ytr = np.concatenate([np.zeros(len(A_tr)), np.ones(len(B_tr))])
            Xte = np.vstack([A_te, B_te])
            sc = StandardScaler().fit(Xtr)
            return sc.transform(Xtr), ytr, sc.transform(Xte)

        # ---- embedding ------------------------------------------------------------------------
        Xtr, ytr, Xte = pair(Xc[c_tr], Xt[fit_rows], Xc[c_te], Xt[te])
        l = LogisticRegression(max_iter=2000, C=1.0, random_state=seed).fit(Xtr, ytr)
        lin.append(auroc(l.decision_function(Xte), lab_te))
        mlp = MLPClassifier(hidden_layer_sizes=(64,), max_iter=400, early_stopping=True,
                            random_state=seed).fit(Xtr, ytr)
        nonlin.append(auroc(mlp.predict_proba(Xte)[:, 1], lab_te))

        # ---- GC, as the cheap-feature reference for the same comparison -----------------------
        Gtr, gytr, Gte = pair(gc_c[c_tr], gc_t[fit_rows], gc_c[c_te], gc_t[te])
        lg = LogisticRegression(max_iter=2000, random_state=seed).fit(Gtr, gytr)
        lin_gc.append(auroc(lg.decision_function(Gte), lab_te))
        mg = MLPClassifier(hidden_layer_sizes=(64,), max_iter=400, early_stopping=True,
                           random_state=seed).fit(Gtr, gytr)
        nonlin_gc.append(auroc(mg.predict_proba(Gte)[:, 1], lab_te))

        if verbose:
            print("  seed %d  emb: linear %.4f -> mlp %.4f (gain %+.4f) | gc: %.4f -> %.4f (%.0fs)"
                  % (seed, lin[-1], nonlin[-1], nonlin[-1] - lin[-1], lin_gc[-1], nonlin_gc[-1],
                     time.time() - t0), flush=True)

    m_lin = float(np.median(lin))
    m_non = float(np.median(nonlin))
    gain = m_non - m_lin
    m_lin_gc = float(np.median(lin_gc))
    m_non_gc = float(np.median(nonlin_gc))

    tight = gain < TIGHT_IF_GAIN_BELOW
    rec = {
        "schema": "promoteratlas-nonlinear-bound-v1",
        "date": "2026-11-01",
        "device": "cpu",
        "question": (
            "Is the control-vs-tile detection signal in this embedding LINEARLY encoded? That bounds "
            "the one caveat the linear removal experiments structurally cannot touch."
        ),
        "reading_rule_fixed_in_advance": (
            "nonlinear gain over linear < %.2f => almost all detectable signal is LINEAR, so the "
            "linear subspace removal deleted almost all of it and the 'nonlinear survives' caveat is "
            "TIGHT. A larger gain => real detection headroom was left behind and the caveat is LOOSE."
            % TIGHT_IF_GAIN_BELOW
        ),
        "n_tiles": len(tiles),
        "n_controls": len(controls),
        "n_seeds": N_SEEDS,
        "fit_subsample_tiles": FIT_SUBSAMPLE,
        "embedding": {
            "linear_auroc_median": round(m_lin, 4),
            "nonlinear_auroc_median": round(m_non, 4),
            "nonlinear_gain": round(gain, 4),
            "per_seed_linear": [round(v, 4) for v in lin],
            "per_seed_nonlinear": [round(v, 4) for v in nonlin],
        },
        "gc_reference": {
            "linear_auroc_median": round(m_lin_gc, 4),
            "nonlinear_auroc_median": round(m_non_gc, 4),
            "nonlinear_gain": round(m_non_gc - m_lin_gc, 4),
            "why": (
                "GC is ONE feature, so a nonlinear probe on it can only exploit a non-monotone "
                "relationship. It is the reference for how much a nonlinear probe buys when there is "
                "almost no structure to find, which keeps the embedding gain in scale."
            ),
        },
        "verdict": "LINEAR_CAVEAT_TIGHT" if tight else "LINEAR_CAVEAT_LOOSE",
        "what_this_does_NOT_license": [
            "It is NOT 'exposure excluded'. A tight bound says the signal a linear projection could "
            "not reach is SMALL on this control cohort and this embedding -- narrower than claiming "
            "no nonlinear exposure exists.",
            "It is a bound on a PROBE, not a removal: nothing nonlinear was projected out.",
            "One MLP architecture (64 hidden units, early stopping). A larger probe might find more, "
            "though it would also overfit 1,406 controls faster.",
            "NOT pre-registered in the frozen bar; corroboration only, and it cannot change the "
            "graded verdict label.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    print("\nembedding: linear %.4f -> nonlinear %.4f  gain %+.4f" % (m_lin, m_non, gain), flush=True)
    print("gc ref   : linear %.4f -> nonlinear %.4f  gain %+.4f"
          % (m_lin_gc, m_non_gc, m_non_gc - m_lin_gc), flush=True)
    print("VERDICT  : %s" % rec["verdict"], flush=True)
    print("-> %s" % OUT, flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    return run(verbose=not ap.parse_args(argv).quiet)


if __name__ == "__main__":
    sys.exit(main())
