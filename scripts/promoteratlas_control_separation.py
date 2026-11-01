"""Experiment A — does PromoterAtlas's 0.5076 survive removing the promoter-DETECTION direction?

The arm's leading alternative explanation, now confirmed as a live mechanism rather than a
conditional (`wiki/promoteratlas_corpus_provenance_2026-11-01.json`): the pretraining corpus is
built by the authors' own `extract_promoter_regions(..., regulatory_region_length=200)`, i.e. the
200 bp immediately upstream of each annotated gene start. So "resembles my pretraining corpus" IS
"looks like a real upstream regulatory region", and that correlates with measured expression. If the
encoder's contribution to the peak-tile rung is mostly that detection signal, the number is
inflated relative to what it claims to measure (how STRONG a promoter is).

**The substrate already ships the control and every arm has discarded it.** 936 `random` + 470
`neg_control` rows, all exactly 150 bp -- the arm's own primary geometry -- thrown away by the
`category == 'tile'` filter in every consumer of `load_peak_tiles`.

Two cells:

  A1  SEPARABILITY. How well does the 128-dim embedding separate controls from tiles, and how well
      does GC alone? The GC comparison is a PRE-COMMITMENT, not decoration: if the controls are not
      separable by a cheap feature either, then "the encoder cannot separate them" is a statement
      about the controls, not about the encoder, and the cell is INDETERMINATE.

  A2  DIRECTION REMOVAL (the decisive cell). Take the direction that best separates controls from
      tiles, project it out of the TILE embeddings, refit the identical ridge on the identical 10
      committed splits, and see what survives.

**The direction is fit WITHOUT the expression target and WITHOUT the test tiles.** It uses train
tiles + controls only, per seed. The control/tile contrast carries no expression term at all, so
this is not label leakage in either direction; holding out the test tiles anyway keeps the geometry
of the removal honest per fold.

**EXPOSURE_EXCLUDED_LINEAR is bounded BY CONSTRUCTION and may never be written up as "exposure
excluded".** Projection is linear; a nonlinear detection signal survives it. That bound is fixed
here in code and in the frozen bar, before any number is read.
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
BAR = REPO / "wiki/promoteratlas_exposure_acceptance_bar.json"
CACHE = REPO / "data/mpra"
EMB_DIR = REPO / "data/promoteratlas"
COMMITTED_PRED = EMB_DIR / "per_tile_predictions_2026-11-01.json"
OUT = REPO / "wiki/promoteratlas_control_separation_2026-11-01.json"

N_SEEDS = 10
FRAC = 0.25
EXPECTED_TILES = 44106


def auroc(scores, labels) -> float:
    """Rank-based AUROC. Ties get mid-ranks -- this repo has measured that `sorted()`-order
    tie-breaking silently shifts medians, so ties are averaged rather than broken."""
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


def project_out(X, w):
    """Remove the component along `w` from every row of X."""
    w = unit(np.asarray(w, dtype=float))
    return X - np.outer(X @ w, w)


def logistic_direction(Xa, Xb, *, seed=0):
    """Direction of a regularised logistic separator of class b (1) from class a (0)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    X = np.vstack([Xa, Xb])
    y = np.concatenate([np.zeros(len(Xa)), np.ones(len(Xb))])
    sc = StandardScaler().fit(X)
    clf = LogisticRegression(max_iter=2000, C=1.0, random_state=seed).fit(sc.transform(X), y)
    # map the coefficient back to the unscaled space so the projection is in embedding coordinates
    return unit(clf.coef_[0] / np.where(sc.scale_ == 0, 1.0, sc.scale_))


def run(verbose=True):
    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.atlas_embed import load_atlas, masked_mean_pool, ridge_spearman
    from dna_decode.glm.tile_controls import control_counts, load_control_tiles

    t0 = time.time()
    arm = json.loads(ARM.read_text(encoding="utf-8"))
    committed_atlas = float(arm["arms"]["native150"]["median"])
    committed_gc = float(arm["gc_anchor"]["median"])
    bar = json.loads(BAR.read_text(encoding="utf-8")) if BAR.is_file() else None
    if bar is None:
        print("REFUSED: the acceptance bar must be frozen BEFORE this runs (%s absent)" % BAR)
        return 3
    margin = float(bar["inherited_margin"])

    tiles_all, _ = gw.load_peak_tiles(CACHE, allow_download=False)
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    if len(tiles) != EXPECTED_TILES:
        print("REFUSED: %d tiles, expected %d" % (len(tiles), EXPECTED_TILES))
        return 3
    controls = load_control_tiles(CACHE)
    cc = control_counts(controls)
    if verbose:
        print("tiles %d, controls %s (total %d)" % (len(tiles), cc, len(controls)), flush=True)

    # ---- embeddings (tiles reuse the reconcile gate's cached pass) -----------------------------
    tile_emb = EMB_DIR / ("tile_embeddings_native150_n%d.npy" % len(tiles))
    if not tile_emb.is_file():
        print("REFUSED: run scripts/promoteratlas_local_reconcile.py first (%s absent)" % tile_emb)
        return 3
    Xt = np.load(tile_emb)
    model = load_atlas()
    ctrl_cache = EMB_DIR / ("control_embeddings_native150_n%d.npy" % len(controls))
    if ctrl_cache.is_file():
        Xc = np.load(ctrl_cache)
    else:
        Xc = masked_mean_pool(model, [c.seq for c in controls], device="cpu", batch=64)
        np.save(ctrl_cache, Xc)
    gc_t = np.asarray(gw.gc_features([t.seq for t in tiles]), dtype=float)
    gc_c = np.asarray(gw.gc_features([c.seq for c in controls]), dtype=float)
    if verbose:
        print("embeddings: tiles %s controls %s" % (Xt.shape, Xc.shape), flush=True)

    idx = {id(t): i for i, t in enumerate(tiles)}
    rng = np.random.default_rng(12345)
    ctrl_perm = rng.permutation(len(controls))

    per_seed = {}
    for seed in range(N_SEEDS):
        train, test, _ = gw.leave_peak_out(tiles, frac=FRAC, seed=seed)
        tr = np.asarray([idx[id(t)] for t in train])
        te = np.asarray([idx[id(t)] for t in test])
        ytr = [tiles[i].expression for i in tr]
        yte = [tiles[i].expression for i in te]

        # controls split the same way in proportion, so A1 is held out too
        k = int(round(len(controls) * (1 - FRAC)))
        c_tr, c_te = ctrl_perm[:k], ctrl_perm[k:]

        # ---- A1 separability, held out ---------------------------------------------------------
        w_emb = logistic_direction(Xc[c_tr], Xt[tr], seed=seed)
        s_emb = np.concatenate([Xc[c_te] @ w_emb, Xt[te] @ w_emb])
        lab = np.concatenate([np.zeros(len(c_te)), np.ones(len(te))])
        a_emb = auroc(s_emb, lab)

        w_gc = logistic_direction(gc_c[c_tr], gc_t[tr], seed=seed)
        s_gc = np.concatenate([gc_c[c_te] @ w_gc, gc_t[te] @ w_gc])
        a_gc = auroc(s_gc, lab)

        # ---- A2 direction removal --------------------------------------------------------------
        mu = unit(Xt[tr].mean(axis=0) - Xc[c_tr].mean(axis=0))
        base, pred_base = ridge_spearman(Xt[tr], ytr, Xt[te], yte)
        rem_log, _ = ridge_spearman(
            project_out(Xt[tr], w_emb), ytr, project_out(Xt[te], w_emb), yte
        )
        rem_mu, _ = ridge_spearman(project_out(Xt[tr], mu), ytr, project_out(Xt[te], mu), yte)

        # control: remove a RANDOM direction. A drop from removing ANY single direction of 128
        # would make the A2 drop meaningless, so this bounds the removal's own cost.
        r = np.random.default_rng(777 + seed).normal(size=Xt.shape[1])
        rem_rand, _ = ridge_spearman(
            project_out(Xt[tr], r), ytr, project_out(Xt[te], r), yte
        )

        gc_base, _ = ridge_spearman(gc_t[tr], ytr, gc_t[te], yte)

        per_seed[seed] = {
            "n_test_tiles": int(len(te)),
            "n_test_controls": int(len(c_te)),
            "sep_auroc_embedding": round(a_emb, 4),
            "sep_auroc_gc_only": round(a_gc, 4),
            "atlas_baseline": round(base, 4),
            "atlas_minus_logistic_dir": round(rem_log, 4),
            "atlas_minus_meandiff_dir": round(rem_mu, 4),
            "atlas_minus_random_dir": round(rem_rand, 4),
            "gc_baseline": round(gc_base, 4),
            "pred_test": pred_base,
        }
        if verbose:
            print(
                "  seed %d  sep(emb)=%.4f sep(gc)=%.4f | atlas %.4f -> logit %.4f / mean %.4f / rand %.4f | gc %.4f"
                % (seed, a_emb, a_gc, base, rem_log, rem_mu, rem_rand, gc_base),
                flush=True,
            )

    def med(key):
        return float(np.median([per_seed[s][key] for s in range(N_SEEDS)]))

    m = {k: round(med(k), 4) for k in (
        "sep_auroc_embedding", "sep_auroc_gc_only", "atlas_baseline",
        "atlas_minus_logistic_dir", "atlas_minus_meandiff_dir",
        "atlas_minus_random_dir", "gc_baseline")}

    # ---- per-tile agreement with the ARM's own committed predictions ---------------------------
    agree = {}
    if COMMITTED_PRED.is_file():
        cp = json.loads(COMMITTED_PRED.read_text(encoding="utf-8"))
        for seed in range(N_SEEDS):
            e = cp.get(str(seed))
            if not e:
                continue
            mine = per_seed[seed]["pred_test"]
            if len(mine) == len(e["pred"]):
                agree[seed] = round(float(np.corrcoef(mine, e["pred"])[0, 1]), 6)
    for s in per_seed.values():
        s.pop("pred_test", None)

    # ---- GRADE against the frozen bar, validity checks FIRST ----------------------------------
    from dna_decode.glm.exposure_verdict import ExposureInput, grade

    rc = EMB_DIR / "local_reconcile_2026-11-01.json"
    reconciled = False
    reconcile_detail = "reconcile report absent -- the gate was never run"
    if rc.is_file():
        rr = json.loads(rc.read_text(encoding="utf-8"))
        reconciled = rr.get("verdict") == "RECONCILED"
        reconcile_detail = "%s (mismatched=%s)" % (rr.get("verdict"), rr.get("mismatched"))

    xin = ExposureInput(
        harness_reconciled=reconciled,
        sep_auroc_gc_only=m["sep_auroc_gc_only"],
        atlas_baseline=m["atlas_baseline"],
        atlas_minus_separation_dir=m["atlas_minus_logistic_dir"],
        atlas_minus_random_dir=m["atlas_minus_random_dir"],
        gc_baseline=m["gc_baseline"],
        margin=margin,
    )
    v = grade(xin)

    drop_log = m["atlas_baseline"] - m["atlas_minus_logistic_dir"]
    drop_rand = m["atlas_baseline"] - m["atlas_minus_random_dir"]
    surviving_over_gc = m["atlas_minus_logistic_dir"] - m["gc_baseline"]

    rec = {
        "schema": "promoteratlas-control-separation-v1",
        "date": "2026-11-01",
        "device": "cpu",
        "bar": "wiki/promoteratlas_exposure_acceptance_bar.json",
        "inherited_margin": margin,
        "n_tiles": len(tiles),
        "control_counts": cc,
        "n_controls": len(controls),
        "n_seeds": N_SEEDS,
        "committed": {"atlas_native150": committed_atlas, "gc_median": committed_gc},
        "per_seed": per_seed,
        "medians": m,
        "derived": {
            "drop_from_removing_separation_direction": round(drop_log, 4),
            "drop_from_removing_a_RANDOM_direction": round(drop_rand, 4),
            "surviving_margin_over_gc": round(surviving_over_gc, 4),
        },
        "per_tile_agreement_with_committed_arm_predictions": agree,
        "harness_reconciled": reconciled,
        "harness_reconcile_detail": reconcile_detail,
        "verdict": v.verdict,
        "verdict_detail": v.detail,
        "verdict_is_valid": v.is_valid,
        "verdict_check_order": list(__import__("dna_decode.glm.exposure_verdict",
                                               fromlist=["CHECK_ORDER"]).CHECK_ORDER),
        "control_expression": {
            "note": "measured from the substrate; the controls are genuinely low-output",
            "tile_median": 0.4338,
            "random_median": 0.3502,
            "neg_control_median": 0.337,
            "tile_max": 158.308,
            "random_max": 10.747,
            "neg_control_max": 10.052,
        },
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    print("\nmedians: %s" % json.dumps(m, sort_keys=True), flush=True)
    print("derived: %s" % json.dumps(rec["derived"], sort_keys=True), flush=True)
    if agree:
        print("per-tile agreement with the committed arm: median r = %.6f"
              % float(np.median(list(agree.values()))), flush=True)
    print("")
    print("VERDICT: %s" % v.verdict, flush=True)
    print("  %s" % v.detail, flush=True)
    print("-> %s" % OUT, flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    return run(verbose=not a.quiet)


if __name__ == "__main__":
    sys.exit(main())
