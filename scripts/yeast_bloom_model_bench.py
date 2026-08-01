"""Layer-2 B1: does a NONLINEAR model (gradient-boosted trees) beat linear ridge on the confound-free
Bloom yeast substrate? Tests the "optimize the model / capture epistasis" hypothesis — Bloom 2013 found
gene-gene interactions contribute up to ~50% of heritability for some traits, which linear ridge misses.

Same segregants + markers as the Layer-1 ridge arm; per trait reports ridge r vs gbm r + the delta.

    uv run python scripts/yeast_bloom_model_bench.py --marker-stride 8
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from dna_decode.eval.genomic_prediction import cv_model_gp, cv_ridge_gp  # noqa: E402
from scripts.yeast_bloom_gp_arm import load_genotype, load_phenotype  # noqa: E402

_TRAITS = ("Cadmium_Chloride,Hydrogen_Peroxide,Copper,Caffeine,Lithium_Chloride,Maltose,YPD,Diamide,"
           "Cobalt_Chloride,Cycloheximide,Neomycin,Zeocin")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Bloom Layer-2 B1: ridge vs gradient-boosted trees.")
    ap.add_argument("--geno", default="D:/dna_decode_cache/bloom/geno_v2.txt")
    ap.add_argument("--pheno", default="D:/dna_decode_cache/bloom/BYxRM_PhenoData.txt")
    ap.add_argument("--traits", default=_TRAITS)
    ap.add_argument("--marker-stride", type=int, default=8)
    ap.add_argument("--date", default="2026-07-31")
    args = ap.parse_args(argv)

    g_ids, markers, G = load_genotype(args.geno, args.marker_stride)
    p_ids, traits, P = load_phenotype(args.pheno)
    g_pos = {sid: i for i, sid in enumerate(g_ids)}
    common = [sid for sid in p_ids if sid in g_pos]
    X = G[[g_pos[s] for s in common]]
    pidx = {s: p_ids.index(s) for s in common}
    print(f"aligned {len(common)} segregants x {X.shape[1]} markers", flush=True)

    wanted = [t.strip() for t in args.traits.split(",") if t.strip() in traits]
    rows = []
    for t in wanted:
        y = P[[pidx[s] for s in common], traits.index(t)]
        rid = cv_ridge_gp(X, y, trait=t, n_perm=20, seed=0)
        gbm = cv_model_gp(X, y, model="gbm", trait=t, n_perm=0, seed=0)
        delta = gbm.predictive_r - rid.predictive_r
        rows.append({"trait": t, "ridge_r": round(rid.predictive_r, 4), "gbm_r": round(gbm.predictive_r, 4),
                     "delta_gbm_minus_ridge": round(delta, 4), "ridge_beats_null": rid.beats_null})
        print(f"  {t:20s} ridge_r={rid.predictive_r:.3f}  gbm_r={gbm.predictive_r:.3f}  "
              f"delta={delta:+.3f}  {'GBM WINS' if delta > 0.01 else ('RIDGE WINS' if delta < -0.01 else 'tie')}",
              flush=True)

    deltas = [r["delta_gbm_minus_ridge"] for r in rows]
    n_gbm_win = sum(1 for d in deltas if d > 0.01)
    n_ridge_win = sum(1 for d in deltas if d < -0.01)
    out = {"substrate": "Bloom_2013_BYxRM_yeast_cross", "n_segregants": len(common),
           "n_markers": X.shape[1], "question": "does nonlinear (gbm) beat linear ridge? = optimize-the-model",
           "n_traits": len(rows), "gbm_wins": n_gbm_win, "ridge_wins": n_ridge_win,
           "mean_delta_gbm_minus_ridge": round(float(np.mean(deltas)), 4),
           "median_delta": round(float(np.median(deltas)), 4), "traits": rows}
    Path(f"wiki/yeast_bloom_model_bench_{args.date}.json").write_text(json.dumps(out, indent=2))
    print(f"\ngbm wins {n_gbm_win}/{len(rows)}, ridge wins {n_ridge_win}/{len(rows)}, "
          f"mean delta {np.mean(deltas):+.4f}  ->  wiki/yeast_bloom_model_bench_{args.date}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
