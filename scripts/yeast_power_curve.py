"""Architecture-vs-power test #1: yeast power curve. Does the nonlinear (gbm) advantage over TUNED ridge
PERSIST as n shrinks to BXD's ~85? If yes, small-n does NOT prevent capturing epistasis when it exists ->
the mouse-BXD Layer-2 non-replication is ARCHITECTURE (additive traits), not POWER.

Uses TUNED ridge (cv_ridge_gp inner-CV alpha) at every n to avoid the fixed-alpha confound of the earlier
quick kill-test. Multiple random subsamples per n, averaged. Local, free.

    uv run python scripts/yeast_power_curve.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from dna_decode.eval.genomic_prediction import cv_model_gp, cv_ridge_gp  # noqa: E402
from scripts.yeast_bloom_gp_arm import load_genotype, load_phenotype  # noqa: E402


def main(argv=None) -> int:
    g_ids, markers, G = load_genotype("D:/dna_decode_cache/bloom/geno_v2.txt", 8)
    p_ids, traits, P = load_phenotype("D:/dna_decode_cache/bloom/BYxRM_PhenoData.txt")
    gpos = {s: i for i, s in enumerate(g_ids)}
    common = [s for s in p_ids if s in gpos]
    X = G[[gpos[s] for s in common]]

    TRAITS = ["Maltose", "Cadmium_Chloride", "Copper", "Indoleacetic_Acid"]  # gbm winners at full n
    NS = [85, 200, 500, 1008]
    N_SUB = 4
    out = {"question": "does gbm-beats-TUNED-ridge persist as n shrinks to BXD's ~85?",
           "traits": {}, "ns": NS}
    for t in TRAITS:
        y = P[[p_ids.index(s) for s in common], traits.index(t)]
        keep = ~np.isnan(y)
        Xk, yk = X[keep], y[keep]
        per_n = {}
        for n in NS:
            n = min(n, len(yk))
            deltas = []
            for seed in range(N_SUB):
                rng = np.random.default_rng(1000 * seed + n)
                idx = rng.choice(len(yk), n, replace=False)
                rr = cv_ridge_gp(Xk[idx], yk[idx], n_perm=0, seed=seed).predictive_r  # TUNED alpha
                gr = cv_model_gp(Xk[idx], yk[idx], model="gbm", n_perm=0, seed=seed).predictive_r
                deltas.append(gr - rr)
            per_n[n] = round(float(np.mean(deltas)), 4)
            print(f"  {t:18s} n={n:4d}  mean gbm-ridge delta = {per_n[n]:+.4f}", flush=True)
        out["traits"][t] = per_n

    # verdict: is the delta at n~85 still clearly positive (gbm wins) like at full n?
    small = np.mean([v[min(NS)] for v in out["traits"].values()])
    large = np.mean([v[max(k for k in v)] for v in out["traits"].values()])
    out["mean_delta_small_n"] = round(float(small), 4)
    out["mean_delta_large_n"] = round(float(large), 4)
    out["verdict"] = ("POWER_RULED_OUT_architecture" if small > 0.03
                      else ("POWER_CONFOUND" if small < 0.0 else "AMBIGUOUS"))
    Path("wiki/yeast_power_curve_2026-08-02.json").write_text(json.dumps(out, indent=2))
    print(f"\nmean gbm-ridge delta: n~85 = {small:+.4f}  vs  n=full = {large:+.4f}")
    print(f"VERDICT: {out['verdict']}  (gbm still beats ridge at n~85 => power ruled out => architecture)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
