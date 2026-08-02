"""Architecture-vs-power test #2: a DENSE PLANT cross — Arabidopsis MAGIC (Gnan et al. 2014, 19-parent
MAGIC; rqtl/qtl2data mirror). 677 lines x 1260 markers x 8 quantitative traits — ~3x denser per-trait than
mouse BXD. If nonlinear (gbm) does NOT beat linear ridge even at this WELL-POWERED n, the Layer-2
non-replication is ARCHITECTURE (additive traits), not power.

Same pipeline as the yeast + BXD arms. Genotype A/B/H (parent-origin collapsed) -> -1/+1/0.

    uv run python scripts/arabmagic_gp_arm.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from dna_decode.eval.genomic_prediction import cv_model_gp, cv_ridge_gp  # noqa: E402

_GENO = "D:/dna_decode_cache/arabmagic/arabmagic_geno.csv"
_PHENO = "D:/dna_decode_cache/arabmagic/arabmagic_pheno.csv"
_CODE = {"A": -1.0, "B": 1.0, "H": 0.0}


def _rows(path):
    return [r for r in csv.reader(open(path)) if r and not r[0].startswith("#")]


def load_genotype(path, stride=1):
    rows = _rows(path)
    lines = rows[0][1:]
    markers, calls = [], []
    for i, r in enumerate(rows[1:]):
        if i % stride:
            continue
        markers.append(r[0])
        calls.append([_CODE.get(c, np.nan) for c in r[1:]])
    G = np.asarray(calls, dtype=float).T
    cm = np.nanmean(G, axis=0)
    inds = np.where(np.isnan(G))
    G[inds] = np.take(cm, inds[1])
    return lines, markers, G


def load_phenotype(path):
    rows = _rows(path)
    traits = rows[0][1:]
    ids, vals = [], []
    for r in rows[1:]:
        ids.append(r[0])
        vals.append([float(x) if x not in ("", "NA", "-") else np.nan for x in r[1:]])
    return ids, traits, np.asarray(vals, dtype=float)


def main(argv=None) -> int:
    g_ids, markers, G = load_genotype(_GENO)
    p_ids, traits, P = load_phenotype(_PHENO)
    gpos = {s: i for i, s in enumerate(g_ids)}
    common = [s for s in p_ids if s in gpos]
    X = G[[gpos[s] for s in common]]
    P_c = P[[p_ids.index(s) for s in common]]
    print(f"aligned {len(common)} MAGIC lines x {X.shape[1]} markers; {len(traits)} traits", flush=True)

    rows_out = []
    for j, t in enumerate(traits):
        y = P_c[:, j]
        if (~np.isnan(y)).sum() < 300:
            continue
        rid = cv_ridge_gp(X, y, trait=t, n_perm=30, seed=0)
        gbm = cv_model_gp(X, y, model="gbm", trait=t, n_perm=0, seed=0)
        d = gbm.predictive_r - rid.predictive_r
        rows_out.append({"trait": t, "n": rid.n, "ridge_r": round(rid.predictive_r, 4),
                         "ridge_beats_null": rid.beats_null, "gbm_r": round(gbm.predictive_r, 4),
                         "delta_gbm_minus_ridge": round(d, 4)})
        print(f"  {t:18s} n={rid.n:3d} ridge_r={rid.predictive_r:+.3f} "
              f"{'BEATS-NULL' if rid.beats_null else 'ns':10s} gbm_r={gbm.predictive_r:+.3f} delta={d:+.3f}",
              flush=True)

    n_decode = sum(1 for r in rows_out if r["ridge_beats_null"])
    n_gbm = sum(1 for r in rows_out if r["delta_gbm_minus_ridge"] > 0.01)
    mean_d = float(np.mean([r["delta_gbm_minus_ridge"] for r in rows_out]))
    out = {"organism": "Arabidopsis thaliana (19-parent MAGIC)", "design": "MAGIC RIL -> confound-free",
           "n_markers": X.shape[1], "n_traits": len(rows_out),
           "layer1_decode_beat_null": n_decode, "layer2_gbm_wins": n_gbm,
           "mean_delta_gbm_minus_ridge": round(mean_d, 4), "traits": rows_out}
    Path("wiki/arabmagic_gp_arm_2026-08-02.json").write_text(json.dumps(out, indent=2))
    print(f"\nLayer1: {n_decode}/{len(rows_out)} decode. Layer2: gbm beats ridge {n_gbm}/{len(rows_out)} "
          f"(mean delta {mean_d:+.3f})  n~{rows_out[0]['n'] if rows_out else 0} -> wiki/arabmagic_gp_arm_2026-08-02.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
