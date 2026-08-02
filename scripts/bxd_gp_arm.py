"""Second confound-free cross — mouse BXD RIL panel (C57BL/6J x DBA/2J). Generality test: does the yeast
finding (genomic prediction decodes quantitative traits + nonlinear beats linear) hold in a MAMMAL?

Data (free, no DUA): R/qtl2 mirror rqtl/qtl2data/BXD -- bxd_geno.csv (7320 markers x 198 strains, B/D/H
parent-origin) + bxd_pheno.csv (198 strains x 5806 traits, sparse). RIL from a two-parent cross ->
structure-free by design (same confound-free property as the Bloom yeast segregants).

Layer 1 (cv_ridge_gp + permutation null): does genomic prediction decode BXD traits?
Layer 2 (cv_model_gp gbm vs ridge): does the nonlinear/epistasis advantage replicate?

    uv run python scripts/bxd_gp_arm.py --marker-stride 3 --n-traits 12
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from dna_decode.eval.genomic_prediction import cv_model_gp, cv_ridge_gp  # noqa: E402

_GENO = "D:/dna_decode_cache/bxd/bxd_geno.csv"
_PHENO = "D:/dna_decode_cache/bxd/bxd_pheno.csv"
_COVAR = "D:/dna_decode_cache/bxd/bxd_phenocovar.csv"
_CODE = {"B": -1.0, "D": 1.0, "H": 0.0}
_EXCLUDE = ("pigment", "coat color", "cofactor", "statistics", "epoch", "phase of produ", "breeding")


def _rows(path):
    return [r for r in csv.reader(open(path)) if r and not r[0].startswith("#")]


def load_genotype(path, stride):
    rows = _rows(path)
    strains = rows[0][1:]
    markers, calls = [], []
    for i, r in enumerate(rows[1:]):
        if i % stride:
            continue
        markers.append(r[0])
        calls.append([_CODE.get(c, np.nan) for c in r[1:]])
    G = np.asarray(calls, dtype=float).T                     # strains x markers
    col_mean = np.nanmean(G, axis=0)
    inds = np.where(np.isnan(G))
    G[inds] = np.take(col_mean, inds[1])
    return strains, markers, G


def load_phenotype(path):
    rows = _rows(path)
    traits = rows[0][1:]
    strains, vals = [], []
    for r in rows[1:]:
        strains.append(r[0])
        vals.append([float(x) if x not in ("", "NA", "-") else np.nan for x in r[1:]])
    return strains, traits, np.asarray(vals, dtype=float)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="BXD mouse confound-free cross genomic-prediction arm.")
    ap.add_argument("--marker-stride", type=int, default=3)
    ap.add_argument("--n-traits", type=int, default=12)
    ap.add_argument("--min-cov", type=int, default=78)
    ap.add_argument("--date", default="2026-08-02")
    args = ap.parse_args(argv)

    g_ids, markers, G = load_genotype(_GENO, args.marker_stride)
    p_ids, traits, P = load_phenotype(_PHENO)
    desc = {r[0]: r[-1] for r in _rows(_COVAR)}
    print(f"genotype {G.shape[0]} strains x {G.shape[1]} markers; phenotype {P.shape[0]} strains x "
          f"{P.shape[1]} traits", flush=True)

    g_pos = {s: i for i, s in enumerate(g_ids)}
    common = [s for s in p_ids if s in g_pos]
    X_all = G[[g_pos[s] for s in common]]
    P_common = P[[p_ids.index(s) for s in common]]
    cov = (~np.isnan(P_common)).sum(0)

    # auto-select distinct quantitative high-coverage traits (dedupe by description prefix)
    order = np.argsort(-cov)
    picked, seen_prefix = [], set()
    for j in order:
        if cov[j] < args.min_cov:
            break
        d = desc.get(traits[j], "")
        if any(k in d.lower() for k in _EXCLUDE):
            continue
        vals = P_common[~np.isnan(P_common[:, j]), j]
        if len(np.unique(vals)) < 10:
            continue
        prefix = d.split(";")[-1].strip()[:22].lower()   # organ/measure prefix
        if prefix in seen_prefix:
            continue
        seen_prefix.add(prefix)
        picked.append(j)
        if len(picked) >= args.n_traits:
            break

    rows_out = []
    for j in picked:
        y = P_common[:, j]
        rid = cv_ridge_gp(X_all, y, trait=traits[j], n_perm=30, seed=0)
        gbm = cv_model_gp(X_all, y, model="gbm", trait=traits[j], n_perm=0, seed=0)
        delta = gbm.predictive_r - rid.predictive_r
        rows_out.append({"trait_id": traits[j], "desc": desc.get(traits[j], "")[:70], "n": rid.n,
                         "ridge_r": round(rid.predictive_r, 4), "ridge_beats_null": rid.beats_null,
                         "null_p95": round(rid.null_r_p95, 4), "gbm_r": round(gbm.predictive_r, 4),
                         "delta_gbm_minus_ridge": round(delta, 4)})
        print(f"  {traits[j]:6s} n={rid.n:3d} ridge_r={rid.predictive_r:+.3f} "
              f"{'BEATS-NULL' if rid.beats_null else 'ns':10s} gbm_r={gbm.predictive_r:+.3f} "
              f"delta={delta:+.3f}  {desc.get(traits[j],'')[:40]}", flush=True)

    n_decode = sum(1 for r in rows_out if r["ridge_beats_null"])
    n_gbm_win = sum(1 for r in rows_out if r["delta_gbm_minus_ridge"] > 0.01)
    out = {"organism": "Mus musculus (BXD RIL, C57BL/6J x DBA/2J)", "substrate": "BXD_rqtl2",
           "design": "recombinant-inbred cross -> confound-free", "n_markers": X_all.shape[1],
           "n_traits": len(rows_out),
           "layer1_traits_decode_beat_null": n_decode, "layer2_gbm_wins": n_gbm_win,
           "mean_delta_gbm_minus_ridge": round(float(np.mean([r["delta_gbm_minus_ridge"] for r in rows_out])), 4),
           "traits": rows_out}
    Path(f"wiki/bxd_gp_arm_{args.date}.json").write_text(json.dumps(out, indent=2))
    print(f"\nLayer1: {n_decode}/{len(rows_out)} traits decode (beat null). "
          f"Layer2: gbm beats ridge {n_gbm_win}/{len(rows_out)} (mean delta {out['mean_delta_gbm_minus_ridge']:+.3f}) "
          f"-> wiki/bxd_gp_arm_{args.date}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
