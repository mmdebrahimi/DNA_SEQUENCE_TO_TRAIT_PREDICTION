"""Bloom-2013 yeast decoding-validation arm — the project's first CONFOUND-FREE genotype->phenotype test.

Loads the Bloom et al. 2013 BYxRM cross (1,008 segregants x 46 quantitative growth traits x 11,623
markers), recodes the B/R parental-origin genotype to -1/+1, aligns on segregant ID, and runs
cross-validated ridge genomic prediction (dna_decode.eval.genomic_prediction.cv_ridge_gp) per trait with
a label-permutation null. The single-cross design removes population structure by construction, so a
held-out predictive r that beats the null is a GENUINE genotype->phenotype signal.

Data (free, no DUA): Princeton BYxRM web supplement (fetched via Wayback; the sandbox IP is 403-blocked
by Princeton directly). Genotype geno_v2.txt (markers x segregants, B/R); phenotype BYxRM_PhenoData.txt
(segregants x 46 traits, colony size). Marker stride subsamples for a fast pilot (genome-wide markers are
LD-redundant); pass --marker-stride 1 for the full set.

    uv run python scripts/yeast_bloom_gp_arm.py --traits Cadmium_Chloride,Hydrogen_Peroxide,Copper --marker-stride 8
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from dna_decode.eval.genomic_prediction import cv_ridge_gp  # noqa: E402

_GENO = "D:/dna_decode_cache/bloom/geno_v2.txt"
_PHENO = "D:/dna_decode_cache/bloom/BYxRM_PhenoData.txt"
_CODE = {"B": -1.0, "R": 1.0}


def load_genotype(path: str, stride: int):
    """Return (segregant_ids, marker_names, X[n_segregants, n_markers]) recoded B/R -> -1/+1."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        seg_ids = header[1:]                              # cols = segregants
        marker_names, rows = [], []
        for i, line in enumerate(f):
            if i % stride:                               # subsample markers (LD-redundant)
                continue
            parts = line.rstrip("\n").split("\t")
            marker_names.append(parts[0])
            rows.append([_CODE.get(c, np.nan) for c in parts[1:]])
    G = np.asarray(rows, dtype=float).T                  # -> segregants x markers
    # impute any missing with marker mean (recombinant panel -> ~0)
    col_mean = np.nanmean(G, axis=0)
    inds = np.where(np.isnan(G))
    G[inds] = np.take(col_mean, inds[1])
    return seg_ids, marker_names, G


def load_phenotype(path: str):
    with open(path) as f:
        traits = f.readline().rstrip("\n").split("\t")[1:]
        seg_ids, rows = [], []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            seg_ids.append(parts[0])
            rows.append([float(x) if x not in ("", "NA") else np.nan for x in parts[1:]])
    return seg_ids, traits, np.asarray(rows, dtype=float)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Bloom-2013 yeast confound-free genomic-prediction arm.")
    ap.add_argument("--geno", default=_GENO)
    ap.add_argument("--pheno", default=_PHENO)
    ap.add_argument("--traits", default="Cadmium_Chloride,Hydrogen_Peroxide,Copper,Caffeine,"
                    "Lithium_Chloride,Maltose,YPD,Diamide", help="comma list or 'all'")
    ap.add_argument("--marker-stride", type=int, default=8)
    ap.add_argument("--n-perm", type=int, default=30)
    ap.add_argument("--date", default="2026-07-31")
    args = ap.parse_args(argv)

    print("loading genotype ...", flush=True)
    g_ids, markers, G = load_genotype(args.geno, args.marker_stride)
    print(f"  genotype: {G.shape[0]} segregants x {G.shape[1]} markers (stride {args.marker_stride})", flush=True)
    p_ids, traits, P = load_phenotype(args.pheno)
    print(f"  phenotype: {P.shape[0]} segregants x {P.shape[1]} traits", flush=True)

    # align on segregant ID
    g_pos = {sid: i for i, sid in enumerate(g_ids)}
    common = [sid for sid in p_ids if sid in g_pos]
    gi = [g_pos[sid] for sid in common]
    pi = [p_ids.index(sid) for sid in common]
    X = G[gi]
    Ptr = {t: P[pi, traits.index(t)] for t in traits}
    print(f"  aligned {len(common)} segregants", flush=True)

    wanted = traits if args.traits == "all" else [t.strip() for t in args.traits.split(",")]
    wanted = [t for t in wanted if t in traits]
    results = []
    for t in wanted:
        res = cv_ridge_gp(X, Ptr[t], trait=t, n_perm=args.n_perm, seed=0)
        results.append(res.as_dict())
        flag = "BEATS NULL" if res.beats_null else "ns"
        print(f"  {t:20s} predictive_r={res.predictive_r:+.3f} r2={res.predictive_r2:.3f} "
              f"nullp95={res.null_r_p95:.3f}  [{flag}]", flush=True)

    n_pos = sum(1 for r in results if r["beats_null"])
    out = {"substrate": "Bloom_2013_BYxRM_yeast_cross", "n_segregants": len(common),
           "n_markers": X.shape[1], "marker_stride": args.marker_stride, "n_perm": args.n_perm,
           "design": "single-cross segregant panel -> NO population structure (confound-free)",
           "n_traits_scored": len(results), "n_traits_beat_null": n_pos, "traits": results}
    Path(f"wiki/yeast_bloom_gp_arm_{args.date}.json").write_text(json.dumps(out, indent=2))
    print(f"\n{n_pos}/{len(results)} traits beat the label-permutation null "
          f"-> wiki/yeast_bloom_gp_arm_{args.date}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
