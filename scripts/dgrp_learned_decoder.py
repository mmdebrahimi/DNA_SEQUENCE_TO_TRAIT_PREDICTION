"""DGRP learned-path test — does a multivariate genotype->trait model beat POPULATION STRUCTURE on a
Drosophila quantitative-trait panel, once de-confounded?

This is the animal-panel instance of the project's recurring de-confounded-embedding test (yeast growth,
Arabidopsis flowering, cipro within-lineage, pathotype) — all of which were NEGATIVE under de-confounding.
Substrate: DGRP2 (205 inbred D. melanogaster lines, dm6 SNPs) x a DGRPool quantitative trait (StarvationRes).

Method (mirrors `scripts/yeast_growth_decoder.py`, reusing the promoted `dna_decode.deconfound` primitives):
1. Build a line x SNP dosage matrix from the VCF (MAF + call-rate filtered, thinned).
2. Derive population-structure GROUPS by genotype-distance clustering (the DGRP analog of yeast clades /
   Arabidopsis population structure — no external Wolbachia/inversion covariate needed; structure is in the
   genotype).
3. Compare a NAIVE ridge (`cv_r2`, whole-genome) against the DE-CONFOUNDED `within_group_r2` (5-fold inside
   each clade, scored on clade-centered residuals so between-clade structure cancels). A permutation null
   (within-clade label shuffle) sets the chance floor.

The naive r2 conflates structure + signal; the within-clade r2 is the honest metric. Expected result:
within-clade r2 ~ 0 (the model learned structure, not the causal trait signal) — a 5th de-confounded negative.
Data on D: (gitignored).
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.deconfound import (  # noqa: E402  (promoted de-confounding primitives)
    cluster_from_distance,
    cv_r2,
    permutation_null,
    r2,
    within_group_r2,
)

_LINE_NUM = re.compile(r"(\d+)$")


def _norm(sample: str) -> str:
    m = _LINE_NUM.search(sample)
    return str(int(m.group(1))) if m else sample        # strip zero-padding: DGRP-021 / DGRP_021 -> '21'


def load_phenotype(pheno_json: Path) -> dict:
    """DGRPool phenotype JSON -> {line_number: mean over F+M replicates}."""
    raw = pheno_json.read_bytes()
    d = json.loads(gzip.decompress(raw).decode("utf-8") if raw[:2] == b"\x1f\x8b" else raw.decode("utf-8"))
    out = {}
    for line, sx in d["original_data"].items():
        vals = []
        for s in ("F", "M"):
            if sx.get(s):
                vals += [float(x) for x in sx[s] if x is not None]
        if vals:
            out[_norm(line)] = float(np.mean(vals))
    return out


def _gt_dosage(field: str) -> float:
    """DGRP inbred GT -> dosage: 0/0->0, 1/1->1, het->0.5, missing->nan."""
    gt = field.split(":", 1)[0]
    a = gt.replace("|", "/").split("/")
    if len(a) != 2 or "." in a:
        return np.nan
    return (int(a[0]) + int(a[1])) / 2.0


def build_snp_matrix(vcf: Path, out_npz: Path, keep_every: int = 140, maf_min: float = 0.05,
                     callrate_min: float = 0.8, max_snps: int = 25000) -> dict:
    """Stream the VCF once; genome-wide sample 1-in-`keep_every` SNP LINES (skip the rest WITHOUT parsing
    genotypes — the pure-Python 205-way GT parse is the bottleneck), MAF+call-rate filter the sampled ones,
    save a line x SNP matrix. Skipping before parsing is the k-mer-LOSO perf lesson applied to VCF streaming."""
    samples = None
    cols = []
    snp_ids = []
    line_i = 0
    sampled = 0
    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line[0] == "#":
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                continue
            line_i += 1
            if line_i % keep_every != 0:                        # skip WITHOUT the 205-way genotype parse
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            sampled += 1
            dos = np.array([_gt_dosage(g) for g in parts[9:]], dtype=float)
            ok = ~np.isnan(dos)
            if ok.mean() < callrate_min:
                continue
            af = np.nanmean(dos[ok]) if ok.any() else 0.0
            if min(af, 1 - af) < maf_min:
                continue
            cols.append(dos)
            snp_ids.append(f"{parts[0]}:{parts[1]}")
            if len(cols) >= max_snps:
                break
    X = np.array(cols, dtype="float32").T                       # lines x SNPs
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_npz, X=X, samples=np.array(samples), snp_ids=np.array(snp_ids))
    return {"n_lines": len(samples), "n_snps_sampled": sampled, "n_snps_kept": X.shape[1]}


def _impute(X: np.ndarray) -> np.ndarray:
    Xi = X.copy()
    col_mean = np.nanmean(Xi, axis=0)
    idx = np.where(np.isnan(Xi))
    Xi[idx] = np.take(col_mean, idx[1])
    return np.nan_to_num(Xi, nan=0.0)


def run(cache_npz: Path, pheno_json: Path, k_clades: int = 6, seed: int = 0) -> dict:
    data = np.load(cache_npz, allow_pickle=True)
    X = _impute(data["X"].astype("float32"))
    samples = list(data["samples"])
    ph = load_phenotype(pheno_json)
    # align lines present in both
    keep = [(i, _norm(s)) for i, s in enumerate(samples) if _norm(s) in ph]
    idx = [i for i, _ in keep]
    y = np.array([ph[n] for _, n in keep], dtype=float)
    Xa = X[idx]
    # DISCRETE-structure diagnostic: euclidean avg-linkage (the yeast/Arabidopsis method). On DGRP (a single
    # Raleigh-derived population) this collapses to ~1 giant clade -> the clade-de-confounding test is
    # DEGENERATE by the substrate's nature, unlike yeast (discrete clades).
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist, squareform
    Deu = squareform(pdist(Xa, metric="euclidean"))
    eu_clades = cluster_from_distance(Deu, k_clades)
    eu_sizes = sorted(np.bincount(eu_clades)[1:], reverse=True)
    discrete_degenerate = bool(eu_sizes and eu_sizes[0] >= 0.9 * len(y))
    # CONTINUOUS-structure de-confounding: cluster on the top genotype PCs (ward) -> balanced groups so the
    # within-group residual test is actually computed across >=3 strata.
    Xc = Xa - Xa.mean(0)
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    pcs = U[:, :5] * S[:5]
    clades = fcluster(linkage(pcs, method="ward"), t=k_clades, criterion="maxclust")
    # naive whole-genome ridge vs de-confounded within-clade
    naive = cv_r2(Xa, y, seed=seed)
    within, n_used = within_group_r2(Xa, y, clades, min_n=15)
    # permutation null on the de-confounded metric (within-clade label shuffle)
    # residualize y by clade to feed permutation_null (matches within_group scoring space)
    yr = y.copy()
    for g in np.unique(clades):
        m = clades == g
        if m.sum() > 1:
            yr[m] = y[m] - np.nanmean(y[m])
    perm = permutation_null(y, yr, clades, n=200)
    perm_p95 = float(np.nanpercentile(np.abs(perm), 95))
    signal = within > 0.05 and within > perm_p95 ** 2
    return {"n_lines": len(y), "n_snps": int(Xa.shape[1]), "k_clades": int(k_clades),
            "structure_method": "pca_ward_5pc",
            "n_clades_used": int(n_used),
            "discrete_structure_degenerate": discrete_degenerate,
            "eu_top_clade_size": int(eu_sizes[0]) if eu_sizes else 0,
            "naive_cv_r2": round(float(naive), 4),
            "within_clade_r2": round(float(within), 4),
            "perm_null_within_p95_abs_corr": round(perm_p95, 4),
            "verdict": ("WITHIN_CLADE_SIGNAL_PRESENT" if signal
                        else "LEARNED_PATH_NEGATIVE_UNDER_DECONFOUNDING")}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", type=Path, default=Path("D:/dna_decode_cache/dgrp/dgrp2_ncsu_snps.vcf.gz"))
    ap.add_argument("--cache", type=Path, default=Path("D:/dna_decode_cache/dgrp/dgrp_snp_matrix.npz"))
    ap.add_argument("--pheno", type=Path, default=Path("D:/dna_decode_cache/dgrp/pheno_2798.json"))
    ap.add_argument("--build", action="store_true", help="build the SNP matrix cache from the VCF first")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "dgrp_learned_scores.json")
    a = ap.parse_args(argv)
    if a.build or not a.cache.exists():
        print("building SNP matrix from VCF (streaming)...", flush=True)
        info = build_snp_matrix(a.vcf, a.cache)
        print("built:", info, flush=True)
    res = run(a.cache, a.pheno)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
