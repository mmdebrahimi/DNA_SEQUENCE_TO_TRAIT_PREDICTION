"""Stage-1b: the POLYGENIC (elastic-net multi-SNP cis) ceiling + the confound falsification.

Tests the row-574 claim that pooled-vs-within inflation is a POLYGENIC phenomenon (not a
single-causal-SNP one): fit elastic-net over the cis-window (PrediXcan-style, the field's
real linear baseline) with 5-fold CV -> OUT-OF-SAMPLE predicted expression per individual,
then pooled vs WITHIN-population. If the polygenic predictor inflates pooled>>within where
the single-SNP did not -> claim SUPPORTED (the confound needs many features to appear).
If it also shows ~0 inflation -> claim FALSIFIED (cis-linear is de-confounded regardless).

Streams plain-gzip VCF once, cis-window dosages (DS field). numpy/scipy/sklearn.
"""
from __future__ import annotations

import argparse, gzip, json
from pathlib import Path
import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import KFold

from dna_decode.organism_multimodal.geuvadis_data import load_expr, parse_sample_population, canon_pop
from dna_decode.organism_multimodal.deconfound_eval import evaluate

CACHE = Path("D:/dna_decode_cache/geuvadis")
GT = CACHE / "genotypes"
VCF = "GEUVADIS.chr{c}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"


def target_genes(chrom, gene_row):
    """eQTL genes on `chrom` (reuse EUR373 best) -> {gene: coord}."""
    out = {}
    with gzip.open(CACHE / "analysis_results/EUR373.gene.cis.FDR5.best.rs137.txt.gz", "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 11 and p[4] == chrom and p[2] in gene_row:
                try: out[p[2]] = int(p[7])       # gene start coord
                except ValueError: pass
    return out


def collect_cis(vcf, genes_coord, half=250_000, maf_min=0.05, cap=400):
    """Stream VCF once; per gene collect up to `cap` cis-window SNP dosage vectors."""
    lo = {g: c - half for g, c in genes_coord.items()}
    hi = {g: c + half for g, c in genes_coord.items()}
    gene_list = sorted(genes_coord, key=lambda g: genes_coord[g])
    coords = np.array([genes_coord[g] for g in gene_list])
    gmin, gmax = min(lo.values()), max(hi.values())        # VCF is position-sorted -> break past gmax
    samples = []
    per_gene = {g: [] for g in genes_coord}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]; continue
            t1 = line.find("\t"); t2 = line.find("\t", t1 + 1)
            pos = int(line[t1 + 1:t2])
            if pos < gmin: continue
            if pos > gmax: break               # sorted VCF: no later variant is in any window
            # genes whose window covers pos: coord in [pos-half, pos+half]
            idx = np.searchsorted(coords, [pos - half, pos + half])
            if idx[0] == idx[1]: continue
            cand = [gene_list[i] for i in range(idx[0], idx[1]) if len(per_gene[gene_list[i]]) < cap]
            if not cand: continue
            p = line.rstrip("\n").split("\t")
            fmt = p[8].split(":")
            di = fmt.index("DS") if "DS" in fmt else None
            vals = np.full(len(samples), np.nan)
            if di is not None:
                for i, s in enumerate(p[9:]):
                    f = s.split(":")
                    if di < len(f) and f[di] not in (".", ""):
                        try: vals[i] = float(f[di])
                        except ValueError: pass
            if np.isnan(vals).mean() > 0.1: continue
            m = np.nanmean(vals) / 2.0
            if min(m, 1 - m) < maf_min: continue
            vals = np.where(np.isnan(vals), np.nanmean(vals), vals)
            for g in cand: per_gene[g].append(vals)
    return samples, per_gene


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--chrom", default="1")
    ap.add_argument("--out", default="wiki/geuvadis_stage1b_elasticnet_2026-07-27")
    ap.add_argument("--max-genes", type=int, default=40)
    ap.add_argument("--half", type=int, default=100_000)
    a = ap.parse_args()
    E = load_expr(CACHE / "GD462.GeneQuantRPKM.txt.gz")
    S = parse_sample_population(CACHE / "E-GEUV-1.sdrf.txt")
    expr_pop = {s: canon_pop(S[s]) for s in E.samples if s in S}
    gene_row = {g: i for i, g in enumerate(E.genes)}
    sidx = {s: i for i, s in enumerate(E.samples)}
    gc_all = target_genes(a.chrom, gene_row)
    gc = dict(sorted(gc_all.items(), key=lambda kv: kv[1])[:a.max_genes])   # first N by coord (tighter window union -> fast stream)
    print(f"chr{a.chrom}: {len(gc)}/{len(gc_all)} eQTL genes (window +-{a.half//1000}kb); streaming cis-windows ...")
    vcf_samples, per_gene = collect_cis(GT / VCF.format(c=a.chrom), gc, half=a.half)
    vcol = {s: i for i, s in enumerate(vcf_samples)}
    common = [s for s in E.samples if s in vcol and s in expr_pop]
    pops = np.array([expr_pop[s] for s in common])
    vpos = [vcol[s] for s in common]

    rows = []
    kf = KFold(n_splits=5, shuffle=True, random_state=0)
    for g, snps in per_gene.items():
        if len(snps) < 3: continue
        X = np.array(snps).T[vpos]                       # individuals x SNPs
        y = np.array([float(E.values[gene_row[g]][sidx[s]]) for s in common])
        if np.std(y) == 0: continue
        pred = np.full(len(y), np.nan)
        try:
            for tr, te in kf.split(X):
                if np.std(y[tr]) == 0: continue
                m = ElasticNetCV(l1_ratio=0.5, n_alphas=8, cv=3, max_iter=1500, n_jobs=1)
                m.fit(X[tr], y[tr]); pred[te] = m.predict(X[te])
        except Exception:
            continue
        ok = ~np.isnan(pred)
        if ok.sum() < 50 or np.std(pred[ok]) == 0: continue
        res = evaluate(y[ok], pred[ok], pops[ok], n_perm=100)
        rows.append({"gene": g, "n_snps": len(snps), "n": int(ok.sum()),
                     "pooled": res.pooled_rho, "within": res.within_rho_mean, "inflation": res.inflation})
    if not rows:
        print("NO GENES"); return 1
    pooled = np.array([abs(r["pooled"]) for r in rows]); within = np.array([abs(r["within"]) for r in rows])
    summ = {"n_genes": len(rows), "chrom": a.chrom, "predictor": "elastic-net multi-SNP cis (5-fold CV, out-of-sample)",
            "mean_abs_pooled": round(float(np.mean(pooled)),4), "mean_abs_within": round(float(np.mean(within)),4),
            "mean_inflation": round(float(np.mean(pooled-within)),4),
            "vs_single_snp": "single-SNP was pooled 0.286 / within 0.291 / inflation -0.005 (in-sample)"}
    Path(a.out + ".json").write_text(json.dumps({"summary": summ, "per_gene": rows}, indent=2), encoding="utf-8")
    print(json.dumps(summ, indent=2)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
