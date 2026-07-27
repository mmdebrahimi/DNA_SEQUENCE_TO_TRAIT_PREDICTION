"""Stage-1: the LINEAR cis-eQTL ceiling on GEUVADIS, pooled vs WITHIN-population.

The project-owned organism-multimodal number. For each gene with a GEUVADIS best
cis-eQTL (EUR373 analysis), extract that SNP's imputed dosage per individual from the
matched VCFs, and measure how well it predicts expression CROSS-INDIVIDUAL -- pooled
across populations vs WITHIN each population (the de-confounding discipline). If the
cis-eQTL is a real causal signal it survives within-population; inflation = pooled-within
is the population-structure confound. This is the ceiling a DNA-encoder arm must beat.

Streams plain-gzip VCFs (no tabix); dosage from the DS FORMAT field. numpy/scipy only.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from dna_decode.organism_multimodal.geuvadis_data import (
    load_expr, parse_sample_population, canon_pop)
from dna_decode.organism_multimodal.deconfound_eval import evaluate

CACHE = Path("D:/dna_decode_cache/geuvadis")
GT = CACHE / "genotypes"
VCF = "GEUVADIS.chr{c}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"


def parse_eqtl_best(path: Path, chroms: set[str]) -> dict[str, tuple[str, int, float]]:
    """gene(ENSG) -> (chrom, snp_pos, eqtl_r) for genes whose best cis-eQTL is on `chroms`."""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11:
                continue
            gene, chrom, pos, r = p[2], p[4], p[6], p[9]
            if chrom in chroms:
                try:
                    out[gene] = (chrom, int(pos), float(r))
                except ValueError:
                    continue
    return out


def extract_dosages(vcf: Path, want_pos: set[int]) -> tuple[list[str], dict[int, np.ndarray]]:
    """Stream a plain-gzip VCF once; return (sample_ids, {pos: dosage array}) for wanted positions."""
    samples: list[str] = []
    ds_idx = None
    out: dict[int, np.ndarray] = {}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            # cheap position check before full split
            t1 = line.find("\t"); t2 = line.find("\t", t1 + 1)
            pos = int(line[t1 + 1:t2])
            if pos not in want_pos:
                continue
            p = line.rstrip("\n").split("\t")
            fmt = p[8].split(":")
            if ds_idx is None or fmt[min(len(fmt)-1,2)] != "DS":
                ds_idx = fmt.index("DS") if "DS" in fmt else None
            vals = np.full(len(samples), np.nan)
            if ds_idx is not None:
                for i, s in enumerate(p[9:]):
                    f = s.split(":")
                    if ds_idx < len(f) and f[ds_idx] not in (".", ""):
                        try: vals[i] = float(f[ds_idx])
                        except ValueError: pass
            else:  # fall back to GT -> alt-allele count
                for i, s in enumerate(p[9:]):
                    gt = s.split(":")[0].replace("|", "/")
                    if gt not in (".", "./."):
                        try: vals[i] = sum(int(x) for x in gt.split("/") if x != ".")
                        except ValueError: pass
            out[pos] = vals
    return samples, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default="1,2,3")
    ap.add_argument("--out", default="wiki/geuvadis_stage1_linear_ceiling_2026-07-27")
    ap.add_argument("--min-samples", type=int, default=50)
    args = ap.parse_args()
    chroms = set(args.chroms.split(","))

    E = load_expr(CACHE / "GD462.GeneQuantRPKM.txt.gz")
    S = parse_sample_population(CACHE / "E-GEUV-1.sdrf.txt")
    expr_pop = {s: canon_pop(S[s]) for s in E.samples if s in S}
    gene_row = {g: i for i, g in enumerate(E.genes)}

    eqtl = parse_eqtl_best(CACHE / "analysis_results/EUR373.gene.cis.FDR5.best.rs137.txt.gz", chroms)
    # keep only eQTL genes present in the expression matrix
    eqtl = {g: v for g, v in eqtl.items() if g in gene_row}
    print(f"eQTL genes on chr{sorted(chroms)} present in expr: {len(eqtl)}")

    per_gene = []
    for c in sorted(chroms, key=int):
        vcf = GT / VCF.format(c=c)
        if not vcf.exists():
            print(f"  chr{c} VCF missing -> skip"); continue
        want = {pos for (ch, pos, _r) in eqtl.values() if ch == c}
        print(f"  chr{c}: streaming for {len(want)} eQTL SNP positions ...")
        vcf_samples, dos = extract_dosages(vcf, want)
        # index VCF sample columns
        vcol = {s: i for i, s in enumerate(vcf_samples)}
        for g, (ch, pos, r) in eqtl.items():
            if ch != c or pos not in dos:
                continue
            d = dos[pos]
            row = E.values[gene_row[g]]
            xs, ys, gs = [], [], []
            for s in E.samples:
                if s in vcol and s in expr_pop:
                    dv = d[vcol[s]]
                    if not np.isnan(dv):
                        xs.append(float(row[E.samples.index(s)])); ys.append(dv); gs.append(expr_pop[s])
            if len(xs) < args.min_samples:
                continue
            res = evaluate(np.array(xs), np.array(ys), np.array(gs), n_perm=100)
            per_gene.append({"gene": g, "chrom": ch, "n": len(xs), "eqtl_r": r,
                             "pooled": res.pooled_rho, "within": res.within_rho_mean,
                             "inflation": res.inflation, "per_group": res.per_group_rho})
        print(f"  chr{c}: scored {sum(1 for x in per_gene if x['chrom']==c)} genes")

    if not per_gene:
        print("NO GENES SCORED"); return 1
    pooled = np.array([abs(x["pooled"]) for x in per_gene])
    within = np.array([abs(x["within"]) for x in per_gene])
    summary = {
        "n_genes": len(per_gene), "chroms": sorted(chroms),
        "mean_abs_pooled_spearman": float(np.mean(pooled)),
        "mean_abs_within_pop_spearman": float(np.mean(within)),
        "median_abs_pooled": float(np.median(pooled)),
        "median_abs_within": float(np.median(within)),
        "mean_inflation_pooled_minus_within": float(np.mean(pooled - within)),
        "frac_genes_within_ge_0.1": float(np.mean(within >= 0.1)),
        "interpretation": ("single best cis-eQTL SNP -> expression, cross-individual; "
                           "pooled vs within-population (CEU/FIN/GBR/TSI/YRI). This is the LINEAR "
                           "ceiling a DNA-encoder arm must beat. In-sample eQTL selection (EUR373); "
                           "the pooled-vs-within CONTRAST is the de-confounding signal."),
    }
    Path(args.out + ".json").write_text(json.dumps(
        {"summary": summary, "per_gene": per_gene[:500]}, indent=2), encoding="utf-8")
    print("\n=== STAGE-1 LINEAR cis-eQTL CEILING ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
