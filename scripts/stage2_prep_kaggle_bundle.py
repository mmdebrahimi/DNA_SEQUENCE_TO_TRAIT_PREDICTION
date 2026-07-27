"""Stage-2 prep: build the Kaggle Enformer bundle for the DNA-encoder arm.

Per top cis-eQTL gene (chr1): TSS-centered 196608bp reference window + the cis-variants
in it (offset, alt base, per-individual dosage) + measured expression + sample/pop. Kaggle
runs Enformer ref vs per-variant-alt at the TSS bin -> effect; aggregate by dosage ->
predicted expression; compare to the Stage-1 linear ceiling (rho~0.19-0.20), within-pop.
"""
from __future__ import annotations
import gzip, json
from pathlib import Path
import numpy as np
from pyfaidx import Fasta
from dna_decode.organism_multimodal.geuvadis_data import load_expr, parse_sample_population, canon_pop

CACHE = Path("D:/dna_decode_cache/geuvadis"); REFS = Path("D:/dna_decode_cache/refs")
FA = Fasta(str(REFS / "GRCh37.chr1.fa"))
VCF = CACHE / "genotypes/GEUVADIS.chr1.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"
WIN = 196608; HALF = WIN // 2
N_GENES = 12; MAX_VARS = 50; VAR_HALF = 25000   # collect cis-variants near TSS (window center, Enformer-accurate)

def gencode_tss(chrom="chr1"):
    tss = {}
    with gzip.open(REFS / "gencode.v19.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.split("\t")
            if p[0] != chrom or p[2] != "gene": continue
            gid = p[8].split('gene_id "')[1].split('"')[0]
            base = gid.split(".")[0]
            tss[base] = int(p[3]) if p[6] == "+" else int(p[4])
    return tss

def top_eqtl_genes(gene_row, tss_base, n):
    cand = []
    with gzip.open(CACHE / "analysis_results/EUR373.gene.cis.FDR5.best.rs137.txt.gz", "rt") as fh:
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) >= 11 and q[4] == "1" and q[2] in gene_row and q[2].split(".")[0] in tss_base:
                cand.append((q[2], abs(float(q[9]))))
    cand.sort(key=lambda x: -x[1])
    return [g for g, _ in cand[:n]]

def main():
    E = load_expr(CACHE / "GD462.GeneQuantRPKM.txt.gz")
    S = parse_sample_population(CACHE / "E-GEUV-1.sdrf.txt")
    gene_row = {g: i for i, g in enumerate(E.genes)}
    sidx = {s: i for i, s in enumerate(E.samples)}
    tss_base = gencode_tss()
    genes = top_eqtl_genes(gene_row, tss_base, N_GENES)
    print(f"target genes (chr1, top {N_GENES} cis-eQTL w/ TSS): {len(genes)}")
    windows = {}   # gene -> (start,end,tss,tss_bin)
    for g in genes:
        tss = tss_base[g.split(".")[0]]
        st = max(0, tss - HALF); en = st + WIN
        windows[g] = (st, en, tss)
    # collect cis-window variants + DS dosage from chr1 VCF (one stream)
    per_gene_vars = {g: [] for g in genes}
    vcf_samples = []
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                vcf_samples = line.rstrip("\n").split("\t")[9:]; continue
            t1 = line.find("\t"); t2 = line.find("\t", t1 + 1)
            pos = int(line[t1 + 1:t2])
            hits = [g for g in genes if abs(pos - windows[g][2]) <= VAR_HALF and len(per_gene_vars[g]) < MAX_VARS]
            if not hits: continue
            p = line.rstrip("\n").split("\t")
            ref, alt = p[3], p[4]
            if len(ref) != 1 or len(alt) != 1: continue          # SNP only
            fmt = p[8].split(":"); di = fmt.index("DS") if "DS" in fmt else None
            if di is None: continue
            ds = np.full(len(vcf_samples), np.nan)
            for i, s in enumerate(p[9:]):
                f = s.split(":")
                if di < len(f) and f[di] not in (".", ""):
                    try: ds[i] = float(f[di])
                    except ValueError: pass
            if np.isnan(ds).mean() > 0.1: continue
            m = np.nanmean(ds) / 2
            if min(m, 1 - m) < 0.05: continue
            ds = np.where(np.isnan(ds), np.nanmean(ds), ds)
            for g in hits: per_gene_vars[g].append((pos, ref, alt, ds))
    vcol = {s: i for i, s in enumerate(vcf_samples)}
    common = [s for s in E.samples if s in vcol and s in S]
    pops = [canon_pop(S[s]) for s in common]
    vpos = [vcol[s] for s in common]
    # build bundle
    bundle = {"samples": common, "pops": pops, "genes": [], "win": WIN}
    arrays = {}
    for g in genes:
        vs = per_gene_vars[g]
        if len(vs) < 5: continue
        st, en, tss = windows[g]
        seq = str(FA["1"][st:en]).upper()
        if len(seq) < WIN: seq = seq + "N" * (WIN - len(seq))
        offs = np.array([v[0] - st for v in vs], dtype=np.int32)
        alts = np.array([ord(v[2]) for v in vs], dtype=np.uint8)
        dosg = np.array([v[3][vpos] for v in vs], dtype=np.float32)     # nvar x ncommon
        expr = np.array([float(E.values[gene_row[g]][sidx[s]]) for s in common], dtype=np.float32)
        arrays[f"{g}__seq"] = np.frombuffer(seq.encode(), dtype=np.uint8)
        arrays[f"{g}__off"] = offs; arrays[f"{g}__alt"] = alts
        arrays[f"{g}__dos"] = dosg; arrays[f"{g}__expr"] = expr
        bundle["genes"].append({"gene": g, "tss": tss, "win_start": st, "n_var": len(vs)})
        print(f"  {g}: {len(vs)} vars, TSS={tss}")
    out = Path("D:/dna_decode_cache/stage2_kaggle"); out.mkdir(exist_ok=True)
    np.savez_compressed(out / "stage2_bundle.npz", **arrays)
    (out / "stage2_meta.json").write_text(json.dumps(bundle), encoding="utf-8")
    print(f"\nbundle: {len(bundle['genes'])} genes, {len(common)} individuals -> {out}/stage2_bundle.npz")
    print("pops:", {p: pops.count(p) for p in set(pops)})

if __name__ == "__main__":
    main()
