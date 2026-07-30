"""VCF input for the pigmentation cell — extract the pigmentation-SNP genotypes from a real genome VCF.

The `--genotypes rsID=GT,...` inline path is fine for a demo, but the real use case is a genome VCF (the same
input shape as `dna-decode pgx`). This parses a VCF by rsID for the requested SNP set (default = the 6 IrisPlex
EYE SNPs; pass the 41-SNP HIrisPlex-S panel for hair/skin) and builds the genotype dict the caller consumes.
A VCF is REFERENCE (+) strand; genotypes are strand-HARMONIZED to each SNP's counted allele (the hair/skin
models count on the webtool strand, so some SNPs are complemented — the eye counted alleles are forward, so
that path is unchanged). DTC-array (23andMe) strand quirks beyond ref-strand remain a documented follow-on.
Pure-python, wheel-only, offline.
"""
from __future__ import annotations

import gzip
from pathlib import Path

from dna_decode.pigment.irisplex import IRISPLEX_SNPS

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
# default SNP set = the 6 IrisPlex eye SNPs with their counted alleles (backward-compatible)
_DEFAULT_SNPS = [(rsid, allele) for rsid, allele, *_ in IRISPLEX_SNPS]


def _harmonize(geno: str, counted: str, site_alleles: set) -> str | None:
    """Return `geno` on the strand where `counted` is a site allele (complement if needed), or None if
    neither `counted` nor its complement is at the site (allele/strand mismatch → omit)."""
    if counted in site_alleles:
        return geno
    if _COMP.get(counted) in site_alleles:
        return "".join(_COMP[b] for b in geno)
    return None


def _open(path: str):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def _gt_to_alleles(ref: str, alts: list[str], gt: str) -> str | None:
    """Map a VCF GT (e.g. '0/1', '1|1') + REF/ALT to a 2-allele genotype string, or None if uncallable."""
    call = gt.split(":", 1)[0]
    sep = "|" if "|" in call else "/"
    idxs = call.split(sep)
    if len(idxs) != 2:
        return None
    alleles = [ref] + alts
    out = []
    for i in idxs:
        if i in (".", ""):
            return None            # missing genotype -> uncallable
        j = int(i)
        if j >= len(alleles):
            return None
        a = alleles[j]
        if len(a) != 1 or a.upper() not in "ACGT":
            return None            # indel / non-SNV at this site -> not an IrisPlex SNV call
        out.append(a.upper())
    return "".join(out)


def genotypes_from_vcf(vcf_path: str, snps=None) -> dict:
    """Return {rsID: genotype-string on the SNP's counted-allele strand} for whichever requested SNPs are
    present + callable in the VCF.

    `snps`: iterable of (rsid, counted_allele); default = the 6 IrisPlex EYE SNPs (backward-compatible).
    Pass the HIrisPlex-S hair/skin panel (41 SNPs) to decode those traits from a genome VCF. Genotypes are
    STRAND-HARMONIZED to each SNP's counted allele: a VCF is reference-strand, but the hair/skin models count
    on the webtool's strand (e.g. rs12913832_T = the reverse-complement of the forward A/G allele), so some
    SNPs are complemented. Eye counted alleles are forward, so the default path is unchanged (no complement).
    Matches the VCF ID column (col 3); FIRST sample column. Absent / uncallable / indel / strand-mismatch
    sites are omitted (the caller's `allow_missing` / required-SNP logic then applies). Raises
    FileNotFoundError if the VCF can't be opened.
    """
    wanted = dict(_DEFAULT_SNPS if snps is None else snps)   # {rsid: counted_allele}
    out: dict = {}
    with _open(vcf_path) as fh:
        for line in fh:
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            vid = cols[2]
            if vid not in wanted:
                continue
            ref, alt_field, gt = cols[3], cols[4], cols[9]
            alts = [] if alt_field in (".", "") else alt_field.split(",")
            geno = _gt_to_alleles(ref, alts, gt)
            if geno is None:
                continue
            site = {ref.upper()} | {a.upper() for a in alts if len(a) == 1}
            h = _harmonize(geno, wanted[vid], site)
            if h is not None:
                out[vid] = h
    return out
