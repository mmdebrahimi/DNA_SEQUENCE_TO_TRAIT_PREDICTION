"""VCF input for the pigmentation cell — extract the 6 IrisPlex SNP genotypes from a real genome VCF.

The `--genotypes rsID=GT,...` inline path is fine for a demo, but the real use case is a genome VCF (the same
input shape as `dna-decode pgx`). This parses a VCF by rsID for the IrisPlex SNPs and builds the genotype dict
`predict_eye_color` consumes. VCF genotypes are on the REFERENCE (+) strand; DTC-array strand harmonization
(23andMe/AncestryDNA report some SNPs on the opposite strand) remains the documented v0.1 follow-on and is NOT
applied here (a VCF is already reference-strand). Pure-python, wheel-only, offline.
"""
from __future__ import annotations

import gzip
from pathlib import Path

from dna_decode.pigment.irisplex import IRISPLEX_SNPS

_WANTED = {rsid for rsid, *_ in IRISPLEX_SNPS}


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


def genotypes_from_vcf(vcf_path: str) -> dict:
    """Return {rsID: genotype-string} for whichever of the 6 IrisPlex SNPs are present + callable in the VCF.

    Matches the VCF ID column (col 3) against the IrisPlex rsIDs; uses the FIRST sample column. Absent /
    uncallable / indel sites are simply omitted (the caller's `allow_missing` / required-SNP logic then
    applies). Raises FileNotFoundError if the VCF can't be opened.
    """
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
            if vid not in _WANTED:
                continue
            ref, alt_field, gt = cols[3], cols[4], cols[9]
            alts = [] if alt_field in (".", "") else alt_field.split(",")
            geno = _gt_to_alleles(ref, alts, gt)
            if geno is not None:
                out[vid] = geno
    return out
