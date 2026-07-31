"""VCF input for the DOG morphology cell — extract the pinned body-size + ear SNP dosages from a real dog
genome VCF (canFam4 / UU_Cfam_GSD_1.0).

The `--dosages LOCUS=n` inline path is fine for a panel/PLINK export, but the real use case is a genome VCF.
Unlike the HUMAN pigment SNPs (matched by rsID in vcf_input.py), the dog body-size/ear causal SNPs have no
stable rsIDs — they are matched by canFam4 COORDINATE (chr:pos) from the pinned catalog (dog_body_size). For
each pinned SNP present + callable, the diploid genotype is converted to the BIG-ALLELE DOSAGE (0/1/2) the
morphology caller consumes. Strand-harmonized to the catalog allele (complement if the VCF site is the other
strand). Chromosome names are normalized (`chr10` == `10`); RefSeq `NC_...` accessions are a documented
follow-on. Pure-python, wheel-only, offline; reuses the general VCF helpers from vcf_input.
"""
from __future__ import annotations

from dna_decode.pigment.dog_body_size import MORPH_LOCI, SIZE_LOCI
from dna_decode.pigment.vcf_input import _COMP, _gt_to_alleles, _open


def _norm_chrom(c: str) -> str:
    """Normalize a chromosome token: strip a leading 'chr', lowercase (so 'chr10' == '10')."""
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def pinned_loci() -> dict:
    """{locus -> (norm_chrom, pos:int, big_allele)} for the 4 size loci + EAR, from the pinned catalog."""
    loci: dict = {}
    for name, L in SIZE_LOCI.items():
        chrom, pos, _ref, _alt = L.canfam4_variant.split(":")
        loci[name] = (_norm_chrom(chrom), int(pos), L.big_allele)
    ear = MORPH_LOCI["EAR"]
    chrom, pos, _ref, _alt = ear.canfam4_variant.split(":")
    loci["EAR"] = (_norm_chrom(chrom), int(pos), ear.high_allele)
    return loci


def dosages_from_vcf(vcf_path: str, loci: dict | None = None) -> dict:
    """Return {locus -> big-allele dosage 0/1/2} for whichever pinned SNPs are present + callable in the VCF.

    `loci`: {locus -> (norm_chrom, pos, big_allele)}; default = the pinned body-size + ear catalog
    (`pinned_loci()`). Matches by (normalized chrom, POS). The big allele is strand-harmonized to the VCF
    site (complemented if the site carries the other strand); the dosage is the count of that allele in the
    diploid genotype. Absent / uncallable (missing GT / indel / '.') / allele-mismatch sites are omitted —
    the caller's partial-panel logic then applies. Raises FileNotFoundError if the VCF can't be opened.
    """
    want = {}  # (norm_chrom, pos) -> (locus, big_allele)
    for name, (chrom, pos, big) in (loci or pinned_loci()).items():
        want[(chrom, pos)] = (name, big.upper())
    out: dict = {}
    with _open(vcf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            key = (_norm_chrom(cols[0]), _safe_int(cols[1]))
            if key is None or key not in want:
                continue
            name, big = want[key]
            ref, alt_field, gt = cols[3], cols[4], cols[9]
            alts = [] if alt_field in (".", "") else alt_field.split(",")
            geno = _gt_to_alleles(ref, alts, gt)     # 2-char diploid on the VCF strand, or None
            if geno is None:
                continue
            site = {ref.upper()} | {a.upper() for a in alts if len(a) == 1}
            big_h = big if big in site else (_COMP.get(big) if _COMP.get(big) in site else None)
            if big_h is None:
                continue                              # allele/strand mismatch -> omit (never a wrong dose)
            out[name] = geno.count(big_h)             # 0/1/2 copies of the big allele
    return out


def _safe_int(s: str):
    try:
        return int(s)
    except (ValueError, TypeError):
        return None
