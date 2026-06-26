"""VKORC1 warfarin-sensitivity caller -- the non-star half of the warfarin pair.

VKORC1 is NOT a star-allele system: the single promoter SNP rs9923231 (c.-1639G>A) drives warfarin
sensitivity (the A allele lowers VKORC1 expression -> less enzyme -> more sensitive -> lower dose).
CPIC's warfarin guideline (Johnson 2017) combines the CYP2C9 diplotype + this VKORC1 genotype for dosing.

STRAND SUBTLETY (load-bearing, grounded at dbSNP): VKORC1 is transcribed on the MINUS strand, so the
cDNA c.-1639G>A is genomic C>T at NC_000016.10:g.31096368 (GRCh38). Therefore in a plus-strand VCF the
ALT allele T == the cDNA 'A' (sensitive) allele. We score genomic ALT-count and translate to the cDNA
genotype + a sensitivity category. (rs9934438, in high LD, is a common clinical surrogate -- not used here.)

Genotype -> sensitivity (Johnson 2017 average daily warfarin dose):
  C/C (cDNA G/G) -> Normal sensitivity  (~6 mg/day)
  C/T (cDNA G/A) -> Intermediate        (~5 mg/day)
  T/T (cDNA A/A) -> High sensitivity    (~3 mg/day, low dose)
NOT a clinical tool.
"""
from __future__ import annotations

from pathlib import Path

GENE = "VKORC1"
ASSEMBLY = "GRCh38"
RSID = "rs9923231"
CHROM = "16"
POS = 31096368
REF = "C"
ALT = "T"          # plus-strand ALT == cDNA 'A' (the sensitive -1639A allele)

# genomic ALT (T) copy count -> (cDNA genotype, sensitivity, approx dose category)
_SENSITIVITY = {
    0: ("G/G", "Normal sensitivity", "normal_dose"),
    1: ("G/A", "Intermediate sensitivity", "intermediate_dose"),
    2: ("A/A", "High sensitivity", "low_dose"),
}


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def call_vkorc1(vcf: str | Path, sample: str | None = None) -> dict:
    """Read rs9923231 from a VCF; return genotype + warfarin-sensitivity category. Absent record -> ref
    (G/G) with an explicit assumed-reference flag (never silent). Raises on a named-but-absent sample."""
    sample_idx = 0
    found = False
    alt_count = 0
    no_call = False
    raw_gt = None
    flags: list[str] = []
    for line in Path(vcf).read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            if sample is not None:
                if sample not in samples:
                    raise ValueError(f"--sample {sample!r} not found in VCF header")
                sample_idx = samples.index(sample)
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 8 or not cols[1].isdigit():
            continue
        if _norm_chrom(cols[0]) != CHROM or int(cols[1]) != POS:
            continue
        found = True
        alts = cols[4].split(",")
        ai = alts.index(ALT) + 1 if ALT in alts else -1
        if len(cols) >= 10:
            fmt = cols[8].split(":")
            col = 9 + sample_idx
            if "GT" in fmt and col < len(cols):
                raw_gt = cols[col].split(":")[fmt.index("GT")]
                no_call = "." in raw_gt
                if ai > 0:
                    nums = [int(a) for a in raw_gt.replace("|", "/").split("/") if a.isdigit()]
                    alt_count = sum(1 for n in nums if n == ai)
        break

    if not found:
        flags.append("assumed_reference_at_uncalled_site")
    if no_call:
        flags.append("no_call")
    cdna_gt, sensitivity, dose = _SENSITIVITY[alt_count]
    return {
        "gene": GENE, "rsid": RSID, "assembly": ASSEMBLY,
        "position": f"chr{CHROM}:{POS}", "genomic_ref_alt": f"{REF}>{ALT}",
        "genomic_gt": raw_gt, "alt_count": alt_count,
        "cdna_genotype": cdna_gt,        # c.-1639 G/A notation (minus-strand corrected)
        "warfarin_sensitivity": sensitivity, "dose_category": dose,
        "status": "ok" if found else "assumed_reference",
        "flags": flags,
        "caveat": ("Single-SNP rs9923231 (c.-1639G>A) warfarin-sensitivity call; VKORC1 is minus-strand so "
                   "genomic C>T == cDNA G>A (the T/A allele = sensitive/low-dose). Combine with the CYP2C9 "
                   "diplotype for CPIC warfarin dosing. NOT a clinical tool."),
    }
