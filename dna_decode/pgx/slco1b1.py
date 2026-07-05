"""SLCO1B1 statin-myopathy caller — single-SNP rs4149056 (c.521T>C, p.Val174Ala), the *5 variant.

SLCO1B1 has a full star system, but its dominant clinical signal is the single c.521T>C variant
(rs4149056, the *5-defining SNP; also carried by *15/*17). The CPIC statin guideline (Cooper-DeHoff 2022)
assigns transporter FUNCTION largely from this 521T>C genotype -> simvastatin-associated myopathy risk. So
a single-SNP rs4149056 -> function cell IS the CPIC-aligned approach for the primary use.

STRAND (grounded): SLCO1B1 is on the PLUS strand, so c.521T>C is genomic T>C directly at
NC_000012.12:g.21178615 (GRCh38) — no strand flip (contrast VKORC1). ALT allele C == the decreased-function
521C allele.

Genotype -> function -> myopathy risk (CPIC simvastatin):
  T/T (521 TT) -> Normal Function       -> typical (low) myopathy risk
  T/C (521 TC) -> Decreased Function     -> intermediate myopathy risk
  C/C (521 CC) -> Poor Function          -> high myopathy risk
NOT a clinical tool.

HONEST TIER: this is a single-SNP genotype->function READOUT (KNOWLEDGE_BASELINE, like VKORC1), NOT a
star-diplotype caller needing independent disambiguation. rs4149056 IS the truth for a 521T>C call, so
"validation" is genotype-readout + trio-Mendelian consistency, never an independent star-concordance number.
"""
from __future__ import annotations

from pathlib import Path

GENE = "SLCO1B1"
ASSEMBLY = "GRCh38"
RSID = "rs4149056"
CHROM = "12"
POS = 21178615
REF = "T"
ALT = "C"          # plus-strand; ALT C == the decreased-function 521C (*5) allele

# genomic ALT (C) copy count -> (521 genotype, star proxy, CPIC function, simvastatin myopathy risk)
_FUNCTION = {
    0: ("T/T", "*1/*1", "Normal Function", "typical_risk"),
    1: ("T/C", "*1/*5", "Decreased Function", "intermediate_risk"),
    2: ("C/C", "*5/*5", "Poor Function", "high_risk"),
}


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def call_slco1b1(vcf: str | Path, sample: str | None = None) -> dict:
    """Read rs4149056 from a VCF; return the 521T>C genotype + CPIC function + statin-myopathy risk.
    Absent record -> ref (T/T) with an explicit assumed-reference flag. Raises on a named-but-absent sample."""
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
    genotype, star_proxy, function, risk = _FUNCTION[alt_count]
    return {
        "gene": GENE, "rsid": RSID, "assembly": ASSEMBLY,
        "position": f"chr{CHROM}:{POS}", "genomic_ref_alt": f"{REF}>{ALT}",
        "genomic_gt": raw_gt, "alt_count": alt_count,
        "variant_genotype": genotype,        # 521T>C genotype (plus-strand; no flip)
        "star_proxy": star_proxy,            # *1/*5 proxy from the single 521 SNP
        "function": function,                # CPIC transporter function
        "myopathy_risk": risk,               # simvastatin-associated myopathy risk band
        "status": "ok" if found else "assumed_reference",
        "flags": flags,
        "caveat": ("Single-SNP rs4149056 (c.521T>C, *5) SLCO1B1 function READOUT -> simvastatin myopathy "
                   "risk (CPIC Cooper-DeHoff 2022). Plus-strand: genomic T>C == cDNA 521T>C. This is a "
                   "single-SNP genotype->function call (KNOWLEDGE_BASELINE), NOT an independently-validated "
                   "star-diplotype number. Full SLCO1B1 star typing needs more variants. NOT a clinical tool."),
    }
