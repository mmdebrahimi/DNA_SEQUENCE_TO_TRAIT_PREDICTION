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

HONEST TIER: this is a single-SNP genotype->function READOUT, NOT a star-diplotype caller.

*** THE "NO STAR-CONCORDANCE NUMBER EXISTS" CLAIM WAS HALF WRONG -- MEASURED 2026-09-22 ***
This file used to assert that validation "is genotype-readout + trio-Mendelian consistency, never an
independent star-concordance number." Measured against GeT-RM CDC multi-lab consensus on 88 real 1000G
samples (wiki/slco1b1_getrm_concordance_2026-09-22.{md,json}; bar + predictions frozen beforehand), the
two halves separate cleanly:

  * RIGHT about the FUNCTION axis. GeT-RM's own star assignment derives from genotypes INCLUDING 521T>C,
    so scoring our 521-derived function call against it is close to self-comparison. It is exactly 1:1
    (88/88, zero off-diagonal) -- which is DEFINITIONAL, a control that bounds plumbing, NOT evidence.
  * WRONG that no star-concordance number exists. One exists, it is measurable, and it is BAD: exact
    star-NAME agreement with the reference is 35/88 = 0.398. The printed `*1/*5`-style label disagrees
    with the reference on the MAJORITY of real samples.

MECHANISM, and it is a COMPLETE explanation rather than a partial one: *1B/*1A/*14/*21 are defined by
388A>G (and more) and *15/*17 are 388+521 -- none of which this caller genotypes. Star-name agreement is
EXACTLY the complement of "carries an allele the caller cannot see", with ZERO exceptions across all 88
samples, and every non-*1 reference allele scores 0/n on naming (*1B 0/23, *15 0/15, *14 0/8, *1A 0/8,
*17 0/4, *21 0/3, *5 0/1).

SAFETY (the axis that matters clinically): ZERO misses -- no truth reduced-function carrier (*5/*15/*17,
all 521C-bearing) was called Normal Function. The function call is right even where the NAME is wrong.

CONSUMER RULE: read the FUNCTION, never the star label. The `*1/*5` string is a function readout wearing
a star-shaped label, not a star call -- a true *1/*15 carrier prints as *1/*5 (same decreased function,
different allele). Full star typing needs 388A>G + promoter variants this cell does not read.
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
        "star_label_disagrees_with_reference_rate": 0.602,   # 53/88 vs GeT-RM, measured 2026-09-22
        "caveat": ("Single-SNP rs4149056 (c.521T>C, *5) SLCO1B1 function READOUT -> simvastatin myopathy "
                   "risk (CPIC Cooper-DeHoff 2022). Plus-strand: genomic T>C == cDNA 521T>C. READ THE "
                   "FUNCTION, NOT THE STAR LABEL: measured vs GeT-RM CDC consensus on 88 real 1000G "
                   "samples (2026-09-22), the printed star name DISAGREES with the reference on 53/88 "
                   "(exact agreement 0.398) because *1B/*1A/*14/*21 (388A>G) and *15/*17 (388+521) are "
                   "invisible to this caller -- a true *1/*15 carrier prints as *1/*5. The FUNCTION call "
                   "was right on all 88 (0 reduced-function carriers called Normal), so the label is "
                   "wrong where the clinical call is right. NOT a clinical tool."),
    }
