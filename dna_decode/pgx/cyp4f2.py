"""CYP4F2 warfarin-dose caller — single-SNP rs2108622 (c.1297G>A, p.V433M), the *3 variant.

CYP4F2 is the THIRD CPIC warfarin gene (completes the triad with VKORC1 + CYP2C9). It is NOT a
metabolizer-phenotype gene — CYP4F2 recycles vitamin K, and the *3 (433M) allele has REDUCED activity, so
*3 carriers accumulate more vitamin K and need a slightly HIGHER warfarin dose (Johnson 2017 CPIC warfarin
guideline adds ~+0.4 mg/day per *3 allele on top of the VKORC1 + CYP2C9 dose).

STRAND (grounded at dbSNP): the cDNA c.1297G>A sits on the CYP4F2 minus strand, so on a plus-strand VCF it
is genomic C>T at NC_000019.10:g.15879621 (GRCh38). ALT allele T == the reduced-function *3 (433M) allele.
So we score genomic ALT (T) copy count directly.

Genotype -> function -> warfarin dose direction (CPIC Johnson 2017):
  C/C (433 Val/Val)  -> Normal Function     -> no CYP4F2 dose adjustment
  C/T (433 Val/Met)  -> Intermediate         -> slightly higher warfarin dose (~+0.4 mg/day)
  T/T (433 Met/Met)  -> Reduced Function      -> higher warfarin dose (~+0.8 mg/day)
NOT a clinical tool.

HONEST TIER: single-SNP genotype->function READOUT. rs2108622 IS the truth for a *3 call.

*** TWO CLAIMS THIS FILE USED TO MAKE WERE MEASURED 2026-09-24 -- ONE WAS WRONG ***
(wiki/cyp4f2_getrm_concordance_2026-09-24.{md,json} + wiki/pgx_af_panel_audit_2026-09-24.{md,json};
acceptance bar AND predictions frozen beforehand.)

  1. "NEVER an independent star-concordance number" -- WRONG, same as the identical sentence in
     slco1b1.py. One exists: exact star-NAME agreement vs GeT-RM CDC multi-lab consensus is 52/64 =
     **0.8125** on real 1000G samples. The 12 disagreements are EXACTLY the 12 samples carrying *2, an
     allele defined by a variant this caller never reads (it genotypes one SNP, rs2108622) -- zero
     exceptions, the same mechanistic identity SLCO1B1 showed. So the printed label is right on the common
     *1//*3 diplotypes and wrong wherever *2 appears.
  2. "*3 T allele ~29% EUR / ~79% EAS" -- the EAS figure was WRONG. Measured over the 3202-sample 1000G
     panel that ships beside this cell: **EUR 0.2773 (the ~29% claim holds) and EAS 0.2128, not 0.79**.
     Note 1 - 0.2128 = 0.7872, so the asserted figure matched the REFERENCE-allele frequency to within 0.3
     percentage points -- consistent with an allele-polarity slip, though the measurement cannot establish
     how the number got there. Full measured spectrum: AFR 0.078 / AMR 0.247 / EAS 0.213 / EUR 0.277 /
     SAS 0.418 / global 0.232.

     It survived because scripts/pgx_af_corroboration.py checks ONE population per variant -- for CYP4F2,
     EUR, the number that is correct. The standing guard is now scripts/pgx_af_panel_audit.py, which
     derives EVERY asserted population frequency from the local panel. On that audit the rest of the PGx
     surface is accurate (UGT1A1 *80 EUR delta 0.0009, UGT1A1 *6 EAS delta -0.0032, CYP2B6 delta -0.0004),
     so this was an isolated wrong number, NOT a pattern.

FUNCTION-axis note: a GeT-RM *3-dosage concordance number is DEFINITIONAL, not validation -- GeT-RM's own
*3 assignment derives from rs2108622, the one SNP this caller reads. It is 64/64 and that is expected.

CONSUMER RULE: read the FUNCTION, not the star label (as in SLCO1B1). *2's FUNCTIONAL status is NOT settled
by any of this -- CPIC's warfarin algorithm dose-adjusts on *3, and whether *2 carries an independent effect
is a separate literature question this repo has not answered. So a *2 carrier printed as *1 is a NAMING
limitation here, not a demonstrated safety gap.
"""
from __future__ import annotations

from pathlib import Path

GENE = "CYP4F2"
ASSEMBLY = "GRCh38"
RSID = "rs2108622"
CHROM = "19"
POS = 15879621
REF = "C"
ALT = "T"          # plus-strand ALT == the reduced-function *3 (433M) allele

# genomic ALT (T) copy count -> (433 genotype, star proxy, CPIC function, warfarin dose direction)
_FUNCTION = {
    0: ("Val/Val", "*1/*1", "Normal Function", "no_adjustment"),
    1: ("Val/Met", "*1/*3", "Intermediate", "slightly_higher_dose"),
    2: ("Met/Met", "*3/*3", "Reduced Function", "higher_dose"),
}


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def call_cyp4f2(vcf: str | Path, sample: str | None = None) -> dict:
    """Read rs2108622 from a VCF; return the *3 (V433M) genotype + CPIC function + warfarin dose direction.
    Absent record -> ref (C/C) with an explicit assumed-reference flag. Raises on a named-but-absent sample."""
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
    genotype, star_proxy, function, dose = _FUNCTION[alt_count]
    return {
        "gene": GENE, "rsid": RSID, "assembly": ASSEMBLY,
        "position": f"chr{CHROM}:{POS}", "genomic_ref_alt": f"{REF}>{ALT}",
        "genomic_gt": raw_gt, "alt_count": alt_count,
        "variant_genotype": genotype,        # 433 Val/Met genotype
        "star_proxy": star_proxy,            # *1/*3 proxy from the single 433 SNP
        "function": function,                # CPIC CYP4F2 function
        "warfarin_dose_direction": dose,     # CPIC warfarin dose adjustment direction
        "status": "ok" if found else "assumed_reference",
        "flags": flags,
        "star_label_disagrees_with_reference_rate": 0.1875,   # 12/64 vs GeT-RM, measured 2026-09-24
        "caveat": ("Single-SNP rs2108622 (c.1297G>A, *3) CYP4F2 function READOUT -> warfarin dose direction "
                   "(CPIC Johnson 2017; the 3rd warfarin gene with VKORC1 + CYP2C9). Minus-strand cDNA: "
                   "genomic C>T == cDNA 433 Val>Met. READ THE FUNCTION, NOT THE STAR LABEL: measured vs "
                   "GeT-RM CDC consensus on 64 real 1000G samples (2026-09-24), the printed star name "
                   "agrees 52/64 = 0.8125, and all 12 disagreements are exactly the samples carrying *2 -- "
                   "an allele defined by a variant this caller never reads, so a *2 carrier prints as *1. "
                   "The FUNCTION call was right on all 64. *2's own functional status is NOT established "
                   "here. NOT a clinical tool."),
    }
