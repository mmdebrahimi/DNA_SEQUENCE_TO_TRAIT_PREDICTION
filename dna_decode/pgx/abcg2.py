"""ABCG2 statin-transporter caller — single-SNP rs2231142 (c.421C>A, p.Q141K), the poor-function variant.

ABCG2 is the CPIC ROSUVASTATIN transporter gene (completes the statin pair with SLCO1B1). ABCG2 effluxes
rosuvastatin; the 141K allele has POOR function, so 141K carriers have INCREASED rosuvastatin systemic
exposure and higher myopathy risk -> CPIC (Cooper-DeHoff 2022) recommends a lower rosuvastatin dose cap for
poor-function patients. Like SLCO1B1, the output is a transporter FUNCTION level, not a metabolizer phenotype.

STRAND (grounded at dbSNP): the cDNA c.421C>A sits on the ABCG2 minus strand, so on a plus-strand VCF it is
genomic G>T at NC_000004.12:g.88131171 (GRCh38). ALT allele T == the poor-function 141K allele. We score
genomic ALT (T) copy count directly.

Genotype -> function -> rosuvastatin exposure (CPIC Cooper-DeHoff 2022):
  G/G (141 Gln/Gln) -> Normal Function     -> typical exposure
  G/T (141 Gln/Lys) -> Decreased Function   -> intermediate exposure
  T/T (141 Lys/Lys) -> Poor Function        -> high exposure (lower rosuvastatin dose cap)
NOT a clinical tool.

HONEST TIER: single-SNP genotype->function READOUT (KNOWLEDGE_BASELINE, like SLCO1B1). rs2231142 IS the
truth for a 141K call, so "validation" is genotype-readout + trio-Mendelian consistency + AF-corroboration
(the 141K T allele ~9% EUR / ~29% EAS), never an independent star-concordance number.
"""
from __future__ import annotations

from pathlib import Path

GENE = "ABCG2"
ASSEMBLY = "GRCh38"
RSID = "rs2231142"
CHROM = "4"
POS = 88131171
REF = "G"
ALT = "T"          # plus-strand ALT == the poor-function 141K allele

# genomic ALT (T) copy count -> (141 genotype, star proxy, CPIC function, rosuvastatin exposure)
_FUNCTION = {
    0: ("Gln/Gln", "141Q/141Q", "Normal Function", "typical_exposure"),
    1: ("Gln/Lys", "141Q/141K", "Decreased Function", "intermediate_exposure"),
    2: ("Lys/Lys", "141K/141K", "Poor Function", "high_exposure"),
}


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def call_abcg2(vcf: str | Path, sample: str | None = None) -> dict:
    """Read rs2231142 from a VCF; return the Q141K genotype + CPIC function + rosuvastatin exposure band.
    Absent record -> ref (G/G) with an explicit assumed-reference flag. Raises on a named-but-absent sample."""
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
    genotype, star_proxy, function, exposure = _FUNCTION[alt_count]
    return {
        "gene": GENE, "rsid": RSID, "assembly": ASSEMBLY,
        "position": f"chr{CHROM}:{POS}", "genomic_ref_alt": f"{REF}>{ALT}",
        "genomic_gt": raw_gt, "alt_count": alt_count,
        "variant_genotype": genotype,        # 141 Gln/Lys genotype
        "star_proxy": star_proxy,            # 141Q/141K proxy from the single SNP
        "function": function,                # CPIC ABCG2 transporter function
        "rosuvastatin_exposure": exposure,   # rosuvastatin systemic exposure band
        "status": "ok" if found else "assumed_reference",
        "flags": flags,
        "caveat": ("Single-SNP rs2231142 (c.421C>A, Q141K) ABCG2 transporter-function READOUT -> "
                   "rosuvastatin exposure (CPIC Cooper-DeHoff 2022; pairs with SLCO1B1 for statins). "
                   "Minus-strand cDNA: genomic G>T == cDNA 141 Gln>Lys. Single-SNP genotype->function call "
                   "(KNOWLEDGE_BASELINE), NOT an independently-validated star number. NOT a clinical tool."),
    }
