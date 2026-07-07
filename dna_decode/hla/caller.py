"""HLA drug-hypersensitivity caller — tag-SNP genotype -> HLA-allele carriage -> CPIC drug action.

Pure-stdlib single-SNP VCF read (mirrors the SLCO1B1/VKORC1 readout cells). Reads the tag SNP for a chosen
HLA allele; >=1 copy of the tag ALT -> CARRIER -> the CPIC drug action. Absence of a record -> reference
(non-carrier) with an explicit assumed-reference flag (never silent). The tag is an LD PROXY -> the record
always carries the proxy tier + the honest "concordance vs real HLA truth is the validation number" caveat.
"""
from __future__ import annotations

import datetime
from pathlib import Path

from dna_decode.hla.catalog import ASSEMBLY, get

SCHEMA = "hla-tag-carriage-v0"


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def call_hla(vcf: str | Path, allele_key: str, sample: str | None = None) -> dict:
    """Read the tag SNP for `allele_key` from a VCF -> HLA-allele carriage + CPIC drug action record.
    Raises on a named-but-absent sample; an absent record -> non-carrier (assumed reference), flagged."""
    a = get(allele_key)
    sample_idx = 0
    found = False
    tag_count = 0
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
        if _norm_chrom(cols[0]) != a.chrom or int(cols[1]) != a.pos:
            continue
        found = True
        alts = cols[4].split(",")
        ai = alts.index(a.tag_alt) + 1 if a.tag_alt in alts else -1
        if len(cols) >= 10:
            fmt = cols[8].split(":")
            col = 9 + sample_idx
            if "GT" in fmt and col < len(cols):
                raw_gt = cols[col].split(":")[fmt.index("GT")]
                no_call = "." in raw_gt
                if ai > 0:
                    nums = [int(x) for x in raw_gt.replace("|", "/").split("/") if x.isdigit()]
                    tag_count = sum(1 for n in nums if n == ai)
        break

    if not found:
        flags.append("assumed_reference_at_uncalled_site")
    if no_call:
        flags.append("no_call")
    carrier = tag_count >= 1
    zygosity = {0: "non-carrier", 1: "heterozygous carrier", 2: "homozygous carrier"}[tag_count]
    return {
        "trait": "hla_drug_hypersensitivity", "allele": a.allele, "allele_key": a.key,
        "organism": "Homo sapiens", "assembly": ASSEMBLY, "schema": SCHEMA,
        "analysis_date": datetime.date.today().isoformat(),
        "tag_rsid": a.rsid, "position": f"chr{a.chrom}:{a.pos}", "tag_ref_alt": f"{a.ref}>{a.tag_alt}",
        "tag_gt": raw_gt, "tag_copies": tag_count,
        "carrier": carrier, "zygosity": zygosity,
        "drug": a.drug, "reaction": a.reaction,
        "risk_call": a.cpic_action if carrier else f"no {a.allele} tag detected — standard {a.drug} risk",
        "proxy_tier": a.proxy_tier,
        "status": "ok" if found else "assumed_reference",
        "flags": flags,
        "caller": {
            "name": f"dna_decode-hla-{a.key}-v0",
            "method": "vcf_tag_snp -> HLA-allele carriage (LD proxy) -> CPIC drug action",
            "is_ld_proxy": True,            # NOT sequence-based typing; carriage inferred via a tag SNP
            "proxy_note": a.proxy_note,
            "reference_tool": "sequence-based HLA typing (HLA*LA / arcasHLA / OptiType); free 1000G HLA truth",
        },
        "caveat": (f"{a.allele} carriage inferred from the TAG SNP {a.rsid} ({a.ref}>{a.tag_alt}) — an LD "
                   f"PROXY, NOT sequence-based HLA typing. {a.proxy_note} VALIDATED vs the free 1000G HLA "
                   "truth (not the literature LD alone). NOT a clinical tool."),
        "source": a.source,
    }
