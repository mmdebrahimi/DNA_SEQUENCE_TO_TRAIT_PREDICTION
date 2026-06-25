"""CYP2C19 PGx runner — assemble the provenance record (the report shape) from a diplotype call."""
from __future__ import annotations

import datetime
from pathlib import Path

from dna_decode.pgx.caller import call_diplotype
from dna_decode.pgx.cyp2c19_catalog import (
    ASSEMBLY,
    CORE_DEFINING,
    GENE,
    PHENOTYPE_ABBREV,
    UNDETECTABLE,
)

SCHEMA = "pgx-diplotype-call-v0"

CAVEAT = (
    "Star-allele CALLING is deterministic + independently validatable vs the GeT-RM consensus panel "
    "(caller_is_independent_baseline=True for the calling step). The metabolizer PHENOTYPE is "
    "FAITHFUL-TO-CPIC -- assigned from the diplotype via CPIC's table, NOT a measured probe-drug PK "
    "phenotype (caller_is_independent_baseline=False for that step; reference tool = PharmCAT). v0 covers "
    "the CORE SNP-defined alleles (*2/*3/*17 + *1); a non-core star allele is mis-called *1 (a flagged "
    "blind spot). Input = a phased VCF, GRCh38. NOT a clinical decision tool."
)


def call_cyp2c19(vcf: str | Path, sample_id: str | None = None,
                 sample_column: str | None = None) -> dict:
    """Run the CYP2C19 caller on a VCF and return the full provenance record."""
    res = call_diplotype(vcf, sample=sample_column)
    sid = sample_id or sample_column or Path(vcf).stem
    rec = {
        "sample_id": sid,
        "trait": "pgx_metabolizer_phenotype",
        "gene": GENE,
        "organism": "Homo sapiens",
        "assembly": ASSEMBLY,
        "analysis_date": datetime.date.today().isoformat(),
        "schema": SCHEMA,
        "status": res.status,
        "diplotype": res.diplotype,
        "allele1": res.allele1,
        "allele2": res.allele2,
        "phenotype": res.phenotype,
        "phenotype_abbrev": PHENOTYPE_ABBREV.get(res.phenotype or "", None),
        "phasing": res.phasing,
        "flags": res.flags,
        "variant_calls": res.variant_calls,
        "caller": {
            "name": "dna_decode-pgx-cyp2c19-v0",
            "method": "vcf_defining_snp -> star_allele -> diplotype -> CPIC_phenotype",
            "calling_is_independent_baseline": True,
            "phenotype_is_independent_baseline": False,
            "reference_tool": "PharmCAT",
        },
        "catalog": {
            "gene": GENE,
            "core_alleles": [d.star for d in CORE_DEFINING] + ["*1"],
            "defining_variants": [
                {"star": d.star, "rsid": d.rsid, "chrom": d.chrom, "pos": d.pos,
                 "ref": d.ref, "alt": d.alt, "cdna": d.cdna}
                for d in CORE_DEFINING
            ],
            "source": "PharmVar CYP2C19 + CPIC/PharmGKB (Caudle 2020; Botton 2021; Gaedigk 2022)",
        },
        "undetectable": UNDETECTABLE,
        "caveat": CAVEAT,
    }
    if res.reason:
        rec["reason"] = res.reason
    return rec
