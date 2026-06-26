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
    "A deterministic CORE-MARKER PROXY caller (the core SNP set *2/*3/*17 + *1 + non-core sentinels), NOT a "
    "full PharmVar star-allele caller. Star-allele CALLING is DESIGNED to be independently validatable vs "
    "the GeT-RM consensus panel; the validation RUN so far is faithful-to-PharmCAT (in-distribution) -- see "
    "caller.independent_validation_status. The metabolizer PHENOTYPE is FAITHFUL-TO-CPIC (assigned from the "
    "diplotype, NOT a measured probe-drug PK phenotype; reference tool = PharmCAT). When a NON-CORE sentinel "
    "(*4 via rs28399504, *35 via rs12769205) proves an allele the core proxy cannot resolve, the phenotype is "
    "WITHHELD (phenotype_status=phenotype_withheld) rather than mis-called. Check phenotype_status, NOT just "
    "status, before consuming the phenotype. Input = a phased VCF, GRCh38. NOT a clinical decision tool."
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
        "status": res.status,                       # PARSE status (did we read a usable VCF?)
        "phenotype_status": res.phenotype_status,   # CONSUMABILITY (ok / phenotype_withheld / phase_ambiguous)
        "phenotype_confidence": res.phenotype_confidence,
        "diplotype": res.diplotype,
        "core_proxy_diplotype": res.core_proxy_diplotype,
        "allele1": res.allele1,
        "allele2": res.allele2,
        "phenotype": res.phenotype,
        "phenotype_abbrev": PHENOTYPE_ABBREV.get(res.phenotype or "", None),
        "alternate_diplotype": res.alternate_diplotype,
        "alternate_phenotype": res.alternate_phenotype,
        "sentinel_hits": res.sentinel_hits,
        "phasing": res.phasing,
        "flags": res.flags,
        "variant_calls": res.variant_calls,
        "caller": {
            "name": "dna_decode-pgx-cyp2c19-v0.1",
            "method": "vcf_core_snp_proxy -> star_allele -> diplotype -> CPIC_phenotype (+ non-core sentinel withhold)",
            # HONESTY: the calling step is DESIGNED to be independently validatable vs GeT-RM, but the only
            # validation RUN so far is faithful-to-PharmCAT (in-distribution). Do not claim achieved independence.
            "calling_independently_validatable": True,
            "independent_validation_status": (
                "GeT-RM consensus: core diplotype 72/72 (100%) on 87 1000G samples (caller independent of "
                "the Astrolabe/Stargazer/Aldy consensus tools); +7 *38==*1 phenotype-equivalent; 2 non-core "
                "correctly withheld; 6/87 (6.9%) non-core silent residual (*8/*13/*15/*39). PharmCAT fixtures 6/6."),
            "phenotype_is_faithful_to_cpic": True,
            "is_core_marker_proxy": True,  # NOT a full PharmVar star-allele caller (core SNP set + sentinels)
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
