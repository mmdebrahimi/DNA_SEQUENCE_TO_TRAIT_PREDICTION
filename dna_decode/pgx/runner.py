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


def call_cyp2c9(vcf: str | Path, sample_id: str | None = None,
                sample_column: str | None = None) -> dict:
    """Run the CYP2C9 caller (activity-score phenotype) on a VCF -> full provenance record."""
    from dna_decode.pgx import cyp2c9_catalog as c9
    res = call_diplotype(vcf, sample=sample_column, defining=c9.CORE_DEFINING, sentinels=c9.SENTINELS,
                         reference_allele=c9.REFERENCE_ALLELE, phenotype_fn=c9.diplotype_phenotype,
                         gene=c9.GENE)
    sid = sample_id or sample_column or Path(vcf).stem
    a_s = c9.activity_score(res.allele1, res.allele2) if res.allele1 else None
    rec = {
        "sample_id": sid, "trait": "pgx_metabolizer_phenotype", "gene": c9.GENE,
        "organism": "Homo sapiens", "assembly": c9.ASSEMBLY,
        "analysis_date": datetime.date.today().isoformat(), "schema": SCHEMA,
        "status": res.status, "phenotype_status": res.phenotype_status,
        "phenotype_confidence": res.phenotype_confidence,
        "diplotype": res.diplotype, "core_proxy_diplotype": res.core_proxy_diplotype,
        "allele1": res.allele1, "allele2": res.allele2,
        "activity_score": a_s,
        "phenotype": res.phenotype,
        "phenotype_abbrev": c9.PHENOTYPE_ABBREV.get(res.phenotype or "", None),
        "alternate_diplotype": res.alternate_diplotype, "alternate_phenotype": res.alternate_phenotype,
        "sentinel_hits": res.sentinel_hits, "phasing": res.phasing, "flags": res.flags,
        "variant_calls": res.variant_calls,
        "caller": {
            "name": "dna_decode-pgx-cyp2c9-v0",
            "method": "vcf_core_snp_proxy -> star_allele -> diplotype -> CPIC_activity_score -> phenotype",
            "calling_independently_validatable": True,
            "independent_validation_status": (
                "GeT-RM consensus: core diplotype concordance reported in "
                "wiki/pgx_getrm_concordance_cyp2c9_2026-06-25 (caller independent of the consensus tools); "
                "non-core *5/*6/*8/*9/*11 mis-called *1 in v0 (sentinel layer = v0.1)."),
            "phenotype_is_faithful_to_cpic": True, "is_core_marker_proxy": True,
            "reference_tool": "PharmCAT",
        },
        "catalog": {
            "gene": c9.GENE, "core_alleles": [d.star for d in c9.CORE_DEFINING] + ["*1"],
            "activity_values": c9.ACTIVITY_VALUE,
            "defining_variants": [
                {"star": d.star, "rsid": d.rsid, "chrom": d.chrom, "pos": d.pos,
                 "ref": d.ref, "alt": d.alt, "cdna": d.cdna} for d in c9.CORE_DEFINING],
            "source": "PharmVar CYP2C9 + CPIC warfarin guideline (Johnson 2017); dbSNP coords",
        },
        "undetectable": c9.UNDETECTABLE,
        "caveat": ("CYP2C9 v0 = core SNP-defined *2/*3 + *1 with CPIC ACTIVITY-SCORE phenotype "
                   "(*1=1.0/*2=0.5/*3=0.0; AS 2=NM, 1-1.5=IM, 0-0.5=PM). NO sentinel layer yet -> a non-core "
                   "allele (*5/*6/*8/*9/*11) is mis-called *1 (documented residual; sentinels = v0.1). "
                   "Phenotype faithful-to-CPIC. NOT a clinical tool."),
    }
    if res.reason:
        rec["reason"] = res.reason
    return rec
