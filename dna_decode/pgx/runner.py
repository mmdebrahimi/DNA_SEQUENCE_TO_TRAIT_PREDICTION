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


def call_cyp2c8(vcf: str | Path, sample_id: str | None = None,
                sample_column: str | None = None) -> dict:
    """Run the CYP2C8 caller (star-allele CALLING only) on a VCF -> full provenance record.

    CYP2C8 has NO CPIC metabolizer-phenotype system (function is substrate-dependent), so this record
    reports the star-allele diplotype + a per-allele PharmVar-clinical FUNCTION annotation, never a
    PM/IM/NM phenotype (has_cpic_phenotype=False). Star-allele calling is independently validatable vs
    the GeT-RM consensus (CYP2C8_getrm_ngs)."""
    from dna_decode.pgx import cyp2c8_catalog as c8
    res = call_diplotype(vcf, sample=sample_column, defining=c8.CORE_DEFINING, sentinels=c8.SENTINELS,
                         reference_allele=c8.REFERENCE_ALLELE, phenotype_fn=c8.diplotype_phenotype,
                         gene=c8.GENE)
    sid = sample_id or sample_column or Path(vcf).stem
    rec = {
        "sample_id": sid, "trait": "pgx_star_allele_diplotype", "gene": c8.GENE,
        "organism": "Homo sapiens", "assembly": c8.ASSEMBLY,
        "analysis_date": datetime.date.today().isoformat(), "schema": SCHEMA,
        "status": res.status, "phenotype_status": res.phenotype_status,
        "phenotype_confidence": res.phenotype_confidence,
        "diplotype": res.diplotype, "core_proxy_diplotype": res.core_proxy_diplotype,
        "allele1": res.allele1, "allele2": res.allele2,
        "has_cpic_phenotype": False,
        "function_annotation": res.phenotype,   # substrate-dependent function descriptor (NOT a phenotype)
        "phenotype": None,                       # explicit: CYP2C8 has no CPIC metabolizer phenotype
        "phenotype_abbrev": None,
        "alternate_diplotype": res.alternate_diplotype, "alternate_phenotype": res.alternate_phenotype,
        "sentinel_hits": res.sentinel_hits, "phasing": res.phasing, "flags": res.flags,
        "variant_calls": res.variant_calls,
        "caller": {
            "name": "dna_decode-pgx-cyp2c8-v0",
            "method": "vcf_core_snp_proxy -> star_allele -> diplotype (CALLING only; no CPIC phenotype)",
            "calling_independently_validatable": True,
            "independent_validation_status": (
                "GeT-RM consensus (CYP2C8_getrm_ngs): core diplotype concordance reported in "
                "wiki/pgx_getrm_concordance_cyp2c8_* (caller independent of the consensus tools); "
                "rare non-core CYP2C8 alleles mis-called *1 in v0 (no sentinel layer yet)."),
            "phenotype_is_faithful_to_cpic": False,   # no CPIC CYP2C8 metabolizer phenotype exists
            "has_cpic_phenotype": False,
            "is_core_marker_proxy": True,
            "reference_tool": "PharmCAT (star-allele calling; PharmCAT likewise emits no CYP2C8 metabolizer phenotype)",
        },
        "catalog": {
            "gene": c8.GENE, "core_alleles": [d.star for d in c8.CORE_DEFINING] + ["*1"],
            "allele_function": c8.ALLELE_FUNCTION,
            "defining_variants": [
                {"star": d.star, "rsid": d.rsid, "chrom": d.chrom, "pos": d.pos,
                 "ref": d.ref, "alt": d.alt, "cdna": d.cdna} for d in c8.CORE_DEFINING],
            "source": ("CYP2C8 *2/*3/*4 defining variants; GRCh38 coords VERIFIED via Ensembl REST; "
                       "star->rsID->function from the CYP2C8 pharmacogenetics review + dbSNP"),
        },
        "undetectable": c8.UNDETECTABLE,
        "caveat": ("CYP2C8 v0 = core SNP-defined *2/*3/*4 + *1, CALLING ONLY. CYP2C8 function is "
                   "substrate-dependent and has NO CPIC metabolizer-phenotype system -> this cell reports "
                   "the star-allele diplotype + per-allele PharmVar-clinical function, never a PM/IM/NM call. "
                   "Star-allele calling is independently validatable vs GeT-RM. NOT a clinical tool."),
    }
    if res.reason:
        rec["reason"] = res.reason
    return rec


def call_cyp3a5(vcf: str | Path, sample_id: str | None = None,
                sample_column: str | None = None) -> dict:
    """Run the CYP3A5 caller (tacrolimus; expressor/non-expressor phenotype) on a VCF -> provenance record.

    Function-pair CPIC phenotype (like CYP2C19): *1 = expressor (normal); *3/*6/*7 = no function. Validated
    against the REAL GeT-RM CDC multi-lab consensus (8/8 on the 1000G-overlapping samples; UNDERPOWERED)."""
    from dna_decode.pgx import cyp3a5_catalog as c3
    res = call_diplotype(vcf, sample=sample_column, defining=c3.CORE_DEFINING, sentinels=c3.SENTINELS,
                         reference_allele=c3.REFERENCE_ALLELE, phenotype_fn=c3.diplotype_phenotype,
                         gene=c3.GENE)
    sid = sample_id or sample_column or Path(vcf).stem
    rec = {
        "sample_id": sid, "trait": "pgx_metabolizer_phenotype", "gene": c3.GENE,
        "organism": "Homo sapiens", "assembly": c3.ASSEMBLY,
        "analysis_date": datetime.date.today().isoformat(), "schema": SCHEMA,
        "status": res.status, "phenotype_status": res.phenotype_status,
        "phenotype_confidence": res.phenotype_confidence,
        "diplotype": res.diplotype, "core_proxy_diplotype": res.core_proxy_diplotype,
        "allele1": res.allele1, "allele2": res.allele2,
        "phenotype": res.phenotype,
        "phenotype_abbrev": c3.PHENOTYPE_ABBREV.get(res.phenotype or "", None),
        "alternate_diplotype": res.alternate_diplotype, "alternate_phenotype": res.alternate_phenotype,
        "sentinel_hits": res.sentinel_hits, "phasing": res.phasing, "flags": res.flags,
        "variant_calls": res.variant_calls,
        "caller": {
            "name": "dna_decode-pgx-cyp3a5-v0",
            "method": "vcf_core_variant_proxy -> star_allele -> diplotype -> CPIC_expressor_phenotype",
            "calling_independently_validatable": True,
            "independent_validation_status": (
                "GeT-RM CDC multi-lab consensus (CYP3A5_getrm_cons): 8/8 core diplotype concordance on the "
                "1000G-overlapping samples covering *1/*3/*6/*7 incl. the *7 insertion + *6/*7 non-expressor "
                "cases (caller independent of the labs). UNDERPOWERED (n=8). Rare non-core alleles mis-called *1."),
            "phenotype_is_faithful_to_cpic": True, "is_core_marker_proxy": True,
            "reference_tool": "PharmCAT",
        },
        "catalog": {
            "gene": c3.GENE, "core_alleles": [d.star for d in c3.CORE_DEFINING] + ["*1"],
            "allele_function": c3.ALLELE_FUNCTION,
            "defining_variants": [
                {"star": d.star, "rsid": d.rsid, "chrom": d.chrom, "pos": d.pos,
                 "ref": d.ref, "alt": d.alt, "cdna": d.cdna} for d in c3.CORE_DEFINING],
            "source": ("CYP3A5 *3/*6/*7 defining variants; GRCh38 coords VERIFIED via Ensembl REST + "
                       "AF-confirmed on 1000G; CPIC tacrolimus guideline (Birdwell 2015)"),
        },
        "undetectable": c3.UNDETECTABLE,
        "caveat": ("CYP3A5 v0 = core *3/*6/*7 (no-function) + *1 (expressor) with CPIC expressor/non-expressor "
                   "phenotype (NM=expressor / IM / PM=non-expressor). Star-allele calling independently "
                   "validatable vs GeT-RM (8/8, UNDERPOWERED). Phenotype faithful-to-CPIC. NO sentinel layer -> "
                   "rare non-core alleles mis-called *1. Tacrolimus-relevant. NOT a clinical tool."),
    }
    if res.reason:
        rec["reason"] = res.reason
    return rec


def call_tpmt(vcf: str | Path, sample_id: str | None = None,
              sample_column: str | None = None) -> dict:
    """Run the TPMT COMPOUND caller (thiopurines) on a VCF -> provenance record.

    First true compound-allele cell: *3A is resolved from BOTH rs1800460 (*3B) + rs1142345 (*3C) on the
    same haplotype (each SNP alone = *3B / *3C). Validated 85/85 vs the GeT-RM CDC consolidated consensus."""
    from dna_decode.pgx import tpmt_catalog as tp
    from dna_decode.pgx.compound_caller import assemble_compound_diplotype
    res = assemble_compound_diplotype(vcf, tp.COMPONENTS, tp.COMPOUND_RULES,
                                      reference_allele=tp.REFERENCE_ALLELE,
                                      phenotype_fn=tp.diplotype_phenotype, gene=tp.GENE, sample=sample_column)
    sid = sample_id or sample_column or Path(vcf).stem
    rec = {
        "sample_id": sid, "trait": "pgx_metabolizer_phenotype", "gene": tp.GENE,
        "organism": "Homo sapiens", "assembly": tp.ASSEMBLY,
        "analysis_date": datetime.date.today().isoformat(), "schema": SCHEMA,
        "status": res.status, "phenotype_status": res.phenotype_status,
        "phenotype_confidence": res.phenotype_confidence,
        "diplotype": res.diplotype, "core_proxy_diplotype": res.core_proxy_diplotype,
        "allele1": res.allele1, "allele2": res.allele2,
        "phenotype": res.phenotype,
        "phenotype_abbrev": tp.PHENOTYPE_ABBREV.get(res.phenotype or "", None),
        "alternate_diplotype": res.alternate_diplotype, "alternate_phenotype": res.alternate_phenotype,
        "sentinel_hits": res.sentinel_hits, "phasing": res.phasing, "flags": res.flags,
        "variant_calls": res.variant_calls,
        "caller": {
            "name": "dna_decode-pgx-tpmt-v0",
            "method": "vcf_compound_haplotype -> star_allele (*3A=*3B+*3C in cis) -> CPIC_phenotype",
            "calling_independently_validatable": True,
            "independent_validation_status": (
                "GeT-RM CDC consolidated consensus (TPMT): 85/85 core-comparable diplotype concordance on "
                "1000G-overlap (truth in *1/*3A/*3B/*3C); the compound *3A path is exercised (6 *3A + 8 *3C "
                "truth samples). Rare non-core alleles (*2/*8/*16...) mis-called *1."),
            "phenotype_is_faithful_to_cpic": True, "is_core_marker_proxy": True,
            "is_compound_allele_caller": True,
            "reference_tool": "PharmCAT",
        },
        "catalog": {
            "gene": tp.GENE, "core_alleles": ["*1", "*3A", "*3B", "*3C"],
            "allele_function": tp.ALLELE_FUNCTION,
            "component_variants": [
                {"tag": d.star, "rsid": d.rsid, "chrom": d.chrom, "pos": d.pos,
                 "ref": d.ref, "alt": d.alt, "cdna": d.cdna} for d in tp.COMPONENTS],
            "compound_rules": {r.star: sorted(r.components) for r in tp.COMPOUND_RULES},
            "source": ("TPMT *3A/*3B/*3C compound; GRCh38 coords VERIFIED via Ensembl + AF-confirmed on "
                       "1000G; CPIC thiopurine guideline (Relling 2019)"),
        },
        "undetectable": tp.UNDETECTABLE,
        "caveat": ("TPMT v0 = COMPOUND *3A (=*3B+*3C in cis) / *3B / *3C + *1, with CPIC thiopurine "
                   "phenotype. First true compound-allele cell (two SNPs on one haplotype -> *3A). Star "
                   "calling independently validatable vs GeT-RM (85/85). Phenotype faithful-to-CPIC. Rare "
                   "non-core alleles mis-called *1 (no sentinel layer). NOT a clinical tool."),
    }
    if res.reason:
        rec["reason"] = res.reason
    return rec


def call_cyp2b6(vcf: str | Path, sample_id: str | None = None,
                sample_column: str | None = None) -> dict:
    """Run the CYP2B6 caller (efavirenz; *6-proxy from the 516G>T signal) on a VCF -> provenance record.

    v0 is a SINGLE-SNP *6-proxy: rs2279343 (785A>G, the 2nd *6 component) is absent from the 1000G 30x
    panel, so this cannot split *6 from *9 (documented). Validated 62/62 on clean *1/*6 truth samples."""
    from dna_decode.pgx import cyp2b6_catalog as c6
    res = call_diplotype(vcf, sample=sample_column, defining=c6.CORE_DEFINING, sentinels=c6.SENTINELS,
                         reference_allele=c6.REFERENCE_ALLELE, phenotype_fn=c6.diplotype_phenotype,
                         gene=c6.GENE)
    sid = sample_id or sample_column or Path(vcf).stem
    rec = {
        "sample_id": sid, "trait": "pgx_metabolizer_phenotype", "gene": c6.GENE,
        "organism": "Homo sapiens", "assembly": c6.ASSEMBLY,
        "analysis_date": datetime.date.today().isoformat(), "schema": SCHEMA,
        "status": res.status, "phenotype_status": res.phenotype_status,
        "phenotype_confidence": res.phenotype_confidence,
        "diplotype": res.diplotype, "core_proxy_diplotype": res.core_proxy_diplotype,
        "allele1": res.allele1, "allele2": res.allele2,
        "phenotype": res.phenotype,
        "phenotype_abbrev": c6.PHENOTYPE_ABBREV.get(res.phenotype or "", None),
        "alternate_diplotype": res.alternate_diplotype, "alternate_phenotype": res.alternate_phenotype,
        "sentinel_hits": res.sentinel_hits, "phasing": res.phasing, "flags": res.flags,
        "variant_calls": res.variant_calls,
        "caller": {
            "name": "dna_decode-pgx-cyp2b6-v0",
            "method": "vcf_single_snp_proxy(516G>T) -> *6-proxy -> CPIC_phenotype",
            "calling_independently_validatable": True,
            "independent_validation_status": (
                "GeT-RM CDC consolidated consensus (CYP2B6): 62/62 on clean *1/*6 truth samples on "
                "1000G-overlap. SINGLE-SNP *6-proxy (516G>T) — cannot split *6 from *9 (rs2279343/785A>G "
                "absent from the 1000G 30x panel); *4/other non-core mis-called. Documented residual."),
            "phenotype_is_faithful_to_cpic": True, "is_core_marker_proxy": True,
            "is_single_snp_proxy": True,
            "reference_tool": "PharmCAT",
        },
        "catalog": {
            "gene": c6.GENE, "core_alleles": ["*1", "*6"],
            "allele_function": c6.ALLELE_FUNCTION,
            "defining_variants": [
                {"star": d.star, "rsid": d.rsid, "chrom": d.chrom, "pos": d.pos,
                 "ref": d.ref, "alt": d.alt, "cdna": d.cdna} for d in c6.CORE_DEFINING],
            "source": ("CYP2B6 516G>T (rs3745274) *6-signal; GRCh38 coord VERIFIED via Ensembl + "
                       "AF-confirmed on 1000G; CPIC efavirenz guideline (Desta 2019)"),
        },
        "undetectable": c6.UNDETECTABLE,
        "caveat": ("CYP2B6 v0 = *6-PROXY from the 516G>T signal (rs3745274) + *1, CPIC efavirenz phenotype. "
                   "SINGLE-SNP proxy: rs2279343 (785A>G, the 2nd *6 component) is absent from the 1000G 30x "
                   "panel, so *6 cannot be split from the rare *9 (documented; a callset with 785 upgrades "
                   "to the compound resolver as in TPMT). Validated 62/62 on clean *1/*6. NOT a clinical tool."),
    }
    if res.reason:
        rec["reason"] = res.reason
    return rec
