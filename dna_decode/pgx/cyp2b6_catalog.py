"""CYP2B6 curated catalog — *6 detection via the 516G>T signal (rs3745274), v0 (efavirenz).

Sixth CYP-cluster-independent PGx gene (chr19q13.2). Efavirenz (and other) metabolism. CYP2B6*6 is the
dominant reduced-function allele.

IMPORTANT DATA-DRIVEN SCOPE (honest, load-bearing):
  CYP2B6*6 is a COMPOUND of two SNPs — rs3745274 (516G>T) + rs2279343 (785A>G) — and 516G>T ALONE is *9,
  785A>G ALONE is *4. The TRUE *6-vs-*9 discrimination needs BOTH SNPs. BUT rs2279343 (785A>G) is ABSENT
  from the 1000G 30x NYGC phased panel (20220422 callset) at chr19:41009358 (empirically verified — the
  callset has records at 41009350/351/368 but not 358). So v0 detects *6 from the 516G>T signal ALONE
  (rs3745274), which is the primary reduced-function/splicing variant. This CANNOT split *6 from the rare
  *9 (documented residual); it is a single-SNP *6-PROXY, not the full 2-SNP compound. When a callset
  carrying rs2279343 is used, the compound_caller path (as in TPMT) upgrades this to a true *6/*9/*4
  resolver (v0.1).

PROVENANCE (grounded, NO fabrication):
  * rs3745274 c.516G>T (p.Gln172His) GRCh38 chr19:41006936 G>T, VERIFIED via Ensembl REST; ALT (T)
    empirically freq 0.320 on 1000G (matches *6+*9 combined global frequency). CYP2B6 is plus-strand.
  * FUNCTION + phenotype: CPIC (Desta 2019 efavirenz guideline) — *1 normal; *6 decreased function ->
    Intermediate/Poor Metabolizer.

HONESTY TIER: star-allele CALLING validatable vs the GeT-RM CDC consolidated consensus on the *6-decodable
subset (truth in *1/*6; validated 62/62 clean *1/*6 samples on 1000G-overlap). caller_is_independent_baseline
=True for the *6-proxy. Phenotype FAITHFUL-TO-CPIC. Reference tool: PharmCAT. This is a SINGLE-SNP proxy
(cannot split *6/*9 without rs2279343) — tiered accordingly.

SCOPE (v0): CORE *6-proxy (516G>T) + *1. *9 (516 alone, would be mis-labelled *6), *4 (785 alone, absent
from callset -> mis-called *1), and *2/*5/*18/*22/*27... are non-core -> documented residual. NOT clinical.
"""
from __future__ import annotations

from dna_decode.pgx.cyp2c19_catalog import DefiningVariant, SentinelVariant

GENE = "CYP2B6"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

# Single-SNP *6-proxy: the 516G>T signal. (rs2279343 785A>G would be the 2nd component, but it is absent
# from the 1000G 30x panel — see module docstring; v0 is 516-only.)
CORE_DEFINING: list[DefiningVariant] = [
    DefiningVariant("*6", "rs3745274", "19", 41006936, "G", "T", "c.516G>T (*6/*9 signal)"),
]
# Non-core-allele SENTINELS: a proven non-core CYP2B6 allele the *1/*6 proxy cannot resolve WITHHOLDS the
# phenotype rather than a silent *1/*6 mis-call. Each row is the allele's DISTINCTIVE SNP (sourced VERBATIM
# from PharmCAT/CLINPGX CYP2B6_translation.json, GRCh38.p14, chr19; each (rsid->pos,ref,alt) Ensembl-verified
# 3/3 OK, 2026-07-28). CRITICAL: the shared *6-haplotype SNP rs2279343 (785A>G) is DELIBERATELY EXCLUDED --
# a sentinel there would false-withhold a valid core *6 call (785G rides on *6). Each row uses the allele's
# OWN distinctive site instead (*7=rs3211371/1459, *18=rs28399499/983, *2=rs8192709/64), none shared with a
# core allele -> accounted_by_core stays None. SCOPE GAP (documented): *4 (785G-alone) and *9 (516T-alone)
# are ABSENCE-defined (distinguished from *6 only by a MISSING partner SNP) -> not presence-sentinel-able;
# they are NOT in the observed GeT-RM leak (which was *18/*7/*2). Exact ALT -> no benign-variant false-withhold.
# VALIDATION NOTE: the distinctive SNP sites span 40991369-41016810, WIDER than the original narrow *6-proxy
# 1000G fetch (41005049-41010995) -> re-fetch a wider region to exercise them:
#   uv run python scripts/fetch_1000g_region.py --chrom chr19 --start 40991000 --end 41017000 --out data/pgx_1000g/cyp2b6_1000g.vcf
# Then `pgx_getrm_concordance --gene cyp2b6` = core 62/62 UNCHANGED + 18 non-core WITHHELD (2026-07-28).
SENTINELS: list[SentinelVariant] = [
    SentinelVariant("rs8192709", "19", 40991369, "C", "T", "*2", "CYP2B6*2 (c.64C>T, R22C) non-core"),
    SentinelVariant("rs3211371", "19", 41016810, "C", "T", "*7", "CYP2B6*7 distinctive SNP (c.1459C>T, R487C); the *6 components 516/785 are shared/excluded"),
    SentinelVariant("rs28399499", "19", 41012316, "T", "C", "*18", "CYP2B6*18 distinctive SNP (c.983T>C, I328T); the shared 785 is excluded"),
]

ALLELE_FUNCTION: dict[str, str] = {
    "*1": "normal",
    "*6": "decreased",
}

_PHENOTYPE_BY_FUNCTION_PAIR: dict[tuple[str, str], str] = {
    ("normal", "normal"): "Normal Metabolizer",
    ("decreased", "normal"): "Intermediate Metabolizer",
    ("decreased", "decreased"): "Poor Metabolizer",
}

PHENOTYPE_ABBREV: dict[str, str] = {
    "Normal Metabolizer": "NM",
    "Intermediate Metabolizer": "IM",
    "Poor Metabolizer": "PM",
    "Indeterminate": "IND",
}

UNDETECTABLE = sorted({
    "star9_vs_star6_needs_785",        # *9 (516 alone) indistinguishable from *6 without rs2279343
    "star4_785_absent_from_callset",   # *4 (785 alone) absent from the 1000G 30x panel -> mis-called *1
    "non_core_star_allele",            # *2/*5/*18/*22/*27... -> mis-called *1
    "cnv_or_gene_deletion",
})


def diplotype_phenotype(allele1: str, allele2: str) -> str:
    """CPIC CYP2B6 metabolizer phenotype (*6-proxy). 'Indeterminate' if either allele unknown."""
    f1 = ALLELE_FUNCTION.get(allele1, "unknown")
    f2 = ALLELE_FUNCTION.get(allele2, "unknown")
    if "unknown" in (f1, f2):
        return "Indeterminate"
    return _PHENOTYPE_BY_FUNCTION_PAIR[tuple(sorted((f1, f2)))]
