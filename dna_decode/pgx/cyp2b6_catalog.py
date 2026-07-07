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

from dna_decode.pgx.cyp2c19_catalog import DefiningVariant

GENE = "CYP2B6"
ASSEMBLY = "GRCh38"
REFERENCE_ALLELE = "*1"

# Single-SNP *6-proxy: the 516G>T signal. (rs2279343 785A>G would be the 2nd component, but it is absent
# from the 1000G 30x panel — see module docstring; v0 is 516-only.)
CORE_DEFINING: list[DefiningVariant] = [
    DefiningVariant("*6", "rs3745274", "19", 41006936, "G", "T", "c.516G>T (*6/*9 signal)"),
]
SENTINELS: list = []

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
