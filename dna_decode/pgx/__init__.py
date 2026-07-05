"""Human pharmacogenomics (PGx) decoders — the first HUMAN cell.

The honest catalog-tractable form of the "higher organism" jump: a deterministic
variant -> star-allele -> diplotype -> CPIC metabolizer-phenotype caller, on human DNA.
v0 = CYP2C19 (SNP-defined core alleles). NOT a clinical tool.

Sibling of the variant->catalog cells (hiv_amr / tb_amr) — NOT a learned model, NOT an embedding
(the frozen-genome-embedding thesis is a closed 0-for-4 negative). The bacterial/viral/fungal AMR
surfaces are untouched (this is a new, non-frozen package).
"""

# Single source of truth for the CLI-routable PGx genes — drives the dna-pgx --gene choices, the
# Evidence-Contract Registry manifest, and the coverage test (so adding a gene can't pass coverage
# vacuously). M1 fix, 2026-06-26.
PGX_GENES: tuple[str, ...] = ("cyp2c19", "cyp2c9", "cyp2c8", "vkorc1")
