"""Human HLA pharmacogenomics — drug-hypersensitivity HLA-allele carriage via validated TAG SNPs.

The new substrate after the CYP-cluster PGx cells (user-ratified 2026-07-06). The clinically-critical HLA
pharmacogenomic alleles (abacavir HSR, allopurinol / carbamazepine SJS-TEN) are each SNP-TAGGABLE via a
near-perfect LD proxy — exactly the deployed clinical screen (e.g. the abacavir rs2395029 B*57:01 test). So
a deterministic tag-SNP -> HLA-allele-carriage -> drug-risk cell IS the CPIC-aligned approach for the primary
use — the SAME curated-catalog deterministic regime as the AMR/PGx cells (NOT full sequence-based 6-digit
HLA typing, which needs read-based tools: HLA*LA / arcasHLA / OptiType / T1K).

LOAD-BEARING HONESTY (distinct from the PGx single-SNP cells): the tag SNP is an LD PROXY, NOT the allele
itself (contrast SLCO1B1 where rs4149056 IS the 521 call). So the caller's validity rests on the LD, and
must be validated against a REAL HLA truth set (the free published 1000G HLA types) — an independent
wrapper-vs-truth concordance, never a literature-asserted LD alone.
"""

# Single source of truth for the CLI-routable HLA alleles (drives dna-hla --allele + the registry + coverage).
# NARROWED to the VALIDATED tag after real 1000G-HLA-truth concordance (2026-07-06): only B*57:01/abacavir
# (rs2395029) cleared deployment (sens 0.979 / spec 0.992 / PPV 0.855). The provisional B*58:01 (rs9263726,
# sens 0.61 / PPV 0.18 — weak) + A*31:01 (rs1061235 — not paneled on 1000G, sens 0.0) FAILED validation and
# are DEMOTED to a documented negative (dna_decode.hla.catalog._UNVALIDATED_TAGS), NOT shipped as routable
# cells. See wiki/hla_validation_2026-07-06.md.
HLA_ALLELES: tuple[str, ...] = ("b5701",)
