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
HLA_ALLELES: tuple[str, ...] = ("b5701", "b5801", "a3101")
