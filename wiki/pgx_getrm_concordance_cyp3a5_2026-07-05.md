# CYP3A5 caller vs GeT-RM consensus on real 1000G (2026-07-05)

**Truth:** GeT-RM NGS consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) via the ursaPGx benchmark star-allele-comparison_common.tsv, column CYP3A5_getrm_cons
**Genotypes:** 1000 Genomes 30x phased panel (CYP3A5 region, pure-Python tabix-over-HTTP (scripts/fetch_1000g_region.py; no Docker); GeT-RM CDC CYP3A4/5 table)

- Overlap samples scored: **8**
- **Core-comparable diplotype concordance: 8/8 (1.0)**  (GeT-RM truth in *1/*3/*6/*7)
- Phenotype-correct incl. *38==*1: **8/8** (+0 *38 phenotype-equivalent samples)
- Correctly WITHHELD by sentinel: **0**
- **Genuine silent mis-call: 0/8 (0.0%)** -- non-core alleles beyond the v0 SNP set (+ sentinels where present); the honest residual blind spot.
- Correct-or-abstains: **8/8**

_GeT-RM CONSENSUS concordance on real 1000G genomes, independent caller. The strongest star-allele-CALLING validation tier available (vs the field's accepted consensus truth set). v0 covers the CORE SNP set; non-core-truth samples are scored separately (the v0.1 sentinel layer should WITHHOLD, not mis-call)._

## Core-comparable samples (GeT-RM truth in the v0 SNP set)

| sample | GeT-RM | predicted | match |
|---|---|---|---|
| HG00436 | *3/*3 | *3/*3 | OK |
| NA10856 | *1/*3 | *1/*3 | OK |
| NA18484 | *1/*7 | *1/*7 | OK |
| NA18518 | *1/*6 | *1/*6 | OK |
| NA18564 | *1/*1 | *1/*1 | OK |
| NA19143 | *6/*7 | *6/*7 | OK |
| NA19819 | *3/*6 | *3/*6 | OK |
| NA19920 | *7/*7 | *7/*7 | OK |
