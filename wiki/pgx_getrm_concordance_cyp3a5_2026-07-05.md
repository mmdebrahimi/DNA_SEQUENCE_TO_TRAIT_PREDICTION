# CYP3A5 caller vs GeT-RM consensus on real 1000G (2026-07-29)

**Truth:** GeT-RM NGS consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) via the ursaPGx benchmark star-allele-comparison_common.tsv, column CYP3A5_getrm_cons
**Genotypes:** 1000 Genomes 30x phased panel (CYP3A5 region, pure-Python tabix-over-HTTP (scripts/fetch_1000g_region.py; no Docker); GeT-RM CDC CYP3A4/5 table)

- Overlap samples scored: **88**
- **Core-comparable diplotype concordance: 88/88 (1.0)**  (GeT-RM truth in *1/*3/*6/*7)
- Phenotype-correct incl. *38==*1: **88/88** (+0 *38 phenotype-equivalent samples)
- Correctly WITHHELD by sentinel: **0**
- **Genuine silent mis-call: 0/88 (0.0%)** -- non-core alleles beyond the v0 SNP set (+ sentinels where present); the honest residual blind spot.
- Correct-or-abstains: **88/88**

_GeT-RM CONSENSUS concordance on real 1000G genomes, independent caller. The strongest star-allele-CALLING validation tier available (vs the field's accepted consensus truth set). v0 covers the CORE SNP set; non-core-truth samples are scored separately (the v0.1 sentinel layer should WITHHOLD, not mis-call)._

## Core-comparable samples (GeT-RM truth in the v0 SNP set)

| sample | GeT-RM | predicted | match |
|---|---|---|---|
| HG00276 | *3/*3 | *3/*3 | OK |
| HG00436 | *3/*3 | *3/*3 | OK |
| HG00589 | *3/*3 | *3/*3 | OK |
| HG01190 | *1/*1 | *1/*1 | OK |
| NA06991 | *3/*3 | *3/*3 | OK |
| NA06993 | *3/*3 | *3/*3 | OK |
| NA07000 | *1/*3 | *1/*3 | OK |
| NA07019 | *3/*3 | *3/*3 | OK |
| NA07029 | *1/*3 | *1/*3 | OK |
| NA07048 | *3/*3 | *3/*3 | OK |
| NA07055 | *3/*3 | *3/*3 | OK |
| NA07056 | *3/*3 | *3/*3 | OK |
| NA07348 | *3/*3 | *3/*3 | OK |
| NA07357 | *3/*3 | *3/*3 | OK |
| NA10831 | *3/*3 | *3/*3 | OK |
| NA10838 | *3/*3 | *3/*3 | OK |
| NA10846 | *3/*3 | *3/*3 | OK |
| NA10847 | *3/*3 | *3/*3 | OK |
| NA10851 | *3/*3 | *3/*3 | OK |
| NA10854 | *1/*3 | *1/*3 | OK |
| NA10855 | *3/*3 | *3/*3 | OK |
| NA10856 | *1/*3 | *1/*3 | OK |
| NA10859 | *3/*3 | *3/*3 | OK |
| NA10865 | *3/*3 | *3/*3 | OK |
| NA11832 | *3/*3 | *3/*3 | OK |
| NA11839 | *1/*3 | *1/*3 | OK |
| NA11881 | *3/*3 | *3/*3 | OK |
| NA11993 | *3/*3 | *3/*3 | OK |
| NA12003 | *1/*3 | *1/*3 | OK |
| NA12006 | *3/*3 | *3/*3 | OK |
| NA12145 | *3/*3 | *3/*3 | OK |
| NA12156 | *3/*3 | *3/*3 | OK |
| NA12236 | *3/*3 | *3/*3 | OK |
| NA12336 | *3/*3 | *3/*3 | OK |
| NA12717 | *1/*3 | *1/*3 | OK |
| NA12753 | *3/*3 | *3/*3 | OK |
| NA12813 | *3/*3 | *3/*3 | OK |
| NA12815 | *3/*3 | *3/*3 | OK |
| NA12873 | *3/*3 | *3/*3 | OK |
| NA12878 | *3/*3 | *3/*3 | OK |
| NA12892 | *3/*3 | *3/*3 | OK |
| NA18484 | *1/*7 | *1/*7 | OK |
| NA18518 | *1/*6 | *1/*6 | OK |
| NA18519 | *1/*6 | *1/*6 | OK |
| NA18526 | *1/*1 | *1/*1 | OK |
| NA18544 | *1/*3 | *1/*3 | OK |
| NA18552 | *3/*3 | *3/*3 | OK |
| NA18563 | *1/*1 | *1/*1 | OK |
| NA18564 | *1/*1 | *1/*1 | OK |
| NA18565 | *1/*3 | *1/*3 | OK |
| NA18572 | *1/*3 | *1/*3 | OK |
| NA18617 | *3/*3 | *3/*3 | OK |
| NA18855 | *3/*6 | *3/*6 | OK |
| NA18861 | *1/*1 | *1/*1 | OK |
| NA18868 | *1/*3 | *1/*3 | OK |
| NA18873 | *1/*1 | *1/*1 | OK |
| NA18942 | *3/*3 | *3/*3 | OK |
| NA18945 | *1/*3 | *1/*3 | OK |
| NA18952 | *3/*3 | *3/*3 | OK |
| NA18959 | *1/*3 | *1/*3 | OK |
| NA18966 | *1/*3 | *1/*3 | OK |
| NA18973 | *1/*3 | *1/*3 | OK |
| NA18980 | *1/*3 | *1/*3 | OK |
| NA18992 | *3/*3 | *3/*3 | OK |
| NA19003 | *3/*3 | *3/*3 | OK |
| NA19007 | *3/*3 | *3/*3 | OK |
| NA19035 | *1/*7 | *1/*7 | OK |
| NA19095 | *1/*3 | *1/*3 | OK |
| NA19109 | *1/*3 | *1/*3 | OK |
| NA19122 | *1/*1 | *1/*1 | OK |
| NA19143 | *6/*7 | *6/*7 | OK |
| NA19147 | *1/*3 | *1/*3 | OK |
| NA19174 | *1/*6 | *1/*6 | OK |
| NA19176 | *1/*3 | *1/*3 | OK |
| NA19207 | *3/*7 | *3/*7 | OK |
| NA19213 | *1/*6 | *1/*6 | OK |
| NA19226 | *1/*6 | *1/*6 | OK |
| NA19238 | *1/*1 | *1/*1 | OK |
| NA19239 | *1/*1 | *1/*1 | OK |
| NA19700 | *1/*3 | *1/*3 | OK |
| NA19785 | *1/*3 | *1/*3 | OK |
| NA19789 | *3/*3 | *3/*3 | OK |
| NA19819 | *3/*6 | *3/*6 | OK |
| NA19908 | *1/*3 | *1/*3 | OK |
| NA19917 | *1/*7 | *1/*7 | OK |
| NA19920 | *7/*7 | *7/*7 | OK |
| NA20296 | *1/*6 | *1/*6 | OK |
| NA20509 | *3/*3 | *3/*3 | OK |
