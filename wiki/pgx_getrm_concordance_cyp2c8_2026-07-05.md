# CYP2C8 caller vs GeT-RM consensus on real 1000G (2026-07-05)

**Truth:** GeT-RM NGS consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) via the ursaPGx benchmark star-allele-comparison_common.tsv, column CYP2C8_getrm_ngs
**Genotypes:** 1000 Genomes 30x phased panel (CYP2C8 region, pure-Python tabix-over-HTTP (scripts/fetch_1000g_region.py; no Docker))

- Overlap samples scored: **87**
- **Core-comparable diplotype concordance: 82/82 (1.0)**  (GeT-RM truth in *1/*2/*3/*4)
- Phenotype-correct incl. *38==*1: **82/87** (+0 *38 phenotype-equivalent samples)
- Correctly WITHHELD by sentinel: **0**
- **Genuine silent mis-call: 5/87 (5.7%)** -- non-core alleles beyond the v0 SNP set (+ sentinels where present); the honest residual blind spot.
- Correct-or-abstains: **82/87**

_GeT-RM CONSENSUS concordance on real 1000G genomes, independent caller. The strongest star-allele-CALLING validation tier available (vs the field's accepted consensus truth set). v0 covers the CORE SNP set; non-core-truth samples are scored separately (the v0.1 sentinel layer should WITHHOLD, not mis-call)._

## Core-comparable samples (GeT-RM truth in the v0 SNP set)

| sample | GeT-RM | predicted | match |
|---|---|---|---|
| HG00276 | *1/*3 | *1/*3 | OK |
| HG00436 | *1/*1 | *1/*1 | OK |
| HG00589 | *1/*1 | *1/*1 | OK |
| HG01190 | *1/*3 | *1/*3 | OK |
| NA06991 | *1/*1 | *1/*1 | OK |
| NA07000 | *1/*1 | *1/*1 | OK |
| NA07019 | *1/*1 | *1/*1 | OK |
| NA07029 | *1/*3 | *1/*3 | OK |
| NA07055 | *1/*1 | *1/*1 | OK |
| NA07056 | *1/*1 | *1/*1 | OK |
| NA07348 | *1/*4 | *1/*4 | OK |
| NA07357 | *1/*4 | *1/*4 | OK |
| NA10831 | *1/*1 | *1/*1 | OK |
| NA10838 | *1/*1 | *1/*1 | OK |
| NA10846 | *1/*1 | *1/*1 | OK |
| NA10847 | *1/*1 | *1/*1 | OK |
| NA10851 | *1/*1 | *1/*1 | OK |
| NA10854 | *3/*3 | *3/*3 | OK |
| NA10855 | *1/*3 | *1/*3 | OK |
| NA10856 | *1/*1 | *1/*1 | OK |
| NA10859 | *1/*1 | *1/*1 | OK |
| NA10865 | *1/*3 | *1/*3 | OK |
| NA11832 | *1/*1 | *1/*1 | OK |
| NA11839 | *1/*3 | *1/*3 | OK |
| NA11881 | *1/*1 | *1/*1 | OK |
| NA11993 | *1/*1 | *1/*1 | OK |
| NA12003 | *1/*3 | *1/*3 | OK |
| NA12006 | *1/*1 | *1/*1 | OK |
| NA12145 | *1/*4 | *1/*4 | OK |
| NA12236 | *1/*1 | *1/*1 | OK |
| NA12336 | *1/*1 | *1/*1 | OK |
| NA12717 | *1/*1 | *1/*1 | OK |
| NA12753 | *1/*4 | *1/*4 | OK |
| NA12813 | *1/*1 | *1/*1 | OK |
| NA12815 | *1/*1 | *1/*1 | OK |
| NA12873 | *1/*1 | *1/*1 | OK |
| NA12878 | *1/*3 | *1/*3 | OK |
| NA12892 | *1/*3 | *1/*3 | OK |
| NA18484 | *1/*4 | *1/*4 | OK |
| NA18518 | *1/*2 | *1/*2 | OK |
| NA18519 | *1/*2 | *1/*2 | OK |
| NA18526 | *1/*1 | *1/*1 | OK |
| NA18544 | *1/*1 | *1/*1 | OK |
| NA18552 | *1/*1 | *1/*1 | OK |
| NA18563 | *1/*1 | *1/*1 | OK |
| NA18564 | *1/*1 | *1/*1 | OK |
| NA18565 | *1/*1 | *1/*1 | OK |
| NA18572 | *1/*1 | *1/*1 | OK |
| NA18617 | *1/*1 | *1/*1 | OK |
| NA18855 | *1/*1 | *1/*1 | OK |
| NA18861 | *1/*1 | *1/*1 | OK |
| NA18868 | *1/*1 | *1/*1 | OK |
| NA18873 | *1/*1 | *1/*1 | OK |
| NA18942 | *1/*1 | *1/*1 | OK |
| NA18945 | *1/*1 | *1/*1 | OK |
| NA18952 | *1/*1 | *1/*1 | OK |
| NA18959 | *1/*1 | *1/*1 | OK |
| NA18966 | *1/*1 | *1/*1 | OK |
| NA18973 | *1/*1 | *1/*1 | OK |
| NA18980 | *1/*1 | *1/*1 | OK |
| NA18992 | *1/*1 | *1/*1 | OK |
| NA19003 | *1/*1 | *1/*1 | OK |
| NA19007 | *1/*1 | *1/*1 | OK |
| NA19035 | *2/*2 | *2/*2 | OK |
| NA19095 | *1/*2 | *1/*2 | OK |
| NA19109 | *2/*2 | *2/*2 | OK |
| NA19122 | *1/*1 | *1/*1 | OK |
| NA19147 | *1/*1 | *1/*1 | OK |
| NA19174 | *1/*1 | *1/*1 | OK |
| NA19176 | *2/*2 | *2/*2 | OK |
| NA19207 | *1/*2 | *1/*2 | OK |
| NA19226 | *1/*1 | *1/*1 | OK |
| NA19238 | *1/*1 | *1/*1 | OK |
| NA19239 | *1/*2 | *1/*2 | OK |
| NA19700 | *1/*1 | *1/*1 | OK |
| NA19785 | *1/*1 | *1/*1 | OK |
| NA19789 | *1/*3 | *1/*3 | OK |
| NA19819 | *1/*2 | *1/*2 | OK |
| NA19908 | *1/*1 | *1/*1 | OK |
| NA19920 | *1/*1 | *1/*1 | OK |
| NA20296 | *1/*1 | *1/*1 | OK |
| NA20509 | *1/*4 | *1/*4 | OK |

## Non-core-truth samples (v0.1 sentinel SHOULD withhold)

| sample | GeT-RM | core-proxy | phenotype_status |
|---|---|---|---|
| NA07048 | *1/*18 | *1/*1 | ok |
| NA12156 | *1/*15 | *1/*1 | ok |
| NA19143 | *1/*17 | *1/*1 | ok |
| NA19213 | *1/*17 | *1/*1 | ok |
| NA19917 | *1/*16 | *1/*1 | ok |
