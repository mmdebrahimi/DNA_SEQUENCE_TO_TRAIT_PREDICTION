# TPMT caller vs GeT-RM consensus on real 1000G (2026-07-05)

**Truth:** GeT-RM NGS consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) via the ursaPGx benchmark star-allele-comparison_common.tsv, column TPMT_getrm_cons
**Genotypes:** 1000 Genomes 30x phased panel (TPMT region, pure-Python tabix-over-HTTP (no Docker); GeT-RM CDC consolidated PGx table; COMPOUND *3A caller)

- Overlap samples scored: **98**
- **Core-comparable diplotype concordance: 85/85 (1.0)**  (GeT-RM truth in *1/*3A/*3B/*3C)
- Phenotype-correct incl. *38==*1: **85/98** (+0 *38 phenotype-equivalent samples)
- Correctly WITHHELD by sentinel: **0**
- **Genuine silent mis-call: 13/98 (13.3%)** -- non-core alleles beyond the v0 SNP set (+ sentinels where present); the honest residual blind spot.
- Correct-or-abstains: **85/98**

_GeT-RM CONSENSUS concordance on real 1000G genomes, independent caller. The strongest star-allele-CALLING validation tier available (vs the field's accepted consensus truth set). v0 covers the CORE SNP set; non-core-truth samples are scored separately (the v0.1 sentinel layer should WITHHOLD, not mis-call)._

## Core-comparable samples (GeT-RM truth in the v0 SNP set)

| sample | GeT-RM | predicted | match |
|---|---|---|---|
| HG00436 | *1/*1 | *1/*1 | OK |
| HG00589 | *1/*3C | *1/*3C | OK |
| HG01190 | *1/*1 | *1/*1 | OK |
| NA06991 | *1/*1 | *1/*1 | OK |
| NA06993 | *1/*1 | *1/*1 | OK |
| NA07000 | *1/*1 | *1/*1 | OK |
| NA07019 | *1/*1 | *1/*1 | OK |
| NA07029 | *1/*1 | *1/*1 | OK |
| NA07048 | *1/*1 | *1/*1 | OK |
| NA07055 | *1/*1 | *1/*1 | OK |
| NA07056 | *1/*1 | *1/*1 | OK |
| NA07348 | *1/*1 | *1/*1 | OK |
| NA07357 | *1/*1 | *1/*1 | OK |
| NA10831 | *1/*1 | *1/*1 | OK |
| NA10838 | *1/*1 | *1/*1 | OK |
| NA10846 | *1/*1 | *1/*1 | OK |
| NA10847 | *1/*1 | *1/*1 | OK |
| NA10851 | *1/*1 | *1/*1 | OK |
| NA10854 | *1/*1 | *1/*1 | OK |
| NA10856 | *1/*1 | *1/*1 | OK |
| NA10859 | *1/*1 | *1/*1 | OK |
| NA10865 | *1/*1 | *1/*1 | OK |
| NA11832 | *1/*1 | *1/*1 | OK |
| NA11839 | *1/*1 | *1/*1 | OK |
| NA11881 | *1/*1 | *1/*1 | OK |
| NA11993 | *1/*1 | *1/*1 | OK |
| NA12003 | *1/*1 | *1/*1 | OK |
| NA12006 | *1/*1 | *1/*1 | OK |
| NA12145 | *1/*1 | *1/*1 | OK |
| NA12156 | *1/*1 | *1/*1 | OK |
| NA12236 | *1/*1 | *1/*1 | OK |
| NA12336 | *1/*1 | *1/*1 | OK |
| NA12717 | *1/*1 | *1/*1 | OK |
| NA12753 | *1/*3A | *1/*3A | OK |
| NA12813 | *1/*1 | *1/*1 | OK |
| NA12815 | *1/*1 | *1/*1 | OK |
| NA12873 | *1/*1 | *1/*1 | OK |
| NA12878 | *1/*1 | *1/*1 | OK |
| NA12892 | *1/*1 | *1/*1 | OK |
| NA18484 | *1/*1 | *1/*1 | OK |
| NA18518 | *1/*1 | *1/*1 | OK |
| NA18519 | *1/*1 | *1/*1 | OK |
| NA18526 | *1/*1 | *1/*1 | OK |
| NA18544 | *1/*1 | *1/*1 | OK |
| NA18552 | *1/*1 | *1/*1 | OK |
| NA18563 | *1/*1 | *1/*1 | OK |
| NA18564 | *1/*1 | *1/*1 | OK |
| NA18565 | *1/*1 | *1/*1 | OK |
| NA18572 | *1/*1 | *1/*1 | OK |
| NA18617 | *1/*1 | *1/*1 | OK |
| NA18855 | *1/*3C | *1/*3C | OK |
| NA18861 | *1/*1 | *1/*1 | OK |
| NA18868 | *1/*1 | *1/*1 | OK |
| NA18873 | *1/*1 | *1/*1 | OK |
| NA18942 | *1/*1 | *1/*1 | OK |
| NA18945 | *1/*1 | *1/*1 | OK |
| NA18952 | *1/*1 | *1/*1 | OK |
| NA18959 | *1/*1 | *1/*1 | OK |
| NA18966 | *1/*3C | *1/*3C | OK |
| NA18973 | *1/*1 | *1/*1 | OK |
| NA18980 | *1/*1 | *1/*1 | OK |
| NA18992 | *1/*1 | *1/*1 | OK |
| NA19003 | *1/*1 | *1/*1 | OK |
| NA19007 | *1/*1 | *1/*1 | OK |
| NA19035 | *3C/*3C | *3C/*3C | OK |
| NA19095 | *1/*1 | *1/*1 | OK |
| NA19109 | *1/*1 | *1/*1 | OK |
| NA19122 | *1/*1 | *1/*1 | OK |
| NA19143 | *1/*1 | *1/*1 | OK |
| NA19147 | *1/*1 | *1/*1 | OK |
| NA19174 | *1/*1 | *1/*1 | OK |
| NA19207 | *1/*1 | *1/*1 | OK |
| NA19213 | *1/*1 | *1/*1 | OK |
| NA19226 | *1/*1 | *1/*1 | OK |
| NA19238 | *1/*1 | *1/*1 | OK |
| NA19239 | *1/*1 | *1/*1 | OK |
| NA19700 | *1/*1 | *1/*1 | OK |
| NA19785 | *1/*1 | *1/*1 | OK |
| NA19789 | *1/*1 | *1/*1 | OK |
| NA19819 | *1/*1 | *1/*1 | OK |
| NA19908 | *1/*1 | *1/*1 | OK |
| NA19917 | *1/*1 | *1/*1 | OK |
| NA19920 | *1/*3C | *1/*3C | OK |
| NA20296 | *1/*3C | *1/*3C | OK |
| NA20509 | *1/*1 | *1/*1 | OK |

## Non-core-truth samples (v0.1 sentinel SHOULD withhold)

| sample | GeT-RM | core-proxy | phenotype_status |
|---|---|---|---|
| HG00133 | *1/*2 | *1/*1 | ok |
| HG00276 | *1/*16 | *1/*1 | ok |
| HG00304 | *3A/*16 | *1/*3A | ok |
| HG01083 | *1/*2 | *1/*1 | ok |
| HG01474 | *1/*40 | *1/*1 | ok |
| HG02496 | *1/*24 | *1/*1 | ok |
| HG03521 | *1/*46 | *1/*1 | ok |
| NA10855 | *1/*32 | *1/*1 | ok |
| NA12044 | *1/*21 | *1/*1 | ok |
| NA12751 | *1/*12 | *1/*1 | ok |
| NA18603 | *1/*6 | *1/*1 | ok |
| NA19026 | *8/*33 or *1/*46 | *1/*1 | ok |
| NA19176 | *1/*8 | *1/*1 | ok |
