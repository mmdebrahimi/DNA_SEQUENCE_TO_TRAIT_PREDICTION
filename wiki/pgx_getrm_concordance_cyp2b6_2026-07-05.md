# CYP2B6 caller vs GeT-RM consensus on real 1000G (2026-07-05)

**Truth:** GeT-RM NGS consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) via the ursaPGx benchmark star-allele-comparison_common.tsv, column CYP2B6_getrm_cons
**Genotypes:** 1000 Genomes 30x phased panel (CYP2B6 region, pure-Python tabix-over-HTTP (no Docker); GeT-RM CDC consolidated PGx table; SINGLE-SNP *6-proxy (785 absent from panel))

- Overlap samples scored: **88**
- **Core-comparable diplotype concordance: 62/62 (1.0)**  (GeT-RM truth in *1/*6)
- Phenotype-correct incl. *38==*1: **62/88** (+0 *38 phenotype-equivalent samples)
- Correctly WITHHELD by sentinel: **0**
- **Genuine silent mis-call: 26/88 (29.5%)** -- non-core alleles beyond the v0 SNP set (+ sentinels where present); the honest residual blind spot.
- Correct-or-abstains: **62/88**

_GeT-RM CONSENSUS concordance on real 1000G genomes, independent caller. The strongest star-allele-CALLING validation tier available (vs the field's accepted consensus truth set). v0 covers the CORE SNP set; non-core-truth samples are scored separately (the v0.1 sentinel layer should WITHHOLD, not mis-call)._

## Core-comparable samples (GeT-RM truth in the v0 SNP set)

| sample | GeT-RM | predicted | match |
|---|---|---|---|
| HG00436 | *1/*6 | *1/*6 | OK |
| HG00589 | *1/*1 | *1/*1 | OK |
| NA06991 | *1/*6 | *1/*6 | OK |
| NA06993 | *1/*1 | *1/*1 | OK |
| NA07000 | *1/*1 | *1/*1 | OK |
| NA07048 | *1/*1 | *1/*1 | OK |
| NA07348 | *1/*1 | *1/*1 | OK |
| NA07357 | *1/*1 | *1/*1 | OK |
| NA10831 | *1/*1 | *1/*1 | OK |
| NA10838 | *6/*6 | *6/*6 | OK |
| NA10846 | *1/*1 | *1/*1 | OK |
| NA10847 | *1/*6 | *1/*6 | OK |
| NA10854 | *1/*1 | *1/*1 | OK |
| NA10855 | *1/*1 | *1/*1 | OK |
| NA10856 | *1/*1 | *1/*1 | OK |
| NA10859 | *1/*6 | *1/*6 | OK |
| NA11832 | *1/*1 | *1/*1 | OK |
| NA11881 | *1/*6 | *1/*6 | OK |
| NA11993 | *1/*1 | *1/*1 | OK |
| NA12003 | *6/*6 | *6/*6 | OK |
| NA12145 | *1/*1 | *1/*1 | OK |
| NA12156 | *1/*1 | *1/*1 | OK |
| NA12336 | *1/*6 | *1/*6 | OK |
| NA12753 | *1/*6 | *1/*6 | OK |
| NA12815 | *1/*1 | *1/*1 | OK |
| NA12873 | *1/*6 | *1/*6 | OK |
| NA12878 | *1/*1 | *1/*1 | OK |
| NA12892 | *1/*1 | *1/*1 | OK |
| NA18518 | *1/*6 | *1/*6 | OK |
| NA18519 | *1/*6 | *1/*6 | OK |
| NA18526 | *1/*1 | *1/*1 | OK |
| NA18544 | *1/*1 | *1/*1 | OK |
| NA18552 | *1/*6 | *1/*6 | OK |
| NA18563 | *1/*1 | *1/*1 | OK |
| NA18564 | *1/*1 | *1/*1 | OK |
| NA18565 | *1/*1 | *1/*1 | OK |
| NA18572 | *1/*1 | *1/*1 | OK |
| NA18617 | *1/*1 | *1/*1 | OK |
| NA18855 | *6/*6 | *6/*6 | OK |
| NA18868 | *1/*6 | *1/*6 | OK |
| NA18952 | *1/*1 | *1/*1 | OK |
| NA18966 | *1/*6 | *1/*6 | OK |
| NA18973 | *1/*6 | *1/*6 | OK |
| NA18980 | *6/*6 | *6/*6 | OK |
| NA18992 | *1/*6 | *1/*6 | OK |
| NA19003 | *1/*6 | *1/*6 | OK |
| NA19007 | *1/*1 | *1/*1 | OK |
| NA19035 | *1/*6 | *1/*6 | OK |
| NA19109 | *1/*6 | *1/*6 | OK |
| NA19122 | *6/*6 | *6/*6 | OK |
| NA19143 | *6/*6 | *6/*6 | OK |
| NA19176 | *6/*6 | *6/*6 | OK |
| NA19207 | *1/*6 | *1/*6 | OK |
| NA19213 | *6/*6 | *6/*6 | OK |
| NA19238 | *1/*6 | *1/*6 | OK |
| NA19239 | *1/*6 | *1/*6 | OK |
| NA19700 | *1/*6 | *1/*6 | OK |
| NA19785 | *1/*1 | *1/*1 | OK |
| NA19819 | *1/*1 | *1/*1 | OK |
| NA19908 | *1/*6 | *1/*6 | OK |
| NA19917 | *6/*6 | *6/*6 | OK |
| NA19920 | *6/*6 | *6/*6 | OK |

## Non-core-truth samples (v0.1 sentinel SHOULD withhold)

| sample | GeT-RM | core-proxy | phenotype_status |
|---|---|---|---|
| HG00276 | *2/(*4) | *1/*1 | ok |
| HG01190 | *1(*27)/*1(*5) | *1/*1 | ok |
| NA07019 | *1(*5) or *1(*22) | *1/*1 | ok |
| NA07029 | *6/(*27) | *1/*6 | ok |
| NA07055 | *6/(*27) | *1/*6 | ok |
| NA07056 | *6/(*22) | *1/*6 | ok |
| NA10851 | *1/*1 (*27) | *1/*1 | ok |
| NA10865 | *1/(*2;*10) | *1/*1 | ok |
| NA11839 | *1/*1 (*15) | *1/*1 | ok |
| NA12006 | *1/*1 (*5) | *1/*1 | ok |
| NA12236 | *1/*1 (*4) | *1/*1 | ok |
| NA12717 | *1/*7 | *1/*6 | ok |
| NA12813 | *1/*2 | *1/*1 | ok |
| NA18484 | *1/*18 | *1/*1 | ok |
| NA18861 | *1/*18 | *1/*1 | ok |
| NA18873 | *1/*18 | *1/*1 | ok |
| NA18942 | *1/*2 | *1/*1 | ok |
| NA18945 | *6/(*2;*10) | *1/*6 | ok |
| NA18959 | *1/*1(*27) | *1/*1 | ok |
| NA19095 | *18/*18 | *1/*1 | ok |
| NA19147 | *1/*18 | *1/*1 | ok |
| NA19174 | *6/*18 | *1/*6 | ok |
| NA19226 | *18/(*20) | *1/*6 | ok |
| NA19789 | *1/*1 (*27) | *1/*1 | ok |
| NA20296 | *1/*2 | *1/*1 | ok |
| NA20509 | *1/*7 | *1/*6 | ok |
