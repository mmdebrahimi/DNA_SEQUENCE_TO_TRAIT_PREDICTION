# CYP2D6 caller vs GeT-RM consensus on real 1000G (2026-07-06)

**Truth:** GeT-RM consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) via ursaPGx star-allele-comparison_common.tsv, column CYP2D6_getrm_cons
**Genotypes:** 1000 Genomes 30x phased panel (CYP2D6 chr22 region, pure-Python tabix-over-HTTP (scripts/fetch_1000g_region.py; no Docker); SNP surface only — structural alleles BAM-required + EXCLUDED)
**Surface:** SNP-decodable star alleles ONLY (structural alleles BAM-required; see note).

- Overlap samples: **87**  (tiered below — no single inflated denominator)
- **Core-SNP diplotype concordance: 46/47 (0.9787)**  (truth in *1/*2/*3/*4/*6/*9/*10/*17/*29/*35/*41)
- Core-SNP PHENOTYPE concordance: 46/47 (0.9787)  (*35==*2==*1 all normal-function)
- Non-core SNP alleles (residual, mis-called; *14/*15/*21/*40/*46): **7**
- Structural EXCLUDED (BAM-required; *5/*xN/*13/*36/*68): **31**
- Ambiguous-truth EXCLUDED (parenthetical alternative annotation): **2**

_GeT-RM CONSENSUS core-diplotype concordance on the SNP-DECODABLE subset, independent caller. Structural + ambiguous truth EXCLUDED (tiered denominators; raw + normalized truth both retained). SNP surface only — NOT full CYP2D6 typing._

_structural alleles (*5 deletion / *xN duplication / *13/*36/*61/*63/*68 CYP2D6-CYP2D7 hybrids) are NOT VCF-decodable -> EXCLUDED from the scored denominator; they are NOT withheld and may be SILENTLY MIS-CALLED (cnv_hybrid_unassessed). Full typing needs a BAM/CRAM + Cyrius-class caller._

## Core-SNP samples (the scored tier)

| sample | raw truth | normalized | predicted | match |
|---|---|---|---|---|
| NA06991 | `*1/*4` | *1/*4 | *1/*4 | OK |
| NA07019 | `*1/*4` | *1/*4 | *1/*4 | OK |
| NA07029 | `*1/*35` | *1/*35 | *1/*35 | OK |
| NA07048 | `*1/*4` | *1/*4 | *1/*4 | OK |
| NA07055 | `*4/*4` | *4/*4 | *4/*4 | OK |
| NA07056 | `*2/*4` | *2/*4 | *2/*4 | OK |
| NA07348 | `*1/*6` | *1/*6 | *1/*6 | OK |
| NA07357 | `*1/*6` | *1/*6 | *1/*6 | OK |
| NA10838 | `*2/*4` | *2/*4 | *2/*4 | OK |
| NA10846 | `*1/*4` | *1/*4 | *1/*4 | OK |
| NA10847 | `*1/*41` | *1/*41 | *1/*41 | OK |
| NA10851 | `*1/*4` | *1/*4 | *1/*4 | OK |
| NA10854 | `*1/*4` | *1/*4 | *1/*4 | OK |
| NA10859 | `*1/*2` | *1/*2 | *1/*2 | OK |
| NA10865 | `*1/*41` | *1/*41 | *1/*41 | OK |
| NA11839 | `*1/*2` | *1/*2 | *1/*2 | OK |
| NA11881 | `*2/*3` | *2/*3 | *2/*3 | OK |
| NA11993 | `*1/*9` | *1/*9 | *1/*9 | OK |
| NA12003 | `*4/*35` | *4/*35 | *4/*35 | OK |
| NA12006 | `*4/*41` | *4/*41 | *4/*41 | OK |
| NA12145 | `*1/*4` | *1/*4 | *1/*4 | OK |
| NA12156 | `*1/*4` | *1/*4 | *4/*4 | X |
| NA12236 | `*1/*4` | *1/*4 | *1/*4 | OK |
| NA12717 | `*1/*1` | *1/*1 | *1/*1 | OK |
| NA12753 | `*2/*3` | *2/*3 | *2/*3 | OK |
| NA12813 | `*2/*4` | *2/*4 | *2/*4 | OK |
| NA12815 | `*2/*41` | *2/*41 | *2/*41 | OK |
| NA12892 | `*2/*3` | *2/*3 | *2/*3 | OK |
| NA18484 | `*1/*17` | *1/*17 | *1/*17 | OK |
| NA18518 | `*17/*29` | *17/*29 | *17/*29 | OK |
| NA18519 | `*1/*29` | *1/*29 | *1/*29 | OK |
| NA18544 | `*10/*41` | *10/*41 | *10/*41 | OK |
| NA18942 | `*2/*2` | *2/*2 | *2/*2 | OK |
| NA18952 | `*2/*2` | *2/*2 | *2/*2 | OK |
| NA18966 | `*1/*2` | *1/*2 | *1/*2 | OK |
| NA19003 | `*1/*1` | *1/*1 | *1/*1 | OK |
| NA19007 | `*1/*1` | *1/*1 | *1/*1 | OK |
| NA19095 | `*1/*29` | *1/*29 | *1/*29 | OK |
| NA19122 | `*2/*17` | *2/*17 | *2/*17 | OK |
| NA19147 | `*17/*29` | *17/*29 | *17/*29 | OK |
| NA19176 | `*1/*2` | *1/*2 | *1/*2 | OK |
| NA19213 | `*1/*1` | *1/*1 | *1/*1 | OK |
| NA19238 | `*1/*17` | *1/*17 | *1/*17 | OK |
| NA19700 | `*4/*29` | *4/*29 | *4/*29 | OK |
| NA19789 | `*1/*1` | *1/*1 | *1/*1 | OK |
| NA20296 | `*1/*2` | *1/*2 | *1/*2 | OK |
| NA20509 | `*4/*35` | *4/*35 | *4/*35 | OK |

### Core mis-call diagnosis
- **NA12156** truth `*1/*4` -> predicted `*4/*4`: likely_structural_confound (predicted homozygous; truth het-with-*1 -> hidden CNV/hybrid; cnv_hybrid_unassessed)

## Non-core SNP (residual — mis-called)

| sample | raw truth | normalized | SNP-proxy predicted |
|---|---|---|---|
| HG00589 | `*1/*21` | *1/*21 | *1/*2 |
| NA18552 | `*1/*14` | *1/*14 | *1/*2 |
| NA18973 | `*1/*21` | *1/*21 | *1/*2 |
| NA19174 | `*4/*40` | *4/*40 | *4/*17 |
| NA19239 | `*15/*17` | *15/*17 | *1/*17 |
| NA19908 | `*1/*46` | *1/*46 | *1/*2 |
| NA19917 | `*1/*40` | *1/*40 | *1/*17 |

## Structural (EXCLUDED — BAM-required; may be silently mis-called)

| sample | raw truth | normalized | SNP-proxy predicted |
|---|---|---|---|
| HG00276 | `*4/*5` | *4/*5 | *4/*4 |
| HG00436 | `*2x2/*71` | *2/*71 | *1/*2 |
| HG01190 | `*68+*4/*5` | *68/*4/*5 | *4/*41 |
| NA10831 | `*4/*5` | *4/*5 | *4/*4 |
| NA10855 | `*1/(*68)+*4` | *1/*68/*4 | *1/*4 |
| NA10856 | `*1/*5` | *1/*5 | *1/*1 |
| NA11832 | `*1/(*68)+*4` | *1/*68/*4 | *1/*4 |
| NA12336 | `*5/*41` | *5/*41 | *41/*41 |
| NA12873 | `*1/*5` | *1/*5 | *1/*1 |
| NA12878 | `*3/(*68)+*4` | *3/*68/*4 | *3/*4 |
| NA18526 | `*1/*36x2+*10` | *1/*36/*10 | *1/*10 |
| NA18563 | `*1/*36+*10` | *1/*36/*10 | *1/*10 |
| NA18564 | `*2A/*36+*10` | *2/*36/*10 | *2/*10 |
| NA18565 | `*10/*36x2` | *10/*36 | *10/*10 |
| NA18572 | `*36+*10/*41` | *36/*10/*41 | *10/*41 |
| NA18617 | `*36+*10/*36+*10` | *36/*10/*36/*10 | *10/*10 |
| NA18855 | `*1/*5` | *1/*5 | *1/*1 |
| NA18861 | `*5/*29` | *5/*29 | *29/*29 |
| NA18868 | `*2/*5` | *2/*5 | *2/*2 |
| NA18873 | `*5/*17` | *5/*17 | *17/*17 |
| NA18945 | `*1/*5` | *1/*5 | *1/*1 |
| NA18959 | `*2/*36+*10` | *2/*36/*10 | *2/*10 |
| NA18980 | `*2/*36+*10` | *2/*36/*10 | *2/*10 |
| NA18992 | `*1/*5` | *1/*5 | *1/*1 |
| NA19035 | `*2/*5` | *2/*5 | *2/*2 |
| NA19109 | `*2x2/*29` | *2/*29 | *2/*29 |
| NA19207 | `*2x2/*10` | *2/*10 | *2/*10 |
| NA19226 | `*2/*2x2` | *2/*2 | *2/*2 |
| NA19785 | `*1/*13+*2` | *1/*13/*2 | *1/*2 |
| NA19819 | `*2/*4x2` | *2/*4 | *2/*4 |
| NA19920 | `*1/*4x2` | *1/*4 | *1/*4 |

## Ambiguous truth (EXCLUDED — uncertain annotation)

| sample | raw truth | normalized | SNP-proxy predicted |
|---|---|---|---|
| NA07000 | `*2 (*35)/*9` | *2/*35/*9 | *9/*35 |
| NA19143 | `*2 (*45)/*10` | *2/*45/*10 | *2/*10 |
