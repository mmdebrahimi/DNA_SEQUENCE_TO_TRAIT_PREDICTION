# CYP2D6 STRUCTURAL surface — read-depth copy-number on real 1000G CRAMs (2026-07-06)

_CYP2D6-body/control read-depth ratio off 1000G 30x CRAMs (samtools + ENA reference auto-fetch, no full-reference download); ratio -> integer copy number -> deletion/normal/duplication._

- Regions: CYP2D6 `chr22:42126000-42130000` / control `chr22:41000000-41050000`
- NORMAL baseline: pinned **1.26** vs measured **1.267** (median of NORMAL-truth ratios)
- Samples measured: **39**  (CN-scored **26**; hybrid/ambiguous excluded **13**)
- **CN-class concordance (deletion/normal/duplication): 26/26 (1.0)**

_Real-CRAM read-depth copy-number caller (the structural surface the SNP cell is blind to). Resolves *5 deletion + *xN duplication; NEVER resolves hybrid IDENTITY (*13/*36/*68 -> excluded from CN scoring, measured only; needs CYP2D6-vs-CYP2D7 PSV analysis, Cyrius-class). Coarse integer CN, not a breakpoint caller._

## CN-scored samples (non-hybrid truth)

| sample | truth | ratio | CN | predicted | truth class | match |
|---|---|---|---|---|---|---|
| NA12873 | `*1/*5` | 0.58 | 1 | deletion | deletion | OK |
| NA18873 | `*5/*17` | 0.59 | 1 | deletion | deletion | OK |
| NA10856 | `*1/*5` | 0.59 | 1 | deletion | deletion | OK |
| NA18855 | `*1/*5` | 0.62 | 1 | deletion | deletion | OK |
| NA18992 | `*1/*5` | 0.62 | 1 | deletion | deletion | OK |
| NA18868 | `*2/*5` | 0.63 | 1 | deletion | deletion | OK |
| NA18945 | `*1/*5` | 0.64 | 1 | deletion | deletion | OK |
| NA18861 | `*5/*29` | 0.65 | 1 | deletion | deletion | OK |
| NA19035 | `*2/*5` | 0.66 | 1 | deletion | deletion | OK |
| HG00276 | `*4/*5` | 0.66 | 1 | deletion | deletion | OK |
| NA12336 | `*5/*41` | 0.67 | 1 | deletion | deletion | OK |
| NA10831 | `*4/*5` | 0.71 | 1 | deletion | deletion | OK |
| HG00589 | `*1/*21` | 1.16 | 2 | normal_copy_number | normal_copy_number | OK |
| NA07000 | `*2(*35)/*9` | 1.23 | 2 | normal_copy_number | normal_copy_number | OK |
| NA07056 | `*2/*4` | 1.24 | 2 | normal_copy_number | normal_copy_number | OK |
| NA07029 | `*1/*35` | 1.26 | 2 | normal_copy_number | normal_copy_number | OK |
| NA07048 | `*1/*4` | 1.27 | 2 | normal_copy_number | normal_copy_number | OK |
| NA07019 | `*1/*4` | 1.29 | 2 | normal_copy_number | normal_copy_number | OK |
| NA06991 | `*1/*4` | 1.35 | 2 | normal_copy_number | normal_copy_number | OK |
| NA07055 | `*4/*4` | 1.40 | 2 | normal_copy_number | normal_copy_number | OK |
| HG00436 | `*2x2/*71` | 1.74 | 3 | duplication | duplication | OK |
| NA19920 | `*1/*4x2` | 1.76 | 3 | duplication | duplication | OK |
| NA19226 | `*2/*2x2` | 1.77 | 3 | duplication | duplication | OK |
| NA19207 | `*2x2/*10` | 1.90 | 3 | duplication | duplication | OK |
| NA19819 | `*2/*4x2` | 1.94 | 3 | duplication | duplication | OK |
| NA19109 | `*2x2/*29` | 1.98 | 3 | duplication | duplication | OK |

## Measured-only (hybrid / ambiguous — NOT CN-scored; identity unresolved)

| sample | truth | ratio | CN | predicted |
|---|---|---|---|---|
| HG01190 | `*68+*4/*5` | 0.81 | 1 | deletion |
| NA11832 | `*1/(*68)+*4` | 1.32 | 2 | normal_copy_number |
| NA12878 | `*3/(*68)+*4` | 1.41 | 2 | normal_copy_number |
| NA10855 | `*1/(*68)+*4` | 1.41 | 2 | normal_copy_number |
| NA18565 | `*10/*36x2` | 1.67 | 3 | duplication |
| NA19785 | `*1/*13+*2` | 1.68 | 3 | duplication |
| NA18572 | `*36+*10/*41` | 1.71 | 3 | duplication |
| NA18959 | `*2/*36+*10` | 1.77 | 3 | duplication |
| NA18563 | `*1/*36+*10` | 1.79 | 3 | duplication |
| NA18980 | `*2/*36+*10` | 1.81 | 3 | duplication |
| NA18564 | `*2A/*36+*10` | 1.87 | 3 | duplication |
| NA18617 | `*36+*10/*36+*10` | 2.15 | 3 | duplication |
| NA18526 | `*1/*36x2+*10` | 2.29 | 4 | duplication |
