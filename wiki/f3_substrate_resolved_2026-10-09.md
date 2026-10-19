# F3 — the positive-control substrate is RESOLVED, and it cost 330 MB

**Authority call #5 asked whether to acquire a natural-population genotype panel so the de-confounder's
positive control could run. The answer is measured, not estimated: the substrate works, it is 330 MB,
and the phenotypes were already on disk.**

This does **not** mean the de-confounder is calibrated. It means the thing that blocked the calibration
is no longer blocking it.

## What was blocked

`wiki/f3_deconfound_calibration_2026-10-09.md`, verdict
`blocked:positive_control_not_fold_compatible`. `dgrp` — the only cached natural population — has **no
rung that is both fold-compatible (largest cluster ≤ 1/k) and not mostly singletons**: p1–p30 give
66–84% singletons, p60 gives a 52% dominant cluster. So grouped CV could not run on it at all.

## What was measured

Same pipeline, same guards, same ladder as the census — only the panel changed.

| rung | clusters | largest | singletons | usable at k=5 |
|---|---|---|---|---|
| p1 | 737 | 0.090 | 0.847 | fail |
| p2 | 389 | 0.095 | 0.609 | fail |
| **p5** | **108** | **0.095** | **0.306** | **PASS** |
| p10 | 62 | 0.208 | 0.226 | fail |
| p15 | 43 | 0.274 | 0.186 | fail |
| p30 | 26 | 0.499 | 0.269 | fail |
| p70 | 7 | 0.862 | 0.000 | fail |

**A fold-compatible grouping exists at percentile 5** (threshold 0.67656): 108 clusters, largest cluster
9.5%, singletons 30.6%. `dgrp` passed at **no** rung; this panel passes cleanly at one.

## The three pieces, each verified

| piece | state |
|---|---|
| genotypes | `1001_SNP_MATRIX.tar.gz`, **330 MB**, HTTP 200, no auth, from `1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/`. **Both MD5s matched the publisher's own `md5sum.txt` exactly.** Contains `imputed_snps_binary.hdf5` (919 MB) **and `kinship_ibs_mac5.hdf5` (10.3 MB), a pre-computed IBS kinship matrix** — so the grouping question was answerable from 10 MB without touching the genotypes. |
| phenotypes | **already on disk.** `data/arabidopsis/FT10_pheno_261.csv` overlaps the genotype accessions on **1,003 (86.3%)**; FT16 on 970 (86.5%); FLC on 112 (67.1%). |
| grouping | measured above: p5, fold-compatible at k=5. |

**The estimate was wrong in the useful direction.** A literature search put this file at "roughly a few
GB, size unverified"; a HEAD request measured **0.33 GB**. The merged VCF really is 19.23 GB, and the full
genome VCF 132 GB — but neither is needed.

## Why FT10 is the right positive control, not merely an available one

FT10 at n=1003 is **exactly the trait and panel of this project's closed Arabidopsis negative**
(`wiki/phase2_arabidopsis_result_2026-06-12.md`): the embedding learned **population structure** rather
than the flowering-time signal, within-group r² −0.13, spearman 0.22 against structure-only 0.48.

That is what makes it an ideal positive control. A control for a *structure* de-confounder needs a panel
where structure is **known** to dominate the phenotype — and this is the one panel here where that has
been measured rather than assumed.

## Honest limits — read before running the calibration

1. **The distance is the 1001G project's OWN kinship (MAC ≥ 5), converted, NOT this repo's
   `ibs_distance_matrix` over raw genotypes.** Closely related, not the same quantity. The grouping
   should be cross-checked against `imputed_snps_binary.hdf5` before the calibration's number is
   published; if the two disagree, that is a finding to diagnose, not a number to pick between.
2. **The genotypes are IMPUTED.** Imputation draws on within-panel information, so it can in principle
   flatten or sharpen relatedness structure — the exact quantity a de-confounding calibration measures.
   Unquantified here.
3. **Fold-compatibility of the GROUPING was measured; the CALIBRATION was not run.** That needs the
   919 MB genotype matrix for the ridge, not just the kinship, and it is the next step — not a result.
4. **86.3% overlap means 159 phenotyped accessions have no genotype in this matrix** and would be
   dropped.
5. A single passing rung (p5) is a narrow window — p2 and p10 both fail, in opposite directions. The
   grouping is therefore sensitive to the ladder, and the `select_threshold` bias (smallest passing rung,
   which under-merges and biases *against* detecting a de-confounding response) applies with full force.
6. This resolves the SUBSTRATE question only. It says nothing about whether the de-confounder will pass
   its own differential criterion, and nothing about the generative north star.

## Reproduce

```
curl -L -o D:/dna_decode_cache/arabidopsis_1001g/1001_SNP_MATRIX.tar.gz \
  https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz
tar -xzf 1001_SNP_MATRIX.tar.gz && md5sum -c md5sum.txt
```

Then convert `kinship_ibs_mac5.hdf5` to a distance (`(max−k)/(max−min)`, symmetrised, zero diagonal)
and call `panel_structure.select_threshold(dm, ids, max_largest=1/5)`.
