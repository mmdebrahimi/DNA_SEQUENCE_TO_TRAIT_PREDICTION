# Stage 1 -- N=38 cipro engineering screen (2026-05-16)

> **Engineering screen, NOT a powered statistical comparison.** Stage 1's role is go/no-go for spending Stage 2 N=150 Databricks burst budget. Stage 2 is the real ship gate (>=5 pp AUROC + biology check on gyrA/parC/parE attribution).

**Cohort:** `data\processed\gate_b_n40_cipro_cohort.parquet` (effective N=38; balance 19R/19S)
**Drug:** ciprofloxacin
**Gate threshold:** >=3 pp (max(NT-XGBoost, NT-logreg) AUROC - k-mer-XGB AUROC)
**Best NT-only head:** NT-logreg (AUROC 0.673)
**Gap vs k-mer-XGB:** +2.5 pp
**Paired bootstrap 95% CI on gap:** [-20.0, +23.3] pp (B=1000 effective 1000, mean +2.5 pp)
**Verdict:** FAIL (gap +2.5 pp < 3 pp threshold)
**Stage 2 action:** `ALTERNATIVE_POOLING_RERUN`

## Per-variant LOSO results

| Variant | AUROC | AUPRC | Gate-bearing? |
|---|---:|---:|:---:|
| NT-XGBoost | 0.615 | 0.682 | yes |
| NT-logreg | 0.673 | 0.676 | yes |
| k-mer-XGB | 0.648 | 0.617 | yes |
| NT+k-mer-fusion-logreg | 0.363 | 0.439 | diagnostic |

## Gate analysis

- best NT-only: **NT-logreg** at **0.673**
- k-mer-XGB AUROC: **0.648**
- Gap: **+2.5 pp** (best NT-only - k-mer-XGB)
- Paired bootstrap 95% CI on gap: **[-20.0, +23.3] pp** (mean +2.5 pp, B=1000 effective 1000)
- Verdict (point-gap function): **FAIL (gap +2.5 pp < 3 pp threshold)**
- Stage 2 action (decision-layer): **`ALTERNATIVE_POOLING_RERUN`**

## Lineage diagnostic

- Unique MLSTs: **38** of 38 strains (uniqueness fraction 100.00%)
- Largest MLST group: `MLST.Escherichia_coli_1.131` with N=1 (0R/1S)
- LOMO note: at N=38 with this MLST cardinality, most LOMO folds are size-1 -> degenerate. Reporting LOSO only; per-strain table below substitutes for LOMO diagnostics.

### Per-strain LOSO predictions

| Strain | MLST | True | NT-best score | k-mer score | NT-best correct? | k-mer correct? |
|---|---|:---:|---:|---:|:---:|:---:|
| 562.28805 | MLST.ecoli_achtman_4.410 | 1 | 0.056 | 0.288 | X | X |
| 562.30362 | MLST.ecoli_achtman_4.4 | 1 | 1.000 | 0.888 | OK | OK |
| 562.28563 | MLST.ecoli_achtman_4.156 | 1 | 0.366 | 0.776 | X | OK |
| 562.17621 | MLST.ecoli_achtman_4.167 | 1 | 0.876 | 0.292 | OK | X |
| 562.17721 | MLST.ecoli_achtman_4.1284 | 1 | 0.997 | 0.914 | OK | OK |
| 1328433.3 | MLST.ecoli_achtman_4.131 | 1 | 0.161 | 0.613 | X | OK |
| 562.22426 | MLST.Escherichia_coli_1.167,MLST.Escherichia_coli_2.2 | 1 | 0.992 | 0.831 | OK | OK |
| 562.12960 | MLST.ecoli_achtman_4.101 | 1 | 0.985 | 0.374 | OK | X |
| 1328434.3 | MLST.ecoli_achtman_4.405 | 1 | 0.357 | 0.547 | X | OK |
| 562.45848 | MLST.ecoli_achtman_4.349 | 1 | 0.618 | 0.133 | OK | X |
| 562.45851 | MLST.ecoli_achtman_4.372 | 1 | 0.032 | 0.190 | X | X |
| 562.50304 | MLST.ecoli_achtman_4.1193 | 1 | 0.932 | 0.296 | OK | X |
| 562.50252 | MLST.ecoli_achtman_4.393 | 1 | 0.274 | 0.344 | X | X |
| 562.50245 | MLST.ecoli_achtman_4.1431 | 1 | 0.984 | 0.095 | OK | X |
| 562.50237 | MLST.ecoli_achtman_4.224 | 1 | 0.577 | 0.673 | OK | OK |
| 562.50250 | MLST.ecoli_achtman_4.44 | 1 | 0.631 | 0.748 | OK | OK |
| 562.7699 | MLST.ecoli_achtman_4.335 | 1 | 0.031 | 0.468 | X | X |
| 562.7572 | MLST.ecoli_achtman_4.301 | 1 | 0.196 | 0.629 | X | OK |
| 562.13502 | MLST.ecoli_achtman_4.10 | 1 | 0.981 | 0.931 | OK | OK |
| 562.28389 | MLST.ecoli_achtman_4.1276 | 0 | 0.055 | 0.237 | OK | OK |
| 562.7627 | MLST.Escherichia_coli_1.4554 | 0 | 0.089 | 0.098 | OK | OK |
| 562.16326 | MLST.ecoli_achtman_4.1809 | 0 | 0.998 | 0.072 | X | OK |
| 562.16325 | MLST.ecoli_achtman_4.1408 | 0 | 0.979 | 0.556 | X | X |
| 562.7641 | MLST.Escherichia_coli_1.5543 | 0 | 0.138 | 0.259 | OK | OK |
| 562.52722 | MLST.ecoli_achtman_4.29 | 0 | 0.962 | 0.971 | X | X |
| 562.28565 | MLST.ecoli_achtman_4.38 | 0 | 0.726 | 0.891 | X | X |
| 562.50301 | MLST.ecoli_achtman_4.1317 | 0 | 0.586 | 0.186 | X | OK |
| 562.7784 | MLST.ecoli_achtman_4.2678 | 0 | 0.062 | 0.123 | OK | OK |
| 1328432.3 | MLST.ecoli_achtman_4.2569 | 0 | 0.227 | 0.146 | OK | OK |
| 562.7789 | MLST.ecoli_achtman_4.123 | 0 | 0.003 | 0.181 | OK | OK |
| 562.45853 | MLST.ecoli_achtman_4.315 | 0 | 0.703 | 0.470 | X | OK |
| 562.50287 | MLST.ecoli_achtman_4.442 | 0 | 0.089 | 0.507 | OK | X |
| 562.50295 | MLST.Escherichia_coli_1.131 | 0 | 0.951 | 0.814 | X | X |
| 562.7690 | MLST.ecoli_achtman_4.2144 | 0 | 0.011 | 0.455 | OK | OK |
| 562.7575 | MLST.ecoli_achtman_4.337 | 0 | 0.070 | 0.622 | OK | X |
| 562.7710 | MLST.ecoli_achtman_4.3 | 0 | 0.017 | 0.264 | OK | OK |
| 562.7695 | MLST.ecoli_achtman_4.2089 | 0 | 0.316 | 0.572 | OK | X |
| 562.7717 | MLST.ecoli_achtman_4.2346 | 0 | 0.041 | 0.170 | OK | OK |

## Notes

- All variants ran with `calibrate=False` (uniform calibration discipline; calibration is small-N footgun per LESSONS_LEARNED 2026-05-14).
- k-mer + fusion use within-fold vocab rebuild from training-set sequences only (no held-out leakage).
- Gene-presence variant NOT included -- RefSeq GFF3 carries `gene=` for ~11% of CDSs -> INDETERMINATE_IDENTIFIER_OOV on this annotation source. See `plans/Gene_Presence_AUROC_Bug_Fix_Plan.md`.
- AMRFinderPlus POINT* SNP-table baseline NOT included -- deferred to Stage 2 per Phase2_Decision_Gate D6. 'Best classical' here is bounded by what was run; gyrA/parC/parE point-mutation features are NOT part of the comparator.
- LOSO at N=38 has +-0.10 noise floor on AUROC; >=3 pp is INSIDE the noise. The bootstrap CI surfaces this honestly.
- Verdict semantics are frozen as a pure function of point gap. The `stage2_action` field is the decision-layer; see `plans/Stage1_N40_Cipro_Engineering_Screen_Plan.md` Verdict-Time Pre-Commitments.

## Next action by stage2_action

- **`BURST_STAGE_2`** -> proceed to Stage 2 Databricks burst with N=150 cohort build.
- **`HOLD_STAGE_2_CI_DEGENERATE`** -> do NOT spend Stage 2 burst budget; next escalation is `ALTERNATIVE_POOLING_RERUN`.
- **`ALTERNATIVE_POOLING_RERUN`** -> Stage 1b with `mean+max` aggregation; if still <3 pp, escalate to `PIVOT_TO_BAKTA`.
- **`PIVOT_TO_BAKTA`** -> Bakta re-annotation + gene-presence comparator pathway per `plans/Gene_Presence_AUROC_Bug_Fix_Plan.md` follow-up.