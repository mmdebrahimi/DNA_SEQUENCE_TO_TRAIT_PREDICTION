# CRyPTIC TB MIC calibration — are the MIC prediction intervals honest? (2026-07-12)

**Verdict: CALIBRATED_MIC_INTERVALS**

Do split-conformal intervals on CRyPTIC BMD-MIC (12 drugs) achieve nominal held-out coverage on resolved MICs (a calibrated quantitative TB decoder)? Substrate: CRyPTIC reuse-table measured BMD-MIC (wet-lab, NOT G1-circular BV-BRC) + WHO grade-1/2 determinant presence from VARIANTS.parquet (frozen tb_amr matcher).

`cover_resolved_90` = held-out coverage on RESOLVED (uncensored) MICs (target 0.90). `interval_fold_factor` = MIC within ×/÷ this factor (2^halfwidth). `consistency_censored` = fraction of censored isolates whose interval respects the censoring bound.

| drug | n | resolved | censored | dets | R2(res) | **cover_90** | MIC ± | consistency | calibrated | informative |
|---|---|---|---|---|---|---|---|---|---|---|
| rifampicin | 12097 | 5012 | 7085 | 17 | 0.42 | **0.9035** | ×/÷3.14 | 0.8388 | YES | YES |
| isoniazid | 12068 | 7669 | 4399 | 5 | 0.7661 | **0.9057** | ×/÷3.56 | 0.9041 | YES | YES |
| ethambutol | 12156 | 11604 | 552 | 13 | 0.4566 | **0.904** | ×/÷2.85 | 0.6793 | YES | YES |
| levofloxacin | 12161 | 11409 | 752 | 11 | 0.5259 | **0.9075** | ×/÷2.25 | 0.3152 | YES | YES |
| moxifloxacin | 12192 | 10635 | 1557 | 10 | 0.3699 | **0.9005** | ×/÷2.18 | 0.2832 | YES | YES |
| amikacin | 12070 | 4317 | 7753 | 2 | -0.7415 | **0.9117** | ×/÷2.83 | 0.9837 | YES | no (interval≈marginal) |
| kanamycin | 12128 | 7847 | 4281 | 6 | -0.2411 | **0.9169** | ×/÷2.16 | 0.9631 | YES | no (interval≈marginal) |
| ethionamide | 12130 | 10890 | 1240 | 8 | 0.2007 | **0.9084** | ×/÷3.07 | 0.5113 | YES | YES |
| bedaquiline | 12066 | 10156 | 1910 | 1 | -0.0797 | **0.9073** | ×/÷3.3 | 0.7733 | YES | no (interval≈marginal) |
| delamanid | 11925 | — | — | — | — | too few / no features / degenerate MIC | — | — | — | — |
| linezolid | 12187 | 12002 | 185 | 1 | 0.0071 | **0.9056** | ×/÷2.56 | 0.0595 | YES | no (interval≈marginal) |
| clofazimine | 12047 | — | — | — | — | too few / no features / degenerate MIC | — | — | — | — |

## Honest caveats
- COVERAGE-VALID != INFORMATIVE: amikacin (R2=-0.7415), kanamycin (R2=-0.2411), bedaquiline (R2=-0.0797), linezolid (R2=0.0071) hit ~0.90 coverage but their determinant model does NOT beat the marginal MIC (R2<=0.05), so the interval is just the marginal spread (few WHO grade-1/2 determinants clear MIN_SUPPORT + resistance is rare). See `informative`.
- CRyPTIC MIC is end-censored (RIF 58% / INH 36% <=/> values) — coverage is on the RESOLVED subset (mid-ladder-enriched), censored isolates scored only for consistency.
- Features = WHO grade-1/2 determinant presence (the catalogue that defines R/S) -> the point model explains R/S better than the exact MIC rung; the conformal WIDTH honestly reflects that.
- Split-conformal gives MARGINAL (not per-genotype) coverage; measured wet-lab MIC, in-distribution vs the WHO catalogue (not independent).

Citation: CRyPTIC Consortium 2022; WHO TB mutation catalogue v2 (2023).