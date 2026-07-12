# CRyPTIC TB MIC calibration — are the MIC prediction intervals honest? (2026-07-12)

**Verdict: CALIBRATED_MIC_INTERVALS**

Do split-conformal intervals on CRyPTIC BMD-MIC (2 drugs) achieve nominal held-out coverage on resolved MICs (a calibrated quantitative TB decoder)? Substrate: CRyPTIC reuse-table measured BMD-MIC (wet-lab, NOT G1-circular BV-BRC) + WHO grade-1/2 determinant presence from VARIANTS.parquet (frozen tb_amr matcher).

`cover_resolved_90` = held-out coverage on RESOLVED (uncensored) MICs (target 0.90). `interval_fold_factor` = MIC within ×/÷ this factor (2^halfwidth). `consistency_censored` = fraction of censored isolates whose interval respects the censoring bound.

| drug | n | resolved | censored | dets | R2(res) | **cover_90** | MIC ± | consistency | calibrated | informative |
|---|---|---|---|---|---|---|---|---|---|---|
| rifampicin | 12097 | 5012 | 7085 | 17 | 0.42 | **0.9035** | ×/÷3.14 | 0.8388 | YES | YES |
| isoniazid | 12068 | 7669 | 4399 | 5 | 0.7661 | **0.9057** | ×/÷3.56 | 0.9041 | YES | YES |

## Honest caveats
- COVERAGE-VALID != INFORMATIVE: none hit ~0.90 coverage but their determinant model does NOT beat the marginal MIC (R2<=0.05), so the interval is just the marginal spread (few WHO grade-1/2 determinants clear MIN_SUPPORT + resistance is rare). See `informative`.
- CRyPTIC MIC is end-censored (RIF 58% / INH 36% <=/> values) — coverage is on the RESOLVED subset (mid-ladder-enriched), censored isolates scored only for consistency.
- Features = WHO grade-1/2 determinant presence (the catalogue that defines R/S) -> the point model explains R/S better than the exact MIC rung; the conformal WIDTH honestly reflects that.
- Split-conformal gives MARGINAL (not per-genotype) coverage; measured wet-lab MIC, in-distribution vs the WHO catalogue (not independent).

Citation: CRyPTIC Consortium 2022; WHO TB mutation catalogue v2 (2023).