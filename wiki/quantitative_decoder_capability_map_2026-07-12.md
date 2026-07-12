# Quantitative-decoder capability map (2026-07-12)

**30/34 powered cells INFORMATIVELY calibrated** (34 coverage-valid, of which 4 are coverage-valid-only — model no better than the marginal) at 90% target coverage across HIV-1, M. tuberculosis; calibrated-cell coverage spans 0.8931–0.9281.
(6 cells under-powered; 40 cells total.)

Sources: HIV `hiv_quantitative_calibration_2026-07-11.json` · TB `tb_mic_calibration_panel_2026-07-12.json`.

The deterministic decoder emits R/S; this is the QUANTITATIVE layer — coverage-valid prediction intervals. `cover_90` = held-out coverage (target 0.90); `fold_factor` = the prediction interval expressed as ×/÷ (fold-change for HIV, MIC dilutions for TB). `informative` = the model beats the marginal label (R²>0.05) — only these intervals carry genotype signal.

| pathogen | class/drug | label | n | R² | cover_90 | interval | calibrated | informative |
|---|---|---|---|---|---|---|---|---|
| HIV-1 | NNRTI/EFV | PhenoSense fold-change | 2168 | 0.8443 | **0.8973** | ×/÷3.52 | YES | YES |
| HIV-1 | PI/NFV | PhenoSense fold-change | 2133 | 0.8676 | **0.902** | ×/÷2.92 | YES | YES |
| HIV-1 | PI/IDV | PhenoSense fold-change | 2098 | 0.8814 | **0.9036** | ×/÷2.67 | YES | YES |
| HIV-1 | PI/SQV | PhenoSense fold-change | 2084 | 0.873 | **0.9038** | ×/÷3.05 | YES | YES |
| HIV-1 | NNRTI/NVP | PhenoSense fold-change | 2052 | 0.7746 | **0.9037** | ×/÷6.1 | YES | YES |
| HIV-1 | PI/FPV | PhenoSense fold-change | 2052 | 0.8628 | **0.8964** | ×/÷2.88 | YES | YES |
| HIV-1 | NRTI/AZT | PhenoSense fold-change | 1853 | 0.7837 | **0.8981** | ×/÷4.84 | YES | YES |
| HIV-1 | NRTI/DDI | PhenoSense fold-change | 1849 | 0.7136 | **0.9021** | ×/÷1.56 | YES | YES |
| HIV-1 | NRTI/D4T | PhenoSense fold-change | 1846 | 0.7296 | **0.9009** | ×/÷1.69 | YES | YES |
| HIV-1 | NRTI/3TC | PhenoSense fold-change | 1839 | 0.9098 | **0.8999** | ×/÷2.13 | YES | YES |
| HIV-1 | PI/LPV | PhenoSense fold-change | 1807 | 0.9086 | **0.903** | ×/÷2.67 | YES | YES |
| HIV-1 | NRTI/ABC | PhenoSense fold-change | 1731 | 0.8181 | **0.9016** | ×/÷1.8 | YES | YES |
| HIV-1 | NRTI/TDF | PhenoSense fold-change | 1548 | 0.6177 | **0.9032** | ×/÷1.96 | YES | YES |
| HIV-1 | PI/ATV | PhenoSense fold-change | 1505 | 0.8825 | **0.8983** | ×/÷2.78 | YES | YES |
| HIV-1 | PI/TPV | PhenoSense fold-change | 1226 | 0.7203 | **0.8997** | ×/÷2.85 | YES | YES |
| HIV-1 | NNRTI/ETR | PhenoSense fold-change | 998 | 0.7273 | **0.8978** | ×/÷3.6 | YES | YES |
| HIV-1 | PI/DRV | PhenoSense fold-change | 993 | 0.8684 | **0.9003** | ×/÷2.67 | YES | YES |
| HIV-1 | INI/EVG | PhenoSense fold-change | 754 | 0.8223 | **0.9086** | ×/÷3.05 | YES | YES |
| HIV-1 | INI/RAL | PhenoSense fold-change | 753 | 0.896 | **0.904** | ×/÷2.26 | YES | YES |
| HIV-1 | INI/DTG | PhenoSense fold-change | 370 | 0.4594 | **0.8941** | ×/÷3.53 | YES | YES |
| HIV-1 | NNRTI/RPV | PhenoSense fold-change | 311 | 0.642 | **0.9087** | ×/÷4.89 | YES | YES |
| HIV-1 | INI/BIC | PhenoSense fold-change | 287 | 0.4655 | **0.9094** | ×/÷2.85 | YES | YES |
| HIV-1 | NNRTI/DOR | PhenoSense fold-change | 130 | 0.8131 | **0.8931** | ×/÷4.0 | YES | YES |
| HIV-1 | INI/CAB | PhenoSense fold-change | 64 | 0.6137 | **0.9281** | ×/÷4.68 | YES | YES |
| HIV-1 | NNRTI/CompMutList | PhenoSense fold-change | 0 | — | — | — | under-powered | — |
| HIV-1 | NRTI/CompMutList | PhenoSense fold-change | 0 | — | — | — | under-powered | — |
| HIV-1 | PI/CompMutList | PhenoSense fold-change | 0 | — | — | — | under-powered | — |
| HIV-1 | INI/CompMutList | PhenoSense fold-change | 0 | — | — | — | under-powered | — |
| M. tuberculosis | moxifloxacin | CRyPTIC BMD-MIC | 12192 | 0.3699 | **0.9005** | ×/÷2.18 | YES | YES |
| M. tuberculosis | levofloxacin | CRyPTIC BMD-MIC | 12161 | 0.5259 | **0.9075** | ×/÷2.25 | YES | YES |
| M. tuberculosis | ethambutol | CRyPTIC BMD-MIC | 12156 | 0.4566 | **0.904** | ×/÷2.85 | YES | YES |
| M. tuberculosis | ethionamide | CRyPTIC BMD-MIC | 12130 | 0.2007 | **0.9084** | ×/÷3.07 | YES | YES |
| M. tuberculosis | rifampicin | CRyPTIC BMD-MIC | 12097 | 0.42 | **0.9035** | ×/÷3.14 | YES | YES |
| M. tuberculosis | isoniazid | CRyPTIC BMD-MIC | 12068 | 0.7661 | **0.9057** | ×/÷3.56 | YES | YES |
| M. tuberculosis | linezolid | CRyPTIC BMD-MIC | 12187 | 0.0071 | **0.9056** | ×/÷2.56 | YES | no (≈marginal) |
| M. tuberculosis | kanamycin | CRyPTIC BMD-MIC | 12128 | -0.2411 | **0.9169** | ×/÷2.16 | YES | no (≈marginal) |
| M. tuberculosis | amikacin | CRyPTIC BMD-MIC | 12070 | -0.7415 | **0.9117** | ×/÷2.83 | YES | no (≈marginal) |
| M. tuberculosis | bedaquiline | CRyPTIC BMD-MIC | 12066 | -0.0797 | **0.9073** | ×/÷3.3 | YES | no (≈marginal) |
| M. tuberculosis | clofazimine | CRyPTIC BMD-MIC | 12047 | — | — | — | under-powered (too few / no features / degenerate MIC) | — |
| M. tuberculosis | delamanid | CRyPTIC BMD-MIC | 11925 | — | — | — | under-powered (too few / no features / degenerate MIC) | — |

## Honest caveats
- CALIBRATED != INFORMATIVE: conformal coverage holds even for a useless model (the interval widens to the marginal label spread). `informative` (R2>0.05) flags the cells whose determinant model actually beats the mean — only those intervals carry genotype signal. The coverage-valid-only cells (some TB second-line drugs with few determinants + rare resistance) hit 0.90 coverage trivially.
- Every cell is IN-DISTRIBUTION vs its own knowledge base (HIV: Stanford catalog features; TB: WHO catalogue determinants) — the interval is a coverage guarantee, NOT an independent-validation claim.
- Split-conformal gives MARGINAL (population-level) coverage, not per-genotype. TB MIC coverage is on the RESOLVED (uncensored) subset; censored isolates are scored only for consistency.
- This map rolls up the QUANTITATIVE (interval) layer only — the deterministic R/S decoder + its provenance-disjoint report card are separate surfaces.