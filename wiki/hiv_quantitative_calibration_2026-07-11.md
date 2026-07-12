# HIV quantitative calibration — are the fold-change prediction intervals honest? (2026-07-11)

**Verdict: CALIBRATED_INTERVALS** — 24/24 powered drug-cells have a 90% interval within 0.05 of nominal coverage (fraction 1.0; PASS bar 0.5).

Does a split-conformal prediction interval on HIV fold-change achieve its nominal held-out coverage (is the quantitative decoder honestly calibrated)? Label = Stanford HIVDB PhenoSense fold-change (free, independent wet-lab). split-conformal on the additive nested-CV OOF residuals (Family A harness); finite-sample quantile; coverage tested on a held-out residual half; averaged over 20 shuffles.

`cover_90` = held-out fraction of isolates whose true fold falls in the 90% interval (target 0.90). `fold_factor_90` = the interval half-width in FOLD units (true fold within x/÷ this factor).

| class:drug | n | R2(oof) | **cover_90** | fold± (90%) | cover_80 | calibrated |
|---|---|---|---|---|---|---|
| NNRTI:EFV | 2168 | 0.8443 | **0.8973** | x/÷3.52 | 0.7925 | YES |
| NNRTI:NVP | 2052 | 0.7746 | **0.9037** | x/÷6.1 | 0.798 | YES |
| NNRTI:ETR | 998 | 0.7273 | **0.8978** | x/÷3.6 | 0.7935 | YES |
| NNRTI:RPV | 311 | 0.642 | **0.9087** | x/÷4.89 | 0.8003 | YES |
| NNRTI:DOR | 130 | 0.8131 | **0.8931** | x/÷4.0 | 0.8246 | YES |
| NNRTI:CompMutList | 0 | — | too few | — | — | — |
| NRTI:3TC | 1839 | 0.9098 | **0.8999** | x/÷2.13 | 0.802 | YES |
| NRTI:ABC | 1731 | 0.8181 | **0.9016** | x/÷1.8 | 0.8006 | YES |
| NRTI:AZT | 1853 | 0.7837 | **0.8981** | x/÷4.84 | 0.7948 | YES |
| NRTI:D4T | 1846 | 0.7296 | **0.9009** | x/÷1.69 | 0.7953 | YES |
| NRTI:DDI | 1849 | 0.7136 | **0.9021** | x/÷1.56 | 0.7986 | YES |
| NRTI:TDF | 1548 | 0.6177 | **0.9032** | x/÷1.96 | 0.7946 | YES |
| NRTI:CompMutList | 0 | — | too few | — | — | — |
| PI:FPV | 2052 | 0.8628 | **0.8964** | x/÷2.88 | 0.8002 | YES |
| PI:ATV | 1505 | 0.8825 | **0.8983** | x/÷2.78 | 0.8042 | YES |
| PI:IDV | 2098 | 0.8814 | **0.9036** | x/÷2.67 | 0.7961 | YES |
| PI:LPV | 1807 | 0.9086 | **0.903** | x/÷2.67 | 0.7955 | YES |
| PI:NFV | 2133 | 0.8676 | **0.902** | x/÷2.92 | 0.7992 | YES |
| PI:SQV | 2084 | 0.873 | **0.9038** | x/÷3.05 | 0.8011 | YES |
| PI:TPV | 1226 | 0.7203 | **0.8997** | x/÷2.85 | 0.8004 | YES |
| PI:DRV | 993 | 0.8684 | **0.9003** | x/÷2.67 | 0.8041 | YES |
| PI:CompMutList | 0 | — | too few | — | — | — |
| INI:RAL | 753 | 0.896 | **0.904** | x/÷2.26 | 0.7993 | YES |
| INI:EVG | 754 | 0.8223 | **0.9086** | x/÷3.05 | 0.8011 | YES |
| INI:DTG | 370 | 0.4594 | **0.8941** | x/÷3.53 | 0.8054 | YES |
| INI:BIC | 287 | 0.4655 | **0.9094** | x/÷2.85 | 0.7951 | YES |
| INI:CAB | 64 | 0.6137 | **0.9281** | x/÷4.68 | 0.8094 | YES |
| INI:CompMutList | 0 | — | too few | — | — | — |

## Honest caveats
- Split-conformal gives MARGINAL coverage (not conditional on the mutation profile).
- Censored folds ('>'/'<') kept at the numeric bound bias the tail slightly.
- Reuses Family A's additive OOF (no interaction terms — A showed they don't help rank; the point estimate is the additive model's).

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; conformal per Vovk/Lei split-conformal.