# HIV NNRTI v0 — validate the catalog vs the underlying tool (OLS) (2026-06-21)

Wrapper = dna_decode.data.hiv_amr (deterministic class-level major-DRM catalog).
Underlying tool = Stanford DRMcv.R OLS (log10 fold ~ binary mutation presence), 5-fold CV, reimplemented in Python/sklearn (R not installed).
Clinical cutoff source = DRMcv.R confusion-matrix lower cutoffs (EFV/NVP/ETR/RPV=3; DOR reused). OLS features (>=10 isolates) = 426.

Both scored on the SAME isolates + the SAME clinical cutoff (fold>=lower). balacc = balanced accuracy ((sens+spec)/2).

| Drug | n | prev R | catalog sens/spec/**balacc** | OLS sens/spec/**balacc** (AUC) | d(OLS-cat) |
|---|---|---|---|---|---|
| efavirenz | 2168 | 0.46 | 0.947/0.904/**0.926** | 0.931/0.951/**0.941** (0.9773) | 0.015 |
| nevirapine | 2052 | 0.513 | 0.906/0.991/**0.948** | 0.923/0.926/**0.925** (0.9729) | -0.023 |
| etravirine | 998 | 0.284 | 0.887/0.571/**0.729** | 0.813/0.885/**0.849** (0.9282) | 0.12 |
| rilpivirine | 311 | 0.379 | 0.822/0.544/**0.683** | 0.737/0.803/**0.77** (0.8671) | 0.087 |
| doravirine | 130 | 0.646 | 0.75/0.348/**0.549** | 0.774/0.761/**0.767** (0.8452) | 0.218 |

**Interpretation:** small delta => the simple deterministic catalog MATCHES the full regression (adds interpretability, not error); large delta on 2nd-gen NNRTIs => per-drug signal the class-level v0 misses (bounds the v0.1 per-drug catalog).

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; method Hedlin/Stanford DRMcv.R.