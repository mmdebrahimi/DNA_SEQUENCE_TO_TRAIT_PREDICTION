# HIV PI v0 — catalog vs underlying tool (OLS) (2026-06-22)

Wrapper = dna_decode.data.hiv_amr (PI catalog) (position-based). Gene = PR.
Underlying tool = Stanford DRMcv.R-style OLS (log10 fold ~ binary mutation presence), 5-fold CV, Python/sklearn reimpl (shared machinery with the NNRTI baseline).
Cutoff: UNIFORM illustrative fold>=3 (NOT per-drug clinical); both models scored at the same cutoff so the DELTA is the fair wrapper-vs-tool signal; absolute calibration needs per-drug clinical cutoffs (v0.1, not sourced here). OLS features (>=10 isolates) = 224.

| Drug | n | prev R | catalog sens/spec/**balacc** | OLS sens/spec/**balacc** (AUC) | d(OLS-cat) |
|---|---|---|---|---|---|
| fosamprenavir | 2052 | 0.431 | 0.998/0.605/**0.802** | 0.922/0.937/**0.929** (0.9801) | 0.127 |
| atazanavir | 1505 | 0.529 | 0.999/0.712/**0.856** | 0.975/0.913/**0.944** (0.9891) | 0.088 |
| indinavir | 2098 | 0.5 | 0.996/0.706/**0.851** | 0.958/0.944/**0.951** (0.99) | 0.1 |
| lopinavir | 1807 | 0.521 | 0.998/0.709/**0.853** | 0.965/0.948/**0.956** (0.9916) | 0.103 |
| nelfinavir | 2133 | 0.572 | 0.987/0.805/**0.896** | 0.962/0.952/**0.957** (0.9893) | 0.061 |
| saquinavir | 2084 | 0.43 | 1.0/0.614/**0.807** | 0.961/0.914/**0.938** (0.9868) | 0.131 |
| tipranavir | 1226 | 0.249 | 0.997/0.469/**0.733** | 0.803/0.911/**0.857** (0.9399) | 0.124 |
| darunavir | 993 | 0.307 | 1.0/0.608/**0.804** | 0.915/0.92/**0.917** (0.9777) | 0.113 |

**Interpretation:** small delta => the deterministic catalog MATCHES the full regression (adds interpretability, not error); large positive delta => per-drug/mutant signal the catalog misses (bounds a v0.1 refinement).

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; method Stanford DRMcv.R.