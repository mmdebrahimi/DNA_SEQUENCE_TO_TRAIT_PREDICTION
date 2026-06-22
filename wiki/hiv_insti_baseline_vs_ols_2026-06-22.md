# HIV INSTI v0 — catalog vs underlying tool (OLS) (2026-06-22)

Wrapper = dna_decode.data.hiv_amr (INSTI catalog) (position-based). Gene = IN.
Underlying tool = Stanford DRMcv.R-style OLS (log10 fold ~ binary mutation presence), 5-fold CV, Python/sklearn reimpl (shared machinery with the NNRTI baseline).
Cutoff: UNIFORM illustrative fold>=3 (NOT per-drug clinical); both models scored at the same cutoff so the DELTA is the fair wrapper-vs-tool signal; absolute calibration needs per-drug clinical cutoffs (v0.1, not sourced here). OLS features (>=10 isolates) = 136.

| Drug | n | prev R | catalog sens/spec/**balacc** | OLS sens/spec/**balacc** (AUC) | d(OLS-cat) |
|---|---|---|---|---|---|
| raltegravir | 753 | 0.413 | 0.939/0.814/**0.877** | 0.942/0.95/**0.946** (0.9799) | 0.069 |
| elvitegravir | 754 | 0.5 | 0.87/0.857/**0.863** | 0.878/0.796/**0.837** (0.9453) | -0.026 |
| dolutegravir | 370 | 0.195 | 1.0/0.534/**0.767** | 0.764/0.923/**0.843** (0.9185) | 0.076 |
| bictegravir | 287 | 0.178 | 1.0/0.419/**0.71** | 0.471/0.928/**0.699** (0.8666) | -0.011 |
| cabotegravir | 64 | 0.812 | 1.0/0.333/**0.667** | 0.962/1.0/**0.981** (0.9712) | 0.314 |

**Interpretation:** small delta => the deterministic catalog MATCHES the full regression (adds interpretability, not error); large positive delta => per-drug/mutant signal the catalog misses (bounds a v0.1 refinement).

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; method Stanford DRMcv.R.