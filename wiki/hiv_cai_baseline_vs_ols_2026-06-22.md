# HIV CAI v0 — catalog vs underlying tool (OLS) (2026-06-22)

Wrapper = dna_decode.data.hiv_amr (CAI catalog) (mutant-level). Gene = CA.
Underlying tool = Stanford DRMcv.R-style OLS (log10 fold ~ binary mutation presence), 5-fold CV, Python/sklearn reimpl (shared machinery with the NNRTI baseline).
Cutoff: UNIFORM illustrative fold>=3 (NOT per-drug clinical); both models scored at the same cutoff so the DELTA is the fair wrapper-vs-tool signal; absolute calibration needs per-drug clinical cutoffs (v0.1, not sourced here). OLS features (>=10 isolates) = 6.

| Drug | n | prev R | catalog sens/spec/**balacc** | OLS sens/spec/**balacc** (AUC) | d(OLS-cat) |
|---|---|---|---|---|---|
| lenacapavir | 140 | 0.793 | 0.982/0.31/**0.646** | 1.0/0.069/**0.534** (0.9105) | -0.112 |

**Interpretation:** small delta => the deterministic catalog MATCHES the full regression (adds interpretability, not error); large positive delta => per-drug/mutant signal the catalog misses (bounds a v0.1 refinement).

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; method Stanford DRMcv.R.