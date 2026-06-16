# External-validation report card — 2026-06-16

External clinical re-validation of the FROZEN decoder on independent measured-MIC cohorts (different country / lab / AST method than the US-NCBI-PD tuning provenance). **Strict tier (HIGH_R/HIGH_S) is the primary metric**; relaxed (+DECISIVE) is secondary. Raw sens/spec is clonality-inflated — the cluster-weighted block (one vote per lineage, Wilson CI) is the honest companion. This is SEPARATE from the frozen decoder report card.

| cohort | drug | metric | sens | spec | n | lineage-wt sens (CI) | lineage-wt spec (CI) | eff-N R/S | scope |
|---|---|---|---|---|---|---|---|---|---|
| oxford | trimethoprim-sulfamethoxazole | binary | 0.922 | 0.977 | 2866 | n/a | n/a | n/a | EXPERIMENTAL_SCORED (scorer_local) |
| sci234 | trimethoprim-sulfamethoxazole | strict | 0.986 | 0.993 | 223 | n/a | n/a | n/a | EXPERIMENTAL_SCORED (scorer_local) |
