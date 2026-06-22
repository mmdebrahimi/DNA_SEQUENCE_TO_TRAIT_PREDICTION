# audit-trail — 2026-06-22-adv1-hiv-ols-baseline

| # | action | gate | result |
|---|---|---|---|
| 1 | assess situation (read ledger Bellman frame + candidate-actions + soraya_runs + git) | auto | AMR banked; forward paths are user decisions; selected the wrapper-vs-tool rigor gap |
| 2 | write scripts/hiv_targetsite_baseline.py (reuse hiv_nnrti_baseline OLS machinery) | auto | created |
| 3 | run baseline for PI/INSTI/CAI (datasets in hand, gitignored) | auto | 3 wiki/*_baseline_vs_ols_*.{md,json} written |
| 4 | verify-in-batch: inspect catalog-vs-OLS deltas | auto | genuine differentiated findings; catalog half consistent w/ shipped validation |
| 5 | wire OLS columns into build_hiv_report_card.py + caveat | auto | 14 PI/INSTI/CAI cells now carry OLS baseline |
| 6 | add report-card regression test | auto | test_pi_insti_cai_cells_carry_ols_baseline |
| 7 | full suite + frozen-surface check | auto | 1568 passed; amr_rules/calibrated_amr_rules byte-unchanged |
| 8 | commit + ledger row 149 | auto | (see recommendation.md) |
