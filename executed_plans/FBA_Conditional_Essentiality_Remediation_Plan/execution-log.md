# Execution Log — FBA_Conditional_Essentiality_Remediation_Plan

Date: 2026-08-13
Waves: 4 (max parallelism 4) — executed sequentially by choice; see below
Files changed: dna_decode/fba/solver_audit.py, dna_decode/fba/conditional_essentiality.py, dna_decode/fba/fitness_browser.py, scripts/fba_infeasibility_probe.py, scripts/fba_conditional_carbon_validate.py, scripts/fba_conditional_essentiality_validate.py, scripts/fba_regulatory_conditional_test.py, scripts/fba_gapfill_carbon_recheck.py, scripts/fba_gapfill_conditional_test.py, tests/test_fba_solver_audit.py, tests/test_fba_fitness_browser_t.py, tests/test_fba_conditional_essentiality.py, tests/test_fba_regulatory_conditional.py, wiki/fba_infeasibility_finding_2026-08-13.md, wiki/fba_conditional_carbon_2026-08-12.md, wiki/fba_gapfill_conditional_answer_2026-08-12.md, wiki/fba_regulatory_conditional_test_2026-08-12.json, CLAUDE.md, LESSONS_LEARNED.md, wiki/decisions-log.md
Sentrux verdict: n/a — sentrux not installed
Commit: 0deb874 (findings) / 3dd088c (docs); waves 0-1 at 3a49024, 6e832a6

## As-built divergence

**The plan's premise was refuted by its own Step 8.** The plan was written to test whether the pFBA
regulatory lift was a solver artifact, with a pre-committed verdict set. The rule fired
`REGULATORY_LIFT_IS_A_SOLVER_ARTIFACT` — and that label was wrong. A re-solve probe (not in the plan)
showed a non-optimal solve here is deterministic ATPM-maintenance infeasibility, i.e. genuine
essentiality: 39/39 deterministic, each the canonical catabolic gene for its carbon source, 38/39
experimentally essential there. Corrected verdict:
`REGULATORY_LIFT_STANDS_ABSTENTION_IS_A_BIASED_LOWER_BOUND`.

Divergences from the plan as written:
- **Added, unplanned:** `scripts/fba_infeasibility_probe.py` + 4 tests — the decisive check, now the
  load-bearing evidence for every other claim in this cluster.
- **Added, unplanned:** `exclude_cells` on `switch_accuracy` / `constant_baselines`. Step 5 required
  abstention but the plan assumed dropping a cell from `predicted` would abstain it. It does not — the
  metric fails OPEN, scoring the dropped cell as dispensable and keeping it in the denominator.
- **Added, unplanned:** the 4-media validator (`fba_conditional_essentiality_validate.py`) was wired to
  the shared audit. The plan's Step 6 named only the two gap-fill scripts, but Verification Signal
  item 2 required all five.
- **Added, unplanned:** repaired `wiki/fba_regulatory_conditional_test_2026-08-12.json`, invalid JSON
  since it was written, plus a guard test that every FBA artifact parses.
- **Step 8 artifact naming:** `fba_regulatory_conditional_recheck_<date>.json` supersedes rather than
  overwrites, as planned. The 08-12 artifact is annotated, not rewritten.
- **Execution mode:** sequential, though the toolkit reported parallel-eligible. Repo convention is
  direct commits to main; waves 2-3 need a live solver and the external D: drive in the parent checkout.
- **Not done (scoped out in the plan, still open):** the label-threshold sweep. Step 3 made
  `GeneFitness.t` reachable (`load_records(..., min_abs_t=)`, `mean_t_matrix`); the sweep itself remains
  its own plan.

## Verification

- Suite: 3,407 passed / 10 skipped (+40 tests over the 3,367 baseline). One pre-existing failure,
  `test_models_foundation.py::test_nt_embed_window_batch_matches_per_sequence` — transformers version
  drift (`find_pruneable_heads_and_indices` no longer exported), unrelated to fba/.
- Carbon headline reproduced exactly: exact-set 23/217, per-cell 0.7368, constant 184/217 = 84.8%.
- Regulatory default-coding arm reproduced exactly: per-cell 0.6157, TP 56, FP 40.
- All five deletion scripts audit solver status; all 21 FBA artifacts parse.
