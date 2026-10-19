# Execution Log — F3_Phylogeny_Deconfounding_Calibrated_Instrument_Plan
Date: 2026-10-09
Waves: 4 (max parallelism 3) — executed SEQUENTIALLY in order 2,3,1,4,5,6; the plan's own Execution Preview records that Step 1 consumes Step 3's module at runtime, so the three Wave-0 steps are not independent in separate worktrees
Files changed: dna_decode/eval/genomic_prediction.py, dna_decode/eval/panel_structure.py (new), dna_decode/eval/regime.py, scripts/f3_substrate_census.py (new), scripts/f3_deconfound_calibration.py (new), tests/test_genomic_prediction_grouped.py (new), tests/test_panel_structure.py (new), tests/test_f3_substrate_census.py (new), tests/test_f3_deconfound_calibration.py (new), tests/test_regime_boundary.py, wiki/f3_substrate_census_2026-10-09.json, wiki/f3_deconfound_calibration_2026-10-09.{md,json}, wiki/plans-index.md, wiki/decisions-log.md, CLAUDE.md, CHANGELOG.md, LESSONS_LEARNED.md, README.md, TODOS.md, NEXT.md, .claude/testable-modules.md
Sentrux verdict: n/a — sentrux not installed
Commit: 81dc31d (+ the archive commit; sequential mode, merged directly to main, no PRs)

## Outcome
Verdict `blocked:positive_control_not_fold_compatible` — a PRE-REGISTERED outcome, not a
failure. The instrument is BUILT and NOT CALIBRATED; no calibration metrics were emitted.
Tests 6238 -> 6365 (+127: 99 F3 + 6 regime + 22 from /test-epilogue), 0 regressions.
Frozen five byte-unchanged throughout; prospective lock re-verified.

## Deviations from the plan as written
1. Step 4's control-selection rule was changed from "highest measured pc1_variance_fraction"
   to selection BY CONSTRUCTION from each panel's `design` field, because pc1 was MEASURED to
   be unfit: under both raw pc1 and excess-over-its-own-null the one natural population ranks
   last-or-third while the top panel is the bi-parental cross that must be the NEGATIVE
   control. pc1 sums LD and structure, and a cross maximises LD by construction.
2. Step 3's fixed low percentile was replaced by a guard-selected sweep, because a fixed rung
   both manufactures clusters on an unstructured panel and shatters a structured one.
3. A verdict state the plan did not anticipate was added:
   `blocked:positive_control_not_fold_compatible`.
4. `wiki/decisions-log.md` was written by /retrospective rather than in Step 6, resolving the
   file-ownership overlap in favour of the post-execution chain's documented owner.
