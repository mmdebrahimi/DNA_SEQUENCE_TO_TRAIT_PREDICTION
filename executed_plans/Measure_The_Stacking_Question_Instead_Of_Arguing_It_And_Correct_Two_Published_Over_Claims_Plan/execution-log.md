# Execution Log — Measure The Stacking Question Instead Of Arguing It, And Correct Two Published Over-Claims
Date: 2026-10-26
Waves: 6 (max parallelism 2) — executed SEQUENTIALLY; the plan's own Execution Preview and its save-time amendment both specify sequential, and `data/promcalc/tile_predictions/` is gitignored so a worktree baseline would be unreconcilable.
Files changed: scripts/glm_tile_hybrid_arm.py, tests/test_glm_tile_hybrid_arm.py, tests/test_tile_chain_substrate_wall_amendment5.py, soraya_runs/arms-ledger.json, project_state/glm-umbrella-2026-10-07.md, wiki/glm_tile_hybrid_arm_2026-10-26.json, wiki/glm_tile_hybrid_arm_2026-10-26.md, wiki/glm_tile_stack_arm_2026-10-26.md (new), wiki/tile_reproduction_scope_correction_2026-10-26.md (new), wiki/decisions-log.md, wiki/plans-index.md, CHANGELOG.md, LESSONS_LEARNED.md, TODOS.md, NEXT.md, FUTURE_FEATURES.md
Sentrux verdict: n/a — sentrux not installed
Commit: n/a — working tree only, nothing committed or pushed by this execution

## Outcome

Tests 73 -> 135 passed, 0 failed, 0 regressions. `tests/test_glm_tile_hybrid_arm.py` 36 -> 63 tests.
Frozen five + `uv.lock` byte-unchanged (verified `git status --porcelain` on exactly those paths: empty).

**The pre-registered prediction was FALSIFIED.** Predicted `NO_STACK` for both new combos; measured
`WEAK_STACK`. `decomposition`-as-third gains +0.0145 over the 2-channel hybrid on 10/10 seeds,
CI95 [0.0108, 0.0172] excluding zero — under the inherited 0.02 margin, which was deliberately NOT
revisited. `gc`+`decomposition` +0.0103 (9/10, CI excludes zero); `gc` alone -0.0039 (4/10, CI
includes zero).

**The orthogonality argument is backwards, not merely insufficient.** `decomposition` is the MOST
redundant channel with the comparator (+0.7812) and is the one that ADDS; `gc` is markedly LESS
redundant (+0.5508 / +0.4593) and SUBTRACTS. Pairwise Spearman between prediction vectors does not
order these channels by whether they stack.

**Regression clean.** Every published 2-channel value byte-identical to `git show HEAD:` —
`_paired_hybrid_vs_comparator` +0.0339 CI95 [0.0238, 0.0409], verdict `HYBRID_ADDS`,
`orthogonality.median` 0.3862.

## Deviations from the plan as written

1. **Step 2's missing-scorer guard was NOT extended to the union over `COMBOS`.** Taken literally
   that would make a conv+max_tx+gc table refuse outright once `decomposition` entered `COMBOS` — so
   the ANSWERED, published question (+0.0339) would become unreproducible from any pre-decomposition
   substrate, and three pre-existing `main()` tests would flip. `required_channels()` is the stacking
   comparator's channels only; the hazard the union reached for (a combo silently scored on a seed
   SUBSET) is closed by three layers instead: a uniform all-or-nothing recorded skip
   (`_skipped_combos`), a per-table defensive guard, and a `_short_arm` refusal. Pinned by
   `test_main_SKIPS_an_optional_combo_and_RECORDS_why_rather_than_refusing_the_whole_arm`.

2. **Step 7's two pre-committed branches did not cover the outcome.** The plan pre-committed
   narrow-and-keep-`planned` on a fail and `measure_arm` on a success; `WEAK_STACK` is neither.
   Recorded as `measured` (metric 0.4404) with the honest verdict in `evidence` — leaving it
   `planned` would hide executed evidence, calling it a success would overstate it. The ceiling
   criterion is still blocked, now through three remaining `planned` siblings, verified by calling
   `blocking()` rather than assumed.

## Plan defect found during execution

**Step 6 carried `Depends on: none` but its content depended on Step 5's result.** It sat in wave 0
and so was written before the run, and its "what survives" section asserted that the CONCLUSION ("no
viable third channel exists in these four scorers") still had support even though its stated reason
did not — which Step 5 then refuted. The wave computation checks file overlap and declared edges; it
cannot see that a memo's PROSE depends on a measurement step. Both artifacts were uncommitted, so
both were corrected in place and the defect recorded in umbrella Action Log row 18 rather than
patched around silently.
