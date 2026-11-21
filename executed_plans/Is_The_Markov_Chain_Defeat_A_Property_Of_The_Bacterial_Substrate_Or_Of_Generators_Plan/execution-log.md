# Execution Log — Is_The_Markov_Chain_Defeat_A_Property_Of_The_Bacterial_Substrate_Or_Of_Generators_Plan
Date: 2026-11-21
Waves: 5 (max parallelism 2) — executed SEQUENTIALLY by deliberate override
Files changed: dna_decode/glm/human_natural.py,scripts/human_adaptor_probe.py,scripts/human_falsifier_margin.py,tests/test_human_falsifier_bar.py,tests/test_human_falsifier_verdict.py,tests/test_human_natural.py,wiki/human_adaptor_probe_2026-11-21.json,wiki/human_falsifier_bar_2026-11-21.json,wiki/human_falsifier_margin_2026-11-21.json,wiki/human_falsifier_margin_2026-11-21.md,wiki/human_falsifier_substrate_2026-11-21.json,wiki/human_falsifier_verdict_2026-11-21.json,wiki/plans-index.md
Sentrux verdict: n/a — sentrux not installed
Commit: 93b89949 (pushed to origin/promoteratlas-exposure-2026-11-01)

## Mode override rationale
The toolkit's wave grouping was eligible for parallel execution. Overridden to sequential
because `data/mpra_human` is gitignored (`/data/*`, verified with `git check-ignore --no-index -v`),
so a worktree would contain NO substrate and Steps 1, 3, 4 and 5 would all fail outright.
Third recorded instance of gitignored-data baseline divergence in this repo.

## Outcome
Primary geometry (frozen in advance): CHAIN_IS_WEAKER_ON_HUMAN, unanimous — and measured to be
an adaptor artifact, not a statement about regulatory structure.
Secondary geometry: CHAIN_IS_COMPARABLE, because the two cells of one panel contradict.
The plan's question is NOT ESTABLISHED either way, which is a legitimate pre-registered outcome.

## Deviations
- Step 6's Files line listed only the memo. Also written: the verdict JSON, a verdict function
  added to the Step 5 script, and the verdict test module — all required by the plan's own test
  strategy, which asks for a test recomputing the verdict from the committed grid.
- Added beyond the plan: `scripts/human_adaptor_probe.py` + `wiki/human_adaptor_probe_2026-11-21.json`,
  because the memo's central correction rested on an UNMEASURED premise (that the flanks are
  constant) and its obvious rival (a 13% shorter-sequence effect) was unaddressed. Both now measured.

## Final suite, reconciled by CAUSE
Baseline: 12 failed / 9054 passed / 191 skipped (738 s)
Final:    12 failed / 9099 passed / 191 skipped (859 s)

Identical failure set, file-for-file and count-for-count: 7 `test_cyp4f2_getrm_concordance.py`,
3 `test_essentiality_report_card_ladder.py`, 1 `test_readme_quickstart_output.py`,
1 `test_pneumo_selection_probe.py`. All four are absent DATA, none a code defect, and all four
predate this plan. **ZERO regressions**; +45 passing.

One accounting discrepancy, chased rather than waved at. 9054 + 49 new tests should be 9103, not
9099. The suite collected 9302 while running; the tree now collects 9298. The 4-test gap is the
plan-hygiene family — the suite began collecting BEFORE the plan folder was archived and
`wiki/plans-index.md` was updated, and this repo has tests that read both, so parametrized collection
legitimately changed underneath the run. Verified benign rather than asserted:
`test_plan_state_hygiene.py` + `test_claude_md_citations.py` are green post-archive (22 passed), so
nothing broke on the archive and no test was lost.

Run state cleaned up (the state file is gitignored, so its removal is not a tracked change).
No worktrees to remove — sequential by deliberate override.
