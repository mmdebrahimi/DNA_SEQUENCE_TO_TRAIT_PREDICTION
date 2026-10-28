# Execution Log — Correct_Two_Wrong_Claims_I_Committed_And_Make_The_Next_Check_Runnable_Plan
Date: 2026-10-28
Waves: 3 (max parallelism 4) — run SEQUENTIAL; Step 5 writes outside the repo
  (`~/.claude/soraya-resume/`), where a worktree gives zero isolation, and the plan forbids
  concurrent execution for exactly that reason.
Scope: FULL (all 8 steps). The narrow option — Steps 1,2,4,5,7 only — was declined on the
  user's standing directive to take the fuller path.
Files changed: dna_decode/glm/lentimpra.py, tests/test_lentimpra.py,
  wiki/lentimpra_metric_resolution_2026-10-28.md
Tests: baseline test_lentimpra 4 passed / 18 skipped -> 8 passed / 18 skipped.
  Final pass 833 passed / 38 skipped across lentimpra, glm, pyramid, tile, frozen-bar,
  arms-ledger. No regressions.
Sentrux verdict: n/a — not installed.
Commit: 1a9a0ef3

## Deviations from the plan, and why

- **Step 3 went further than its letter.** The plan specified extracting a CLASSIFIER taking
  `(first, cells, cur, methods)`. A classifier alone leaves each caller carrying its own copy
  of the state transitions, which is the same duplication one level down. The whole WALK was
  extracted instead (`_walk_rows`), so the carried state lives once and both accessors are
  pure consumers. Satisfies the plan's intent (one rule, shared) more completely than its text.
- **Step 3's source-shape assertion had to change form.** The plan said to assert
  `published_incumbents` contains no `in CELL_TYPES` / `== "test fold"` test. As a substring
  check that can never pass: Step 2's corrected prose must QUOTE the header
  (`test fold | <method names>`) to explain why reading it cannot settle the metric. Asserted
  on the AST instead — only a literal in an executable `Compare` counts; a literal in a
  docstring or a returned dict is prose and is exempt structurally.
- **Step 4 found a third instance.** The header-row claim appeared in the findings list as
  well, calling the header "the one place that could state the metric". Caught by an assertion
  after the first two were fixed, and corrected in place with the superseded text labelled.
- **Step 5 was largely already discharged.** An earlier checkpoint rewrite this session had
  already removed the provenance framing and the header-row instruction; only the
  consumption-rule `facts[]` entry remained.
- **Step 7 ended in a report, as designed.** The arms ledger does NOT assert the arm is
  metric-blocked — its rationale already scopes the metric to comparison. Its md5
  (`09674340d04df6e8679b5997d064bdb2`) is asserted byte-unchanged before and after. Surfaced
  and left alone: that rationale carries the OTHER over-claim ("the sheet NEVER STATES THE
  METRIC"), whose correction is a user authority call.
- **Step 8 was a one-file no-op by half.** Step 7 named exactly one prose surface, inside the
  declared `Files:`. `wiki/lentimpra_acquisition_2026-10-27.md` needed no edit — its wording
  was already comparison-scoped — and per the plan a no-op was not padded into one.
- **Step 6 did not settle table identity.** Three routes; only one result is evidence about
  the source. Recorded attempted-and-open, with the Europe PMC 503-parsed-as-zero-hits trap
  written down.
