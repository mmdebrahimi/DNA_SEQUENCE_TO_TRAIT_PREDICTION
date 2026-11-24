# Execution Log — Score The Checkpoint On Both Eval Sets, And Measure What The Orthology Filter Costs

Date: 2026-11-24
Waves: 3 declared (max parallelism 2) — EXECUTED SEQUENTIALLY. The toolkit reported
parallel-eligible and was overridden on two facts it cannot see: `data/` is gitignored
(`.gitignore:4` MG1655 cache, `:184` checkpoint) so a worktree would carry neither and Steps 2-3
would be impossible there; and Step 2 declared `Depends on: none` while its body calls Step 1's
gate, a real import dependency invisible to a parser that reads only the `Depends on:` line.
Files changed: .claude/testable-modules.md, CLAUDE.md, FUTURE_FEATURES.md, LESSONS_LEARNED.md,
README.md, TODOS.md, dna_decode/glm/eval_admissibility.py, scripts/glm_beyond_kmer.py,
scripts/glm_eval_set_inventory.py, scripts/glm_two_set_contrast.py, soraya_runs/arms-ledger.json,
tests/test_eval_admissibility.py, tests/test_glm_beyond_kmer_eval_sets.py,
tests/test_glm_eval_set_inventory.py, wiki/decisions-log.md,
wiki/glm_beyond_kmer_two_sets_2026-11-23.json, wiki/glm_eval_set_inventory_2026-11-23.json,
wiki/glm_two_set_reading_rule_2026-11-24.json, wiki/plans-index.md (plus this plan folder)
Sentrux verdict: clean — "No degradation detected", quality 4651 -> 4650, coupling 0.45 -> 0.45,
cycles 3 -> 3, god files 0 -> 0. Figures come from the `sentrux gate` run's stdout;
`.sentrux/.last-compare.txt` is 0 bytes and dated 2026-11-04, so this build does not write it. The
tracked `.sentrux/baseline.json` was restored after `--save` overwrote it.
Commit: 7bce1a4f (branch promoteratlas-exposure-2026-11-01; committed directly, no PR)

## Outcome, which is not the one the plan set out to produce

The plan was scoped to measuring the near-ortholog leak. It found the premise behind the committed
2026-11-22 result was wrong, in two stages, and the plan's own Step 4 outcome set had to be
re-frozen before the scoring run because all three of its named outcomes presupposed a refusal that
measurement had already undercut.

* |C| = 648 list entries / 634 DISTINCT sequences, not 6. The scorer refuses only when
  `len(of) < --n-windows`, so its own default of 120 would have run.
* Q1 ADMISSIBLE on all four sets at 3600 predictions: the committed
  EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM does NOT stand, and the accessory-gene out-of-distribution
  explanation recorded beside it does not describe a set its own null predicts above chance.
* The replacement is a RESOLVED NEGATIVE: the order-3 chain beats the model on 3 of 4 sets. The
  committed delta was +0.0192 (model ahead); at power it is -0.0191 (model behind). The sign flipped.
* Q2 UNRESOLVED_CONTRAST_SPANS_THE_MARGIN: -0.015372, CI [-0.030708, -0.001625] against a derived
  0.01-nat margin. The interval excludes zero but not the margin; only the first is supported.
* The re-pretrain arm is DECLARED and blocking under `generator-ceiling:`. Its rationale is weakened,
  not refuted; retiring it is a user authority call.

Tests: baseline 15 failed / 9579 passed / 172 skipped -> final 14 failed / 9708 passed / 172
skipped. Zero regressions; failure SET matches baseline minus one pre-existing failure this run
fixed. +129 tests. One 125-test regression was introduced mid-run (a second `ceiling:` criterion in
the shared arms ledger) and fixed; only the full suite could see it.
