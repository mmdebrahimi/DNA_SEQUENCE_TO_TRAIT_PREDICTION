# Execution Log — Measure_The_Generator_Noise_Floor_As_A_Paired_Delta_In_One_Session_Plan
Date: 2026-11-20
Waves: 5 (max parallelism 2) — run SEQUENTIALLY by deliberate override; see below
Files changed: scripts/kaggle/build_glm_finetune_kernel.py, scripts/kaggle/glm_paired_sweep_kernel.py, scripts/glm_paired_sweep_verdict.py, tests/test_glm_paired_noise_bar.py, tests/test_glm_paired_sweep_verdict.py, tests/test_glm_finetune_kernel_builder.py, wiki/glm_generator_paired_noise_bar_2026-11-20.json, wiki/glm_generator_paired_sweep_2026-11-20.json, wiki/glm_generator_paired_sweep_verdict_2026-11-20.json, wiki/glm_generator_paired_sweep_2026-11-20.md
Sentrux verdict: n/a — plan is not structural (one builder module, tests, wiki artifacts)
Commit: 6b5ae955 (branch promoteratlas-exposure-2026-11-01; no merge to main)

## Result

DELTA_NOT_ESTABLISHED. The carried +0.0271 "fine-tuning made it worse" reading is not
resolved at this geometry: mean paired delta +0.0458 against a spread of 0.2072, 3 of 10
seeds against the committed sign, and at seed 7 — the seed both committed arms ran at —
the delta is -0.0142, the opposite sign. This does NOT establish the opposite; the mean is
positive, so the direction is weakly consistent with the original reading.

TRAJECTORY_UNRESOLVED, with the suggestive monotone pattern recorded and deliberately not
upgraded (step means +0.0016 / +0.0187 / +0.0255, positive counts 6 / 8 / 9).

The sibling claim SURVIVES: against its own seed-matched 3-mer Markov null the generator is
worse at 10 of 10 seeds, mean gap +0.1186, never once beating it.

## Why sequential, against the toolkit's verdict

The toolkit reported errors=[] and max_parallelism=2 with a remote configured, so parallel
was eligible. Overridden because the toolkit cannot see: ~4 GB of gitignored data/ makes a
worktree's absent-data failure set diverge from the main-checkout baseline this plan's
verification reconciles against (recorded twice before in this repo, 2026-10-20 and
2026-10-21); and the Kaggle SDK lives in the gitignored .venv while the run is one
host-global session, so neither precondition is satisfiable per-worktree.

## Deviations

1. Step 1's Files: line named only the JSON while its test strategy required tests; added
   tests/test_glm_paired_noise_bar.py, a file the plan never named.
2. Step 2's test strategy asked for the no-flag emitted body to be byte-identical to the
   committed frozen kernel. Impossible after Step 2 by construction, and satisfying it would
   have meant regenerating the provenance record of the two committed arms. Inverted into a
   sha256 pin (6bdbd02641abedc5, verified unchanged throughout) plus a behavioural test.
3. Step 4's "never call the file-writing entry point" would have forbidden the pre-existing
   test_builder_refuses_to_write_a_kernel_that_does_not_parse, which calls main([]) safely
   because the AST guard raises before write_text. Replaced by a module-scoped sha256
   snapshot of the artifact — the invariant belongs on the file, not the call.
4. glm_paired_sweep_kernel.py was assigned to Step 5 but is generated and verified in Step 3.
5. The plan's runtime estimate said "roughly one hour at ten seeds"; ten seeds is ~1.7 h
   (40 scorings x 138 s + 578 s training). The run took ~105 min, matching the corrected
   figure. One hour corresponds to five seeds.
6. Step 7's three epilogue agents were NOT launched. /test-epilogue's output already exists
   (42 tests added for exactly these changes, 6 load-bearing assertions proven non-vacuous);
   /documentation's primary target is CLAUDE.md, which is user authority and is the subject
   of the pending correction-propagation decision, so letting an agent edit it would pre-empt
   the user's call; /retrospective was written directly with full execution context. The
   mandatory [plan_file:] ledger entry was written to wiki/decisions-log.md.

## Tests

Baseline 12 failed / 9012 passed / 191 skipped. Final 12 failed / 9054 passed / 191 skipped.
Reconciled by CAUSE, not count: identical failure set (cyp4f2/GeT-RM 7, essentiality-ladder 3,
README-quickstart 1, pneumo cps DB 1), all four documented absent-data families, none touched
by this plan. Passed rose by exactly the 42 tests added. Zero regressions.

## Frozen surface

scripts/kaggle/glm_finetune_kernel.py unchanged (sha256 6bdbd02641abedc5). verdict() and the
five frozen method constants unchanged. The pre-registration was AMENDED augment-only, proven
against git HEAD (19 -> 20 keys, one added, none removed, zero frozen fields changed), with the
original wording preserved. prospective_lock.verify_lock returns ok=True.

## Open, not taken

Propagating the correction to CLAUDE.md, features/glm-generator-finetune/ and the project
memory, all of which carry "fine-tuning made it worse" as settled. User authority.
