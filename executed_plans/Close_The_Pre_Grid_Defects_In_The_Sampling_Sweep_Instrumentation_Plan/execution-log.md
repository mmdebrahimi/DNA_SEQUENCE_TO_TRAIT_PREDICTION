# Execution Log — Close the Pre-Grid Defects in the Sampling-Sweep Instrumentation
Date: 2026-10-21
Waves: 3 (max parallelism 3) — executed SEQUENTIALLY by override: `data/` is gitignored, so a fresh worktree lacks the MG1655 reference, and Step 6's only signal is reconciling the full-suite count against a baseline measured in this checkout.
Files changed: scripts/glm_generator_falsifier.py, dna_decode/glm/sweep_gates.py, scripts/glm_run_lock.py, run_glm_anchor_detached.bat, tests/test_glm_geometry_contract.py, tests/test_glm_sweep_gates.py, tests/test_glm_run_lock.py, tests/test_glm_knob_probe_path.py, plans/Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan/technical-plan.md, CHANGELOG.md, LESSONS_LEARNED.md, TODOS.md, NEXT.md, FUTURE_FEATURES.md, wiki/decisions-log.md, wiki/plans-index.md, .claude/testable-modules.md
Sentrux verdict: n/a — sentrux not installed
Commit: 937c8ae4 (main; six plan commits af1f978, 2d4aba9, 5bbbeed, 4b79b25, 0add19a, af9d878 plus the epilogue)

## What the execution itself found

Three defects in work committed during this same plan, each found by running it rather than reading it:

1. **The lock could not start a job.** The `run` form the launcher uses was rejected by argparse
   (`unrecognized arguments: -- uv run python ...`, exit 2), so the guarded command never ran and
   the sentinel would have said `EXIT=2` — "the job failed" rather than "never launched". All 21
   lock tests were green and structurally had to be: every one crosses `run_guarded` in-process,
   never the argparse layer the launcher goes through.
2. **That defect is INTERPRETER-DEPENDENT**, and the first write-up of it was incomplete in a way
   that made it look false. It reproduces under the venv (3.11.17, what `uv run` and the suite
   use) and not under the 3.14.7 on PATH. Five verification attempts under bare `python` all
   passed. The launcher now pins `uv run python` so it no longer depends on PATH order.
3. **`acquire` + `release` is the wrong shape for a launcher** — `acquire` records the pid of the
   process that calls it, which exits at once, so the lock was reclaimed as stale and a second
   job started (measured: exit 0, job ran). Found by the manual two-launch check; the in-process
   tests could not see it. Fixed with a `run` wrapper holding the lock for the child's lifetime.

Five of this execution's own guards were initially vacuous or matched the wrong text: two matched
the launcher's REM comments documenting the traps, one matched the `import` of `verify_knob_live`
instead of its call (reporting an ordering backwards), one asserted a plan string that already
existed before the amendment, and one mutation-based non-vacuity check used a `sed` that matched
nothing and exited 0 — briefly making a good test look vacuous. Each is fixed and each fix is
proven by a mutation confirmed to have landed.

## Deviations from the plan

- The sentinel fix (`echo EXIT=%RC%>>` redirecting handle 0 and leaving the file empty) was not in
  the plan. Same completion-signalling defect in the same file as Step 5, so it landed there.
- The plan's own Pre-Change Baseline (15/6726/185) did not match the measured baseline
  (14/6759/186). Explained, not adopted silently: the 25 geometry cases executing for the first
  time, torch now installed, and one data-gated test moving failed→skipped with `E:` attached.
- Step 2's file briefly carried Step 4's assertion; deferred rather than committing a known-red test.
