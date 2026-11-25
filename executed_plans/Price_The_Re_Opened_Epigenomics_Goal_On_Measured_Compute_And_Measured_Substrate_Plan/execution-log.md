# Execution Log — Price_The_Re_Opened_Epigenomics_Goal_On_Measured_Compute_And_Measured_Substrate_Plan

Date: 2026-11-24
Waves: 4 (max parallelism 2 — reported parallel-eligible by the toolkit and deliberately run sequentially)
Files changed: scripts/encode_substrate_census.py, tests/test_encode_substrate_census.py, tests/test_encode_roadmap_budget_sources.py, tests/test_encode_roadmap_gate_screen.py, wiki/encode_substrate_census_2026-11-24.json, wiki/encode_roadmap_training_budget_sources.json, wiki/encode_roadmap_cloud_sizing.md, wiki/encode_roadmap_gate_screen_2026-11-24.json, wiki/baseline_junit_2026-11-24/chunk{0..3}.xml, soraya_runs/constraint-expiry.json, project_state/glm-umbrella-2026-10-07.md, wiki/decisions-log.md, LESSONS_LEARNED.md, README.md, TODOS.md, .claude/testable-modules.md
Sentrux verdict: n/a — gate not run; artifact-and-data plan, and `.sentrux/.last-compare.txt` is 0 bytes with an mtime of 2026-11-04
Commit: n/a — branch only (`promoteratlas-exposure-2026-11-01`); 9 commits `568e6a9e`, `07d77269`, `615c5467`, `15998565`, `847a9e34`, `1e8063d9`, `876ac0e0`, `0bda0bfd`, `c7c2028c`

---

## What the plan was for, and what it concluded

It priced a goal abandoned on 2026-05-12 — a tissue-specific gene-regulatory-element activity
predictor, Enformer-shaped, on ENCODE / Roadmap multi-tissue data — without building or training
anything. **Verdict: compute cleared, data open.** The reference's full 196,608 bp context fits on a
free T4 at 1,024 channels and 134.63M parameters in 10.86 of 14.56 GiB, so the original cut (made on
a 4 GiB local card) no longer describes the world. The one remaining blocker is a storage figure that
needs reading in a browser.

## Deviations from the plan as written

Four, all recorded rather than silently absorbed:

1. **The plan's premise was false.** Step 5 was to flip a registered constraint's verdict. That
   constraint was never in the register — commit `cfdf4452`, whose message says it registered it,
   changed one line of an unrelated ledger. Cause: `ConstraintRegister.register()` mutates and
   returns a populated object without persisting. It was registered properly here instead.
2. **The requested derived ratio is not computable.** "Reference training compute / one measured T4
   session" needs a wall-clock or step count, and neither is published (a scan of all 78,591
   readable characters found zero duration-shaped expressions). Recorded `unsourced` with no number;
   a batch-size normaliser ships in its place.
3. **Three of the plan's own factual claims were corrected** from the artifacts: a parameter-count
   comparison with no sourced right-hand side, a single AMP speedup figure that varies 14.0x to 4.2x
   by shape, and a regime cell assigning an organism endpoint to a molecular readout.
4. **The frozen Verification Signal was not satisfiable.** It asked for a full-suite failure-SET
   comparison; the baseline crashed this host at ~28 minutes on an unmodified tree. Splitting into
   four chunks with per-chunk junit XML worked (9,924 unique test ids, zero overlap) and the XMLs are
   archived at `wiki/baseline_junit_2026-11-24/`.

## Result

Tests: 9,924 collected / 9,738 passed / **14 failed** / 0 errors / 172 skipped. All 14 pre-existing;
**zero** in the three test files this plan added. 54 new tests across those three files.

Register: 16 → 18 rows. `encode-roadmap-regulatory-goal-unsized-on-cloud` → **EXPIRED**;
`regulatory-training-tensor-vs-unverified-cloud-storage` → **still-binding**, with an expiry test
that fires the moment `wiki/kaggle_storage_limits.md` is written.

One published claim was retracted mid-execution: commit `615c5467` asserted a 64,046.7 GiB data
wall; `15998565` retracted it as a 212x overstatement with both the numerator and the denominator the
wrong quantity.
