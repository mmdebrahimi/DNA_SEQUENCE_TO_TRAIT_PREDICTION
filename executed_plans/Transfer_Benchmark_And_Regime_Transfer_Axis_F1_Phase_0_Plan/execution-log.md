# Execution Log — Transfer_Benchmark_And_Regime_Transfer_Axis_F1_Phase_0_Plan

Date: 2026-10-03
Waves: 6 (max parallelism 1 — sequential mode; no git remote configured for worktree PRs)
Files changed: dna_decode/eval/regime.py, dna_decode/eval/transfer.py, dna_decode/eval/transfer_gauntlet.py, scripts/regime_map.py, scripts/transfer_benchmark_falsifier.py, wiki/phase2_arabidopsis_result_2026-06-12.json, tests/test_regime_boundary.py, tests/test_eval_transfer.py, tests/test_transfer_gauntlet.py, tests/test_transfer_benchmark_falsifier.py, tests/test_g2_replay_fixture.py
Sentrux verdict: n/a — sentrux not installed
Commit: 61a4b12 (archive); defect fixes + test rewrites in the preceding commit on `main`

## Status: as-built

Execution diverged from the plan in four recorded ways. Each is a deviation from the
saved text, not a shortfall against the plan's purpose.

### 1. The substrate file the plan named is CORRUPT

The plan pinned `BYxRM_GenoData.txt` for Arm A1. It is truncated — exactly 1 MiB,
511 of 11,623 rows, and its last row carries 626 fields against a 1009-field header.
`geno_direct.txt` beside it is a 199-byte HTTP 403 page under a genotype filename.
`geno_v2.txt` (23.7 MB) is intact and is what shipped as `BLOOM_GENO`.

The plan's own "511 markers" figure was read off the truncated file and is wrong.
`assert_genotype_file_intact` was added for exactly this: a ragged-row check that
refuses rather than scoring a partial panel, since a truncated genotype file scores
cleanly and silently.

### 2. The cross-panel transfer control was DROPPED, not deferred

The plan proposed holding out a whole species panel as a transfer positive. The result
is knowable without running it: the panels share no feature ontology, so this is a
representation-incompatibility negative, not a transfer measurement. Running it would
have produced a number that reads as a transfer failure and is actually a schema
mismatch. Dropped outright — the honest form is that no transfer-positive arm exists on
disk, which is why `TRANSFER_POSITIVE_UNMEASURED` is a first-class verdict.

I mis-parsed the panel headers three times while chasing this (comment lines,
transposition, individual-vs-marker IDs) before concluding the string-overlap
measurement was the wrong instrument.

### 3. `structure_only_baseline` was demoted from rival to diagnostic, on a measurement

The plan listed it among the verdict's rival controls. Under leave-one-group-out the
held-out group's one-hot column is all-zero in every training row, so the model cannot
represent that group's offset at all — it scores worse the more structure explains
(−0.777 where group explains ~everything vs −0.001 where it explains nothing),
presenting its easiest bar exactly where it was added to bite. The Arabidopsis
comparison that motivated it was on POOLED metrics, not the within-group frame the
verdict reads, where any purely structural predictor is identically 0 and therefore
redundant with `kshot_group_offset_baseline`. Still computed and reported; not a rival.

### 4. `WITHIN_GROUP_UNSCORABLE` became a required verdict mid-execution

The plan assumed `MIN_SCORED = 10` was the binding floor. `within_group_r2`'s own
`min_n = 30` is higher, so the within-group cell can be unscorable on a partition that
clears every other gate. The k-shot grid is now DERIVED from
`min_group_size − max(MIN_SCORED, WITHIN_GROUP_MIN_N)` and raises `KShotGridError`
rather than emitting a grid whose largest k leaves the cell unscorable.

## Post-execution defects found by `/test-epilogue`

The epilogue wrote seven tests named after live defects in my own just-shipped code.
All seven are fixed and the tests rewritten to assert the fix: the unreachable
shuffled-label null (gating nothing, and on mismatched units underneath), a NaN control
reporting as run, `power_check` conflating NOTHING_COMPARED with a measured zero,
`regime_map --screen` exiting 0 on a typo, `fold_report`'s four-vs-six key shape gap,
and `_shot_ids` keyed on PYTHONHASHSEED-salted `hash()` (reproducible within a process,
not across them). Detail is in that commit's message.

## Not actioned — user authority

- **Is F1 worth building at all?** It can currently only measure an ABSENCE: all 8
  regimes read `organism_transfer=unmeasured` and no transfer-positive arm exists on
  disk. The harness is honest about that, and the headline says so loudly, but a reader
  could reasonably ask whether a benchmark with no positive control earns its keep yet.
- **3 stranded commits** on `lane/refseq-cache-integrity` (2) and `lane/fba-structure`
  (1), pre-existing and not from this run.
- **CLAUDE.md is 5 words from its ceiling** (35,395 against ~35,400). The next entry
  must compress an existing one.
