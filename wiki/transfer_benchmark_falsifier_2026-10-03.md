# Transfer-benchmark falsifier — 2026-10-03

**Status:** `DISCRIMINATED_BOTH_DIRECTIONS` (exit 0)

**transfer_positive:** `TRANSFER_POSITIVE_UNMEASURED`

every genotype-phenotype panel on disk is a single constructed population, so a held-out-organism fold has nothing to hold out. This is a measured absence, not an omission.

## Arms

- **A1 — A1_scoring_core_control**
  - `positive`: True
  - `is_transfer_evidence`: False
- **B — B_synthetic_structure_negative**
  - `negative`: True
  - `verdict`: WITHIN_NOISE
  - `machinery_exercised`: True
- **C — C_verdict_replay_regression**
  - `machinery_exercised`: False
  - `counts_toward_pass`: False
  - `all_seeds_non_positive`: True
  - `won_global_r2_on_every_seed`: True

## Pre-registered bar

```json
{
  "pass_requires": [
    "A1_positive",
    "B_negative"
  ],
  "A1_positive": "cv_ridge_gp predictive r beats its own label-permutation null p95",
  "B_negative": "gauntlet verdict is NOT BEATS_ALL_CONTROLS on a zero-within-group-signal input",
  "C_counts_toward_pass": false,
  "C_requirement": "every replayed seed must be NON-positive DESPITE the embedding winning global r2",
  "transfer_positive": "TRANSFER_POSITIVE_UNMEASURED unless a live held-out-organism substrate is supplied",
  "frozen_at": "plan v3, before any arm was run"
}
```
