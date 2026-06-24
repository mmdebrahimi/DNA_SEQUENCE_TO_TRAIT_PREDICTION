# Workhorse Convergence Recovery Checklist (2026-06-23)

## Purpose

Turn the current source-control diagnosis into an exact recovery sequence for the next session or the other machine.

## Current truth

This Databricks-bundle checkout is **not** the tracked clone.

Measured here:

- local `.git` exists
- branch has no commits
- no remote exists
- `git pull` is impossible from this checkout

See:

- `reports/workhorse_repo_convergence_status_2026-06-22.md`

## Strongest provenance clue already in-repo

`PORT_INSTRUCTIONS.md` explicitly says:

- preferred channel was git, not the zip
- intended tracked source was:
  - `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION.git`
- the Databricks-bundle copy was the fallback path for a non-git bundle checkout
- the cited historical Path B port point was:
  - `HEAD 83aaa14`

This is enough to define the recovery target, even though this checkout itself cannot reach it by `git pull`.

## Recovery goal

Reattach the authoritative workhorse artifacts to the intended tracked repo lineage, instead of letting them remain stranded in the orphan bundle checkout.

## Best currently known lineage anchor

Found locally on this machine:

- `C:\Users\b0652085\PycharmProjects\PythonProject\skills_manager\.tmp\mmdebrahimi_dna_sequence_to_trait_prediction`

What is true there:

- it is a real tracked clone
- `.git/config` points to:
  - `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION.git`
- local branch state is:
  - `main...origin/main [behind 7]`
- commit history exists

What is **not** true there:

- it does not already contain these newer workhorse packets:
  - `reports/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.md`
  - `reports/arabidopsis_af002_ft16_bounded_execution_2026-06-21.md`
  - `reports/nf001_branch_state_decision_2026-06-20.md`
  - `reports/dna_decoder_integrated_handoff_2026-06-22.md`

So it is the best local recovery anchor, but not yet the converged destination.

## Exact artifacts worth preserving

### Load-bearing branch-state packets

- `reports/nf001_branch_state_decision_2026-06-20.md`
- `reports/nf001_branch_state_decision_2026-06-20.json`
- `reports/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.md`
- `reports/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.json`
- `reports/arabidopsis_af002_ft16_bounded_execution_2026-06-21.md`
- `reports/arabidopsis_af002_ft16_bounded_execution_2026-06-21.json`

### Cross-session / portfolio state

- `reports/dna_decoder_cross_session_handoff_2026-06-22.md`
- `reports/dna_decoder_integrated_handoff_2026-06-22.md`
- `reports/workhorse_repo_convergence_status_2026-06-22.md`
- `reports/project_forward_value_map_2026-06-23.md`

### Supporting frontier docs

- `reports/eukaryotic_frontier_recommendation_2026-06-15.md`
- `reports/eukaryotic_active_frontier_registry_2026-06-15.md`
- `reports/eukaryotic_active_frontier_registry_2026-06-15.json`
- `reports/arabidopsis_af002_ft16_linkage_feasibility_card_2026-06-16.md`
- `reports/arabidopsis_af002_ft16_linkage_feasibility_card_2026-06-16.json`

### Load-bearing code / test delta from the last Path B cycle

- `dna_decode/models/foundation.py`
- `tests/test_models_foundation.py`

These two matter because the DNABERT2 safe-attention patch is part of the reproducible FT16 support-check path.

## Best recovery sequence

### 1. Find the real tracked clone

On the other machine or any known tracked checkout, verify:

1. `git remote -v`
2. `git rev-parse --is-inside-work-tree`
3. `git log --oneline -n 5`

Wanted outcome:

- real commit history exists
- intended remote lineage is visible

Current best local candidate already found:

- `C:\Users\b0652085\PycharmProjects\PythonProject\skills_manager\.tmp\mmdebrahimi_dna_sequence_to_trait_prediction`

### 2. Verify whether the tracked clone already contains the workhorse packets

Check for:

- `reports/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.md`
- `reports/arabidopsis_af002_ft16_bounded_execution_2026-06-21.md`
- `reports/nf001_branch_state_decision_2026-06-20.md`

If they already exist there:

- compare contents before doing anything else

If they do not exist there:

- this bundle copy is the only current source of those packets
- current local tracked-clone anchor is in exactly that state

### 3. Port only the authoritative files, not the whole orphan `.git`

Do:

- copy the exact report / json / code / test files listed above into the real tracked clone

Do not:

- copy this orphan `.git`
- treat this bundle checkout as the canonical git source

### 4. Reconcile before push

In the tracked clone:

1. diff the incoming files
2. verify no newer tracked versions already supersede them
3. only then commit on a branch

### 5. Only after tracked recovery, resume portfolio movement

Once the workhorse artifacts are attached to real history:

- keep `AF-001` as the active decision packet
- keep `AF-002` as support evidence only
- keep `NF-001` frozen under the old binary objective

## What not to do

- do **not** add a blind remote to this orphan checkout
- do **not** push from this checkout as if it were the real repo
- do **not** spend more compute on FT16 or ERG11-only fungal work before convergence is repaired

## Short port sentence

`reports/workhorse_convergence_recovery_checklist_2026-06-23.md` is the exact next-step recovery guide: use the real tracked clone to recover commit lineage first, then port the listed authoritative workhorse packets and the DNABERT2 patch there before doing any new science work.
