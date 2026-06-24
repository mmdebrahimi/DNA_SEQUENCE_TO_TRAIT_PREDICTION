# AF-002 FT16 Bounded Execution Summary (2026-06-21)

## Purpose

Capture the first real `AF-002` execution packet after the FT16 dry-manifest went green, and decide whether `AF-002` stays the active branch or collapses back into `AF-001` stop-vs-escalate triage.

## Starting state

- dry-manifest:
  - `wiki/g2_dry_manifest_2026-06-21_ft16_remat.md`
  - verdict = `GREEN`
  - `analysis_n = 970`
  - `n_groups = 9`
  - `pseudogenome_template = pseudo{id}.fasta.gz`
- phenotype:
  - `data/arabidopsis/FT16_pheno_262.csv`
- support assets recovered from:
  - `C:\Users\b0652085\OneDrive - Bombardier\Laptop Backup\B0652085\downloaded_files\arabidopsis_pathb`
- execution venue:
  - repo `.venv`
  - local CUDA torch available
  - local cached Hugging Face snapshots available for:
    - `gena_lm`
    - `nucleotide_transformer`
    - `dnabert2`

## Structure-only baseline on the FT16 retained slice

This baseline stayed constant across all bounded runs:

- `r2 = -0.2210`
- `spearman = 0.4605`
- `within_group_r2 = 0.1393`

## GENA-LM packet

### Smoke (`window_budget = 4`)

Artifact:

- `wiki/pathb_ft16_bounded_genalm_smoke_2026-06-21.json`

Metrics:

- `r2 = -0.0951`
- `spearman = -0.3010`
- `within_group_r2 = -0.0064`

### Bounded main (`window_budget = 32`, seed `42`)

Artifact:

- `wiki/pathb_ft16_bounded_genalm_mvp970_2026-06-21.json`

Metrics:

- `r2 = -0.0650`
- `spearman = 0.0004`
- `within_group_r2 = -0.0133`

### Seed stability (`window_budget = 32`)

Artifacts:

- `wiki/pathb_ft16_bounded_genalm_seed7_mvp970_2026-06-21.json`
- `wiki/pathb_ft16_bounded_genalm_seed99_mvp970_2026-06-21.json`

| seed | embedding r2 | embedding spearman | embedding within_group_r2 |
|---|---:|---:|---:|
| `42` | `-0.0650` | `0.0004` | `-0.0133` |
| `7` | `-0.1195` | `-0.2212` | `-0.0337` |
| `99` | `-0.1287` | `-0.2646` | `-0.0426` |

Interpretation:

- `GENA-LM` is consistently less bad than the structure-only baseline on `r2`
- but it loses decisively on `spearman`
- and it fails the within-group signal test on all seeds

## Nucleotide Transformer packet

### Smoke (`window_budget = 4`)

Artifact:

- `wiki/pathb_ft16_bounded_nt_smoke_2026-06-21.json`

Metrics:

- `r2 = -0.0971`
- `spearman = 0.0939`
- `within_group_r2 = -0.0849`

### Bounded main (`window_budget = 32`, seed `42`)

Artifact:

- `wiki/pathb_ft16_bounded_nt_mvp970_2026-06-21.json`

Metrics:

- `r2 = -0.0411`
- `spearman = 0.1423`
- `within_group_r2 = 0.0085`

### Seed stability (`window_budget = 32`)

Artifacts:

- `wiki/pathb_ft16_bounded_nt_seed7_mvp970_2026-06-21.json`
- `wiki/pathb_ft16_bounded_nt_seed99_mvp970_2026-06-21.json`

| seed | embedding r2 | embedding spearman | embedding within_group_r2 |
|---|---:|---:|---:|
| `42` | `-0.0411` | `0.1423` | `0.0085` |
| `7` | `-0.1114` | `0.0325` | `-0.0273` |
| `99` | `-0.0688` | `0.0785` | `-0.0148` |

Interpretation:

- `nucleotide_transformer` is the strongest `AF-002` model tried so far on `r2`
- it keeps weak positive `spearman` on all bounded seeds
- but it still loses clearly to the structure-only baseline on `spearman`
- within-group signal is near-zero to slightly negative, not convincingly positive

## Comparison to the old FT10 bounded fallback

Reference packet:

- `wiki/pathb_bounded_genalm_mvp114_2026-06-09.json`
- `wiki/pathb_bounded_genalm_seed_stability_2026-06-09.md`

Old FT10 bounded `GENA-LM` state:

- all seeds negative on `r2`
- all seeds negative on `spearman`
- structure-only baseline directionally stronger

What changed on FT16:

- the retained slice is much larger:
  - `970` vs old bounded `114`
- both tested embedding models now beat the structure-only baseline on `r2`
- `nucleotide_transformer` keeps weakly positive rank correlation instead of collapsing below zero
- but no tested model beats the structure-only baseline on the full metric set required for a clean Path B pass

## DNABERT2 packet

### Smoke (`window_budget = 4`)

Artifact:

- `wiki/pathb_ft16_bounded_dnabert2_smoke_2026-06-21.json`

Metrics:

- `r2 = -0.1761`
- `spearman = 0.0595`
- `within_group_r2 = -0.1580`

### Bounded main (`window_budget = 32`, seed `42`)

Artifact:

- `wiki/pathb_ft16_bounded_dnabert2_mvp970_2026-06-21.json`

Metrics:

- `r2 = -0.0966`
- `spearman = 0.1231`
- `within_group_r2 = -0.0886`

Interpretation:

- `dnabert2` is now executable locally after forcing its wrapper onto the safe PyTorch attention path
- it improves over the structure-only baseline on `r2`
- it remains weaker than `nucleotide_transformer`
- it still loses decisively to the structure-only baseline on `spearman`
- and it stays negative on within-group signal

## Decision

- branch_status: `active_mixed_signal`
- decision_label: `support_check_complete`

Why:

- `AF-002` is still the best active no-money branch
- it is no longer blocked on setup or missing assets
- it now provides stronger and more informative evidence than the old FT10 fallback packet
- but it has **not** crossed into a clean `PASS` story
- and `dnabert2` does not change that conclusion

## Best next move

Highest-VOI next action is no longer another cached FT16 model run.

1. use this packet to force the `AF-001` stop-vs-escalate decision

Reason:

- the remaining open queue item is portfolio choice, not another cheap support check
- three cached supported models have now been executed on the same FT16 retained slice
- no model has converted FT16 into a clean Path B pass

Fallback if `AF-002` is deliberately deprioritized:

1. none; this packet already feeds that decision directly

## Not the right next move

- do **not** reopen `NF-001`
- do **not** spend cycles re-proving FT16 dry-manifest green
- do **not** claim a Path B `PASS` from the current FT16 packet
- do **not** collapse `AF-002` back to `prep_next`; it is past that stage
