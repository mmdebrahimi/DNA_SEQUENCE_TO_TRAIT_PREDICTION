# AF-001 Stop-vs-Escalate Decision (2026-06-21)

## Purpose

Force the next honest `AF-001` decision after the `AF-002` FT16 support check completed.

## Inputs used

- `wiki/g2_dry_manifest_2026-06-09_ft10_full1003.json`
- `wiki/pathb_bounded_genalm_mvp114_2026-06-09.json`
- `wiki/pathb_bounded_genalm_seed_stability_2026-06-09.json`
- `reports/arabidopsis_af002_ft16_bounded_execution_2026-06-21.md`
- `wiki/pathb_ft16_bounded_dnabert2_mvp970_2026-06-21.json`

## What is already true

### FT10 (`AF-001`)

- full dry-manifest was green at `analysis_n = 1003`
- bounded `GENA-LM` packet on the old `114`-accession slice stayed negative across all tested seeds
- structure baseline on that FT10 bounded packet was directionally stronger than the embedding packet on every reported metric

### FT16 (`AF-002`)

- dry-manifest was green at `analysis_n = 970`
- three cached supported models were run on the same frozen retained slice:
  - `GENA-LM`
  - `nucleotide_transformer`
  - `dnabert2`
- all three improve over the FT16 structure baseline on `r2`
- none beats the FT16 structure baseline on `spearman`
- none beats the FT16 structure baseline on `within_group_r2`
- strongest tested FT16 model is still `nucleotide_transformer`, but its best bounded packet remains mixed rather than clean

## Forced decision

- decision_label: `stop_default_no_money_pathb_execution`
- branch_state: `mixed_signal_not_escalation_ready`

## Why

`AF-002` did the job it was supposed to do: it checked whether a second Arabidopsis phenotype could turn Path B into a cleaner story than the old FT10 bounded packet.

It did not.

What changed:

- Path B is no longer only an `AF-001` single-slice failure
- the FT16 check is materially better than the old FT10 bounded packet on `r2`
- but the signal still does not survive the rank / within-group tests that matter for a clean claim

What did **not** change:

- there is still no clean no-money Path B `PASS`
- current local evidence still does not justify treating Arabidopsis embedding as the default forward branch

## Escalate only if

Escalation is still allowed, but only under a stricter gate:

1. a larger-GPU or rewritten inference path is explicitly available
2. the next packet is defined in advance as an escalation packet, not another cheap support check
3. success is judged against the same full metric surface, not `r2` alone

## Recommended next move

1. treat `AF-001` as decided for the current no-money local objective:
   - stop default Path B execution
2. keep `AF-002` as support evidence showing Path B is mixed rather than purely dead
3. move the active portfolio away from default Arabidopsis embedding execution unless a deliberate compute-backed escalation is chosen

## Short readout

`AF-002` improved the evidence quality but did not rescue Path B. Best judgment is to stop default no-money Arabidopsis embedding execution and treat further Path B work as an explicit escalation-only choice.
