# DNA Decoder Cross-Session Handoff (2026-06-22)

## Purpose

Carry one current source of truth into the next session.

This file answers:

- what is decided
- what is still active
- what is worth doing next
- what is not worth doing next

## Source-control reality

Current measured state of this Databricks-bundle checkout:

- local `.git` exists
- current branch has no commits
- no remote is configured

Meaning:

- this workhorse side is still not converged with the decoder-side tracked repo
- do not assume `origin/main` already contains the workhorse packets from this session
- source-control recovery is now a higher-value operational task than more default model execution

## Current portfolio state

### `NF-001`

Status:

- `freeze_negative_result`

Authoritative packet:

- `reports/nf001_branch_state_decision_2026-06-20.md`
- `reports/nf001_branch_state_decision_2026-06-20.json`

Meaning:

- the current binary objective is finished
- naive cross-lineage `ERG11`-only deterministic calling failed
- do not reopen by default

Reopen only if all are true:

1. new objective is lineage-aware interpretation, not binary calling
2. a downstream consumer is named
3. one precise beyond-`ERG11` question is named

### `AF-002`

Status:

- `support_check_complete`
- branch state = `active_mixed_signal`

Authoritative packet:

- `reports/arabidopsis_af002_ft16_bounded_execution_2026-06-21.md`
- `reports/arabidopsis_af002_ft16_bounded_execution_2026-06-21.json`

What is now true:

- FT16 dry-manifest is green at `analysis_n = 970`
- three cached supported models ran on the same retained slice:
  - `GENA-LM`
  - `nucleotide_transformer`
  - `dnabert2`
- all three improve over structure on `r2`
- none beats structure on `spearman`
- none beats structure on `within_group_r2`

Meaning:

- `AF-002` added real decision value
- it did **not** rescue Path B into a clean no-money pass
- do not keep spending default cycles on more like-for-like FT16 support runs

### `AF-001`

Status:

- current live portfolio decision point

Authoritative packet:

- `reports/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.md`
- `reports/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.json`

Forced readout:

- `decision_label = stop_default_no_money_pathb_execution`
- `branch_state = mixed_signal_not_escalation_ready`

Meaning:

- default no-money Path B execution should stop
- Path B is only worth reopening as an explicit compute-backed escalation

### `NF-002`

Status:

- `park`

Meaning:

- still too broad
- only worth touching if it changes a real next-cycle decision

## Current ranking

1. workhorse-side source-control convergence hygiene
2. `AF-001` - decision / closure / escalation policy
3. `AF-002` - support evidence for that decision
4. `NF-001` - archival / reopen-trigger only
5. `NF-002` - parked

## Best next value-add lanes

### Lane 1 - converge the workhorse side into real source control

Highest-value operational move now.

Why:

- decoder-side state is already banked on `origin/main`
- workhorse-side state is still local-only in practice
- that is the biggest remaining drift risk

Read first:

- `reports/workhorse_repo_convergence_status_2026-06-22.md`

### Lane 2 - close the current cycle cleanly

Highest-value non-research move:

- treat `AF-001` as the current portfolio decision
- stop default no-money Path B execution
- make sure future sessions do not drift back into "one more cheap model run"

This is mostly governance / clarity value, not more model execution.

### Lane 3 - define a real Path B escalation contract

Only do this if there is appetite to continue Arabidopsis embedding work.

Value-add artifact:

- a narrow escalation memo or technical plan that answers:
  - what bigger GPU / memory tier is required
  - what exact rerun would be performed
  - what success bar must be met
  - what would count as failure even after escalation

Important:

- do not reopen Path B without this contract
- do not use `r2` alone as the success criterion

### Lane 4 - identify a genuinely new high-VOI frontier

The current branches have already paid out most of their local no-money value.

So the next big value likely comes from:

- a new dataset / phenotype family with a much cleaner readout surface than current Path B
- or a tightly scoped new branch whose smallest credible experiment is green and decision-bearing

Use this standard:

1. exact source pinned
2. smallest experiment checkable
3. time-to-first-packet bounded
4. real portfolio decision unblocked

### Lane 5 - only reopen fungal work under a new contract

`NF-001` can still add value, but not under the old binary objective.

Only valuable fungal follow-on:

- a new lineage-aware interpretation branch

That branch should start with:

- named consumer
- explicit non-binary output contract
- one precise beyond-`ERG11` question

Without that, fungal work will drift back into already-falsified framing.

## What not to do next

- do **not** reopen `NF-001` under the old binary objective
- do **not** run more default FT16 bounded model checks just because the stack is working
- do **not** treat `AF-002` mixed signal as a hidden pass
- do **not** keep Path B alive without an explicit escalation gate
- do **not** start a broad new scouting branch unless it clears the same decision-bearing bar

## Exact files to read first next session

Read in this order:

1. `reports/dna_decoder_cross_session_handoff_2026-06-22.md`
2. `reports/workhorse_repo_convergence_status_2026-06-22.md`
3. `reports/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.md`
4. `reports/arabidopsis_af002_ft16_bounded_execution_2026-06-21.md`
5. `reports/eukaryotic_frontier_recommendation_2026-06-15.md`
6. `reports/nf001_branch_state_decision_2026-06-20.md`

## Short port sentence

`reports/dna_decoder_cross_session_handoff_2026-06-22.md` is the current source of truth: workhorse-side source control is still unconverged despite local `.git`, `NF-001` is frozen as a negative result, `AF-002` completed a mixed-signal FT16 support check across three cached models, and the live portfolio decision is `AF-001` - stop default no-money Path B execution unless a stricter compute-backed escalation is explicitly chartered.
