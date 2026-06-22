# DNA Decoder Cross-Session Handoff (2026-06-22)

## Purpose

Carry one current source of truth into the next session.

This file answers:

- what is decided
- what is still active
- what is worth doing next
- what is not worth doing next

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

1. `AF-001` — decision / closure / escalation policy
2. `AF-002` — support evidence for that decision
3. `NF-001` — archival / reopen-trigger only
4. `NF-002` — parked

## Best next value-add lanes

### Lane 1 — close the current cycle cleanly

Highest-value non-research move:

- treat `AF-001` as the current portfolio decision
- stop default no-money Path B execution
- make sure future sessions do not drift back into “one more cheap model run”

This is mostly governance / clarity value, not more model execution.

### Lane 2 — define a real Path B escalation contract

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

### Lane 3 — identify a genuinely new high-VOI frontier

The current branches have already paid out most of their local no-money value.

So the next big value likely comes from:

- a new dataset / phenotype family with a much cleaner readout surface than current Path B
- or a tightly scoped new branch whose smallest credible experiment is green and decision-bearing

Use this standard:

1. exact source pinned
2. smallest experiment checkable
3. time-to-first-packet bounded
4. real portfolio decision unblocked

### Lane 4 — only reopen fungal work under a new contract

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
2. `reports/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.md`
3. `reports/arabidopsis_af002_ft16_bounded_execution_2026-06-21.md`
4. `reports/eukaryotic_frontier_recommendation_2026-06-15.md`
5. `reports/nf001_branch_state_decision_2026-06-20.md`

## Short port sentence

`reports/dna_decoder_cross_session_handoff_2026-06-22.md` is the current source of truth: `NF-001` is frozen as a negative result, `AF-002` completed a mixed-signal FT16 support check across three cached models, and the live portfolio decision is `AF-001` — stop default no-money Path B execution unless a stricter compute-backed escalation is explicitly chartered.

---

## Laptop-side reconciliation (2026-06-22, Claude on the GTX 860M laptop)

Ingested this handoff into the laptop repo (the AF/NF packets themselves live in the Databricks bundle and are NOT synced here — cross-machine drift; this laptop is also ahead of origin with unrelated HIV + label-survey work). Two reconciliations against the laptop repo's standing position:

1. **`AF-001` "compute-backed escalation" door vs the repo's 2026-06-12 closure — RECONCILE TOWARD CLOSED.** CLAUDE.md (the laptop's authoritative position, 2026-06-12) already closed Path B as a de-confounded NEGATIVE with an explicit principled rule: "do NOT scale embeddings on a bigger/paid GPU — a negative de-confounded metric is a signal-vs-structure problem, not a window-budget one." `AF-002` (3 more cached models; within_group_r2 still does not beat structure) is the 4th de-confounded embedding failure — it STRENGTHENS that closure. So Lane-2's escalation door should be treated as effectively closed: a future escalation contract is legitimate ONLY if it addresses the signal-vs-structure problem itself (a fundamentally different label or de-confounding design), NOT merely a bigger GPU/compute tier. Bigger compute on the same readout surface is the exact move the repo forbids. (User ratifies if they disagree.)

2. **Three-way convergence (the real signal).** Independently of this handoff, the laptop's 2026-06-22 label-acquisition survey (`research_outputs/noncircular-label-sources-2026-06-22.md`) reached the SAME conclusion from the data-acquisition side: no new free, downloadable, non-circular, genome-linked label source exists for the executor to gather; the binding constraint is labels/branch-selection, not models or local runs. So repo closure + workhorse AF-001 stop + laptop label-survey ceiling all agree: stop low-VOI local execution; the next real value is a USER branch/label decision, not another packet.

**Forward fork (user decision — same on both machines):**
- Free + executor-eligible: prospective-lock validation of the frozen deterministic decoder (no new label).
- External walls (money / academic-MTA / human curation / contact): ARESdb (verified proprietary), a hand-curated post-2023 TB independent gold set, or a contact-gated reference collection.
- Do NOT: rerun FT16/Path B without a signal-vs-structure-addressing contract; reopen NF-001 under the old binary objective; start a broad new branch that doesn't clear the pinned-source/green-smallest-experiment/real-decision bar.

**Sync action recommended (not auto-done):** a bidirectional sync is needed to converge — this laptop's HIV + label-survey commits aren't on origin, and the workhorse's AF/NF reports aren't pulled here. Neither machine sees the other's latest until then.
