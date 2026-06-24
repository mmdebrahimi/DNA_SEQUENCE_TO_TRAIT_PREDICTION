# NF-001 Branch State Decision (2026-06-20)

## Decision

- decision: `freeze_negative_result`
- confidence: `HIGH`
- summary: Freeze NF-001 as a successful negative result against naive cross-lineage ERG11-only deterministic calling.

## Current packet

- scored_rows: `13`
- cross_lineage_status: `FAIL`
- scoped_safe_claim: `none`
- overall_verdict: `FAIL`
- accuracy: `0.6153846153846154`
- sensitivity: `0.8571428571428571`
- specificity: `0.3333333333333333`
- tp/tn/fp/fn: `6/2/4/1`

## Why

- All tested scope slices now fail under the current observed packet.
- No safe positive claim survives under the current ERG11-only catalog.
- Lineage III contains replicated determinant-positive susceptible rows through MIC8 and MIC16.
- Lineage IV is fully sampled in the current packet and remains uncatalogued/non-discriminative.

## Per-lineage status

### Lineage I

- status: `positive_control_replicated`
- n_scored: `3`
- summary: Observed lineage I remains a clean resistant positive-control slice.

### Lineage III

- status: `conflicted_false_positive`
- n_scored: `8`
- summary: Observed lineage III contains determinant-positive susceptible rows, so current ERG11-only calling is conflicted inside lineage III.
- false_positive_isolate_ids: `2566, 3345, 3561_77, 6277`

### Lineage IV

- status: `uncatalogued_non_discriminative`
- n_scored: `2`
- summary: Observed lineage IV is fully sampled in the current packet and remains non-discriminative under the current catalog.

## Reopen conditions

- A new follow-on branch defines non-binary lineage-aware interpretation as the objective rather than binary resistance calling.
- A downstream consumer is named for interpretation-only outputs.
- One precise beyond-ERG11 mechanism question is identified that current outputs cannot represent.

## Do not do next

- Do not run more like-for-like ERG11-only observation inside this packet by default.
- Do not start a broad mammoth/de-extinction or generic comparative-genomics research sweep.
- Do not add more determinant loci onto the current binary scoring surface without redefining outputs first.
