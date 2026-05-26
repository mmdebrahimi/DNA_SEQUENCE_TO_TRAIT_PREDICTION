# DNA Decoder v0.1 genome-input cipro plan

## Summary

Move from the closed v0 cached-strain cipro predictor to a v0.1 **real genome-input cipro decode** path. Keep scope tight: one previously unseen E. coli genome, one drug (ciprofloxacin), same audit-aware output shape, no multi-drug expansion yet.

## Why this is the next phase

v0 is now closed around the cached-strain cipro predictor:

- leakage-safe model exists
- canonical audited `predict` path exists
- release packet and real example exist
- cipro interpretability is explicitly scope-limited for v0

So the biggest remaining product gap is no longer stabilization. It is the gap between:

- current input surface: `--strain-id` for already-cached strains
- intended decoder intuition: user-supplied genome input

## Goal

Define and ship the smallest credible v0.1 path that can accept a real genome input and produce the same class of output as v0:

- `prediction`
- `calibrated_probability`
- `confidence_tier`
- `audit_verdict`
- `provenance`
- optional exploratory attribution where feasible

## Scope

In scope:

- ciprofloxacin only
- one previously unseen genome input path
- CLI only
- audit-aware output retained
- reuse existing NT model and trained cipro classifier

Out of scope:

- cef or multi-drug expansion
- web UI / REST
- new interpretability research
- broad model retraining
- re-opening v0 semantics

## Key decisions

- D1: v0.1 starts with **real genome-input cipro decode**, not cef cached-strain expansion.
- D2: Keep the output contract close to v0. Input changes first; output shape stays stable.
- D3: Reuse existing genome parsing / CDS extraction / embedding machinery where possible; do not invent a second parallel embedding path.
- D4: Keep audit framing even for real genome-input decode. The training cohort remains the same noisy substrate; changing input type does not remove that obligation.
- D5: Treat attribution as optional/exploratory in v0.1 as well unless implementation cost is low and output remains honest.

## Current reusable building blocks

- `dna_decode/data/annotations.py`
  - `parse_gff3`
  - `parse_genbank`
  - `extract_cds_sequences`
- `scripts/populate_cache.py`
  - `resolve_strain_assets`
  - model build + cache populate pattern
- `dna_decode/models/cache.py`
  - `EmbeddingCache`
- `scripts/pipeline.py`
  - current audited predict/report path

## Main product question for v0.1

What is the smallest acceptable input contract for real genome decode?

Candidate options:

1. FASTA + GFF3 required
   - lowest engineering risk
   - fastest to ship
   - less friendly for users

2. FASTA + GenBank required
   - still structured
   - can reuse `parse_genbank`
   - less common for some workflows

3. FASTA only
   - best product feel
   - highest engineering complexity because annotation/CDS discovery must be solved in-flow

Recommended default:

- start with **FASTA + annotations required**
- expand toward FASTA-only only after the first v0.1 decode path works

## Proposed v0.1 slice

### Slice 1: input-contract upgrade

Add a new CLI path that accepts real genome assets rather than `--strain-id`.

Minimal acceptable interface:

- `--genome-fasta`
- `--annotations`
- `--model-path`
- `--audit-merge-json`
- `--output`

Behavior:

- parse annotations
- extract CDS sequences
- embed CDS sequences with the configured foundation model
- aggregate features the same way the cached-strain path does
- score with the same trained cipro classifier
- emit the same JSON + markdown surface

### Slice 2: temporary in-memory vs on-disk embedding decision

Two viable implementation strategies:

1. in-memory embed-and-predict
   - simplest product semantics
   - no cache write required
   - likely best first slice

2. temporary cache-backed flow
   - reuses more of the existing cache machinery
   - may be easier if current embedding path assumes cache semantics
   - noisier operationally

Recommended default:

- prefer **in-memory embed-and-predict** first
- fall back to temporary cache-backed flow only if the current embedding APIs make in-memory use awkward

## First success gate

v0.1 Slice 1 is successful when all of these are true:

1. One previously unseen cipro genome can be scored from local genome assets.
2. Output still contains audit-aware provenance and the same core prediction fields.
3. The new path does not break the current cached-strain v0 path.
4. A focused end-to-end test covers the new CLI contract with synthetic or tiny real fixtures.

## Risks

- annotation dependency may make the first real-input UX feel clunky
- embedding runtime for new inputs may be slow on this machine
- current model was trained on cached cohort flows; real-input path may expose shape assumptions not visible in v0
- attribution may be expensive or awkward for first-pass genome-input inference

## Recommended implementation order

1. Decide the exact v0.1 input contract.
2. Build the narrowest genome-input predict path around that contract.
3. Preserve the existing v0 output schema and audit framing.
4. Add focused tests.
5. Only then decide whether to expand to FASTA-only or cef follow-on support.

## Not the next step

Do not do these before Slice 1 lands:

- cef expansion
- Mash/clade research
- attribution-method redesign
- UI work

## Handoff references

- `reports/dna_decoder_v0_closeout_handoff_2026-05-24.md`
- `reports/dna_decoder_v0_release_candidate_2026-05-24.md`
- `wiki/decoder_v0_ux_and_success_criterion.md`
- `scripts/pipeline.py`
- `scripts/populate_cache.py`
- `dna_decode/data/annotations.py`
