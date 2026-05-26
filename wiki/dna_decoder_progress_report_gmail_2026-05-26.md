# DNA Decoder progress report for transfer - 2026-05-26

## Executive summary

Current state:

- `v0` is closed
- `v0.1` cipro genome-input slice is real and validated
- `v0.1` cef expansion is now real and strongly validated
- next highest-value step is **cef audit-aware packet design**

What is now true:

- cached-strain cipro prediction works
- genome-input cipro prediction from `FASTA + annotations` works
- cached-strain cef prediction works
- genome-input cef prediction matches cached cef on the full usable real panel

What is not closed yet:

- cef is still debug-mode only because no cef audit sidecar exists yet
- non-AMR phenotype scoping has **not** started yet

## Phase status

Completed:

- `v0` closeout
- leakage-safe cipro retrain
- cipro release packet
- cipro genome-input `v0.1` slice
- cef cached-strain substrate/cache/model build
- cef release-candidate packet
- cef cached-vs-genome-input full-panel validation

Current phase:

- `v0.1` expansion and packaging

Recommended immediate next phase step:

- **cef audit-aware packet design**

## Canonical current artifacts

### v0 closeout

- `reports/dna_decoder_v0_release_candidate_2026-05-24.md`
- `reports/dna_decoder_v0_closeout_handoff_2026-05-24.md`
- `reports/cipro_v0_scope_limit_decision_2026-05-23.md`

### v0.1 cipro genome-input

- `reports/dna_decoder_v0_1_genome_input_release_candidate_2026-05-25.md`
- `reports/dna_decoder_v0_1_genome_input_example_2026-05-25.md`
- `reports/dna_decoder_v0_1_genome_input_validation_2026-05-25.md`

### v0.1 cef expansion

- `reports/dna_decoder_v0_1_cef_cached_release_candidate_2026-05-26.md`
- `reports/dna_decoder_v0_1_cef_overnight_handoff_2026-05-26.md`
- `reports/current_cef_duplicate_accession_audit_2026-05-25.md`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_validation_2026-05-26.md`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_full_validation_2026-05-26.md`

## Current cipro state

Model posture:

- leakage-safe
- primary CV strategy: `leave_one_accession_out`
- primary CV AUROC: about `0.870`

Product posture:

- canonical `v0` is cached-strain cipro prediction
- canonical `v0.1` slice 1 is genome-input cipro prediction

Important scope note:

- cipro interpretability is still explicitly scope-limited for product use
- predictive output is available and audit-aware
- interpretability is exploratory, not a hard product promise

## Current cef state

Substrate:

- cohort: `data/processed/gate_b_cohort.parquet`
- cef pool size: `50`
- label balance: `26R / 24S`

Cache:

- dedicated NT cache built for gate-B cohort
- cache file:
  `C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_cache/embeddings/nt_gate_b_cohort_67.h5`
- `64` complete
- `3` absent with `expected=0` annotation-gap cases

Model:

- `data/processed/models/ceftriaxone_nucleotide_transformer.pkl`
- usable train/eval set: `49`
- class balance: `25R / 24S`
- CV strategy: `loso`
- CV grouping: `strain_id`
- primary CV AUROC: `0.895`
- AUPRC: `0.838`

Duplicate-accession audit:

- verdict: `PASS`
- no duplicated non-empty `assembly_accession` values in the cef pool

## Overnight cef result

Real full-panel validation completed across the full usable cef panel.

Headline:

- requested cef-pool strains: `50`
- completed usable strains: `49`
- cached-strain vs genome-input prediction concordance: `49 / 49`
- label alignment on both paths: `47 / 49`
- max absolute probability delta: `0.063148`
- mean absolute probability delta: `0.002305`

Two shared model misses:

- `562.28389` expected `S`, both paths predicted `R`
- `562.7695` expected `R`, both paths predicted `S`

Interpretation:

- cef genome-input and cef cached-strain are now strongly aligned
- remaining issues are model-side, not path-consistency issues

## Important caveat

Current cef reporting is still **debug-mode only**.

Why:

- no cef audit merge sidecar equivalent to the cipro merge packet exists yet

So current cef commands rely on:

- `--allow-missing-audit`
- usually `--no-attribution`

That is the main reason cef is not yet packaged at the same discipline level as cipro.

## Best next step

Do **cef audit-aware packet design** next.

Reason:

- highest-value remaining gap is reporting/product discipline
- predictive viability is already proven
- path consistency is already proven
- more cef validation now would mostly repeat the same evidence

## What should not come next yet

Not yet:

- more cipro work
- more cef path-consistency validation
- non-AMR phenotype scoping

Reason:

- there is still one obvious unfinished AMR packaging task: cef audit-aware closeout

## After cef audit-aware closeout

Then the next strategic decision becomes reasonable:

1. non-AMR phenotype scoping
2. another AMR expansion
3. more input-surface work

Current recommendation:

- after cef audit-aware closeout, **non-AMR phenotype scoping becomes a strong next candidate**

## Suggested prompt for the other computer

```text
We have closed v0, landed v0.1 genome-input cipro, and now validated cef cached vs genome-input on the full usable real panel (49/49 prediction concordance, 47/49 label alignment, two shared model misses). Please analyze the repo and propose the smallest credible design for a cef audit-aware packet that matches the cipro release discipline without reopening broader scope.
```

## Best files to attach or open first

- `reports/dna_decoder_v0_1_cef_overnight_handoff_2026-05-26.md`
- `reports/dna_decoder_v0_1_cef_cached_release_candidate_2026-05-26.md`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_full_validation_2026-05-26.md`
- `reports/current_cef_duplicate_accession_audit_2026-05-25.md`
- `reports/dna_decoder_v0_1_genome_input_release_candidate_2026-05-25.md`
