# DNA Decoder v0.1 cef overnight handoff - 2026-05-26

## What completed overnight

1. Formal cef cached-strain release packet written:
   - `reports/dna_decoder_v0_1_cef_cached_release_candidate_2026-05-26.md`

2. Full real-panel cef cached-vs-genome validation completed:
   - `reports/dna_decoder_v0_1_cef_cached_vs_genome_full_validation_2026-05-26.json`
   - `reports/dna_decoder_v0_1_cef_cached_vs_genome_full_validation_2026-05-26.md`

## Current cef state

- drug: `ceftriaxone`
- cohort file: `data/processed/gate_b_cohort.parquet`
- cef pool size: `50`
- dedicated NT cache:
  - `C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_cache/embeddings/nt_gate_b_cohort_67.h5`
- trained model:
  - `data/processed/models/ceftriaxone_nucleotide_transformer.pkl`
- usable train/eval set: `49`
- primary CV AUROC: `0.895`
- duplicate-accession audit: `PASS`

## Real validation headline

- completed usable samples: `49 / 50`
- cached-strain vs genome-input prediction concordance: `49 / 49`
- label alignment on both paths: `47 / 49`
- max absolute probability delta: `0.063148`
- mean absolute probability delta: `0.002305`

Shared misses:

- `562.28389` expected `S`, both paths predicted `R`
- `562.7695` expected `R`, both paths predicted `S`

Interpretation:

- cef genome-input and cef cached-strain now agree across the full usable real panel
- remaining errors are model-side, not input-path-side

## Important caveat

Current cef reporting is still debug-mode only:

- no cef audit merge sidecar exists yet
- current cef examples and validation runs therefore use:
  - `--allow-missing-audit`
  - `--no-attribution`

So cef is now strongly viable, but not yet at the same audit-aware packaging level as cipro.

## Recommended next move

Do **cef audit-aware packet design** next.

Reason:

- biggest remaining gap is product/reporting discipline, not predictive viability
- lower risk than widening scope to another new input mode or another new drug
- current overnight evidence is already strong enough to stop spending on more cef path-consistency validation

## Best files to open first

- `reports/dna_decoder_v0_1_cef_cached_release_candidate_2026-05-26.md`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_full_validation_2026-05-26.md`
- `reports/current_cef_duplicate_accession_audit_2026-05-25.md`

## Suggested prompt for the other machine

```text
We now have a formal cef cached-strain release candidate plus full real-panel cached-vs-genome validation (49/49 path concordance, 47/49 label alignment, two shared model misses). Please propose the smallest credible design for a cef audit-aware packet that matches the cipro release discipline without reopening broader scope.
```
