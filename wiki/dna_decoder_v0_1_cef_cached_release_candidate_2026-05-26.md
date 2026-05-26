# DNA Decoder v0.1 cef cached-strain release candidate - 2026-05-26

## Status

- Release state: `candidate`
- Scope: cached-strain ceftriaxone predictor
- Canonical command surface: `python -m scripts.pipeline predict`
- Predictive posture: viable on real gate-B substrate
- Reporting posture: debug-only until a cef audit sidecar exists
- Interpretability posture: exploratory and out of scope for this packet

## What this slice adds

This v0.1 follow-on adds a second drug on the existing cached-strain surface:

- v0: cached-strain ciprofloxacin
- v0.1 slice 1: genome-input ciprofloxacin
- v0.1 follow-on: cached-strain ceftriaxone

The cef path keeps the simpler cached-strain contract:

- read one strain already present in the dedicated cef embedding cache
- score with a trained nucleotide-transformer cef classifier
- emit the same JSON + markdown output shape as the cipro path

## What this slice still is not

- not genome-input cef
- not audit-aware canonical cef reporting yet
- not a mechanism-first interpretability product
- not a clinical decision support tool

## Gold-path command

```bash
uv run python -m scripts.pipeline predict \
  --model-path data/processed/models/ceftriaxone_nucleotide_transformer.pkl \
  --strain-id 562.12960 \
  --cache "C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_cache/embeddings/nt_gate_b_cohort_67.h5" \
  --allow-missing-audit \
  --no-attribution \
  --output reports/dna_decoder_v0_1_cef_cached_example_R_2026-05-25.json
```

## Required inputs

- trained cef model pickle
- cached strain ID already present in the cef HDF5 cache
- cef cache path

Current recommended local cache:

- `C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_cache/embeddings/nt_gate_b_cohort_67.h5`

Current required flags for honest local use:

- `--allow-missing-audit`
- `--no-attribution`

Why:

- no cef audit merge sidecar equivalent to the cipro merge packet is wired yet
- this packet is about predictive viability on the second drug, not audit-closeout

## Expected outputs

- `result.json`
- `result.md`

Required output fields remain:

- `prediction`
- `calibrated_probability`
- `confidence_tier`
- `audit_verdict`
- `provenance`

Current expected provenance highlights:

- `training_cohort`
- `reporting_mode = non_canonical_missing_audit`
- `cv_strategy`
- `cv_auroc`
- `trained_on`

## Real substrate and model

Working cef substrate:

- cohort file: `data/processed/gate_b_cohort.parquet`
- cef pool size: `50`
- cef label balance: `26R / 24S`

Dedicated cef cache:

- `C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_cache/embeddings/nt_gate_b_cohort_67.h5`

Cache probe result:

- `64` complete strains
- `3` absent strains with `expected=0`
- those `3` are annotation-gap cases, not partial cache corruption

Trained model:

- `data/processed/models/ceftriaxone_nucleotide_transformer.pkl`
- usable train/eval set: `49`
- class balance: `25R / 24S`
- CV strategy: `loso`
- CV grouping: `strain_id`
- primary CV AUROC: `0.895`
- AUPRC: `0.838`

Duplicate-accession audit:

- `reports/current_cef_duplicate_accession_audit_2026-05-25.md`
- verdict: `PASS`

## Real reference runs

Debug-mode example artifacts:

- resistant example:
  - `reports/dna_decoder_v0_1_cef_cached_example_R_2026-05-25.json`
  - `reports/dna_decoder_v0_1_cef_cached_example_R_2026-05-25.md`
- susceptible example:
  - `reports/dna_decoder_v0_1_cef_cached_example_S_2026-05-25.json`
  - `reports/dna_decoder_v0_1_cef_cached_example_S_2026-05-25.md`

Headline:

- `562.12960` -> `R`, `p(R)=0.753`
- `562.7572` -> `S`, `p(R)=0.204`

## Cross-path validation

Validation artifact:

- `reports/dna_decoder_v0_1_cef_cached_vs_genome_validation_2026-05-26.json`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_validation_2026-05-26.md`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_full_validation_2026-05-26.json`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_full_validation_2026-05-26.md`

Current mixed-panel result:

- `8` real samples (`4R / 4S`)
- cached-strain path vs genome-input path on the same trained cef model
- prediction concordance: `8 / 8`
- label alignment on both paths: `7 / 8`
- max absolute probability delta: `0.063148`
- mean absolute probability delta: `0.008866`

Important interpretation:

- cef cached and cef genome-input agree on every tested sample
- one borderline susceptible sample (`562.28389`) is miscalled by both paths
- that points to model/label behavior on that sample, not an input-path inconsistency

Overnight full-panel result:

- requested cef-pool strains: `50`
- completed usable strains: `49`
- cached-strain path vs genome-input path prediction concordance: `49 / 49`
- label alignment on both paths: `47 / 49`
- max absolute probability delta: `0.063148`
- mean absolute probability delta: `0.002305`
- only two shared model misses:
  - `562.28389` expected `S`, both paths predicted `R`
  - `562.7695` expected `R`, both paths predicted `S`

Interpretation:

- the cef genome-input path is now strongly validated against the cached cef path on essentially the full usable panel
- remaining misses are shared model errors, not path-consistency errors

## Current release assertion

Use this wording consistently:

> Cef cached-strain prediction is now product-viable on the current gate-B substrate and matches the genome-input path on a real mixed validation panel, but it is still a debug-mode surface until a cef audit sidecar is defined.

## Recommended next decision

Best next move is no longer "can cef work at all?" It can.

Next question:

1. add a cef audit-aware packet
2. add cef genome-input release packaging
3. broaden cef validation further

Recommended next step:

- do **cef audit-aware packet design first**

Reason:

- it closes the main gap between cef viability and cef release discipline
- lower risk than changing drug and input mode at the same time
- keeps the current strong predictive substrate intact

## Handoff references

- `reports/dna_decoder_v0_1_cef_cached_handoff_2026-05-25.md`
- `reports/current_cef_duplicate_accession_audit_2026-05-25.md`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_validation_2026-05-26.md`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_full_validation_2026-05-26.md`
- `reports/dna_decoder_v0_1_genome_input_release_candidate_2026-05-25.md`
- `reports/dna_decoder_v0_1_parallel_handoff_2026-05-25.md`
