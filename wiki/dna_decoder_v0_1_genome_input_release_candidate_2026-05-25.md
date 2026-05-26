# DNA Decoder v0.1 genome-input release candidate - 2026-05-25

## Status

- Release state: `candidate`
- Scope: genome-input ciprofloxacin predictor
- Canonical command surface: `python -m scripts.pipeline predict`
- Input contract: `FASTA + annotations`
- Predictive posture: leakage-safe inherited cipro model
- Interpretability posture: exploratory support with explicit scope limit

## What v0.1 adds

v0.1 slice 1 keeps the current cipro classifier and output contract, but adds a second `predict` input mode:

- v0: `--strain-id` for cached strains
- v0.1: `--genome-fasta` + `--annotations` for one local genome

The genome-input path:

- parses local annotations
- extracts CDS sequences
- embeds them in memory with the same foundation model used in training
- mean-pools gene embeddings the same way as the cached path
- scores with the same trained leakage-safe cipro classifier
- emits the same JSON + markdown shape as v0

## What v0.1 still is not

- not FASTA-only decode
- not multi-drug
- not a mechanism-first interpretability product
- not a clinical decision support tool

## Gold-path command

```bash
uv run python -m scripts.pipeline predict \
  --model-path data/processed/models/ciprofloxacin_nucleotide_transformer.pkl \
  --genome-fasta path/to/sample.fna \
  --annotations path/to/sample.gff3 \
  --sample-id external_sample_001 \
  --audit-merge-json wiki/cipro_mechanism_phenotype_merge_2026-05-17.json \
  --output reports/genome_input_predict_example.json
```

Supported annotation suffixes:

- `.gff`
- `.gff3`
- `.gbk`
- `.gbff`
- `.genbank`

## Required inputs

- trained cipro model pickle
- local genome FASTA
- local annotation file
- audit merge JSON sidecar for canonical reporting

Optional:

- `--sample-id`
  Overrides the default output label derived from the FASTA filename stem.

Internal/debug only:

- `--allow-missing-audit`
  Use only for non-canonical local runs. Output will carry `audit_verdict = null` and `provenance.reporting_mode = non_canonical_missing_audit`.

## Expected outputs

- `result.json`
- `result.md`

Required output fields remain:

- `prediction`
- `calibrated_probability`
- `confidence_tier`
- `top_k_attribution`
- `audit_verdict`
- `provenance`

Expected provenance highlights:

- `training_cohort`
- `reporting_mode`
- `cv_strategy`
- `cv_auroc`
- `trained_on`
- `input_mode`

## Real reference run

Concrete example artifact from the current leakage-safe model:

- `reports/dna_decoder_v0_1_genome_input_example_2026-05-25.json`
- `reports/dna_decoder_v0_1_genome_input_example_2026-05-25.md`

Real run headline:

- sample: `GCA_002201835.1`
- drug: `ciprofloxacin`
- prediction: `R`
- calibrated probability: `0.8619441390037537`
- confidence tier: `MEDIUM`
- reporting mode: `canonical_audit_aware`
- input mode: `genome_input`
- CV strategy: `leave_one_accession_out`
- primary CV AUROC: `0.8697607851155182`

Cross-path validation artifact:

- `reports/dna_decoder_v0_1_genome_input_validation_2026-05-25.json`
- `reports/dna_decoder_v0_1_genome_input_validation_2026-05-25.md`

Validation headline:

- 4 real cohort genomes
- 2 resistant + 2 susceptible
- cached-strain path vs genome-input path on the same trained model
- prediction concordance: `4 / 4`
- max absolute probability delta: `0.011599`

Audit behavior on this run:

- the external genome is not present in the training-cohort merge JSON
- canonical output therefore falls back to **cohort-level** audit framing
- `SUSPEND_CONDITION_4` still propagates honestly
- no fake per-strain audit row is invented

## Current cipro scope limit

Use this wording consistently:

> Predictive output is available and audit-aware. Mechanistic attribution is exploratory and sometimes recovers known QRDR loci, but it does not yet meet the locked v0 consistency target for held-out resistant strains.

Reference artifact:

- `reports/cipro_v0_scope_limit_decision_2026-05-23.md`

## Canonical release assertions

- genome-input `predict` works on a real local genome package
- genome-input and cached-strain paths are concordant across a 4-sample mixed panel
- cached-strain v0 path remains intact
- canonical reporting still requires `--audit-merge-json`
- audit-free prediction is explicitly debug-only
- focused genome-input + cached-strain predict tests are green
- the real runtime path required batched live embedding to avoid GPU OOM

## Recommended next decision

The next highest-value choice is no longer "does genome-input decode work?" It does.

The next question is:

1. FASTA-only decode next
2. cef follow-on next
3. more real genome-input validation examples first

Recommended next step:

- do **more real genome-input validation examples first**

Reason:

- it strengthens confidence in the new input mode before broadening scope again
- cheaper than FASTA-only contract expansion
- keeps cef expansion from racing ahead of the new decoder path

## Handoff references

- `plans/DNA_Decoder_v0_1_Genome_Input_Cipro_Plan.md`
- `reports/dna_decoder_v0_closeout_handoff_2026-05-24.md`
- `reports/dna_decoder_v0_release_candidate_2026-05-24.md`
- `reports/dna_decoder_v0_1_genome_input_example_2026-05-25.md`
- `wiki/decoder_v0_ux_and_success_criterion.md`
