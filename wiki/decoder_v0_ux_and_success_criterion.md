# DNA Decoder v0 - UX + Success Criterion

**Status:** RELOCKED 2026-05-23. Canonical v0 surface now matches the implemented cached-strain predictor (`pipeline.py predict --strain-id ...`) rather than the earlier genome-input decoder draft.

This note anchors downstream Phase 2 decisions. The repo had drifted into two competing v0 definitions:

1. the original genome-input decoder concept
2. the implemented cached-strain cipro predictor

This relock makes the second one authoritative for v0. Real previously unseen genome-input decode remains a v0.1+ decision unless explicitly re-opened.

## Vision

A tool that scores an E. coli strain for drug resistance and emits provenance plus audit framing around that score. In the current v0 implementation, the strain must already exist in the embedding cache. Real user-supplied genome decode remains the intended direction, but is not the shipped v0 contract.

## v0 User Story

> "I have a strain that is already present in the current cipro embedding cache. I want to query: for ciprofloxacin, is this strain predicted R or S? Which gene-level regions look most influential? And what training-cohort audit framing should accompany the result?"

The user runs a single command against a cached strain and gets back a structured answer. Real previously unseen genome input is deferred to v0.1+.

## v0 Interface

CLI command via `scripts/pipeline.py predict` subcommand. Web UI and REST remain out of scope for v0.

## v0 Inputs

- `--strain-id` - cached BV-BRC strain ID already present in the embedding HDF5
- `--cache` - path to the NT embedding HDF5 cache (`nt_n147_cipro.h5`)
- `--model-path` - path to the trained classifier pickle produced by `pipeline.py train`
- `--annotations` - GFF3 path for attribution (recommended)
- `--audit-merge-json` - merge-gate JSON sidecar for explicit training-cohort verdict propagation (required for canonical reporting)
- `--output` - JSON output path; markdown sidecar auto-generated alongside it

Drug scope for the current v0 surface is cipro-only.

## v0 Outputs

Single JSON object plus a human-readable markdown sidecar.

```json
{
  "strain_id": "562.12345",
  "drug": "ciprofloxacin",
  "prediction": "R",
  "calibrated_probability": 0.87,
  "confidence_tier": "HIGH",
  "top_k_attribution": [
    {"gene_id": "gene-gyrA", "gene_symbol": "gyrA", "score": 0.42, "tier": "Tier 1"},
    {"gene_id": "gene-parC", "gene_symbol": "parC", "score": 0.18, "tier": "Tier 1"}
  ],
  "audit_verdict": {
    "suspend_gate_fired": true,
    "verdict_explanation": "Training cohort had noisy phenotype evidence; prediction is informational only."
  },
  "provenance": {
    "model": "nucleotide_transformer + XGBoost",
    "training_cohort": "stage2_n150_cipro_cohort",
    "reporting_mode": "canonical_audit_aware",
    "cv_strategy": "leave_one_accession_out",
    "cv_auroc": 0.87,
    "trained_on": "2026-05-22"
  }
}
```

Canonical v0 reporting requires `--audit-merge-json`. When that sidecar is omitted, `predict` only runs with `--allow-missing-audit`, marks the result as non-canonical internal/debug output, and leaves `audit_verdict` as `null`.

## v0 Success Criteria

A v0 ship is achieved when all of these pass on cipro:

1. **Functional:** `pipeline.py predict --strain-id X --cache Y.h5 --model-path Z.pkl --audit-merge-json A.json --output result.json` runs end-to-end on a cached held-out strain without crashing and emits both JSON and markdown outputs. Hard gate.
2. **Predictive:** primary CV AUROC on the current cipro cohort is at least `0.70`. The leakage-safe reference path is accession-grouped CV when duplicate accessions exist.
3. **Honest (audit-aware):** when `--audit-merge-json` is supplied, the `audit_verdict` field propagates from the merge gate. If `SUSPEND_CONDITION_4` fired on the training cohort, the prediction output says so explicitly. Hard gate for canonical reporting.
4. **Interpretability stance:** gene-level attribution is included and may recover known cipro loci, but it does **not** need to meet the earlier top-K=10 recovery threshold for v0 ship. Current stance is exploratory attribution with explicit scope limits, not a mechanism-first ship gate.
5. **Documentation:** `README.md` has a quick start showing the cached-strain predict command and points to the current cipro decision packet.

## Explicit non-criteria for v0

- multi-drug support
- LOMO-clade-out validation
- 4-bit Evo / DNABERT-2 / GENA-LM expansion
- per-gene NT windows
- strict-MIC training labels
- real previously unseen FASTA plus GFF decode path
- top-K=10 known-mechanism recovery on at least 50% of held-out R strains

## v0 Non-Goals

- not a clinical decision support tool
- not a classifier benchmark paper
- not multi-organism
- not a polished UI

## Current cipro scope limit

Current cipro interpretability packet says:

- predictive path is healthy on the leakage-safe model
- attribution is biologically plausible but inconsistent
- ranking-only rescue failed on the 12-strain falsifier
- real Mash/clade diagnosis remains tooling-blocked on this machine

So v0 should describe attribution as exploratory support for the prediction, not as a consistent mechanism decoder. See `reports/cipro_v0_scope_limit_decision_2026-05-23.md`.

## What the Databricks cache unlocked

The N=147 cipro NT cache made the current cached-strain v0 possible:

1. `pipeline.py train` for the cipro classifier
2. `pipeline.py predict --strain-id ...` for JSON plus markdown output
3. merge-gate `audit_verdict` propagation

Remaining work is not basic pipeline wiring. It is contract clarity plus any future v0.1 decision about real genome-input inference.

## Decisions

1. **Interface:** CLI via `scripts/pipeline.py predict`
2. **Input surface:** cached-strain inference via `--strain-id` for v0; real genome-input decode deferred to v0.1+
3. **AUROC bar:** primary CV AUROC at least `0.70` on current cipro categorical labels
4. **Interpretability stance:** exploratory and audit-aware for v0; not a hard ship gate at the earlier top-K=10 threshold
5. **Drug scope:** cipro-only for v0; cef is v0.1
6. **Output format:** JSON primary plus markdown sidecar

## What this unblocks

With these criteria written, the next work items have clearer shapes:

- **Current v0 closeout:** use the leakage-safe cached-strain predictor and current scope-limit artifact honestly
- **v0.1 planning:** decide whether to build real genome-input decode and/or cef follow-on support

Canonical release/handoff packet for the current v0 surface:

- `reports/dna_decoder_v0_release_candidate_2026-05-24.md`
- `reports/dna_decoder_v0_closeout_handoff_2026-05-24.md`

Current v0.1 genome-input follow-on packet:

- `reports/dna_decoder_v0_1_genome_input_release_candidate_2026-05-25.md`

Without this relock, code, tests, and docs continue to justify incompatible claims about what v0 is.
