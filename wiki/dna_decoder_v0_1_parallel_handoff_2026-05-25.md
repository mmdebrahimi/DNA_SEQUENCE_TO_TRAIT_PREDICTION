# DNA Decoder v0.1 parallel handoff - 2026-05-25

## Executive state

- v0 is closed.
- v0.1 slice 1 is landed and real.
- Current repo now supports:
  - cached-strain cipro predict (`--strain-id`)
  - genome-input cipro predict (`--genome-fasta` + `--annotations`)

## What landed in this session

### 1. New genome-input predict path

Implemented in:

- `scripts/pipeline.py`
- `dna_decode/data/annotations.py`

New input mode:

- `--genome-fasta`
- `--annotations`
- optional `--sample-id`

Behavior:

- loads annotations from GFF3 or GenBank
- extracts CDS sequences
- embeds CDS sequences in memory with the trained foundation-model family
- mean-pools the same way as the cached-strain path
- scores with the existing leakage-safe cipro classifier
- emits the same JSON + markdown schema as v0

### 2. Runtime hardening

Real first smoke exposed a root-cause bug:

- genome-input path embedded all CDS in one giant batch
- NT model OOMed on GPU

Fix:

- live embedding now batches in small chunks, matching the existing cache-populate behavior

### 3. Audit behavior cleanup

Canonical audited prediction now falls back to **cohort-level audit framing** when:

- `--audit-merge-json` is supplied
- but the specific sample is missing from `per_strain`

This now behaves consistently for both:

- cached-strain path
- genome-input path

## Verification status

### Focused tests

- `uv run pytest tests/test_data_annotations.py -q` -> `26 passed`
- `uv run pytest tests/test_pipeline_cli.py tests/test_pipeline_predict_v0.py tests/test_pipeline_predict_genome_input_e2e.py tests/test_pipeline_predict_e2e.py -q` -> `43 passed`

### Real genome-input smoke

Artifacts:

- `reports/dna_decoder_v0_1_genome_input_example_2026-05-25.json`
- `reports/dna_decoder_v0_1_genome_input_example_2026-05-25.md`

Headline:

- sample: `GCA_002201835.1`
- drug: `ciprofloxacin`
- prediction: `R`
- probability: `0.8619441390037537`
- confidence: `MEDIUM`
- reporting mode: `canonical_audit_aware`
- input mode: `genome_input`

### Cross-path consistency validation

Artifacts:

- `reports/dna_decoder_v0_1_genome_input_validation_2026-05-25.json`
- `reports/dna_decoder_v0_1_genome_input_validation_2026-05-25.md`

Panel:

- 4 real cohort genomes
- 2 resistant, 2 susceptible
- each scored both ways:
  - cached-strain path
  - genome-input path

Headline:

- prediction concordance: `4 / 4`
- max absolute probability delta: `0.011599`

Interpretation:

- genome-input path is no longer just synthetically correct
- it is now functionally consistent with the existing cached-strain path on a mixed real panel

## Current product truth

We can now honestly say:

> We can predict ciprofloxacin resistance for cached E. coli strains and for local genome inputs with audit-aware output.

We still should **not** say:

> We can identify whether a strain is resistant to an antibiotic

That broader claim is still too strong because:

- current scope is cipro only
- current genome-input contract still requires annotations
- interpretability remains exploratory
- multi-drug support is not closed

## Current gold-path artifacts

- v0 contract: `wiki/decoder_v0_ux_and_success_criterion.md`
- v0 release packet: `reports/dna_decoder_v0_release_candidate_2026-05-24.md`
- v0.1 release packet: `reports/dna_decoder_v0_1_genome_input_release_candidate_2026-05-25.md`
- v0.1 real example: `reports/dna_decoder_v0_1_genome_input_example_2026-05-25.md`
- v0.1 validation packet: `reports/dna_decoder_v0_1_genome_input_validation_2026-05-25.md`
- cipro scope limit: `reports/cipro_v0_scope_limit_decision_2026-05-23.md`

## Recommended next path

### Recommendation

Do **not** spend the next cycle on:

- more cipro validation of the same kind
- FASTA-only UX
- more cipro interpretability work

Recommended next path:

> **cef expansion next**

### Why

1. Genome-input cipro is now proven enough for slice 1:
   - one real smoke
   - one 4-sample concordance check

2. FASTA-only is a bigger engineering jump with less immediate product value:
   - requires annotation / CDS discovery in-flow
   - higher runtime/tooling complexity
   - does not broaden drug usefulness

3. Cef adds more product value than polishing cipro input UX further:
   - broadens drug scope
   - keeps current CLI/output contract largely stable
   - lets us isolate the **drug axis** next instead of changing both drug and annotation contract together

### Best implementation shape for the next phase

If you take cef next, prefer:

- **cached-strain cef first**

not:

- genome-input cef first

Reason:

- only change one axis at a time
- input-mode expansion is already validated on cipro
- next clean question should be: does the same product surface hold for a second drug?

## What the other Claude session can do in parallel

Highest-value parallel work:

1. Review current cef substrate readiness under the now-relied-on v0.1 product framing.
   - Is cached-strain cef viable enough for the next slice?
   - What exact cohort / model / audit gaps remain?

2. Produce a short plan for:
   - `cef cached-strain v0.1 slice`
   - not genome-input cef
   - not FASTA-only cipro

3. Check whether any existing cef artifacts in repo should be promoted into the same release-packet style used for cipro.

## Suggested ask for the other session

Use this prompt:

> We have closed v0 and landed v0.1 slice 1 for genome-input ciprofloxacin decode. Please analyze the repo and recommend the smallest credible next slice for **cef cached-strain expansion**, including what already exists, what is missing, and the minimum success gate. Use the current v0/v0.1 release packets and validation artifacts as the source of truth.

## Bottom line

- v0.1 cipro genome-input slice is real and validated.
- The highest-value next move is no longer more cipro proof.
- The next clean product step is **cef expansion**, ideally **cached-strain cef first**.
