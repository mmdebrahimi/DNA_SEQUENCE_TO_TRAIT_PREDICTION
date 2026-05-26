# DNA Decoder v0.1 cef audit-aware release packet — TEMPLATE (pending audit results)

> Fill-in-the-blank template for the audit-aware cef release packet that will replace
> `wiki/dna_decoder_v0_1_cef_cached_release_candidate_2026-05-26.md` once Codex's
> Artifacts 2-4 land (cef MIC audit + mechanism × MIC merge + 4 canonical predicts).
>
> This template pre-commits the verdict-conditional release wording so the
> post-Codex-push update is mechanical: pick the branch matching the actual
> `gate_verdict` from `wiki/cef_mechanism_phenotype_merge_2026-05-26.json`, fill in
> the 5 numerics (`{{signal_quality}}` / `{{clean_count}}` / `{{suspect_count}}` /
> `{{opacity_count}}` / `{{total_merged}}`), and remove the other two verdict
> branches.
>
> See `plans/Cef_Audit_Aware_Packet_Design.md` (commit `e551808`) edit #5 for the
> pre-commitment rationale.

---

# DNA Decoder v0.1 cef audit-aware release - 2026-05-26 (or {{run_date}})

## Status

- Release state: `audit-aware release`
- Scope: cached-strain ceftriaxone predictor with audit-aware honest output
- Canonical command surface: `python -m scripts.pipeline predict`
- Predictive posture: viable on real gate-B substrate (AUROC 0.895 / AUPRC 0.838 on N=49 usable strains)
- Reporting posture: **canonical_audit_aware** (was: debug-only)
- Interpretability posture: exploratory; gene-level ISM available via predict CLI without `--no-attribution`

## What changed from the 2026-05-26 debug-mode candidate

- Removed `--allow-missing-audit` requirement
- Added `--audit-merge-json wiki/cef_mechanism_phenotype_merge_2026-05-26.json` to gold-path command
- Added Audit framework section
- Added 4 canonical example artifacts (was 2)
- Updated current release assertion to verdict-conditional wording (below)

---

## ===== VERDICT BRANCH SELECTION (delete two of three) =====

### Branch A — `RUN_FULL_AND_CLEAN` (signal_quality ≥ 0.70)

```
> Cef cached-strain prediction is now product-viable + audit-aware on the gate-B
> substrate. The cef cohort cleared the SUSPEND gate at signal_quality={{signal_quality}}
> (≥0.70). Predictions ship audit-aware; categorical labels are interpretable.
> AUROC 0.895 / AUPRC 0.838 stand without scope-limit qualification beyond the
> standard v0 disclaimer (cached-strain only; not a clinical decision support tool).
```

### Branch B — `MIXED` (0.40 ≤ signal_quality < 0.70)

```
> Cef cached-strain prediction is now audit-aware on the gate-B substrate. The
> cef cohort fired MIXED at signal_quality={{signal_quality}} (in [0.40, 0.70)).
> Predictions ship audit-aware with the framing that label noise is a real
> confounder: of {{total_merged}} cohort strains, {{clean_count}} are CLEAN,
> {{suspect_count}} are SUSPECT, {{opacity_count}} have mechanism opacity. Treat
> AUROC 0.895 as descriptive of the noisy-label set, not a deployable performance
> claim. Clean-subset rerun recommended if `clean_count >= 20`.
```

### Branch C — `SUSPEND_CONDITION_4` (signal_quality < 0.40)

```
> Cef cached-strain prediction is now audit-aware on the gate-B substrate, but
> the cef cohort fired SUSPEND_CONDITION_4 at signal_quality={{signal_quality}}
> (<0.40): of {{total_merged}} cohort strains, only {{clean_count}} are CLEAN,
> {{suspect_count}} are SUSPECT, {{opacity_count}} have mechanism opacity. The
> AUROC 0.895 figure is uninterpretable against this label-noise profile;
> predictions ship with explicit informational-only framing + `suspend_gate_fired=True`
> propagation in every per-prediction JSON. Do not deploy as clinical decision
> support. Two follow-ups recommended:
>   1. if `opacity_count >= 5`, debug AMRFinder coverage before declaring labels unusable
>   2. otherwise, expand cef pool via PATRIC/NARMS/EuSCAPE backfill (the BV-BRC
>      4-drug feasibility census 2026-05-18 already established cef-S labels are scarce
>      at strict-MIC; relaxed-MIC cohort + audit framework is the canonical path).
```

## ===== END VERDICT BRANCH SELECTION =====

---

## Audit framework

The cef cohort was audited against three orthogonal signals per the cipro audit
discipline:

- **Mechanism audit** (`wiki/ceftriaxone_mechanism_audit_2026-05-26.json`):
  AMRFinderPlus per-strain mechanism detection across 50 cef-pool strains.
  Verdict: `PRIMARY_DOMINANT`. 26/26 R have at least one primary cef mechanism
  (acquired β-lactamase / extended-spectrum β-lactamase / AmpC hyperproduction).
  23/24 S have a silent primary mechanism hit — disambiguation needed via MIC
  tier.
- **MIC-tier audit** (`wiki/ceftriaxone_mic_audit_2026-05-26.json`):
  per-strain MIC classification under CLSI 2024 + EUCAST 14.0 E. coli
  breakpoints. Tier distribution: {{mic_tier_distribution}}.
- **Mechanism × MIC merge** (`wiki/cef_mechanism_phenotype_merge_2026-05-26.json`):
  joint per-strain `noise_class` + cohort-level `gate_verdict`. Gate verdict:
  `{{gate_verdict}}` at signal_quality={{signal_quality}}.

Per-prediction `audit_verdict` field now propagates:

- `noise_class` (per-strain bucket: CLEAN/SUSPECT/OPAQUE/NOISY/OTHER × R/S)
- `mechanism_opacity_flag` (True if R has no known primary mechanism — explains
  the "model thinks it's S" miss pattern)
- `mic_tier` (HIGH_R / HIGH_S / DECISIVE / BORDERLINE / AMBIGUOUS / CONFLICT /
  NO_MIC under CLSI + EUCAST joint call)
- `primary_mechanisms` (list of cef-relevant mechanisms detected; plural even
  when one)
- `co_resistance_modifiers` (efflux / regulatory / porin_loss — modifiers, not
  primary cef drivers)
- `cohort_gate_verdict` (cohort-level — the gate this strain trained against)
- `suspend_gate_fired` (True iff cohort_gate_verdict contains "SUSPEND")
- `verdict_explanation` (full text if suspend_gate_fired=True)

Schema-contract regression tests pin these consumer fields at
`tests/test_drug_mechanism_phenotype_merge_contract.py` — any future merge-output
drift fails CI.

## Gold-path command (audit-aware)

```bash
uv run python -m scripts.pipeline predict \
  --drug ceftriaxone \
  --strain-id 562.12960 \
  --model-path data/processed/models/ceftriaxone_nucleotide_transformer.pkl \
  --cache "C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_cache/embeddings/nt_gate_b_cohort_67.h5" \
  --audit-merge-json wiki/cef_mechanism_phenotype_merge_2026-05-26.json \
  --no-attribution \
  --output result.json
```

Note: `--allow-missing-audit` is no longer needed; `--no-attribution` still
optional (drop it to get gene-level ISM at the cost of ~3-5 min/strain).

## Required inputs

- trained cef model pickle (`data/processed/models/ceftriaxone_nucleotide_transformer.pkl`)
- cached strain ID present in the cef HDF5 cache
- cef cache path
- audit merge JSON (`wiki/cef_mechanism_phenotype_merge_2026-05-26.json`)

## Expected outputs

- `result.json` + `result.md` (markdown sidecar)

Required output fields:

- `prediction`
- `calibrated_probability`
- `confidence_tier`
- `audit_verdict` (NEW — full audit context per the 9 fields above)
- `provenance` with `reporting_mode = canonical_audit_aware`
- `attribution_scope_confidence` (HIGH/PARTIAL/INDETERMINATE)

## Real substrate and model

Working cef substrate:

- cohort file: `data/processed/gate_b_cohort.parquet`
- cef pool size: `50`
- cef label balance: `26R / 24S`
- usable train/eval set: `49`

Dedicated cef cache:

- `C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_cache/embeddings/nt_gate_b_cohort_67.h5`

Trained model:

- `data/processed/models/ceftriaxone_nucleotide_transformer.pkl`
- CV strategy: `loso`
- CV grouping: `strain_id`
- primary CV AUROC: `0.895`
- AUPRC: `0.838`

Duplicate-accession audit:

- `reports/current_cef_duplicate_accession_audit_2026-05-25.md`
- verdict: `PASS` (AUROC 0.895 is NOT inflated by same-genome leakage)

## 4 canonical example artifacts

Per the revised design memo (4 examples, not 3 — both model misses surfaced for
release honesty):

| Strain ID | Label | Predicted | Notes |
|---|---|---|---|
| `562.12960` | R | R | CLEAN_R + primary mechanism |
| `562.7572` | S | S | CLEAN_S + no primary mechanism |
| `562.28389` | S | R | FP miss; surfaces OPAQUE_R_no_mechanism / borderline behavior |
| `562.7695` | R | S | FN miss; surfaces SUSPECT_S_silent_primary_mechanism or BORDERLINE_R |

Artifacts:

- `reports/dna_decoder_v0_1_cef_canonical_example_R_2026-05-26.{json,md}`
- `reports/dna_decoder_v0_1_cef_canonical_example_S_2026-05-26.{json,md}`
- `reports/dna_decoder_v0_1_cef_canonical_example_FP_miss_2026-05-26.{json,md}`
- `reports/dna_decoder_v0_1_cef_canonical_example_FN_miss_2026-05-26.{json,md}`

## Cross-path validation (preserved from debug-mode candidate)

- requested cef-pool strains: `50`
- completed usable strains: `49`
- cached-strain path vs genome-input path prediction concordance: `49 / 49`
- label alignment on both paths: `47 / 49`
- max absolute probability delta: `0.063148`
- mean absolute probability delta: `0.002305`
- only two shared model misses:
  - `562.28389` expected `S`, both paths predicted `R`
  - `562.7695` expected `R`, both paths predicted `S`

Both misses are at decision-boundary probabilities and both are flagged as the
canonical examples (FP_miss + FN_miss) above.

## Current release assertion (selected verdict branch text from above)

> {{insert selected branch text from "VERDICT BRANCH SELECTION" section}}

## Handoff references

- `reports/dna_decoder_v0_1_cef_cached_handoff_2026-05-25.md`
- `reports/current_cef_duplicate_accession_audit_2026-05-25.md`
- `reports/dna_decoder_v0_1_cef_cached_vs_genome_full_validation_2026-05-26.md`
- `wiki/ceftriaxone_mechanism_audit_2026-05-26.md`
- `wiki/ceftriaxone_mic_audit_2026-05-26.md` (NEW — Artifact 2)
- `wiki/cef_mechanism_phenotype_merge_2026-05-26.md` (NEW — Artifact 3)
- `reports/dna_decoder_v0_1_cef_audit_progress_2026-05-26.md`

---

## Template usage instructions

1. Confirm `wiki/cef_mechanism_phenotype_merge_2026-05-26.json` exists.
2. Open the JSON; read `gate_verdict`, `signal_quality`, `clean_count`,
   `suspect_count`, `opacity_count`, and `len(per_strain)` (= total_merged).
3. Open `wiki/ceftriaxone_mic_audit_2026-05-26.json`; tally tier distribution
   into a short string (e.g., "HIGH_R=12, HIGH_S=18, BORDERLINE=15, NO_MIC=5").
4. Save this template as `wiki/dna_decoder_v0_1_cef_audit_aware_release_2026-05-26.md`.
5. In the saved file:
   - delete the two verdict branches that don't match
   - replace `{{signal_quality}}` etc. with actual numerics from step 2
   - replace `{{mic_tier_distribution}}` with the string from step 3
   - replace `{{run_date}}` with `_date.today().isoformat()`
   - replace the "Current release assertion" placeholder with the selected branch text
6. `git add wiki/dna_decoder_v0_1_cef_audit_aware_release_2026-05-26.md`
7. `git commit -m "feat(cef-audit): audit-aware release packet (verdict={{gate_verdict}})"`
8. `git push origin main`

Mechanical update path: ~10 min once the merge JSON exists.
