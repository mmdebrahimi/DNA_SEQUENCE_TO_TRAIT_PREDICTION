# Cef Cached-Strain v0.1 Slice — Plan

> Smallest-credible cef cached-strain v0.1 slice per Codex 2026-05-25 handoff ("cached-strain cef first; one axis at a time"). Maps onto EP-2 cef sub-track in `plans/Post_V0_EP_Ladder_Plan.md`.

**Status:** DRAFT 2026-05-26. Pre-execution; gated on Codex pushing v0.1 cipro slice 1 artifacts (retrained model + runtime `pipeline.py` changes + reports/ release packets) to origin first.

---

## What this slice IS

Train + ship a cef cached-strain decoder on the existing N=12 mini cohort. ONE axis change from v0.1 cipro: drug = ceftriaxone instead of ciprofloxacin. Everything else (cached-strain input mode, v0 JSON schema, audit-aware output, locus-tag-prefix proxy for attribution_scope_confidence) stays.

## What this slice is NOT

- NOT a real N≥75 cef cohort — that needs cef-S label backfill from an alternative source (PATRIC / NARMS / EuSCAPE) per `FUTURE_FEATURES.md` "v0.1 cef decoder substrate via cef-S label backfill — 18/05/2026". Out of scope; v0.2 territory.
- NOT genome-input cef — per Codex's framing, "one axis at a time"; cef genome-input is v0.2+.
- NOT a new architecture — uses the same NT-frozen-pool + XGBoost as cipro v0/v0.1.
- NOT a multi-drug refactor — drug-specific paths stay independent; mic_tiers.py already handles cef.

---

## Existing assets (what already exists in repo on origin)

| Asset | Path | Status |
|---|---|---|
| Mini cef cohort (N=12, 6R/6S) | `data/processed/gate_b_mini_cef_cohort.parquet` | ✓ Built 2026-05-17; reuses cipro N=38 strain pool |
| Cef smoke verdict | `wiki/smoke_gate_12strain_ceftriaxone_2026-05-17.md` + `wiki/EP2_cef_tet_verdict_2026-05-17.md` | ✓ PASS (NT-XGBoost 0.833 = k-mer 0.833 on smoke) |
| Shared NT cache w/ cef strain embeddings | `D:/dna_decode_cache/embeddings/nt_n40_cipro.h5` (Codex's Precision 7780; needs pull-side equivalent) | ✓ Contains all 12 cef mini-cohort strains |
| Cef per-drug catalogs | `dna_decode/data/mic_tiers.py` (`DRUG_BREAKPOINTS["ceftriaxone"]`, `DRUG_AMRFINDER_CLASSES["ceftriaxone"]`, `DRUG_LOCI_BY_MECHANISM["ceftriaxone"]`) | ✓ Configured: CLSI R≥4 S≤1; AMRFinder classes = BETA-LACTAM + CARBAPENEM + CEPHALOSPORIN + MULTIDRUG; loci = β-lactamase families (CTX-M, SHV, TEM, OXA, CMY, AmpC) + porin (ompC/F) |
| Drug-agnostic mechanism audit | `scripts/drug_mechanism_audit.py` | ✓ Shipped 2026-05-24 |
| EP-2.5 cohort build mode (drug-agnostic) | `scripts/bvbrc_strict_mic_4drug_census.py --build` | ✓ Shipped 2026-05-25; works for cef if categorical-MIC labels are sufficient |
| v0.1 cipro `pipeline.py predict --strain-id` path | `scripts/pipeline.py` | ✓ on origin (cached-strain mode; drug-parameterized via `--drug`) |
| v0.1 cipro `pipeline.py predict --genome-fasta` path | Codex's Precision 7780 ONLY (NOT yet pushed) | ✗ Awaiting push |

---

## Gaps (what's missing for this slice to ship)

| Gap | Effort | Blocker |
|---|---|---|
| **Trained leakage-safe cef classifier on N=12 mini cohort** | ~30 min (training on 12 strains is fast) | Needs cef cohort + cef NT cache subset; both exist |
| **Cef merge-gate audit JSON** | ~1 hr (run AMRFinder via `drug_mechanism_audit.py --drug ceftriaxone` on the 12 cef strains; emit per-strain mechanism + opacity flags) | Needs Docker + AMRFinder DB; Codex has on Precision 7780 |
| **Cef predict end-to-end smoke** | ~10 min (`pipeline.py predict --drug ceftriaxone --strain-id X --model-path <cef.pkl> --audit-merge-json <cef-merge.json>`) | Needs above two |
| **Cef release packet** | ~30 min (markdown template matching `reports/dna_decoder_v0_release_candidate_2026-05-24.md` pattern) | Needs above three |
| **Cef cross-path consistency validation** (cached-strain only here, so NOT applicable — Codex's 2026-05-25 cross-path validation was cipro-cached-vs-cipro-genome-input) | n/a | n/a |

Total slice effort: **~2-3 hr on Precision 7780** (Codex has the model + Docker + cache). Zero blocking work on Claude side.

---

## Implementation steps (Codex executes on Precision 7780 after v0.1 push)

### Step 1: Confirm cef mini cohort + NT cache integrity

```bash
uv run python scripts/probe_nt_cache.py \
  --cohort data/processed/gate_b_mini_cef_cohort.parquet \
  --cache D:/dna_decode_cache/embeddings/nt_n40_cipro.h5 \
  --refseq-cache D:/dna_decode_cache/refseq
# Expected: ALL_COMPLETE on 12 cef strains
```

### Step 2: Train cef classifier (leakage-safe `leave_one_accession_out` CV)

```bash
uv run python -m scripts.pipeline train \
  --drug ceftriaxone \
  --model nucleotide_transformer \
  --cohort data/processed/gate_b_mini_cef_cohort.parquet \
  --cache D:/dna_decode_cache/embeddings/nt_n40_cipro.h5 \
  --output data/processed/models/ceftriaxone_nucleotide_transformer_mini.pkl
# Expected: CV AUROC ≥ 0.70 (smoke verdict 2026-05-17 was 0.833)
```

### Step 3: Build cef merge-gate audit JSON

```bash
uv run python -m scripts.drug_mechanism_audit \
  --drug ceftriaxone \
  --cohort data/processed/gate_b_mini_cef_cohort.parquet \
  --refseq-cache D:/dna_decode_cache/refseq \
  --out-root data/amrfinder_runs/ceftriaxone_mini \
  --output wiki/ceftriaxone_mechanism_audit_2026-05-26.md
# Then build merge gate JSON (per-strain mechanism + opacity flag) - 
# reuse cipro_mechanism_phenotype_merge.py logic OR write a thin
# drug-parameterized merge if cipro_*.py refuses to generalize cleanly.
```

### Step 4: Cef end-to-end predict smoke

```bash
# Hold one cef strain out; predict against the trained classifier.
uv run python -m scripts.pipeline predict \
  --drug ceftriaxone \
  --strain-id <held-out-cef-strain-id> \
  --model-path data/processed/models/ceftriaxone_nucleotide_transformer_mini.pkl \
  --cache D:/dna_decode_cache/embeddings/nt_n40_cipro.h5 \
  --annotations <held-out-cef.gff3> \
  --audit-merge-json wiki/ceftriaxone_mechanism_phenotype_merge_2026-05-26.json \
  --output reports/cef_v0_1_predict_example_2026-05-26.json
# Verify: JSON has prediction + calibrated_probability + confidence_tier +
# attribution_scope_confidence + top_k_attribution + audit_verdict + provenance.
```

### Step 5: Cef release packet (matches cipro pattern)

Write `reports/dna_decoder_cef_v0_1_release_candidate_2026-05-26.md` with:
- Sample command
- Real example output (from Step 4)
- HONEST scope caveats (N=12 mini cohort; statistical floor; v0.2 needs N≥75 with cef-S label backfill)
- Comparison to cipro v0/v0.1: same predict CLI, same JSON schema, same honest-output discipline, smaller cohort

### Step 6: Tests for cef path

- Verify `tests/test_pipeline_predict_e2e.py` synthetic-fixture tests pass with `--drug ceftriaxone` (drug-agnostic path)
- Add 1-2 tests pinning that cef per-drug catalogs in mic_tiers.py + cef predict round-trip emit correct provenance

---

## Minimum success gate

A cef v0.1 cached-strain slice ships when ALL of:

1. **Functional (HARD gate):** `pipeline.py predict --drug ceftriaxone --strain-id X --cache Y.h5 --model-path Z.pkl --audit-merge-json A.json --output result.json` runs end-to-end on a held-out cef mini-cohort strain WITHOUT crashing AND emits v0-schema JSON + markdown.
2. **Predictive:** primary CV AUROC ≥ 0.70 on the N=12 mini cef cohort. (Note: variance is high at N=12; this is a smoke-tier gate, NOT clinical-grade.)
3. **Honest (audit-aware, HARD gate):** when `--audit-merge-json` is supplied with the cef merge-gate JSON, the `audit_verdict` field propagates. SUSPEND framing carries if the cef cohort's signal quality is below threshold.
4. **Documentation:** `reports/dna_decoder_cef_v0_1_release_candidate_2026-05-26.md` exists + names the explicit scope: "N=12 mini cohort; cef-S label backfill is the next axis to expand."
5. **Reproducibility:** the exact commands above reproduce the result on Precision 7780 from a fresh checkout.

---

## Honest scope caveats (must be in the release packet)

- **N=12 cef mini cohort is statistically thin.** AUROC noise band at this N is ±0.20. A measured 0.833 on smoke is consistent with anything from ~0.65 to ~1.0 at 95% CI. Don't overclaim.
- **Mini cohort is FILTERED FROM cipro N=38 strain pool, not built fresh from cef AST data.** Strain selection may reflect cipro biology, not cef biology. Mechanism distribution may be skewed.
- **BV-BRC strict-MIC is INFEASIBLE for cef at scale** per the 2026-05-18 census (66 HIGH_R / 2 HIGH_S of 4,567 cef AST rows). At categorical-MIC bar (relaxed; includes DECISIVE), counts may be better — but not yet tested.
- **No external benchmark yet** vs AMRFinderPlus / RGI for cef phenotype prediction. That's EP-1B territory once cipro v0.1 finishes its EP-1B benchmark.

---

## Why this slice (vs alternatives)

| Alternative | Why not |
|---|---|
| Genome-input cef (cef + novel input simultaneously) | Codex's framing: "one axis at a time." Two-axis changes raise the conflation risk. |
| Build real N≥75 cef cohort first | Cef-S label backfill is the actual blocker; that's v0.2 work. Smallest-credible slice reuses what exists. |
| Skip cef; do tet instead | Tet is the distributed-mechanism FAIL drug (AUROC 0.400 anti-predictive); needs EP-1.5 architecture decision first. NOT a "smallest" slice. |
| Tighten cipro interpretability before cef | Per Codex's 2026-05-25 recommendation: "do NOT spend the next cycle on more cipro validation of the same kind / more cipro interpretability work." Cef expansion adds more product value. |

---

## Risk flags

- **R1 (HIGH):** Codex's v0.1 slice 1 (genome-input cipro) is NOT yet on origin per 2026-05-26 sync check. This plan is written against the post-v0.1-push state; if Codex never pushes, this side's view of `pipeline.py` will lack the `--allow-missing-audit` flag + `cv_strategy` train-time persistence Codex added. Mitigation: this plan's Step 2 (cef classifier train) depends on Codex's updated `pipeline.py train` to write `cv_strategy` provenance; if those changes aren't on origin, Codex needs to either push OR Step 2 happens on Precision 7780 anyway (so the unpushed code is in scope).
- **R2 (MEDIUM):** N=12 may not produce a stable CV AUROC. If cef CV at the leakage-safe `leave_one_accession_out` setting drops below 0.70, the slice ships with the scope-limit doc (matching v0 cipro FAIL ship pattern); do NOT chase calibration tweaks at this cohort scale.
- **R3 (MEDIUM):** Cef merge-gate audit JSON requires either generalizing `scripts/cipro_mechanism_phenotype_merge.py` to drug-agnostic OR writing a small cef-specific merge alongside. Per CLAUDE.md gotcha "scripts/cipro_*.py left as cipro-specific by design," the cleaner path is to write `scripts/drug_mechanism_phenotype_merge.py` as a parallel to the existing `drug_mechanism_audit.py`. ~1 hr extra work; flag explicitly.
- **R4 (LOW):** ε=0.01 same-strain-parity tolerance in EP-1A memo was based on the cipro v0.1 finding (Codex reports max delta 0.011599 — just barely above ε). For cef, expect similar drift. The slice does NOT require cross-path parity (cef-cached-strain only); R4 is informational only.

---

## What this plan deliberately does NOT cover

- The actual held-out cef strain pick for Step 4 — Codex picks on Precision 7780.
- Cef cross-organism extension — out of scope.
- Cef external benchmark vs AMRFinderPlus / RGI — that's EP-1B-style work for cef, separate slice.
- Multi-drug refactor of `scripts/pipeline.py` predict — not needed; cef path already drug-agnostic via `--drug`.
- Cef v0.2 (real cohort N≥75 with cef-S backfill) — separate plan; needs cef-S label source first.

---

## Open questions for Codex / user

1. Which held-out cef strain serves as the predict-smoke target? Recommend reserving 1 of the 6 cef-R + 1 of the 6 cef-S for held-out; train on the other 10. OR run LOSO across all 12 + pick the best example for the release packet.
2. Should cef merge-gate include the same SUSPEND_CONDITION_4 threshold as cipro? Cef cohort signal quality is untested. Recommend running the merge first, then evaluating SUSPEND threshold post-hoc.
3. Should this cef slice tag a release (`v0.1-cef`) or share `v0.1` with cipro genome-input? Codex picks; recommend `v0.1-cef` as a sub-tag.
4. Is the "cef-S label backfill from PATRIC/NARMS/EuSCAPE" v0.2 work the right substrate expansion, OR should v0.2 pivot to cef strict-MIC at categorical-MIC bar (untested feasibility)?

---

## Recommended sequencing

1. **First:** Codex pushes v0.1 slice 1 artifacts to origin (retrained cipro pickle, runtime `pipeline.py` changes, scope-limit doc, release packets, genome-input + validation packets). This unblocks cross-machine sync — same blocker as 2026-05-24.
2. **Then:** Codex executes Steps 1-6 above on Precision 7780.
3. **Then:** push cef v0.1 artifacts to origin.
4. **Then on this side:** sync check + update ledger.

Total elapsed: ~3-4 hours of Codex work + the push roundtrip.

## Bottom line

Smallest credible cef cached-strain v0.1 slice = reuse the existing N=12 cef mini cohort + train a cef classifier + emit a cef release packet matching the cipro pattern, all on Precision 7780. ~2-3 hr of compute work, zero new architecture. Honest scope: smoke-tier N=12 cef; v0.2 cef requires cohort expansion via cef-S label backfill (separate plan).
