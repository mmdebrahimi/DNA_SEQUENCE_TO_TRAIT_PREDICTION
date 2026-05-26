# Cef Audit-Aware Packet — Design Memo

> Smallest credible design for the cef audit-aware closeout that matches the cipro release discipline **without reopening broader scope**. Per Codex's 2026-05-26 ask ("propose the smallest credible design for a cef audit-aware packet").

**Status:** DRAFT 2026-05-26.
**Anchors on:** `wiki/dna_decoder_v0_1_cef_overnight_handoff_2026-05-26.md` (Codex's overnight ship: 49/49 concordance + 47/49 label alignment + AUROC 0.895 + leakage audit PASS); `plans/Cef_V0_1_Promotion_Slice_Plan.md` (Steps 2-5 of which this memo refines); cipro audit infrastructure 2026-05-17 (the analog this matches).
**Supersedes:** `plans/Cef_V0_1_Promotion_Slice_Plan.md` Steps 2-5 (this memo provides the design at higher fidelity).

---

## What "matches cipro release discipline" means concretely

The cipro release packet (`reports/dna_decoder_v0_release_candidate_2026-05-24.md`) carries these audit-aware features that cef currently lacks:

1. **`audit_verdict` field populated** in the v0 JSON output, sourced from a per-strain merge-gate JSON.
2. **`suspend_gate_fired` boolean** + explanation text when training cohort signal quality is below threshold.
3. **`noise_class` per strain** (CLEAN_R_primary_mechanism / OPAQUE_R_no_mechanism / SUSPECT_S_silent_mechanism / etc.).
4. **`mechanism_opacity_flag`** distinguishing tool-failure vs label-failure.
5. **`mic_tier`** (HIGH_R / HIGH_S / DECISIVE_R / DECISIVE_S / BORDERLINE / AMBIGUOUS / CONFLICT / NO_MIC) per strain.
6. **`primary_mechanisms`** list + **`co_resistance_modifiers`** list per strain.
7. **`reporting_mode = canonical_audit_aware`** (vs the current cef `--allow-missing-audit` debug mode).

The cipro discipline came from 4 audit tiers: AMRFinder mechanism audit + raw BV-BRC MIC rejoin + mechanism×MIC merge with structurally-enforced SUSPEND gate + cohort-quality classification.

Cef needs the same 4 tiers — adapted to cef AMRFinder Class filter + cef breakpoints + cef mechanism catalog. All per-drug catalogs already exist in `dna_decode/data/mic_tiers.py`.

---

## The smallest credible design

Five concrete artifacts; ~3-4 hr of Codex compute on Precision 7780.

### Artifact 1 — Cef AMRFinder mechanism audit

**What:** Run AMRFinder on all 50 cef-pool strains; emit per-strain mechanism calls filtered to cef-relevant AMRFinder classes (BETA-LACTAM + CARBAPENEM + CEPHALOSPORIN + MULTIDRUG per `mic_tiers.amrfinder_classes_for("ceftriaxone")`).

**How:** Reuse `scripts/drug_mechanism_audit.py` (shipped 2026-05-24; drug-agnostic; in the prior Claude→Codex bundle):

```bash
uv run python -m scripts.drug_mechanism_audit \
  --drug ceftriaxone \
  --cohort data/processed/gate_b_cohort.parquet \
  --refseq-cache "C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_cache/refseq" \
  --out-root data/amrfinder_runs/ceftriaxone_gate_b \
  --output wiki/ceftriaxone_mechanism_audit_2026-05-26.md
```

**Cost:** ~50 strains × ~95 s = ~80 min Docker on Precision 7780.

**Output:** `wiki/ceftriaxone_mechanism_audit_2026-05-26.{md,json}` with per-strain `mechanisms_present`, `primary_mechanism_class`, `mech_hits` per mechanism class.

**Cef-specific mechanisms to expect (per `mic_tiers.DRUG_LOCI_BY_MECHANISM["ceftriaxone"]`):**
- ESBL_class_C (CMY-2, AmpC variants)
- ESBL_extended_spectrum (CTX-M family — primary driver)
- carbapenemase (KPC, NDM, OXA-48 — rare but informative)
- intrinsic_beta_lactamase (TEM, SHV, OXA — common)
- porin_loss (ompC, ompF — co-resistance modifier)

### Artifact 2 — Cef MIC tier audit

**What:** Per-strain cef MIC re-join + tier classification using `mic_tiers.classify_tier`. Output: per-strain `mic_tier` field for the merge.

**How:** New thin script `scripts/cef_mic_audit.py` OR generalize `scripts/cipro_mic_audit.py` to take `--drug` (per CLAUDE.md gotcha "leave cipro_*.py cipro-specific by design"; recommend new file).

**Inputs:**
- `data/processed/gate_b_cohort.parquet` (cef pool strain_ids)
- Raw BV-BRC AST CSV (filtered to ceftriaxone rows)

**Tier classification per strain:** Aggregate MIC values per strain → call `mic_tiers.classify_tier(mics, distinct_calls, mic_tiers.breakpoints_for("ceftriaxone"))` → emit `wiki/ceftriaxone_mic_audit_2026-05-26.{md,json}` with per-strain MIC tier.

**Cef breakpoints (per `mic_tiers.DRUG_BREAKPOINTS["ceftriaxone"]`):** CLSI R ≥ 4 / S ≤ 1; EUCAST R ≥ 2 / S ≤ 1.

**Cost:** ~15 min (pandas-only; no GPU; no Docker).

### Artifact 3 — Cef mechanism × MIC merge with SUSPEND gate

**What:** Per-strain merge of Artifacts 1 + 2; computes `noise_class` + `mechanism_opacity_flag` + cohort-level signal quality + `gate_verdict`.

**How:** Write NEW `scripts/drug_mechanism_phenotype_merge.py` parallel to `scripts/cipro_mechanism_phenotype_merge.py` (per CLAUDE.md gotcha). Drug-parameterized; uses `mic_tiers.classify_gene_symbol(drug, ...)` + `mic_tiers.primary_mechanisms_for(drug)`.

**Cef-specific SUSPEND threshold (calibrated below):**
- `clean_count` = strains with HIGH_R MIC + primary mechanism (CTX-M / CMY / SHV / etc.) present
- `opacity_count` = strains with HIGH_R MIC but NO primary mechanism (tool failed to detect biology that should be there)
- `signal_quality` = clean_count / (clean_count + opacity_count + suspect_count)
- **Cef SUSPEND threshold:** `signal_quality < 0.50` → SUSPEND_CONDITION_4. (Cipro used 0.40; cef can be tighter because plasmid β-lactamases are easier to detect than QRDR point mutations.)

**Output:** `wiki/cef_mechanism_phenotype_merge_2026-05-26.{md,json}` with:
- Per-strain `noise_class` (CLEAN_R_primary / OPAQUE_R_no_mechanism / SUSPECT_S_silent_mechanism / etc.)
- Per-strain `mechanism_opacity_flag` (True if HIGH_R + no primary mechanism)
- Per-strain `primary_mechanisms` list
- Per-strain `co_resistance_modifiers` (porin_loss, efflux)
- Cohort-level `gate_verdict` (RUN_FULL_AND_CLEAN / MIXED / SUSPEND_CONDITION_4)
- Cohort-level `signal_quality` numeric

**Cost:** ~30 min (code new + run).

### Artifact 4 — Re-run cef predict examples in canonical_audit_aware mode

**What:** Run `pipeline.py predict` on the same 2 strains Codex used overnight (562.12960 R + 562.7572 S) AND a third strain that exposes the merge gate's behavior — recommended `562.28389` (a shared model miss: S labeled, but predicted R; expected to surface `OPAQUE_R_no_mechanism` IF the model is calling R without mechanism support; informative for the scope-limit doc).

**How:** Drop `--allow-missing-audit` + `--no-attribution`; supply the new cef merge JSON:

```bash
uv run python -m scripts.pipeline predict \
  --drug ceftriaxone \
  --strain-id 562.12960 \
  --model-path data/processed/models/ceftriaxone_nucleotide_transformer.pkl \
  --cache <cef-nt-cache> \
  --annotations <refseq-cache>/<acc>.gff3 \
  --audit-merge-json wiki/cef_mechanism_phenotype_merge_2026-05-26.json \
  --output reports/dna_decoder_v0_1_cef_canonical_example_R_2026-05-26.json
# Same for 562.7572 (S) and 562.28389 (S labeled, R predicted -- diagnostic)
```

**Cost:** ~30 min.

**Verify outputs include:**
- `audit_verdict.noise_class` populated
- `audit_verdict.mechanism_opacity_flag` populated
- `audit_verdict.cohort_gate_verdict` populated
- `audit_verdict.primary_mechanisms` list (e.g., CTX-M-15)
- `provenance.reporting_mode = canonical_audit_aware`
- NOT `provenance.reporting_mode = debug_internal`

### Artifact 5 — Updated cef release candidate (audit-aware)

**What:** Replace `reports/dna_decoder_v0_1_cef_cached_release_candidate_2026-05-26.md` (or write a new dated version) that:
- Updates the "Important caveat" section: cef is NOW canonical_audit_aware (not debug-mode).
- Includes the 3 canonical example outputs from Artifact 4.
- Adds an "Audit framework" section explaining the cef merge gate.
- Names the cohort-level `gate_verdict` (likely RUN_FULL_AND_CLEAN given the 47/49 label alignment + AUROC 0.895; cef's signal quality should be HIGH).

**Cost:** ~30 min.

---

## Total effort

~3.5 hr Codex compute on Precision 7780:
- Artifact 1 (AMRFinder audit): ~80 min Docker
- Artifact 2 (MIC audit): ~15 min CPU
- Artifact 3 (merge + new script): ~30 min code + run
- Artifact 4 (canonical examples): ~30 min compute
- Artifact 5 (release packet update): ~30 min writing

---

## Success gate (the cef audit-aware closeout ships when ALL of these pass)

1. `wiki/ceftriaxone_mechanism_audit_2026-05-26.{md,json}` exists with per-strain mechanism calls for 50 cef strains.
2. `wiki/ceftriaxone_mic_audit_2026-05-26.{md,json}` exists with per-strain MIC tier.
3. `wiki/cef_mechanism_phenotype_merge_2026-05-26.{md,json}` exists with per-strain `noise_class` + cohort-level `gate_verdict`.
4. 3 cef predict examples in `reports/dna_decoder_v0_1_cef_canonical_example_*.json` with `provenance.reporting_mode = canonical_audit_aware` + populated `audit_verdict`.
5. Updated `reports/dna_decoder_v0_1_cef_cached_release_candidate_<date>.md` references the audit artifacts + drops the "debug-mode only" caveat.
6. Pushed to origin.

---

## Honest scope limits (must be in the updated release packet)

- **Cef cohort N=49 usable.** Smaller than cipro N=146. AUROC noise band ~ ±0.06 at N=49 (vs ~±0.04 at N=146). Don't overclaim.
- **Cef cohort is FILTERED from `gate_b_cohort.parquet`**, not built fresh from cef AST. Strain selection inherits any cipro biology bias.
- **2 shared model misses (562.28389, 562.7695)** at decision-boundary probabilities (0.49-0.59). May reflect label noise OR architectural mismatch on edge cases. Surface as known limitations.
- **Cef strict-MIC infeasible at scale** per 2026-05-18 census (66 HIGH_R / 2 HIGH_S of 4,567 AST rows). v0.2 cef requires cef-S backfill from PATRIC / NARMS / EuSCAPE.
- **Cef interpretability NOT promoted** beyond exploratory. Same posture as cipro v0/v0.1: predictive output is canonical; interpretability is informational.

---

## What this design deliberately does NOT cover

- **Cef external benchmark vs AMRFinderPlus + RGI** — that's EP-1B-equivalent for cef; separate slice. Codex's recommendation: "non-AMR phenotype scoping becomes a strong next candidate" AFTER this audit-aware closeout.
- **Cef genome-input** — already shipped + validated (49/49 concordance); not in this slice.
- **Cef multi-organism extension** — Phase 3 territory; not this slice.
- **Tet / gent / other AMR drugs** — EP-1.5 architecture-decision-gated; not this slice.
- **Architecture rework** — none. Same NT mean-pool + XGBoost as cipro.

---

## Risk flags

- **R1 (LOW):** Cef SUSPEND threshold (0.50) is calibrated heuristically; if cohort signal quality is near 0.50, it'll be borderline. Mitigation: report the raw `signal_quality` numeric; threshold is a guideline, not a hard ship gate.
- **R2 (LOW):** Cef MIC tier audit needs raw BV-BRC cef AST CSV; if not on Precision 7780, needs to be downloaded. Mitigation: probably already present (Codex used it for the 4-drug census 2026-05-18).
- **R3 (MEDIUM):** Writing `scripts/drug_mechanism_phenotype_merge.py` is new code (~1 hr). Should match the cipro version's API + per-strain JSON shape so the existing `_load_audit_verdict` helper in `pipeline.py` works unchanged. Mitigation: reference `scripts/cipro_mechanism_phenotype_merge.py` (already in repo).
- **R4 (LOW):** 2 shared model misses (562.28389, 562.7695) — could be label noise OR architectural failure. Mitigation: surface in the release packet's "known limitations" section; defer investigation to post-ship.

---

## Recommended order (Codex executes on Precision 7780)

1. Artifact 1 (AMRFinder audit) — kick off in background (~80 min).
2. While running: write `scripts/drug_mechanism_phenotype_merge.py` (the new file).
3. Artifact 2 (MIC audit) — run after Artifact 1 (or in parallel; they don't overlap).
4. Artifact 3 (merge) — runs after 1 + 2.
5. Artifact 4 (canonical predict examples) — runs after 3.
6. Artifact 5 (updated release packet) — write last.
7. Commit + push to origin.

---

## After cef audit-aware closeout

Per Codex's progress report 2026-05-26 + this side's `plans/Trait_Decoding_Roadmap.md`:

> "after cef audit-aware closeout, non-AMR phenotype scoping becomes a strong next candidate"

Concretely: invoke `/idea-anchor` with the EP-4 recommended candidate from `plans/EP_4_Non_AMR_Phenotype_Candidates.md` → **E. coli pathotype prediction** (EPEC / EHEC / ETEC / UPEC / EAEC / commensal multiclass) using EnteroBase substrate.

---

## Bottom line

5 artifacts; ~3.5 hr Codex compute; zero new architecture. Closes the last "debug-mode" gap in the cef release surface + matches the cipro discipline that v0 + v0.1-cipro already have. After this ships, non-AMR phenotype scoping (Phase 4) is the cleanest next jump per the roadmap.
