# Cef Mechanism Audit Plan

> Smaller alternative to the rejected 11-step EP2-Step-3 technical plan. Determines whether NT-XGBoost's cef AUROC 0.833 is genuine β-lactamase attribution vs lineage-tracking. Plan revised 2026-05-17 after /brainstorm round 2 caught 9 grounded issues — corrections absorbed in-place per the 2026-05-14 HIGH-salience LESSON.

---

## Problem Statement

EP2 cef smoke PASSED 2026-05-17 with NT-XGBoost AUROC 0.833 = k-mer 0.833 (gap 0 pp). The load-bearing follow-up question is: **is NT's cef AUROC genuine β-lactamase attribution, or is it lineage-tracking on a 12-strain mini-cohort?**

The /technical-plan originally generated for EP2 source-plan Step 3 was rejected by /review (smoke-runner variant expansion at N=12 with ±0.19 CI cannot answer mechanism attribution). The first version of this plan (3 steps) was then /brainstorm-reviewed in 2 rounds and surfaced 9 grounded issues:

1. **Step 3 gate inverted** — skip-attribution-on-mechanism-dominance is the opposite of what the cef question needs.
2. **Verdict bands brittle at N=6** — cipro hit 18/20 (90%) well above its 70% cliff; cef would be one-strain-flip-away. Plus the CI-lower-bound rule is mathematically unreachable for n=6.
3. **Plan confused cohort biology with model behavior** — Step 2 couldn't answer "is NT 0.833 mechanism-driven?" without joining NT predictions to AMRFinder calls.
4. **Cross-tab diagnostic missing** — the cheap direct comparison (`true × NT-pred × AMRFinder-hit`) was absent.
5. **Cohort provenance under-treated** — cef strains filtered from cipro N=38 cohort, not BV-BRC-wide cef-balanced.
6. **Per-strain prediction artifact unavailable** — today's cef smoke markdown has only aggregate variant rows; cross-tab can't be produced from existing artifacts.
7. **Wrong NT head** — the saved plan tested NT-LR but the cef PASS was NT-XGBoost.
8. **Attribution port doesn't include per-prediction join** — cipro preflight pattern aggregates across R strains; cef needs per-strain join to answer "were the LOSO PASS calls mechanism-driven?"
9. **Mash-clade lineage diagnostic absent** — if the question is explicitly lineage-vs-mechanism, the lineage baseline isn't optional.

This revision (2026-05-17 post-brainstorm) restructures into 8 implementation steps (0-5 plus 2a + 2b parallel diagnostics) across 6 waves.

## Design Decisions

### D1: Drop the 11-step EP2-Step-3 technical plan

**Decision:** Do NOT execute the rejected `EP2 Step 3 — Mechanism-Aware Baselines` technical plan from earlier this session.

**Rationale:** /review concluded smoke-runner variant expansion at N=12 with ±0.19 AUROC CI cannot distinguish lineage from mechanism via AUROC ranking. The technical plan built ~10 new files for a single binary question that AMRFinder + a per-prediction join answer directly.

### D2: Cut Bakta entirely from cef-PASS-substrate validation

**Decision:** No Bakta integration. AMRFinder's `main.tsv` calls β-lactamases via gene-presence — same signal as a Bakta gene-presence variant.

**Rationale:** Bakta wallclock 1-6 hr unverified; `--skip-cds` flag in the rejected plan would have produced empty annotations. No NEW resistance-mechanism information beyond AMRFinder.

### D3: Port `cipro_mechanism_audit.py` to `cef_mechanism_audit.py`; emit evidence object NOT verdict label

**Decision:** Create `scripts/cef_mechanism_audit.py` as a port of the existing cipro mechanism audit. Inline `CEF_LOCI_BY_MECHANISM` + `AMRFINDER_CEF_RELEVANT_CLASSES` parallel to the cipro version. The output is an **evidence object** (exact counts + Wilson 95% CI + posterior `Pr(p ≥ 0.70 | data)` under beta(1,1) prior + per-strain mechanism table + R/S × mechanism cross-tab), NOT a single `DOMINANT/MIXED/UNKNOWN` verdict label.

**Rationale (revised post-brainstorm Issue 2):** verdict bands `≥70% / 40-70% / <40%` are brittle at N=6 (one-strain flip changes the verdict; cipro got away with it because data landed at 90%). CI-based threshold rules are unreachable for n=6 (even 6/6 has Wilson CI lower bound 54.1%, well below 0.70). Honest output is the full evidence, NOT a derived label that hides the uncertainty.

**Trade-off:** Two copies of the audit pipeline (cipro + cef) until consolidation. Acceptable; the cipro version emits a verdict label (preserved for backward compat); the cef version emits an evidence object. Consolidation deferred to a separate refactor.

### D4: Fix smoke-runner silent-variant-drop bug + add JSON sidecar in ONE commit

**Decision:** Before any cef audit work, patch `scripts/smoke_gate_12strain_cipro.py` to (a) fix the exception-swallowing bug at `:374-382` AND (b) emit a per-strain predictions JSON sidecar (analogous to the Cipro Decision Bundle Technical Plan's Step 8 design). Then re-fire today's cef smoke to produce the sidecar.

**Rationale (revised post-brainstorm Issue 6):** The cross-tab diagnostic (Step 2a) requires per-strain NT-XGBoost predictions which don't exist in today's markdown artifact. The bug fix and the sidecar are both small surgical changes to the same file; bundling them is minimum-diff.

### D5: Attribution preflight is NON-OPTIONAL; adds per-prediction join

**Decision (CHANGED from prior version):** Step 4 (cef attribution preflight) always fires after Step 3. Plus a join layer: per-strain `(NT-XGB LOSO score, AMRFinder β-lactamase hit, mutagenesis prediction delta on β-lactamase loci)`. The join is what answers "were the LOSO PASS calls mechanism-driven?" — not just "is mechanism present in the cohort?"

**Rationale (revised post-brainstorm Issues 1 + 8):** The originally-proposed gate was inverted. `scripts/cipro_mechanism_audit.py:397` framed `QRDR_DOMINANT` as "If verdict is QRDR_DOMINANT: NT's preflight INCONCLUSIVE_MISS is a MODEL failure" — for interpreting model FAILURE on known mechanism. The cef question is opposite: was the cef PASS a model SUCCESS or lineage artifact? Attribution is MOST informative when mechanism dominates the cohort, not least. Plus the cipro preflight aggregates across R strains without joining to per-strain predictions; the cef version needs the join to be question-answering.

**Trade-off:** Step 4 always runs (5-10 min cached). Gating saved little and risked skipping the load-bearing test.

### D6: External publication deferred until BOTH cipro + cef have full audit chains

**Decision (unchanged):** Defer arXiv / blog until cipro EP1 closeout + cef mechanism audit + cef attribution audit all complete, AND the cross-drug architectural finding ("NT-frozen-pooling works on concentrated-signal mechanisms; fails on distributed mobile-element mechanisms") is solidly grounded.

### D7 (NEW): Cohort-provenance scope — PC1_cef gates BV-BRC rebuild

**Decision:** Conclusions from this plan are scoped to "within this reused mini-cohort" UNLESS PC1_cef = `publication_facing_cross_drug_claim`, which triggers a BV-BRC-wide cef cohort rebuild (separate plan; ~5-7 hr NT cache populate on GTX 860M + cohort build). All output packets carry the scope caveat in their headers.

**Rationale (post-brainstorm Issue 5):** The 12 cef strains were filtered from the N=38 cipro cohort, not built BV-BRC-wide. Mechanism distribution may reflect cipro-cohort selection artifacts, not cef biology. Cross-drug architectural claims require a cef-native cohort.

**Trade-off:** Internal-triage path is much cheaper; publication path multiplies the work 5-10×. Most likely PC1_cef = `internal_triage_only` per the EP1 closeout's discipline.

### D8 (NEW): Pre-conditions discipline (PC1_cef + PC2_cef)

**Decision:** Analogous to the cipro decision bundle's PC1/PC2 lock. User declares PC1_cef + PC2_cef in `wiki/cef_decision_bundle_pre_conditions_<date>.md` BEFORE any runtime step fires.

- **PC1_cef:** `internal_triage_only` vs `publication_facing_cross_drug_claim`. Determines whether D7's BV-BRC rebuild fires.
- **PC2_cef:** estimand list (NOT a single threshold). Required estimands:
  - Binary agreement at score=0.5 (cross-tab of `(true × NT-XGB-pred × AMRFinder-hit)`).
  - AUROC rank concordance (rank-correlation between NT-XGB LOSO scores and AMRFinder β-lactamase count per strain).
  - Discordant-case classification (NT-R-without-mechanism / NT-S-with-mechanism / etc.).
  - Mash-clade-only baseline AUROC for direct lineage-tracking-hypothesis test.

**Rationale (post-brainstorm Medium-issue M4):** Pre-locked pre-conditions prevent post-hoc rationalization at runtime decision points. Matches the cipro decision bundle's discipline.

## Implementation Plan

### Step 0: Pre-conditions artifact freeze (PC1_cef + PC2_cef)
Files: wiki/cef_decision_bundle_pre_conditions_<date>.md (new)
Depends on: none

**What changes:**
- User writes `wiki/cef_decision_bundle_pre_conditions_<date>.md` declaring PC1_cef + PC2_cef per D8.
- This is a user-action step (not a code step). Subsequent Steps 4 and 5 STRUCTURALLY refuse to write their output artifacts unless this file exists.
- Recommended default (per EP1 closeout discipline): PC1_cef = `internal_triage_only`; PC2_cef = the 4-estimand list from D8.

**Verification:** file exists; contains PC1_cef + PC2_cef sections; ledger row appended.

### Step 1: Fix smoke-runner silent-variant-drop bug + emit per-strain JSON sidecar
Files: scripts/smoke_gate_12strain_cipro.py, tests/test_smoke_gate_12strain.py (new)
Depends on: Step 0

**What changes:**
- `scripts/smoke_gate_12strain_cipro.py:374-382` — replace `try/except Exception` with named-exception handling:
  1. Re-raise on unexpected errors. Only swallow `ClassifierTrainingError`, `FileNotFoundError`, `IndeterminateOOVError`.
  2. Render `INDETERMINATE_<reason>` rows in the markdown packet instead of silently dropping.
  3. Stderr surface for the dropped variant + error type.
- Emit a per-strain predictions JSON sidecar at `wiki/smoke_gate_12strain_<drug>_<date>.predictions.json`. Schema:
  ```json
  {
    "drug": str,
    "cohort_path": str,
    "n_total": int,
    "per_strain": [
      {
        "strain_id": str, "mlst": str, "y_true": int,
        "nt_xgb_score": float | null,
        "nt_lr_score": float | null,
        "kmer_xgb_score": float | null,
        "gene_presence_xgb_score": float | null
      }
    ]
  }
  ```
- Sidecar emission is unconditional (every smoke run produces both .md and .json).

**Test strategy:**
- NEW file `tests/test_smoke_gate_12strain.py`:
  - `test_variant_with_classifiertrainingerror_renders_indeterminate` — assert INDETERMINATE row in markdown, not missing row.
  - `test_variant_with_filenotfounderror_renders_indeterminate_with_strain_id` — same with strain_id in error.
  - `test_unexpected_exception_propagates` — assert `main()` does NOT swallow `RuntimeError("unexpected")`.
  - `test_predictions_json_sidecar_emitted` — assert sidecar exists alongside markdown.
  - `test_predictions_json_schema_pins_required_fields` — assert per_strain entries have required keys.

### Step 2: Re-fire cef smoke under patched runner
Files: (runtime only — no source changes)
Depends on: Step 1

**What runs:**
```bash
HF_HOME=D:/hf_cache uv run python scripts/smoke_gate_12strain_cipro.py \
  --cohort data/processed/gate_b_mini_cef_cohort.parquet \
  --nt-cache D:/dna_decode_cache/embeddings/nt_n40_cipro.h5 \
  --refseq-cache D:/dna_decode_cache/refseq \
  --drug ceftriaxone
```
- Produces `wiki/smoke_gate_12strain_ceftriaxone_<date>.md` AND `wiki/smoke_gate_12strain_ceftriaxone_<date>.predictions.json`.
- Expected wallclock: ~9 min (per today's empirics).
- Both Steps 2a + 2b + 3 + 4 consume the JSON sidecar.

**Verification:** sidecar JSON has 12 per_strain entries with `nt_xgb_score` populated for 11 (one strain dropped per today's empirics) + `kmer_xgb_score` populated for all 12.

### Step 2a: Cef NT-vs-mechanism cross-tab diagnostic
Files: scripts/cef_nt_vs_mechanism_crosstab.py (new), tests/test_cef_nt_vs_mechanism_crosstab.py (new)
Depends on: Step 2

**What changes:**
- NEW file `scripts/cef_nt_vs_mechanism_crosstab.py`.
- Input: `--predictions-json` (Step 2's sidecar), `--amrfinder-output-root` (default `data/amrfinder_runs/` — will run AMRFinder via `tools.docker_runner` on cohort strains if not cached).
- Primary head: **NT-XGBoost** (the cef-PASS head). NT-LR included as **sensitivity row** for cross-head consistency check.
- Per-strain table: `(strain_id, true_R, nt_xgb_score, nt_xgb_pred_R, nt_lr_score, nt_lr_pred_R, amrfinder_beta_lactamase_hit, amrfinder_ampc_hit, amrfinder_porin_hit, amrfinder_efflux_hit, mash_clade_id)` — mash_clade_id populated by Step 2b.
- Cross-tab statistics:
  - Binary-agreement count at `score >= 0.5` threshold.
  - AUROC rank concordance (Spearman) between `nt_xgb_score` and `amrfinder_beta_lactamase_count`.
  - Discordant case classes (NT-R-without-any-AMRFinder-cef-mech, NT-S-with-AMRFinder-mech, true-R-with-no-mech, true-S-with-mech).
- Output: `wiki/cef_nt_vs_mechanism_crosstab_<date>.{md,json}`.
- AMRFinder invocation reuses `scripts/cipro_mechanism_audit.py`'s `_run_amrfinder` pattern (lifted as a helper; not full module refactor yet).

**Test strategy:**
- NEW file `tests/test_cef_nt_vs_mechanism_crosstab.py`:
  - `test_crosstab_uses_nt_xgboost_as_primary_head` — synthetic sidecar with both NT-XGB + NT-LR; assert primary cross-tab built on NT-XGB.
  - `test_crosstab_emits_rank_concordance` — assert Spearman computed.
  - `test_crosstab_classifies_discordant_cases` — synthetic edge cases (NT-R-no-mech, NT-S-with-mech, etc.); assert each is bucketed correctly.
  - `test_crosstab_requires_pre_conditions_artifact` — missing `--pre-conditions-md` → RuntimeError.

### Step 2b: Cef Mash-clade lineage baseline
Files: scripts/cef_mash_clade_baseline.py (new), tests/test_cef_mash_clade_baseline.py (new)
Depends on: Step 2

**What changes:**
- NEW file `scripts/cef_mash_clade_baseline.py`.
- Uses existing `dna_decode/eval/phylogeny.py:compute_mash_distances` (Mash-batched 2-call refactor validated 2026-05-15) on the 12 cef strain FASTAs.
- Cluster the 12 strains into Mash clades (use existing hierarchical clustering helper if present; else single-linkage at distance threshold ~0.05).
- Train a **clade-only XGBoost** under LOSO (features = one-hot clade ID; binary label = cef R/S).
- Reports clade-only AUROC. If clade-only AUROC ≥ 0.70, lineage-tracking hypothesis is corroborated.
- Annotates Step 2a's cross-tab table with `mash_clade_id` per strain.
- Output: `wiki/cef_mash_clade_baseline_<date>.{md,json}`.

**Test strategy:**
- NEW file `tests/test_cef_mash_clade_baseline.py`:
  - `test_mash_clade_clustering_returns_valid_assignment` — synthetic strain distances; assert clusters non-trivial.
  - `test_clade_only_xgboost_loso_auroc_computed` — synthetic clades + labels; assert AUROC reported.
  - `test_clade_id_annotated_per_strain` — assert each strain has a clade_id in the output.

### Step 3: Cef mechanism audit (evidence object output)
Files: scripts/cef_mechanism_audit.py (new), tests/test_cef_mechanism_audit.py (new)
Depends on: Step 2

**What changes:**
- NEW file `scripts/cef_mechanism_audit.py`. Port of `scripts/cipro_mechanism_audit.py` (~411 lines).
- Replace `CIPRO_LOCI_BY_MECHANISM` with `CEF_LOCI_BY_MECHANISM`:
  - `beta_lactamase_acquired = {blaCTX-M, blaSHV, blaTEM, blaOXA, blaCMY, blaKPC, blaNDM}`
  - `ampC_chromosomal = {ampC}`
  - `porin_loss = {ompC, ompF}`
  - `efflux = {acrA, acrB, tolC}`
- Replace AMR class filter: `AMRFINDER_CEF_RELEVANT_CLASSES = {BETA-LACTAM, CARBAPENEM, CEPHALOSPORIN}` (speculative until first run; Step 3 first-run will surface mismatch).
- **Output = evidence object, NOT verdict label** (per D3 + brainstorm Issue 2):
  ```json
  {
    "n_R": 6,
    "n_S": 6,
    "r_strains_with_primary_mech_count": <int>,
    "r_strains_with_primary_mech_fraction": <float>,
    "wilson_ci_95_lower": <float>,
    "wilson_ci_95_upper": <float>,
    "posterior_pr_p_geq_0.70_beta_1_1": <float>,
    "per_strain_mechanism_table": [...],
    "rs_x_mechanism_crosstab": {...},
    "descriptive_label": "high_evidence_for_dominance" | "moderate_evidence" | "low_evidence" | "absent",
    "cohort_provenance_caveat": "filtered from cipro N=38 cohort"
  }
  ```
- `descriptive_label` is descriptive only — NOT decision-bearing. Decision logic lives in Step 5.
- Default `--cohort data/processed/gate_b_mini_cef_cohort.parquet`. `--out-root data/amrfinder_runs/` (shared with cipro cache — different accession subdirs, no collision).
- Output: `wiki/cef_mechanism_audit_<date>.{md,json}`.

**Test strategy:**
- NEW file `tests/test_cef_mechanism_audit.py`:
  - `test_cef_loci_includes_beta_lactamase_families` — assert dict keys.
  - `test_classify_symbol_blactxm_variant_prefix_match` — `classify_symbol("blaCTX-M-15")` → `"beta_lactamase_acquired"`.
  - `test_evidence_object_emits_wilson_ci` — synthetic 5/6 R-with-mech; assert Wilson CI bounds correct ([41.8%, 99.0%]).
  - `test_evidence_object_emits_posterior_prob` — synthetic data; assert posterior `Pr(p≥0.70 | data)` computed under beta(1,1).
  - `test_no_decision_bearing_verdict_label_emitted` — assert output does NOT contain a single `DOMINANT/MIXED/UNKNOWN` field; ONLY descriptive_label.

### Step 4: Cef attribution preflight + per-prediction join (NON-OPTIONAL)
Files: scripts/cef_attribution_preflight.py (new), tests/test_cef_attribution_preflight.py (new)
Depends on: Step 2, Step 3

**What changes:**
- NEW file `scripts/cef_attribution_preflight.py`. Port of `scripts/cipro_attribution_preflight.py` (~250 lines).
- Replace cipro locus set with `CEF_LOCI_BY_MECHANISM`.
- **Per-prediction join (NEW vs cipro version, per brainstorm Issue 8):** for each cef cohort strain, emit:
  - `nt_xgb_loso_score` (from Step 2 sidecar)
  - `nt_xgb_pred_R_at_0.5` (binary)
  - `nt_xgb_correct` (whether the LOSO prediction matched the true label)
  - `top_k_attribution_loci` (top-K=20 genes by signed-positive mutagenesis delta)
  - `n_top_k_loci_in_beta_lactamase_set` (intersection size with `CEF_LOCI_BY_MECHANISM["beta_lactamase_acquired"]` + `ampC_chromosomal`)
  - `mutagenesis_delta_sum_on_beta_lactamase_loci`
- Cohort-level aggregates: % of correctly-predicted-R strains where top-K attribution intersects β-lactamase set; same for incorrectly-predicted strains. Comparison is the load-bearing diagnostic.
- Refuses to write output unless `wiki/cef_decision_bundle_pre_conditions_<date>.md` exists (Step 0).
- Output: `wiki/cef_attribution_preflight_<date>.{md,json}`.
- Cache-only; reuses `D:/dna_decode_cache/embeddings/nt_n40_cipro.h5`; no Docker.
- Wallclock: ~5-10 min.

**Test strategy:**
- NEW file `tests/test_cef_attribution_preflight.py`:
  - `test_per_prediction_join_emits_all_fields` — synthetic; assert each per-strain row has all 6 join fields.
  - `test_aggregates_compare_correct_vs_incorrect_predictions` — synthetic with mixed correct/incorrect; assert both aggregates computed.
  - `test_refuses_without_pre_conditions_artifact` — missing PC artifact → RuntimeError.

### Step 5: Combined interpretation packet
Files: scripts/cef_decision_packet.py (new), tests/test_cef_decision_packet.py (new), wiki/cef_decision_packet_<date>.md (output), wiki/cef_decision_packet_<date>.json (output)
Depends on: Step 2a, Step 2b, Step 3, Step 4

**What changes:**
- NEW file `scripts/cef_decision_packet.py`. Consumes:
  - `wiki/cef_nt_vs_mechanism_crosstab_<date>.json` (Step 2a)
  - `wiki/cef_mash_clade_baseline_<date>.json` (Step 2b)
  - `wiki/cef_mechanism_audit_<date>.json` (Step 3)
  - `wiki/cef_attribution_preflight_<date>.json` (Step 4)
  - `wiki/cef_decision_bundle_pre_conditions_<date>.md` (Step 0)
- Synthesizes into a single structured decision artifact + markdown narrative.
- JSON schema:
  ```json
  {
    "pc1_cef_framing": str,
    "pc2_cef_estimands": [...],
    "cross_tab_summary": {...},
    "mash_clade_only_auroc": float,
    "mechanism_audit_evidence_object": {...},
    "attribution_per_prediction_aggregates": {...},
    "lineage_vs_mechanism_finding": "mechanism_driven_strong" | "mechanism_driven_weak" | "lineage_explained" | "indeterminate",
    "rationale": str,
    "recommended_next_step": str
  }
  ```
- `lineage_vs_mechanism_finding` derivation rule (pre-declared, NOT post-hoc):
  - `mechanism_driven_strong`: AMRFinder β-lactamase agreement ≥75% AND Mash-clade-only AUROC < 0.70 AND ≥75% of correct-R NT predictions have top-K attribution in β-lactamase set.
  - `mechanism_driven_weak`: ≥2-of-3 above, but not all 3.
  - `lineage_explained`: Mash-clade-only AUROC ≥ 0.75.
  - `indeterminate`: none of the above.
- Refuses to write output unless `wiki/cef_decision_bundle_pre_conditions_<date>.md` exists.

**Test strategy:**
- NEW file `tests/test_cef_decision_packet.py`:
  - 4 tests for the 4 `lineage_vs_mechanism_finding` cells.
  - `test_refuses_without_pre_conditions_artifact` — missing PC artifact → RuntimeError.

### Step 6 (deferred): BV-BRC-wide cef cohort rebuild
ONLY fires if PC1_cef = `publication_facing_cross_drug_claim`. NOT in this plan's scope. Separate plan if triggered.

## Execution Preview

```
Wave 0 (1 sequential): Step 0 — user-locked PC1_cef + PC2_cef artifact
Wave 1 (1 sequential): Step 1 — smoke-runner bug fix + JSON sidecar
Wave 2 (1 sequential, manual runtime): Step 2 — re-fire cef smoke
Wave 3 (2 parallel): Step 2a (cross-tab) + Step 2b (Mash-clade baseline) + Step 3 (mechanism audit)
                     Actually 3 parallel — all consume Step 2's sidecar, no file overlap
Wave 4 (1 sequential): Step 4 — attribution preflight (consumes Step 2 + Step 3 outputs)
Wave 5 (1 sequential): Step 5 — combined interpretation packet

Critical path: Step 0 → Step 1 → Step 2 → Step 3 → Step 4 → Step 5 (6 waves)
Max parallelism: 3 agents (Wave 3)
```

## Verification (end-to-end)

After all Steps 0-5 land:

- `uv run pytest tests/ -v` passes with new tests in Steps 1, 2a, 2b, 3, 4, 5 added cleanly with 0 regressions.
- `wiki/cef_decision_bundle_pre_conditions_<date>.md` exists with PC1_cef + PC2_cef declarations.
- All 5 wiki artifacts (cef_nt_vs_mechanism_crosstab + cef_mash_clade_baseline + cef_mechanism_audit + cef_attribution_preflight + cef_decision_packet) exist with the required schemas.
- `cef_decision_packet_<date>.json` has `lineage_vs_mechanism_finding` ∈ {mechanism_driven_strong, mechanism_driven_weak, lineage_explained, indeterminate} per the pre-declared derivation rule.
- The user can answer "is cef NT-XGB 0.833 mechanism-driven or lineage?" with grounded yes/no.

**Decision gate:** if `lineage_vs_mechanism_finding == "mechanism_driven_strong"`, commit to writing the EP1 + EP2 + cef-attribution cross-drug architectural-finding narrative (still NOT external publication per D6 + PC1_cef = internal_triage_only). If `lineage_explained`, the cef PASS becomes weaker evidence — re-frame to "NT-pooling at N=12 lineage-tracks even for concentrated-signal mechanisms; need bigger cohort for real mechanism claim." If `indeterminate`, BV-BRC-wide cef rebuild becomes the next experiment (per D7).

## Out of scope

- Bakta integration (cut per D2).
- 5-variant smoke-runner expansion (cut by /review).
- Consolidation of `CIPRO_LOCI_BY_MECHANISM` duplicates (deferred per D3).
- Rename `smoke_gate_12strain_cipro.py` (pre-existing debt; not blocking).
- Tet mechanism audit (H17 falsified; architectural-bottleneck verdict already established).
- BV-BRC-wide cef cohort rebuild (only fires if PC1_cef = publication_facing).
- External publication (deferred per D6).
