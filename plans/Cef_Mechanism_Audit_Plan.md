# Cef Mechanism Audit Plan

> Smoke-runner bug fix + cheap cef NT-vs-mechanism cross-tab diagnostic. Reduced 2026-05-17 (round-3 framing critique) from 8 steps to 3. The cef audit is an **appendix-tier corroborating diagnostic**, NOT a decision-bearing framework — the cross-drug architectural finding is already captured in EP1+EP2 verdicts. Real next deliverable is `plans/EP1_EP2_Cross_Drug_Synthesis_Plan.md`.

---

## Problem Statement

The EP2 cef PASS (NT-XGBoost 0.833 = k-mer 0.833) raises the question "is this mechanism-driven or lineage-tracking?" But /brainstorm round 3 surfaced a load-bearing framing issue: **the EP2 verdict packet (`wiki/EP2_cef_tet_verdict_2026-05-17.md`) already states the cross-drug architectural pattern**:

> "NT-frozen-pooling architecture works on concentrated-signal mechanisms (cipro QRDR + cef plasmid β-lactamases) and fails on distributed mobile-element mechanisms (tet)."

That cross-drug claim is the load-bearing scientific finding from EP1+EP2. The cef mechanism audit, even at its most thorough, is a CORROBORATING-NOT-DECISION-BEARING diagnostic at N=12 with ±0.19 AUROC CI. The plan revision history is telling:

- v1 (rejected): 11-step technical plan with Bakta + 4 new modules → /review verdict: SCOPE REDUCTION
- v2 (post-/review): 3-step plan (bug fix + mechanism audit + optional attribution)
- v3 (post-round-2): 8-step plan (PC artifacts + cross-tab + Mash-clade + evidence object + per-prediction join + decision packet)
- **v4 (post-round-3 framing critique, this revision):** 3 steps. Treats cef audit as appendix; primary deliverable is the EP1+EP2 cross-drug synthesis (separate plan).

The round-3 critique identified scope re-inflation (v3's "smaller alternative" language was preserving /review's authority while authorizing the same expansion pattern), threshold substitution (75%/0.70/75% thresholds at N=6 are as brittle as the verdict bands they replaced), and Mash-clade-at-N=12 contradicting prior project lessons about singleton-clade degeneracy. The plan now strips back to its load-bearing core.

## Design Decisions

### D1: Cef audit is appendix-tier, NOT decision-bearing

**Decision (round-3 reframing):** The cef NT-vs-mechanism cross-tab (Step 3 below) is a corroborating diagnostic that may inform the EP1+EP2 synthesis writeup's residual-uncertainty section. It is NOT a yes/no adjudication on the architectural claim. The cross-drug story is already grounded in the 3-experiment chain (cipro EP1 closeout + cef PASS + tet FAIL); the cef cross-tab CORROBORATES it or it doesn't, but neither outcome changes the synthesis verdict.

**Rationale:** N=12 with ±0.19 AUROC CI cannot adjudicate mechanism vs lineage with statistical confidence. The EP2 verdict has already claimed the architectural pattern. Round-3 framing critique: "the plan asks 'is cef 0.833 mechanism-driven?' but the project decision now is 'is the EP1+EP2 packet already writeup-worthy?' — those are different questions."

### D2: Smoke-runner bug fix is the only load-bearing step

**Decision:** Step 1 (smoke-runner silent-variant-drop bug + per-strain JSON sidecar) is the only step that has decision-bearing impact beyond cef. It fixes a regression in today's cef + tet reports (both missing the gene-presence variant row) AND adds a reusable JSON sidecar that the Cipro Decision Bundle Technical Plan's Step 8 already specifies. Both serve future EPs.

### D3: Drop framework discipline disproportionate to N=12

**Decision (round-3 reductions):** Drop PC1_cef + PC2_cef pre-conditions discipline, Mash-clade lineage baseline (Step 2b), mechanism audit evidence object schema (Wilson CI + posterior), attribution preflight + per-prediction join, and the 4-cell `lineage_vs_mechanism_finding` decision matrix. **All deferred** to a hypothetical future N≥150 cef cohort if such ever fires.

**Rationale:** The 2026-05-15 project state already records that mini-cohorts with 12 unique MLSTs make clade-only baselines degenerate. Wilson CI + posterior machinery is "rigor theater" when no downstream consumer uses the uncertainty fields. PC artifact + 4-cell decision matrix are disproportionate to the cohort scale.

### D4: External publication still deferred per EP1 closeout

**Decision (unchanged):** Defer arXiv / blog until BOTH cipro + cef + tet have completed audit chains AND the cross-drug architectural finding is solidly grounded. The synthesis writeup (separate plan) is INTERNAL closeout, not external publication.

## Implementation Plan

### Step 1: Fix smoke-runner silent-variant-drop bug + emit per-strain JSON sidecar
Files: scripts/smoke_gate_12strain_cipro.py, tests/test_smoke_gate_12strain.py (new)
Depends on: none

**What changes:**
- `scripts/smoke_gate_12strain_cipro.py:374-382` — replace `try/except Exception` with named-exception handling:
  1. Re-raise on unexpected errors. Only swallow `ClassifierTrainingError`, `FileNotFoundError`, `IndeterminateOOVError`.
  2. Render `INDETERMINATE_<reason>` rows in the markdown packet instead of silently dropping.
  3. Stderr surface for dropped variant + error type.
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
- Sidecar emission unconditional (every smoke run produces both .md and .json).

**Test strategy:**
- NEW file `tests/test_smoke_gate_12strain.py`:
  - `test_variant_with_classifiertrainingerror_renders_indeterminate`
  - `test_variant_with_filenotfounderror_renders_indeterminate_with_strain_id`
  - `test_unexpected_exception_propagates`
  - `test_predictions_json_sidecar_emitted`
  - `test_predictions_json_schema_pins_required_fields`

### Step 2: Re-fire cef smoke under patched runner
Files: (runtime only)
Depends on: Step 1

**What runs:**
```bash
HF_HOME=D:/hf_cache uv run python scripts/smoke_gate_12strain_cipro.py \
  --cohort data/processed/gate_b_mini_cef_cohort.parquet \
  --nt-cache D:/dna_decode_cache/embeddings/nt_n40_cipro.h5 \
  --refseq-cache D:/dna_decode_cache/refseq \
  --drug ceftriaxone
```

- Produces `wiki/smoke_gate_12strain_ceftriaxone_<date>.md` + `.predictions.json`.
- Expected wallclock: ~9 min.

**Verification:** sidecar JSON has 12 per_strain entries; gene-presence row is now visible (either AUROC value or `INDETERMINATE_MISSING_GFF3`).

### Step 3: Cef NT-vs-mechanism cross-tab diagnostic (appendix-tier)
Files: scripts/cef_nt_vs_mechanism_crosstab.py (new), tests/test_cef_nt_vs_mechanism_crosstab.py (new)
Depends on: Step 2

**What changes:**
- NEW file `scripts/cef_nt_vs_mechanism_crosstab.py`.
- Input: `--predictions-json` (Step 2's sidecar), `--amrfinder-output-root` (default `data/amrfinder_runs/`; will run AMRFinder via `tools.docker_runner` on cohort strains if not cached, mirroring `scripts/cipro_mechanism_audit.py:_run_amrfinder`).
- Per-strain table: `(strain_id, true_R, nt_xgb_score, nt_xgb_pred_R, nt_lr_score, nt_lr_pred_R, amrfinder_beta_lactamase_hit, amrfinder_ampc_hit, amrfinder_porin_hit)`.
- Primary head: **NT-XGBoost** (the cef-PASS head). NT-LR included only as sensitivity row.
- Output statistics:
  - Binary agreement count at `score >= 0.5`.
  - Spearman rank correlation between `nt_xgb_score` and `amrfinder_beta_lactamase_count`.
  - Discordant case classes (NT-R-no-mech / NT-S-with-mech / true-R-with-no-mech / true-S-with-mech).
- Output: `wiki/cef_nt_vs_mechanism_crosstab_<date>.{md,json}`.
- **Header explicitly notes appendix-tier scope: N=12, ±0.19 CI, no decision-bearing claim.**

**Test strategy:**
- NEW file `tests/test_cef_nt_vs_mechanism_crosstab.py`:
  - `test_crosstab_uses_nt_xgboost_as_primary_head`
  - `test_crosstab_emits_rank_concordance`
  - `test_crosstab_classifies_discordant_cases`

**Interpretation (no decision rule — narrative only):**
- High NT-XGBoost agreement with AMRFinder β-lactamase calls + high rank concordance → CORROBORATES the EP2 architectural claim. Synthesis writeup may cite as supporting evidence.
- Low agreement → introduces residual uncertainty. Synthesis writeup notes "cef mechanism question requires a larger cef-native cohort to adjudicate" in the residual-uncertainty section.
- Either way, the synthesis writeup is the primary deliverable, NOT this cross-tab.

## Execution Preview

```
Wave 0 (1 sequential): Step 1 — smoke-runner bug fix + JSON sidecar
Wave 1 (1 sequential, runtime): Step 2 — re-fire cef smoke
Wave 2 (1 sequential): Step 3 — cef NT-vs-mechanism cross-tab

Critical path: Step 1 → Step 2 → Step 3 (3 waves)
Max parallelism: 1 agent (each step depends on the prior)
```

This plan is now genuinely smaller than the rejected 11-step plan, in spirit AND structure.

## Verification (end-to-end)

After Steps 1 + 2 + 3 land:

- `uv run pytest tests/ -v` passes with the new tests in Steps 1 + 3.
- Re-fired cef smoke packet shows 3 variant rows (NT-XGBoost + k-mer + gene-presence-or-INDETERMINATE) — variant-drop bug is fixed.
- `wiki/cef_nt_vs_mechanism_crosstab_<date>.md` has the 12-strain per-strain table + the 3 cohort statistics + the explicit "appendix-tier, no decision-bearing claim" header.
- The cef cross-tab artifact is then available to cite (or NOT cite) in the EP1+EP2 synthesis writeup.

**The synthesis writeup is the primary deliverable.** See `plans/EP1_EP2_Cross_Drug_Synthesis_Plan.md`.

## Out of scope (dropped post-round-3)

- Bakta integration (cut per v2 D2).
- 5-variant smoke-runner expansion (cut by /review v1).
- PC1_cef + PC2_cef pre-conditions discipline (disproportionate to N=12; was v3 Step 0).
- Mash-clade lineage baseline (degenerate at N=12 per project state; was v3 Step 2b).
- Mechanism audit evidence object with Wilson CI + posterior (rigor theater; was v3 Step 3).
- Attribution preflight + per-prediction join (was v3 Step 4; defer to N≥150 cohort if ever).
- 4-cell `lineage_vs_mechanism_finding` decision matrix (was v3 Step 5).
- BV-BRC-wide cef cohort rebuild (was v3 D7 conditional; deferred unconditionally now).
- External publication (deferred per D4).
- Tet mechanism audit (H17 falsified; architectural-bottleneck verdict already established).
