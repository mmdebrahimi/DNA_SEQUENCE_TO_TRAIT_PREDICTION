# Stage 1 Refactor & Test Hardening — Technical Plan

> Convert the /review synthesis into a 3-step refactor: pre-commit decision rules in the Stage 1 plan, reduce `scripts/stage1_n40_cipro.py` to thin orchestration over existing infrastructure, and pin the two critical untested behaviors (fusion-exclusion from gate + `calibrate=False` discipline).

**Revision history:**
- 2026-05-14 — initial save (post-/review)
- 2026-05-14 PM — incorporated /brainstorm Round 1+2 corrections (4 critical + 3 medium issues). Specifically: shared k-mer API made order-explicit; verdict semantics split from `stage2_action`; `CONTIG_SEPARATOR` made a module-level constant; `_train_baseline_logreg` calibration-skip branch spelled out; `CVResult.strain_ids` property added; alignment validation in `compute_gate_outcome`; 3rd calibration test added (k-mer-XGB call site); smoke gate regression tolerance loosened to ±0.005 AUROC.

---

## Problem Statement

The /review (2026-05-14) of `plans/Stage1_N40_Cipro_Engineering_Screen_Plan.md` surfaced three accepted concerns about Stage 1 readiness:

1. **Decision discipline gap (product):** Stage 2 actions per verdict bucket are soft ("review diagnostics," "options include"). At the moment Stage 1 results land, motivated reasoning is at maximum. Pre-commit deterministic Stage 2 actions per bucket now.
2. **Implementation duplication (engineering):** `scripts/stage1_n40_cipro.py` reimplements `leave_one_strain_out_cv` (from `dna_decode/eval/cv.py:68-110`), duplicates ~80% of `run_kmer_xgboost` from the smoke gate, writes a fresh `_train_logreg` that omits the existing `_train_baseline_logreg`'s minority-class safety net, and silently falls back to `train_y.mean()` on `ClassifierTrainingError` — which masks the exact bug modes Stage 1 is supposed to surface.
3. **Test gap (engineering):** The 11 existing tests cover helpers but miss the two gate-bearing behaviors that the /brainstorm specifically flagged — fusion correctly excluded from the gate, and `calibrate=False` actually passed at each call site (regression guard for the N=11 isotonic-collapse bug).

Defensible defaults for the two clarifying questions: **Stage 2 burst is treated as atomic** (one-shot commit, no Stage 1b middle path) so NOISY PASS must hard-resolve via the CI lower bound; **duplication with the smoke gate was accidental** under time pressure, refactor is welcome. User can override either default by editing Step 1's text.

## Codebase Context

Files researched at C:/Users/Farshad/PythonProjects/dna_decode:

- `dna_decode/eval/cv.py:68-110` — `leave_one_strain_out_cv(features, labels, strain_ids, train_fn, predict_fn, drug="")` returns a `CVResult` with `all_y_true` / `all_y_score` concatenated in strain order. Single-class training folds emit `np.nan` in `y_score`. Stage 1 LOSO for NT variants reuses this directly (fixed feature matrix).
- `dna_decode/eval/metrics.py:43-94` — `compute_metrics` uses `_nan_safe_filter` to drop NaN-scored entries cleanly. The NaN-emit pattern from `leave_one_strain_out_cv` round-trips correctly through `compute_metrics`.
- `dna_decode/models/classical_baselines.py:133-181` — `_train_baseline_logreg(X, y, drug_name)` provides logreg training with `MIN_TRAINING_SAMPLES` guard + minority-class-bounded calibration; current Stage 1 `_train_logreg` ignores this and omits the safety net.
- `dna_decode/models/classical_baselines.py:230` — `contig_separator: str = "N" * 100` is the canonical separator default; Stage 1 currently re-creates `sep = "N" * 100` at multiple sites.
- `scripts/smoke_gate_12strain_cipro.py:116-162` — `run_kmer_xgboost(cohort, refseq_root, drug, k=8, top_n=10000)` is the canonical within-fold-vocab-rebuild k-mer LOSO. Stage 1's `loso_kmer_xgb` is a near-duplicate with cosmetic differences. Factor into a shared module so both runners import from one source.
- `scripts/smoke_gate_12strain_cipro.py:67-113` — `run_nt_xgboost(cohort, cache_path, drug)` is the canonical NT-XGBoost LOSO. Reuses `leave_one_strain_out_cv` from `cv.py`. Stage 1 should match.
- `tests/test_stage1_n40_cipro.py` — 11 tests for `verdict_label`, `paired_bootstrap_ci`, `per_mlst_breakdown`. Imports from `scripts.stage1_n40_cipro`. Confirmed 2026-05-14 passing.
- `dna_decode/models/cache.py` — `EmbeddingCache.bulk_get` already exists; Stage 1 uses it correctly.
- `dna_decode/data/annotations.py` — `parse_gff3` is imported by `scripts/stage1_n40_cipro.py:37` but never called. Dead import.
- `wiki/decisions-log.md` HIGH-salience 2026-05-14: "AUROC≈0.000 with symmetric scores = calibration/label/row-order/class-mapping bug signature" → reinforces Step 3's calibration-discipline regression test.

## Design Decisions

### D1: Treat Stage 2 burst as atomic; pre-commit deterministic `stage2_action` separately from verdict label (per /brainstorm Round 1)

**Decision:** The Stage 2 N=150 Databricks burst is committed in one shot or not at all — no Stage 1b middle path. The CI-lower-bound rule drives a NEW field `stage2_action` (literal values: `BURST_STAGE_2` / `HOLD_STAGE_2_CI_DEGENERATE` / `ALTERNATIVE_POOLING_RERUN` / `PIVOT_TO_BAKTA`), NOT the verdict label. Verdict semantics stay frozen as a pure function of point gap (CLEAN ≥5 pp, NOISY 3-5 pp, FAIL <3 pp). Result packet renders both.

**Rationale:** The /review's dominant risk was "motivated reasoning at NOISY PASS" — soft language at the moment of truth lets the verdict slide. Pre-commit removes the slide. The /brainstorm Round 1 critique flagged that the original phrasing "convert NOISY → FAIL if CI ≤ 0" changes verdict semantics, conflicting with the stated non-goal. Splitting verdict from action preserves both: verdict reports point-gap result honestly; action codifies the budget decision.

**Trade-off:** Loses optionality for a "Stage 1b tightening run before committing the rest" path. Acceptable because the user has stated a preference for atomic decisions over staged ambiguity, and Stage 2 itself is the real ship gate (Option-C threshold + biology check). The two-field design (verdict + action) adds one packet line but preserves audit trail of the underlying point gap.

### D2: Refactor `scripts/stage1_n40_cipro.py` to reuse existing infrastructure rather than re-implement; order-explicit shared API (per /brainstorm Round 1)

**Decision:** Stage 1's runner becomes thin orchestration over existing modules: `leave_one_strain_out_cv` for NT variants, factored-out `dna_decode/eval/loso_kmer.py` for k-mer + fusion, `_train_baseline_logreg(..., calibrate=False)` for logreg path. The shared k-mer API uses ORDER-EXPLICIT parameters (`seqs_by_strain` dict + `labels_by_strain` dict + explicit `strain_ids` list) — NOT the original `cohort` + `refseq_root` signature which would have silently re-ordered or re-subsetted. `_train_baseline_logreg` gets an explicit calibration-skip branch (not just a kwarg). Adds `@property strain_ids` to `CVResult` so all LOSO producers share one alignment contract. Eliminates ~200 LOC of duplicated logic.

**Rationale:** The /review found the runner reimplements `leave_one_strain_out_cv` from `cv.py:68-110`, duplicates `run_kmer_xgboost` from the smoke gate, and writes a fresh `_train_logreg` that omits the existing minority-class safety net. /brainstorm Round 1 caught that a naive factoring with `cohort` parameter would let Stage 1 evaluate k-mer on a different strain subset/order than NT (Stage 1 skips strains without NT cache; smoke runner sorts strain IDs from all drug-labeled FASTAs). Order-explicit API prevents that silent divergence. /brainstorm Round 2 caught that `CVResult` lacks `strain_ids` metadata at the top level — adding the property gives every LOSO producer the same alignment contract for free.

**Trade-off:** Adds one new module file (`dna_decode/eval/loso_kmer.py`), a `calibrate: bool` kwarg + explicit branch on `_train_baseline_logreg`, a `CONTIG_SEPARATOR` module constant, and a one-line `CVResult.strain_ids` property. All additions are minimal; the smoke runner also benefits from the factored module.

### D3: Replace silent mean-fallback on `ClassifierTrainingError` with re-raise

**Decision:** When `train_xgboost_classifier` raises `ClassifierTrainingError` in the k-mer / fusion / NT-XGBoost paths, re-raise rather than fall back to `train_y.mean()`. The NT LOSO variants get the same discipline implicitly by switching to `leave_one_strain_out_cv` (which emits `np.nan` on single-class folds, handled by `compute_metrics._nan_safe_filter`).

**Rationale:** At N=38 with balanced 19R/19S, the exception only fires on real bugs (xgboost import failure, NaN inputs). Silent fallback to the training-set mean produces a plausible-looking AUROC and hides the failure Stage 1 is supposed to catch.

**Trade-off:** Stage 1 now fails loudly on classifier-training bugs instead of degrading gracefully. Correct behavior for an engineering screen; rejected the soft-failure alternative.

### D4: Add regression tests for fusion-exclusion + `calibrate=False` discipline at ALL gate-bearing call sites (per /brainstorm Round 2 expansion)

**Decision:** Step 3 adds dedicated tests pinning (a) fusion correctly excluded from the gate (synthetic test where fusion AUROC=0.99 while all gate-bearing variants are 0.50 → asserted FAIL verdict), (b) `calibrate=False` actually passed at EACH gate-bearing call site (NT-XGBoost + NT-logreg + k-mer-XGB — three tests, not two; the /brainstorm Round 2 critique caught that the k-mer-XGB call site is also gate-bearing and was missed in the original plan), (c) `stage2_action` decision-layer mapping (4 tests pinning the bucket × CI-lo decision table), (d) strain_ids alignment validation (raise on NT-vs-k-mer mismatch; suppress fusion note on fusion mismatch).

**Rationale:** The /brainstorm flagged fusion-passing-while-NT-fails as the "k-mer carrying the result" failure mode; the calibration discipline is one careless edit away from re-triggering the N=11 isotonic-collapse bug. Both behaviors are gate-bearing glue currently untested. The Round 2 critique caught that the k-mer-XGB call site is structurally identical and would have remained an untested regression surface. The `stage2_action` tests pin Step 1's plan-file rules as code-tested decision logic rather than prose. Alignment tests pin the Round 1 order-explicit-API contract.

**Trade-off:** Adds 13 tests total (vs 6 in the original plan) + 1 `monkeypatch` import; small implementation cost. The test count expanded specifically because /brainstorm Round 2 caught two missing coverage areas (k-mer calibration call site + alignment validation). Considered combining with the existing 11 tests under loose naming — rejected for clarity (these tests are regression guards, not unit-helpers).

### D5: Loud MLST handling instead of silent "unknown" fallback

**Decision:** Replace `getattr(s, "mlst", "unknown")` in `load_features` with `getattr(s, "mlst", None)` + `raise ValueError(...)` on `None`. Cohorts with missing MLST must be fixed upstream before Stage 1 runs.

**Rationale:** The MLST diagnostic appendix is the load-bearing lineage check. Silent fallback to `"unknown"` collapses missing-MLST strains into a same-string bucket, falsely inflating `largest_mlst_group` and degrading the lineage signal.

**Trade-off:** Strictness could break edge cohorts with incomplete metadata. Verified: 12-strain + N=40 current cohorts all have MLST. Edge case is theoretical; explicit error message tells the user how to fix.

### D6: Bootstrap-skip-count surfaces in the result packet

**Decision:** `paired_bootstrap_ci` returns 4 elements `(mean, lo, hi, n_effective)` where `n_effective` is the count of non-degenerate resamples. Result packet's bootstrap line reads `"B=1000 (effective N_eff)"`; if `n_effective < 800`, append "CI honesty degraded; investigate cohort imbalance."

**Rationale:** Silent degeneracy in bootstrap resamples quietly degrades CI honesty. Surfacing the effective count keeps the diagnostic transparent.

**Trade-off:** Changes the function signature from 3-tuple to 4-tuple; the 4 existing bootstrap tests need single-character updates. Trivial.

## Implementation Plan

### Step 1: Add verdict-time pre-commitments to the saved Stage 1 plan (verdict semantics UNCHANGED)
Files: plans/Stage1_N40_Cipro_Engineering_Screen_Plan.md
Depends on: none

**What changes:**
- `plans/Stage1_N40_Cipro_Engineering_Screen_Plan.md` — append a new top-level section `## Verdict-Time Pre-Commitments (locked 2026-05-14)` immediately before `## Verification`. **The section is a decision-layer ON TOP OF the unchanged verdict semantics — it never redefines the 3-bucket verdict labels.** Verdict is purely a function of point gap (≥5 pp CLEAN PASS, 3-5 pp NOISY PASS, <3 pp FAIL) per the locked Stage 1 plan; the pre-commitments determine `stage2_action` based on (verdict + CI lower bound + fusion behavior).

  Four deterministic rules:
  1. **CI-lower-bound budget-decision rule** (drives `stage2_action`, NOT the verdict label).
     - If CI lower bound on `(NT-best − k-mer-XGB)` is > 0 → `stage2_action` follows the verdict bucket directly.
     - If CI lower bound ≤ 0 AND verdict is NOISY PASS (gap ∈ [3, 5) pp) → `stage2_action = HOLD_STAGE_2_CI_DEGENERATE` (effective FAIL handling for budget purposes; verdict label STAYS "NOISY PASS").
     - If CI lower bound ≤ 0 AND verdict is CLEAN PASS (gap ≥ 5 pp) → `stage2_action = BURST_STAGE_2` (proceed) but annotate the packet "wide CI; lineage check load-bearing." Verdict label STAYS "CLEAN PASS."
  2. **Asymmetric-cost statement.** "We are more willing to make a false-FAIL than a false-PASS at Stage 1 because false-PASS only costs Stage 2 budget (atomic, recoverable via Stage 2 verdict) while false-FAIL silently kills a potentially useful track. The ≥3 pp threshold reflects this trade."
  3. **Fusion-outperforms-primary rule.** If `fusion_AUROC − max(NT-XGBoost, NT-logreg) ≥ 3 pp` AND fusion alignment is valid (see Step 2.6), log it in the result packet as "Stage 2 architecture note: fusion outperformed both NT-only heads, revisit at Stage 2." If fusion alignment is invalid (subset/order mismatch with NT-best), suppress the note entirely with "fusion alignment mismatch — diagnostic suppressed." Fusion NEVER alters the gate-bearing gap computation or the verdict label.
  4. **Per-bucket `stage2_action` mapping** (deterministic, no "review and decide"; literal values: `BURST_STAGE_2` / `HOLD_STAGE_2_CI_DEGENERATE` / `ALTERNATIVE_POOLING_RERUN` / `PIVOT_TO_BAKTA`):
     - **CLEAN PASS + CI lo > 0** → `BURST_STAGE_2`. Proceed to Stage 2 Databricks burst with N=150 cohort build. No further deliberation.
     - **CLEAN PASS + CI lo ≤ 0** → `BURST_STAGE_2` with packet annotation "wide CI."
     - **NOISY PASS + CI lo > 0** → `BURST_STAGE_2`. Annotate "diagnostics checked, CI clears 0."
     - **NOISY PASS + CI lo ≤ 0** → `HOLD_STAGE_2_CI_DEGENERATE`. Do NOT spend Stage 2 burst budget. Next: `ALTERNATIVE_POOLING_RERUN` (Stage 1b with `mean+max` aggregation instead of `mean`); if still NOISY+degenerate, `PIVOT_TO_BAKTA`.
     - **FAIL** → `ALTERNATIVE_POOLING_RERUN`. Run Stage 1 once with `mean+max` aggregation; if still <3 pp, `PIVOT_TO_BAKTA` (re-annotation + gene-presence comparator pathway per `plans/Gene_Presence_AUROC_Bug_Fix_Plan.md`). Do NOT spend Stage 2 burst budget on NT-only.

**Key details:**
- Verdict labels (CLEAN PASS / NOISY PASS / FAIL) remain pure functions of point gap — Step 2.5's `compute_gate_outcome` returns `verdict_bucket` (point-gap function) AND `stage2_action` (the decision-layer output of the rules above). The result packet renders BOTH.
- No code changes in Step 1 itself; pure plan-file edit. Step 2.5 implements the `stage2_action` computation.
- All defaults are defensible; user can edit any rule before Stage 1 executes.

**Test strategy:**
- None for the plan edit. Step 3 adds tests that pin `stage2_action` semantics.
- Spot-check the file renders correctly and the new section sits between `## Design Decisions` (D4) and `## Verification`.

---

### Step 2: Refactor `scripts/stage1_n40_cipro.py` to reuse existing infrastructure + eliminate silent failure modes
Files: scripts/stage1_n40_cipro.py, scripts/smoke_gate_12strain_cipro.py, dna_decode/eval/loso_kmer.py, dna_decode/eval/cv.py, dna_decode/models/classical_baselines.py
Depends on: none

**What changes:**

1. **Factor out the k-mer LOSO into a shared module with order-explicit API.**
   - `dna_decode/eval/loso_kmer.py` (NEW) — extract two pure functions with NO `cohort` / `refseq_root` parameters (caller supplies strain order explicitly):
     ```python
     def run_kmer_xgboost_loso(
         seqs_by_strain: dict[str, str],
         labels_by_strain: dict[str, int],
         strain_ids: list[str],
         drug: str,
         k: int = 8,
         top_n: int = 10_000,
     ) -> CVResult: ...

     def run_fusion_loso(
         X_nt: np.ndarray,                    # (len(strain_ids), nt_dim)
         seqs_by_strain: dict[str, str],
         labels_by_strain: dict[str, int],
         strain_ids: list[str],
         drug: str,
         k: int = 8,
         top_n: int = 10_000,
     ) -> CVResult: ...
     ```
     Both functions return a `CVResult` (consistent with `leave_one_strain_out_cv`). Fold order = input `strain_ids` order — NO internal sorting. Within-fold vocab rebuild using training-set sequences only. Use `CONTIG_SEPARATOR` from `classical_baselines` (added in sub-task 8).
   - `scripts/smoke_gate_12strain_cipro.py` — replace `run_kmer_xgboost` body with a wrapper that builds `seqs_by_strain` / `labels_by_strain` / `strain_ids = sorted(strain_contigs.keys())` from cohort.strains, then calls `dna_decode.eval.loso_kmer.run_kmer_xgboost_loso(...)`. Smoke's broader subset behavior is preserved (sorted strain IDs from all drug-labeled FASTAs).
   - Stage 1 will call the same factored functions with its NT-cache-present `strain_ids` (the alignment-preserving subset).

2. **Reuse `leave_one_strain_out_cv` for NT variants.**
   - `scripts/stage1_n40_cipro.py` — delete the local `loso_scores` function. Replace its two call sites (NT-XGBoost, NT-logreg) with `leave_one_strain_out_cv(X_nt, y, strain_ids, train_fn, predict_fn, drug="cipro_stage1")` calls. Use `result.all_y_true` and `result.all_y_score` to feed into `compute_metrics` and the per-strain table. Single-class training folds will emit `np.nan` per `cv.py:83-95`; `compute_metrics` strips them via `_nan_safe_filter`. This eliminates the silent mean-fallback.
   - **Strain-order extraction:** use the new `CVResult.strain_ids` property (added in sub-task 9) rather than reconstructing from `[f.held_out_id for f in result.folds]` at each call site.

3. **Replace `_train_logreg` with `_train_baseline_logreg(calibrate=False)` — and add the explicit calibration-skip branch.**
   - `dna_decode/models/classical_baselines.py` — extend `_train_baseline_logreg` signature with `calibrate: bool = True` (backward-compat default). Add an EXPLICIT calibration-skip branch BEFORE the existing CalibratedClassifierCV path:
     ```python
     def _train_baseline_logreg(X, y, drug_name, *, calibrate: bool = True) -> TrainedClassifier:
         # ... existing MIN_TRAINING_SAMPLES + single-class validation ...
         base = LogisticRegression(max_iter=1000, solver="liblinear", random_state=42)
         if not calibrate:
             base.fit(X, y)
             return TrainedClassifier(
                 model=base, drug_name=drug_name, feature_dim=X.shape[1], calibrated=False,
             )
         # ... existing CalibratedClassifierCV(base, ...) path ...
     ```
     The signature change is backward-compat: only one existing internal caller (`train_kmer_baseline` at line 275, positional args). Default preserved.
   - `scripts/stage1_n40_cipro.py` — delete the local `_train_logreg`. Import `_train_baseline_logreg` from `dna_decode.models.classical_baselines`. Wrap as:
     ```python
     def _nt_logreg_train(X, y):
         return _train_baseline_logreg(X, y, drug_name="cipro_stage1_nt_logreg", calibrate=False)
     def _nt_logreg_predict(clf, X):
         return clf.model.predict_proba(X)[:, 1].astype(np.float32)
     ```

4. **Eliminate the silent mean-fallback in k-mer + fusion paths.**
   - In `dna_decode/eval/loso_kmer.py`, when `train_xgboost_classifier` raises `ClassifierTrainingError`, **re-raise** rather than fall back to `train_y.mean()`. At Stage 1 N=38 with balanced 19R/19S, this exception only fires on a real bug (xgboost import, NaN inputs); silent fallback masks exactly the failure Stage 1 is supposed to catch.
   - Smoke gate's outer `main` already wraps each variant in `try/except` (smoke_gate_12strain_cipro.py:303-313) and surfaces failures in the packet — no behavioral regression to smoke at the orchestration layer; only the inner mean-fallback is removed.

5. **Split `write_packet`; new `compute_gate_outcome` returns BOTH `verdict_bucket` (unchanged) AND `stage2_action` (new decision-layer).**
   - `scripts/stage1_n40_cipro.py` — extract a new pure function:
     ```python
     def compute_gate_outcome(results: list[VariantResult]) -> dict:
         """Pure-function: takes variant results, returns
            {nt_best_name, nt_best_auroc, nt_best_scores, kmer_auroc,
             gap_pp, ci_mean_pp, ci_lo_pp, ci_hi_pp, ci_n_effective,
             verdict_bucket,       # pure function of gap_pp; unchanged 3-bucket logic
             stage2_action,        # decision-layer: BURST_STAGE_2 / HOLD_STAGE_2_CI_DEGENERATE
                                   # / ALTERNATIVE_POOLING_RERUN / PIVOT_TO_BAKTA
             fusion_note,
             fusion_outperforms_primary,
             fusion_alignment_valid}.
            DOES NOT write files. DOES compute paired bootstrap CI."""
     ```
   - **Alignment validation inside `compute_gate_outcome`:** before computing the gap and bootstrap CI, assert `nt_best.strain_ids == kmer.strain_ids` (element-wise, ordered). If mismatched, raise `ValueError("Stage 1 alignment: NT-best and k-mer strain_ids diverge; refactor regression")` — this is gate-bearing and must fail loudly.
   - **Fusion alignment is permissive:** check `fusion.strain_ids == nt_best.strain_ids`. If mismatched, set `fusion_alignment_valid = False` and `fusion_note = "fusion alignment mismatch — diagnostic suppressed"`. Do NOT raise — fusion isn't gate-bearing.
   - `verdict_bucket = verdict_label(gap_pp)` — `verdict_label` UNCHANGED from current 3-bucket logic (a pure function of `gap_pp`). Verdict semantics frozen.
   - `stage2_action` computed by a new helper `decide_stage2_action(verdict_bucket, ci_lo_pp, fusion_outperforms_primary) -> Literal["BURST_STAGE_2", "HOLD_STAGE_2_CI_DEGENERATE", "ALTERNATIVE_POOLING_RERUN", "PIVOT_TO_BAKTA"]` implementing the 4-rule mapping from Step 1.
   - `write_packet` becomes a thin wrapper: call `compute_gate_outcome(results)` → render markdown showing BOTH `verdict` and `stage2_action` lines → write file → return summary.

6. **Add `strain_ids` to `VariantResult` for alignment metadata.**
   - `scripts/stage1_n40_cipro.py` — extend the `VariantResult` dataclass with `strain_ids: list[str]`. Populated at construction from the producing variant's CV output (`CVResult.strain_ids` for NT variants; pass through input `strain_ids` for k-mer/fusion).

7. **Loud MLST handling.**
   - `scripts/stage1_n40_cipro.py:load_features` — replace `getattr(s, "mlst", "unknown")` with `getattr(s, "mlst", None)`. If any strain returns `None`, raise `ValueError(f"strain {s.strain_id} missing MLST — cohort metadata incomplete; fix before Stage 1")` rather than collapsing into a "unknown" bucket that would inflate `largest_mlst_group`.
   - `load_features` signature change: return `(X_nt, seqs_by_strain, labels_by_strain, strain_ids, mlsts)` so all variants consume the same `strain_ids` list (alignment guarantee from a single source of truth).

8. **Bootstrap-skip-count reporting.**
   - `scripts/stage1_n40_cipro.py:paired_bootstrap_ci` — return `(mean_gap, lo, hi, n_effective)` where `n_effective` is the count of non-degenerate resamples. `compute_gate_outcome` surfaces `n_effective` in its output dict; result packet renders `"B=1000 (effective {n_effective})"`. If `n_effective < 800`, append "CI honesty degraded; investigate cohort imbalance."

9. **Add `CVResult.strain_ids` property + `CONTIG_SEPARATOR` module constant.**
   - `dna_decode/eval/cv.py` — add a one-line property to `CVResult`:
     ```python
     @property
     def strain_ids(self) -> list[str]:
         return [f.held_out_id for f in self.folds]
     ```
     Single-line addition, no behavior change. Used by Step 2.2 to extract NT-variant strain order.
   - `dna_decode/models/classical_baselines.py` — add module-level constant at the top (just after the imports section, before `DEFAULT_KMER_K`):
     ```python
     CONTIG_SEPARATOR: str = "N" * 100
     ```
     Update `train_kmer_baseline` parameter default to reference it: `contig_separator: str = CONTIG_SEPARATOR`. The new `dna_decode/eval/loso_kmer.py` imports `from dna_decode.models.classical_baselines import CONTIG_SEPARATOR`. Stage 1 runner does the same in `load_features`. Removes the three magic-string `"N" * 100` copies.

10. **Trivial cleanups.**
    - Delete unused `from dna_decode.data.annotations import parse_gff3` import at `scripts/stage1_n40_cipro.py:37`.

**Key details:**
- New module path: `dna_decode/eval/loso_kmer.py` — placed under `eval/` because it's a CV-strategy variant tied to k-mer features; matches `cv.py`'s home.
- Backward-compat: `_train_baseline_logreg` gets an optional `calibrate: bool = True` kwarg; existing internal caller unaffected.
- Error semantics: re-raise `ClassifierTrainingError` from k-mer/fusion paths; the runner's `main` catches it and exits non-zero with a clear message rather than swallowing.
- The 11 existing tests must continue to pass after the refactor. If `compute_gate_outcome` extraction changes function-import locations, tests are updated accordingly (single file, minimal change). The existing 4 paired_bootstrap_ci tests need to unpack 4 elements instead of 3 — single-character changes each.
- Grep verifications before each sub-task (per `feedback_verify_stale_content_claims`):
  - Verify `loso_scores` exists at `scripts/stage1_n40_cipro.py` before deleting
  - Verify `_train_logreg` exists at `scripts/stage1_n40_cipro.py:75-80` before deleting
  - Verify `parse_gff3` import is at `scripts/stage1_n40_cipro.py:37`
  - Verify `contig_separator` is currently a PARAMETER (line 230) not a module constant in `classical_baselines.py`

**Test strategy:**
- Re-run all 11 existing tests in `tests/test_stage1_n40_cipro.py`; they must pass post-refactor (after 4-tuple unpacking adjustment).
- Manual smoke: run the helper imports (`from scripts.stage1_n40_cipro import compute_gate_outcome, decide_stage2_action, paired_bootstrap_ci, verdict_label, per_mlst_breakdown`) to confirm the new module structure is importable.
- Manual smoke: run `from dna_decode.eval.loso_kmer import run_kmer_xgboost_loso, run_fusion_loso` to confirm the new module loads.
- Manual smoke: run `from dna_decode.models.classical_baselines import CONTIG_SEPARATOR; from dna_decode.eval.cv import CVResult; r = CVResult(strategy="loso", drug="x"); assert hasattr(r, "strain_ids")` to confirm the constant + property addition.
- **Smoke gate runner regression: `HF_HOME=D:/hf_cache uv run python scripts/smoke_gate_12strain_cipro.py` should produce the same verdict + AUROCs within ±0.005 of the 2026-05-14 PASS run** (NT 0.750 / k-mer 0.694 / gene-presence INDETERMINATE / PASS). The ±0.005 tolerance accommodates dtype/ordering implementation details while still catching behavioral drift.
- Full test suite green: `uv run pytest tests/ -m 'not slow' -q` — must remain at 369 passed / 1 skipped post-refactor (Step 3 then bumps to 376 passed / 1 skipped).

---

### Step 3: Add fusion-exclusion + calibration-discipline + bootstrap-skip-count tests
Files: tests/test_stage1_n40_cipro.py
Depends on: Step 2

**What changes:**

1. **`TestGateOutcomeFusionExcluded` test class** — 2 tests against the new `compute_gate_outcome` function:
   - `test_fusion_excluded_from_gate_when_fusion_wins`: build synthetic `VariantResult` list where `NT+k-mer-fusion-logreg.auroc = 0.99`, `NT-XGBoost.auroc = 0.50`, `NT-logreg.auroc = 0.50`, `k-mer-XGB.auroc = 0.50`. Assert `compute_gate_outcome(results)["verdict_bucket"]` matches the FAIL label. Assert `compute_gate_outcome(results)["fusion_outperforms_primary"] is True`. Pins /brainstorm's exact failure mode.
   - `test_fusion_ignored_when_nt_xgboost_wins_alone`: NT-XGBoost = 0.80, NT-logreg = 0.70, k-mer-XGB = 0.70, fusion = 0.50. Assert verdict is CLEAN PASS (gap = 10 pp ≥ 5 pp) and `nt_best_name == "NT-XGBoost"`.

2. **`TestCalibrationDiscipline` test class** — **3 tests** using `monkeypatch`:
   - `test_nt_xgboost_call_site_passes_calibrate_false`: monkeypatch `train_xgboost_classifier` to record `(args, kwargs)`. Run the Stage 1 NT-XGBoost path with a tiny synthetic 4-strain cohort. Assert the recorded `kwargs["calibrate"]` is `False`.
   - `test_nt_logreg_call_site_passes_calibrate_false`: monkeypatch `_train_baseline_logreg`. Assert kwargs contain `calibrate=False` (post-Step-2, the helper accepts this kwarg).
   - **(NEW per /brainstorm Round 2)** `test_kmer_xgb_call_site_passes_calibrate_false`: monkeypatch `train_xgboost_classifier`. Run `run_kmer_xgboost_loso(...)` with a tiny synthetic cohort. Assert kwargs contain `calibrate=False`. Pins the k-mer-XGB call site that the original test plan missed.
   - All three are regression guards against the N=11 isotonic-collapse bug. The "easy revert" PR that adds calibration back would trip these tests.

3. **`TestStage2Action` test class** — **(NEW per /brainstorm)** 4 tests pinning the decision-layer:
   - `test_clean_pass_ci_clears_zero_returns_burst`: verdict=CLEAN PASS, ci_lo=2.0 → `stage2_action == "BURST_STAGE_2"`.
   - `test_clean_pass_ci_negative_still_returns_burst`: verdict=CLEAN PASS, ci_lo=-0.5 → `stage2_action == "BURST_STAGE_2"` (wide-CI annotation but still proceed).
   - `test_noisy_pass_ci_negative_returns_hold`: verdict=NOISY PASS, ci_lo=-1.2 → `stage2_action == "HOLD_STAGE_2_CI_DEGENERATE"`.
   - `test_fail_returns_alternative_pooling_rerun`: verdict=FAIL, any ci_lo → `stage2_action == "ALTERNATIVE_POOLING_RERUN"`.
   - These pin the Step 1 pre-commitment rules as code-tested decision logic, not just plan-file prose.

4. **`TestBootstrapSkipCount` test class** — 1 test:
   - `test_paired_bootstrap_returns_n_effective`: with all-same-class labels (degenerate), assert `paired_bootstrap_ci(...)` returns a 4-tuple `(mean, lo, hi, n_effective)` and that `n_effective < n_iterations` (since all resamples are degenerate and skipped). This pins the Step 2.8 signature change.

5. **`TestMlstLoudHandling` test class** — 1 test:
   - `test_load_features_raises_on_missing_mlst`: build a mock cohort with a strain whose `mlst` attribute is `None`. Assert `load_features(...)` raises `ValueError` mentioning the strain ID. (May require minimal mocking of the EmbeddingCache + FASTA loading — keep it small; the test is about the MLST guard, not real LOSO.)

6. **`TestStrainIdAlignment` test class** — **(NEW per /brainstorm Round 2)** 2 tests:
   - `test_compute_gate_outcome_raises_on_nt_vs_kmer_mismatch`: build `VariantResult` list where NT-XGBoost and k-mer-XGB have different `strain_ids` lists. Assert `compute_gate_outcome(results)` raises `ValueError` mentioning "alignment" or "strain_ids."
   - `test_compute_gate_outcome_suppresses_fusion_note_on_fusion_mismatch`: NT-XGBoost / NT-logreg / k-mer-XGB all aligned; fusion has different `strain_ids`. Assert `compute_gate_outcome(results)` does NOT raise; `fusion_alignment_valid is False`; `fusion_note` contains "alignment mismatch."

7. **Adjust existing test imports + 4-tuple unpacking.**
   - Add `from scripts.stage1_n40_cipro import compute_gate_outcome, decide_stage2_action` to test imports.
   - The 4 existing `paired_bootstrap_ci` tests use 3-tuple unpacking (`mean, lo, hi = ...`). Update to 4-tuple (`mean, lo, hi, _ = ...` or `mean, lo, hi, n_eff = ...`). Single-character diff per test.

**Key details:**
- Total new tests: **13** (2 fusion-exclusion + 3 calibration-discipline + 4 stage2_action + 1 bootstrap-skip + 1 MLST-loud + 2 strain_id-alignment). Total test count post-Step-3: **24 in `tests/test_stage1_n40_cipro.py`** (was 11), **382 in the full suite** (was 369).
- `monkeypatch` is a pytest fixture; already available.
- Mocking strategy for the MLST test: use a minimal `dataclass` mock for the cohort strain rather than instantiating the real `CandidateStrain` (which has many required fields). If mocking proves disproportionate, accept the test as a smoke-only assertion that uses a real fixture cohort entry with `mlst=None` injected post-load.

**Test strategy:**
- All 24 tests pass: `uv run pytest tests/test_stage1_n40_cipro.py -v`.
- Full suite remains green: `uv run pytest tests/ -m 'not slow' -q` → 382 passed / 1 skipped.
- Stretch: a parametrized version of fusion-exclusion covering 3 verdict buckets × fusion-wins-or-not = 6 micro-tests is acceptable but optional.

## Execution Preview

```
Wave 0 (2 parallel):  Step 1 — Add pre-commitment language to plan,  Step 2 — Refactor runner + factor k-mer LOSO
Wave 1 (1 sequential): Step 3 — Add fusion-exclusion + calibration-discipline + bootstrap + MLST tests
```

Critical path: Step 2 → Step 3 (2 waves)
Max parallelism: 2 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- **Step 2 is wide-ranging within a single file** — 8 sub-tasks all modify `scripts/stage1_n40_cipro.py`. Cannot be parallelized further without extracting helpers to new modules (would add structure for a small project). Mitigation: sub-task list is enumerated in Step 2's "What changes"; each sub-task is small and can be committed individually if needed.
- **Step 2 also touches `scripts/smoke_gate_12strain_cipro.py`** — factoring `run_kmer_xgboost` out alters the smoke runner. Mitigation: explicit regression check is "smoke gate re-run produces 2026-05-14 result packet PASS verdict + same AUROCs" in the test strategy.
- **`_train_baseline_logreg` gets a new kwarg `calibrate: bool = True`** — backward-compat default preserves all existing callers. Verified no other tests pass `calibrate` to it (no current signature collision).
- **Step 1's pre-commitments are defaults** — user can edit any of the 4 rules before Stage 1 executes. The Step doesn't force consensus, it forces explicitness.
- **Bootstrap-skip-count signature change** — `paired_bootstrap_ci` returns 4 elements instead of 3. Existing 4 bootstrap tests need ≤4-character changes each (e.g., `mean, lo, hi = ...` → `mean, lo, hi, _ = ...`). Single-file impact.
- **MLST loud-handling could break edge cohorts** — if any historical cohort parquet was saved with a strain having `mlst=None`, Stage 1 will now refuse to load it. Mitigation: explicit error message says "fix cohort metadata before Stage 1"; the 12-strain + N=40 current cohorts have MLST per the earlier diagnostic check; smoke gate's 12-strain run on the new code path will verify before populate finishes.
- **No file overlaps in Wave 0** — Step 1 touches plans/, Step 2 touches scripts/ + dna_decode/eval/. Step 3 (Wave 1) touches tests/.
- **Restructuring applied:** none required; original 3-step shape held after dependency analysis. The runner-refactor monolithic step (Step 2) cannot be split into parallel sub-steps without extracting helpers to new modules, which the /review explicitly noted should be done in the spirit of "minimal diff over large rewrites."

## Verification

- All 24 tests pass in `tests/test_stage1_n40_cipro.py` after Step 3 lands: `uv run pytest tests/test_stage1_n40_cipro.py -v`.
- Full test suite: `uv run pytest tests/ -m 'not slow' -q` returns 382 passed / 1 skipped / 0 failed.
- **Smoke gate regression (tolerance-loosened per /brainstorm Round 2):** `HF_HOME=D:/hf_cache uv run python scripts/smoke_gate_12strain_cipro.py` produces verdict PASS, NT-XGBoost AUROC within ±0.005 of 0.750, k-mer AUROC within ±0.005 of 0.694, gene-presence INDETERMINATE_IDENTIFIER_OOV. The ±0.005 tolerance accommodates dtype/ordering implementation details from the factoring; tighter than the visible-rounding (±0.001) but looser than exact-match (catches behavioral drift without flagging on float-ordering noise).
- Plan file renders: `plans/Stage1_N40_Cipro_Engineering_Screen_Plan.md` shows the new `## Verdict-Time Pre-Commitments` section between Design Decisions and Verification; all 4 rules are present AND structured as a decision layer over the unchanged verdict.
- Post-populate end-to-end (only verifiable after the N=40 NT cache finishes populating in the background): `HF_HOME=D:/hf_cache uv run python scripts/stage1_n40_cipro.py` produces `wiki/stage1_n40_cipro_<date>.md` with: 4 variant AUROCs, gate analysis line, paired bootstrap CI with `B=1000 (effective N_eff)`, per-MLST table, per-strain LOSO predictions table, AND TWO lines: `verdict` (point-gap function only) + `stage2_action` (decision-layer output). Exit code 0 for PASS-like `stage2_action`, 1 for FAIL-like.
- Honest failure surfacing: a synthetic `ClassifierTrainingError` injected into NT-XGBoost training (via local edit) causes Stage 1 to exit non-zero with a clear traceback rather than silently producing AUROC≈0.5 from mean-fallback.
- Alignment validation triggers loudly: if `compute_gate_outcome` receives NT-best and k-mer variants with mismatched `strain_ids`, raises `ValueError` (not silent comparison on misaligned arrays).
