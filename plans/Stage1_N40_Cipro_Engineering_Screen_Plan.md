# Stage 1 N=40 Cipro Engineering Screen — Plan

> Run a 4-experiment matrix (NT-XGBoost gate + NT-logreg sanity + k-mer-XGB classical + NT+k-mer-fusion-logreg diagnostic) under LOSO on the N=40 cipro cohort (effective N=38) with paired bootstrap CI, MLST diagnostic appendix, and a 3-bucket verdict to decide whether to spend Stage 2 N=150 Databricks burst budget.

---

## Problem Statement

Phase 2 entry HARD gate (12-strain cipro smoke gate) PASSED 2026-05-14 with NT-XGBoost AUROC=0.750 vs k-mer-XGB 0.694 (+5.6 pp NT lift). Stage 1 of Phase 2 is now active: a larger-N engineering screen that decides go/no-go on spending Stage 2 N=150 Databricks burst budget.

Locked criterion (2026-05-14 user lock):
> **Stage 1 PASSES iff max(NT-XGBoost, NT-logreg) AUROC ≥ k-mer-XGB AUROC + 3 pp under LOSO at N=40 cipro (effective N=38 — 2 GCA strains missing GFF3, split 1R+1S, new balance 19R/19S).**

Stage 2 then tightens to Option-C (≥5 pp AUROC + top-K attribution includes ≥1 of {gyrA, parC, parE}) under leave-one-Mash-clade-out CV.

The /brainstorm review (2026-05-14) flagged three critical issues with the initial proposed 3-experiment matrix:
1. The proposed matrix had drifted from the locked spec — NT-logreg was the primary head, but the locked criterion is NT-XGBoost.
2. Calibration was an uncontrolled variable across variants (smoke runner already disabled calibration after the N=11 isotonic collapse bug).
3. LOSO alone is too weak as a lineage-confounding screen at high MLST cardinality.

This plan captures the revised matrix and adds the diagnostic appendix that survives those critiques.

## Design Decisions

### D1: Restore NT-XGBoost as the primary gate-bearing head; add NT-logreg as sanity-check baseline

**Decision:** Stage 1 runs a 4-experiment matrix:
- **NT-XGBoost** (primary, gate-bearing) — matches the locked criterion + the 12-strain smoke gate.
- **NT-logreg** (sanity-check baseline, gate-bearing) — confirms XGBoost isn't winning by overfitting; if NT-logreg materially beats NT-XGBoost, plumbing-bug check triggers per H13 in the project ledger.
- **k-mer-XGB** (classical comparator, gate-bearing) — the "would a simple k-mer baseline equally predict?" control.
- **NT+k-mer-fusion-logreg** (diagnostic only, NOT gate-bearing) — surfaces complementary-signal question; flags k-mer-carrying-the-result scenarios.

Gate formula: `max(NT-XGBoost, NT-logreg) − k-mer-XGB ≥ 3 pp`. Fusion is excluded from the gate. Fusion passing while both NT-only heads fail = Stage 1 FAIL (k-mer carrying it).

**Rationale:** The locked Stage 1 spec consistently says "NT-XGBoost ≥3 pp" (project ledger Mid-term Milestone #3, Evidence row 28, Phase2_Decision_Gate_Plan.md). The initial proposed matrix that swapped NT-logreg as the primary head was justified by a /research finding (GBDT beats deep MLPs at small N) — but that finding is about MLPs vs GBDT, not logreg vs GBDT, and on dense 512-dim continuous features the question is different. Adding NT-logreg alongside NT-XGBoost gives a head-comparison sanity check without changing the gate-bearing variant.

**Trade-off:** Two NT-head experiments cost ~2× the NT-classifier wallclock vs running one. Acceptable — LOSO at N=38 with cached embeddings runs in minutes per variant. Considered keeping NT-logreg only (cheaper, more interpretable) — rejected because it drifts from the locked spec.

### D2: All variants run with `calibrate=False` for primary AUROC

**Decision:** Every gate-bearing AND diagnostic variant uses `calibrate=False` in `train_xgboost_classifier` and a raw `LogisticRegression` (no `CalibratedClassifierCV` wrapper) for the logreg paths.

**Rationale:** `CalibratedClassifierCV` runs an internal CV that can re-order scores → AUROC measurement changes despite AUROC being rank-based in theory. Mixing calibrated and uncalibrated variants makes Stage 1 measure calibration-wrapper behavior rather than representation quality. The smoke runner already chose `calibrate=False` after the N=11 isotonic collapse bug; Stage 1 inherits that discipline.

**Trade-off:** If calibration metrics (Brier, ECE) are desired in the result packet, they need to be computed as a secondary measurement, never used as primary gate. Not implementing them at Stage 1 — defer to Stage 2 where N supports it. Considered leaving smoke-default calibration on — rejected because the calibration wrapper behavior is the load-bearing risk.

### D3: Add diagnostic appendix to Stage 1 result packet (MLST + per-strain + paired bootstrap CI)

**Decision:** The result packet includes a "Lineage diagnostic" section reporting:
- Unique MLST count, uniqueness fraction, largest MLST group with R/S split.
- Per-strain LOSO predictions table (strain ID, MLST, true label, NT-best score, k-mer score, per-classifier correctness check).
- Paired bootstrap 95% CI on the gap `NT-best − k-mer-XGB` (B=1000, paired resampling of strain indices).
- 3-bucket verdict label: **CLEAN PASS (≥5 pp)** / **NOISY PASS (3-5 pp, diagnostics flagged)** / **FAIL (<3 pp)**.

LOMO is NOT computed at Stage 1 — at N=38 with high MLST cardinality (most MLST groups will be size 1), LOMO degenerates to LOSO. The per-MLST table + per-strain predictions substitute as the lineage diagnostic.

**Rationale:** A +3 pp gap inside the ±0.10 LOSO noise floor at N=38 can plausibly be cohort-structure artifact rather than biology. The bootstrap CI surfaces the clean-vs-noisy distinction honestly. The per-strain MLST overlay lets a reviewer detect "won by lineage memorization." The 3-bucket verdict gives an explicit tie-breaker label without re-litigating the locked ≥3 pp threshold.

**Trade-off:** Paired bootstrap with B=1000 adds ~30 seconds of compute; negligible. Considered DeLong's analytic test instead — rejected because paired bootstrap is more flexible and matches what most ML papers use for AUROC differences.

### D4: Gene-presence + AMRFinderPlus baselines explicitly out of scope for Stage 1

**Decision:** Stage 1 reports "best classical baseline = k-mer-XGB" with an explicit caveat that gene-presence and AMRFinderPlus POINT* SNP-table baselines are absent from the comparison. The packet notes "'Best classical' here is bounded by what was run; gyrA/parC/parE point-mutation features are NOT part of the comparator."

**Rationale:**
- Gene-presence is INDETERMINATE_IDENTIFIER_OOV on this annotation source — RefSeq GFF3 carries `gene=` for ~11% of CDSs (see `plans/Gene_Presence_AUROC_Bug_Fix_Plan.md`). Including it would require Bakta re-annotation; that's a separate deferred decision.
- AMRFinderPlus POINT* is deferred to Stage 2 per `plans/Phase2_Decision_Gate_Plan.md` D6.

**Trade-off:** Stage 1 may declare PASS while a SNP-table baseline (if added) would close the gap to <3 pp. Acceptable risk because Stage 2 (N=150 + Mash-clade-out CV + Option-C threshold ≥5 pp + biology check) is the real ship gate — Stage 1 is the engineering screen that gates spending Stage 2 budget. The explicit caveat in the packet prevents over-claiming.

## Implementation Plan

1. **Pre-flight check: confirm R/S split of the 2 missing GCA strains.**
   - Files: query against `data/processed/gate_b_n40_cipro_cohort.parquet`
   - Result confirmed 2026-05-14: GCA_902807115.1 = R, GCA_008727135.1 = S → 1R+1S → effective N=38 stays 19R/19S balanced.

2. **Write Stage 1 runner script.**
   - File: `scripts/stage1_n40_cipro.py` (new)
   - Loads NT embeddings from `D:/dna_decode_cache/embeddings/nt_n40_cipro.h5` + FASTA contigs from `D:/dna_decode_cache/refseq/<accession>/genome.fna`.
   - Functions: `load_features`, `loso_scores` (generic fixed-X LOSO), `loso_kmer_xgb` (within-fold vocab rebuild), `loso_fusion_logreg`, `paired_bootstrap_ci`, `verdict_label`, `per_mlst_breakdown`, `write_packet`.
   - Defaults: `GATE_THRESHOLD_PP = 3.0`, `BOOTSTRAP_ITERATIONS = 1000`, `BOOTSTRAP_SEED = 42`, `kmer_k = 8`, `kmer_top_n = 10_000`.
   - Output: `wiki/stage1_n40_cipro_<date>.md`. Exit 0 if "PASS" in verdict, else 1.

3. **Add unit tests for the verdict + bootstrap helpers.**
   - File: `tests/test_stage1_n40_cipro.py` (new)
   - Coverage: 5 verdict-label tests (boundary at 3 pp, boundary at 5 pp, fail case, negative gap, pin GATE_THRESHOLD_PP=3.0), 4 paired-bootstrap-CI tests (identical scores → zero-centered CI; A strictly better → positive CI; return-shape; degenerate-resample-handling), 2 per-MLST-breakdown tests.

4. **Update CLAUDE.md + TODOS.md.**
   - File: `CLAUDE.md` — add Stage 1 runner to "Common commands" block; add populate command for completeness.
   - File: `TODOS.md` — close the Stage 1 criterion entry; record the runner shipped.

5. **Wait for N=40 NT populate to finish** (running in background since 2026-05-14 PM; ~5-7 hr GTX 860M).

6. **Execute Stage 1 runner against the populated cache.**
   - Command: `HF_HOME=D:/hf_cache uv run python scripts/stage1_n40_cipro.py`
   - Expected wallclock: ~1-2 hr (LOSO over 38 strains, 4 variants).
   - Output: `wiki/stage1_n40_cipro_<date>.md` result packet + console verdict.

7. **Synthesize verdict + next-action recommendation.**
   - **CLEAN PASS (≥5 pp)** → proceed to Stage 2 Databricks burst with N=150 cohort build.
   - **NOISY PASS (3-5 pp)** → review diagnostic appendix (MLST overlap, per-strain table, CI lower bound). If CI lower bound < 0, treat as borderline; tighten 1-2 experiments before Stage 2 commitment.
   - **FAIL (<3 pp)** → NT track doesn't justify Stage 2 spend at this cohort scale. Options: rerun with alternative pooling / dimensionality reduction / drop to k-mer-only Stage 2 / pivot to Bakta annotation + gene-presence comparator.

## Verdict-Time Pre-Commitments (locked 2026-05-14)

This section is a **decision-layer ON TOP OF the unchanged verdict semantics**. Verdict labels (`CLEAN PASS` ≥5 pp / `NOISY PASS` 3-5 pp / `FAIL` <3 pp) remain pure functions of the point gap `(NT-best − k-mer-XGB)` per the locked criterion above. The pre-commitments determine a separate `stage2_action` field based on (verdict + CI lower bound + fusion behavior). The result packet renders BOTH `verdict` and `stage2_action` so the underlying point gap stays auditable while the budget decision is unambiguous.

### Rule 1: CI-lower-bound budget-decision rule
Drives `stage2_action`, **NOT** the verdict label.
- **CI lower bound > 0** → `stage2_action` follows the verdict bucket directly (no degradation).
- **CI lower bound ≤ 0 AND verdict = NOISY PASS** (gap ∈ [3, 5) pp) → `stage2_action = HOLD_STAGE_2_CI_DEGENERATE` (effective FAIL handling for budget purposes). Verdict label STAYS `NOISY PASS`.
- **CI lower bound ≤ 0 AND verdict = CLEAN PASS** (gap ≥ 5 pp) → `stage2_action = BURST_STAGE_2` (proceed) but result packet annotates "wide CI; lineage check load-bearing." Verdict label STAYS `CLEAN PASS`.

### Rule 2: Asymmetric-cost statement
We are more willing to make a false-FAIL than a false-PASS at Stage 1. False-PASS only costs Stage 2 burst budget (atomic, recoverable via Stage 2 verdict — Option-C threshold + biology check); false-FAIL silently kills a potentially useful track. The ≥3 pp threshold reflects this trade.

### Rule 3: Fusion-outperforms-primary rule
If `fusion_AUROC − max(NT-XGBoost, NT-logreg) ≥ 3 pp` AND fusion alignment is valid:
- Log "Stage 2 architecture note: fusion outperformed both NT-only heads, revisit at Stage 2" in the result packet.

If fusion alignment is invalid (strain_ids mismatch with NT-best):
- Suppress the note entirely with "fusion alignment mismatch — diagnostic suppressed."

Fusion NEVER alters the gate-bearing gap computation or the verdict label.

### Rule 4: Per-bucket `stage2_action` mapping (deterministic; no "review and decide")
Literal values: `BURST_STAGE_2` / `HOLD_STAGE_2_CI_DEGENERATE` / `ALTERNATIVE_POOLING_RERUN` / `PIVOT_TO_BAKTA`.

| verdict_bucket | CI lo > 0 | CI lo ≤ 0 |
|---|---|---|
| CLEAN PASS (≥5 pp) | `BURST_STAGE_2` | `BURST_STAGE_2` (annotate "wide CI") |
| NOISY PASS (3-5 pp) | `BURST_STAGE_2` (annotate "CI clears 0") | `HOLD_STAGE_2_CI_DEGENERATE` |
| FAIL (<3 pp) | `ALTERNATIVE_POOLING_RERUN` | `ALTERNATIVE_POOLING_RERUN` |

Action semantics:
- **`BURST_STAGE_2`** → Stage 2 Databricks burst with N=150 cohort build. No further deliberation.
- **`HOLD_STAGE_2_CI_DEGENERATE`** → DO NOT spend Stage 2 burst budget. Next: `ALTERNATIVE_POOLING_RERUN` (Stage 1b with `mean+max` aggregation instead of `mean`); if still NOISY+degenerate, escalate to `PIVOT_TO_BAKTA`.
- **`ALTERNATIVE_POOLING_RERUN`** → Run Stage 1 once with `mean+max` aggregation; if still <3 pp, escalate to `PIVOT_TO_BAKTA`. Do NOT spend Stage 2 burst on NT-only.
- **`PIVOT_TO_BAKTA`** → re-annotation + gene-presence comparator pathway per `plans/Gene_Presence_AUROC_Bug_Fix_Plan.md` follow-up.

## Verification

- `uv run pytest tests/test_stage1_n40_cipro.py -v` — 11 tests pass (verdict bucketing + paired bootstrap CI + per-MLST breakdown). Confirmed 2026-05-14.
- Full test suite remains green: `uv run pytest tests/ -m 'not slow' -q` — currently 369 passed / 1 skipped after the gene-presence fix shipped earlier this session.
- Post-populate end-to-end: Stage 1 runner produces a valid result packet at `wiki/stage1_n40_cipro_<date>.md` with all 4 variant AUROCs, gate analysis, paired bootstrap CI, per-MLST table, per-strain LOSO predictions, and a single verdict line. Exit code matches verdict (0 for PASS, 1 for FAIL).
- Quick correctness check: if NT-XGBoost AUROC ≥ NT-logreg + 3 pp, the per-strain MLST table should NOT show NT-XGBoost only winning on MLST groups absent from training folds (lineage memorization signal). If lineage signal dominates, the verdict label should annotate this even on a numerical pass.
