# Phase 2 Decision Gate Plan

> Split the originally-conflated "12-strain decision gate" into TWO gates: a 12-strain **smoke/falsification** gate that runs immediately, and a 50-150-strain **real decision** gate that adds the /research-suggested classifier alternatives. Based on the /research output 2026-05-13 + /brainstorm reframe of the same date.

---

## Problem Statement

The /research run on 2026-05-13 surfaced 5 candidate additions to the dna_decode decision gate: Random Forest as downstream classifier, TabPFN for small-sample regime, SNP-table feature variant for cipro, classical-baseline >90% accuracy bar to beat, and SNP > gene-presence framing for fluoroquinolones. The natural impulse was to expand the originally-planned 5-variant 12-strain decision gate into an 8-variant comparison.

The follow-up /brainstorm flagged that this expansion is **mis-sized for N=12**:

1. **Statistical power problem.** LOSO at N=12 means each misclassification shifts accuracy ±8.3%. An apparent gap of one strain between two classifier variants is noise, not signal. Adding RF + TabPFN + SNP-table to a 5-variant comparison turns the gate into an 8-way multiple-comparison problem with no power to rank.

2. **TabPFN dimension mismatch.** TabPFN-2 supports up to 500 features; NT outputs 512-dim. Adding PCA changes the experimental variable from a clean head swap to `NT + reducer + TabPFN`. Reproducibility requires pinning TabPFN version explicitly.

3. **Classical side not credible at N=12.** The Talamantes-Becerra >90% bar is on 256 genomes × 5-fold CV; we'd be comparing against 12 strains × LOSO. SNP-table support requires a variant-calling/alignment pipeline (not yet present in `classical_baselines.py`). And the clade-only baseline at `scripts/pipeline.py:285` uses `hash(s.mlst) % 10` placeholder — comparing NT lift against a synthetic baseline is misleading.

The reframe: 12-strain run as **smoke/falsification gate** (does NT show ANY signal vs. an honest clade-only baseline?), real **decision gate** at 50-150 strains with the /research-suggested additions implemented.

## Design Decisions

### D1: 12-strain = smoke gate, 50-150-strain = decision gate

**Decision:** Split into two distinct gates with different acceptance bars.

- 12-strain smoke: acceptance bar = "NT-XGBoost is not obviously worse than k-mer + clade-only-fixed." Sufficient to falsify "NT is broken on real data." No promotion or demotion based on smoke result.
- 50-150-strain decision: acceptance bar = "NT path beats classical baseline by a margin larger than single-fold noise (~3-5%)." This is where promotion vs demotion of NT happens.

**Rationale:** N=12 has 8.3%/strain noise floor — can reject obvious failures but cannot rank classifier variants. N=50 brings the per-fold noise floor to 2%; N=150 to 0.67%. The literature bar (>90% on Talamantes-Becerra) is comparable to N=150 LOSO.

**Trade-off:** Considered keeping a single decision-gate at 12 strains with conservative thresholds (rejected — multiple-comparison risk dominates), and considered defining decision gate at N=12 with 80% bar (rejected — too easy; not falsifying anything meaningful).

### D2: Fix clade-only placeholder BEFORE running smoke gate

**Decision:** Replace `scripts/pipeline.py:285` `hash(s.mlst) % 10` with real MLST groups. This is a prerequisite for the smoke gate to produce meaningful clade-only baseline numbers.

**Rationale:** Comparing NT lift against a `hash(mlst) % 10` baseline is comparing against a synthetic random partition, not actual lineage structure. The /brainstorm flagged this as load-bearing. MLST values are already in `CandidateStrain`; the fix is small.

**Trade-off:** Considered running smoke with placeholder + noting the caveat (rejected — risk of misreading the result), and considered jumping directly to Mash CLI install (rejected — Mash CLI is a Windows-install pain; MLST grouping is sufficient for smoke).

### D3: Keep 12-strain smoke at existing 5 variants — no RF, TabPFN, SNP-table additions yet

**Decision:** Smoke gate runs exactly: AMRFinder / k-mer / gene-presence / clade-only-fixed / NT-XGBoost. The /research-suggested additions (RF, TabPFN, SNP-table) wait until the 50-150-strain decision gate.

**Rationale:** Adding 3 variants at N=12 is multiple-comparison noise without ranking power. Better to keep smoke narrow + ship fast, then expand at decision gate where the differences are statistically detectable.

**Trade-off:** Considered running RF as a cheap-to-add additional variant (rejected — even a single additional variant adds noise and gives false confidence that the gate "compared classifiers"; better to be honest that smoke ≠ comparison).

### D4: Document the smoke result as smoke, not decision

**Decision:** Result packet from the 12-strain gate will explicitly say "smoke gate — multiple-comparison statistical power at N=12 forbids classifier ranking. Real decision gate scheduled at 50-150 strains."

**Rationale:** Future-self (or any reader of `wiki/`) might re-read the 12-strain result and treat it as a head-to-head. The honest framing prevents that.

**Trade-off:** Considered just reporting AUROC numbers without framing (rejected — invites misinterpretation).

### D5: Decision gate as a tiered N=50 → N=150 staged path

**Decision:** Two-stage decision gate, not a single fixed-N gate.

- **Stage 1 — N=50 local screen** (~4 hours of GTX 860M time + train + attribute):
  - Cohort: 50 strains balanced cipro (extend `gate_b_mini_cohort.parquet` via `build_mini_cohort.py --per-class 25`)
  - Acceptance criterion: NT-XGBoost shows **ANY positive AUROC lift over clade-only-fixed baseline**
  - Purpose: cheap screen that catches broken NT before burning Databricks credit
  - PASS → proceed to Stage 2
  - FAIL → demote NT track; classical baselines become primary project spine; no Databricks needed
- **Stage 2 — N=150 Databricks burst** (~one bounded cloud job, only if Stage 1 passes):
  - Cohort: 150 strains balanced cipro (build via expanded `build_mini_cohort.py`)
  - Acceptance criterion: BOTH must hold:
    - (a) NT AUROC ≥ best classical baseline + **5 pp**
    - (b) Top-10 attributed genes for NT include at least one of `gyrA`, `parC`, or `parE` (known cipro resistance loci — biological-plausibility check)
  - Purpose: rigorous decision with biological validity, apples-to-apples-ish vs. Talamantes-Becerra (256 genomes)
  - PASS → NT confirmed as primary path; proceed with full cohort + Phase 2 backlog
  - FAIL → NT demoted to "interesting but not primary"; classical baselines become spine

**Rationale:** N=50 local has ~2% noise floor — sufficient to falsify broken NT cheaply. N=150 Databricks brings noise floor to 0.67% — sufficient to detect a 5 pp signal cleanly. The staged path means we only spend cloud credit when local evidence says the path is worth pursuing. The biological-plausibility check (b) catches clade-leakage failures that pure AUROC numbers can miss — for cipro specifically, gyrA/parC/parE are textbook resistance loci, so a model that doesn't surface them is finding the wrong signal regardless of AUROC.

**Trade-off:** Considered direct-to-N=150 (rejected — wastes Databricks credit if NT is obviously broken at N=50); considered N=50-only with 5 pp + biology threshold (rejected — N=50 noise floor 2% can't cleanly detect 5 pp deltas); considered numeric-only threshold at Stage 2 without biology check (rejected — cipro clade-leakage is the dominant failure mode the gate is supposed to catch).

### D6: TODOS additions for the deferred /research items

**Decision:** Add four `[OPEN]` TODOS entries gated on the smoke gate's outcome:

1. **`[OPEN] Replace clade-only hash(mlst) % 10 placeholder with real MLST groups`** — small fix, gates the smoke gate.
2. **`[OPEN] Add RF wrapper to classifiers.py`** — straightforward, for Stage 2 decision gate.
3. **`[OPEN] Add TabPFN wrapper to classifiers.py with PCA-to-≤2000-feature step; pin package version`** — moderate, for Stage 2 decision gate. NT outputs 512-dim → fits TabPFN-2.6's 2000-feature envelope directly, but PCA still recommended for stability + reproducibility.
4. **`[OPEN] Add SNP-table feature variant — parse AMRFinderPlus output (--organism Escherichia --mutation_all) for POINT* rows`** — `gyrA_S83I` / `parC_S80I`-style point mutations are already in AMRFinderPlus's standard output (Method column = `POINT`, `POINTP`, `POINTX`, or `POINTN`). Treat as binary features alongside existing gene-presence; gyrA / parC / parE rows feed cipro-specific classical baseline. Scope: hours, not days (no variant-calling pipeline needed).

**Rationale:** Tracks the work without forcing it now. SNP-table scope correction per Codex 2026-05-13: AMRFinderPlus mutation rows are first-class output, not a separate feature extraction project.

**Trade-off:** Considered filing as `[BLOCKED]` (rejected — they're deferred for size reasons, not architecturally blocked).

## Implementation Plan

### Step 1: Configure smoke gate to skip clade-only (B-B path)

**Files:** `scripts/pipeline.py` invocation in the smoke gate runner.

**Updated 2026-05-14 per B-B decision lock** (`Sidework_Sequence_Ship_Path_Plan.md` D1). The originally-planned "fix clade-only placeholder before smoke" step is replaced with "drop clade-only from smoke entirely":

- Smoke gate invocation does NOT pass `--include-clade-baseline`. Clade-only baseline is meaningless at N=12 (12 unique MLST → singleton clades → AUROC ≈ 0.5 regardless of foundation model). Confirmed empirically via post-save /brainstorm.
- The `--include-clade-baseline` flag remains in `pipeline.py` for Stage 1 N=50 use.
- The `mlst_to_clade_id` helper + deferred per_clade_baseline strain-keying fix are scheduled for Stage 1 prep, not for the 12-strain smoke gate.

### Step 2: Run 12-strain smoke gate

**Cohort:** existing `data/processed/gate_b_mini_cohort.parquet` (12 strains, 6R/6S cipro).

Variants to run (**3 total** — clade-only dropped per B-B lock 2026-05-14; AMRFinder deferred 2026-05-14 — cohort's persisted `plasmid_resistance_genes + chromosome_resistance_genes` fields are empty for the 12-strain mini cohort, and no per-strain AMRFinderPlus CLI infrastructure exists yet; defer to Stage 1 N=50 prep when AMRFinder driver lands as a separate TODO):
- k-mer + XGBoost (existing in `classical_baselines.py`)
- gene-presence + XGBoost (existing)
- NT-XGBoost (existing; uses `mini_cipro_nt_cache.h5`)

Capture result packet per variant:
- LOSO AUROC, AUPRC
- Label balance check
- Top-K attribution genes for NT-XGBoost (Tier 1-5 classification; check if gyrA / parC / parE appear)

Acceptance bar (smoke): **NT-XGBoost is not obviously worse than the best classical baseline (max of k-mer / gene-presence / AMRFinder).** "Obviously worse" = NT-XGBoost AUROC < (best_classical - 0.15) — that's a 15-percentage-point gap, well outside N=12 noise. Clade-only is intentionally dropped (B-B lock 2026-05-14): 12 unique MLST strings make it degenerate; deferred to Stage 1 N=50.

### Step 3: Document smoke result as smoke

**Files:** `wiki/smoke_gate_12strain_cipro_2026-05-13.md` (new) or appended to `wiki/phase1_ship_report.md`.

Header: "12-strain cipro smoke gate — multiple-comparison statistical power at N=12 forbids classifier ranking. Real decision gate scheduled at 50-150 strains."

Sections:
- Per-variant AUROC + AUPRC table
- Label balance + cohort composition
- Top-K attribution genes for NT-XGBoost
- Smoke acceptance: PASS or FAIL (binary)
- Open questions surfaced (typically: which classical baseline performed best, how big the NT gap was, any anomalies)

### Step 4: File TODOS for decision gate

**Files:** `TODOS.md`

Add four `[OPEN]` entries per D6.

### Step 5: Plan the real decision gate (separate plan document, not in this one)

**Files:** `plans/Phase2_Decision_Gate_Real_Plan.md` (future plan to be saved when smoke passes; NOT written here).

This plan defines the smoke gate and the staging for the real gate; the real gate's plan will be written after the smoke gate runs. Expected contents:
- Cohort size choice (50 vs 150) based on compute availability
- Variant list (5 smoke variants + RF + TabPFN-with-PCA + SNP-table)
- Compute path (local-only vs Databricks burst)
- Acceptance bar definition (5%+ AUROC gap, calibration check, attribution sanity)

## Verification

- `scripts/pipeline.py:285` no longer uses `hash(s.mlst) % 10`; clade-only baseline produces real MLST-group-based predictions.
- 5-variant smoke gate runs end-to-end on `gate_b_mini_cohort.parquet`; result table generated.
- Smoke acceptance evaluated: PASS or FAIL recorded.
- `wiki/smoke_gate_12strain_cipro_*.md` exists with explicit "this is smoke, not decision" framing.
- `TODOS.md` contains the 4 new `[OPEN]` entries for the deferred /research-suggested work.

## Resolved Decisions (2026-05-13)

The three "Open Decisions for User" from the initial save have been resolved via /brainstorm critique + Codex web verification:

1. **Minimum N tiers:** N=12 = smoke, N=50 = engineering screen, N=150 = decision-grade. Locked into D5's two-stage gate.
2. **AMRFinderPlus point mutations:** Yes — emitted as first-class output rows with Method=`POINT*` when run with `--organism Escherichia --mutation_all`. SNP-table TODO scope collapsed from "days" to "hours" (parsing task, not pipeline build). Locked into D6.
3. **Cohort size 50 vs 150:** Both, staged. N=50 local first; N=150 Databricks only if N=50 passes the "any positive lift" screen. Locked into D5.

Actionable threshold (locked into D5):
- Stage 1 (N=50, local): NT-XGBoost shows any positive lift over clade-only-fixed → pass
- Stage 2 (N=150, Databricks): NT ≥ best classical + 5 pp AUROC AND top-10 NT attribution includes ≥1 of {gyrA, parC, parE} → pass
