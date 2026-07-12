# Genome world-model — creative data-reuse plan (2026-07-11)

**Mode:** `/soraya decompose & --plan` (strategy + plan-first; STOPS for user ratification — no build started).
**Author:** Soraya-J2. **Status:** Family **A EXECUTED** 2026-07-11 (user ratified option 1, flagship-first) →
verdict `FAIL_ADDITIVE_SUFFICES` (`wiki/hiv_epistasis_result_2026-07-11.{json,md}`; `scripts/hiv_epistasis.py`;
`tests/test_hiv_epistasis.py` 13/13; frozen surface `verify_lock` OK). Family **C EXECUTED** 2026-07-11
(user ratified option 2, depth+breadth) → verdict `PASS_LINKAGE_STRUCTURE` (raw 143/144; clonality-corrected
114/121 = 94.2%; `scripts/determinant_cooccurrence.py`; `tests/test_determinant_cooccurrence.py` 10/10;
`wiki/determinant_cooccurrence_result_2026-07-11{,_dedup}.{json,md}`).
Family **B EXECUTED** 2026-07-11 (deep dive) → `CALIBRATED_INTERVALS` (24/24 HIV drug-cells;
`scripts/hiv_quantitative_calibration.py`; 7 tests). Family **C-deep EXECUTED** → `PASS_CORESISTANCE_IMPUTABLE`
(class-level imputation AUC 0.9–0.97; virulence cross-axis unavailable from AMR-only runs — Bakta sweep stays
the Databricks follow-on; `scripts/coresistance_imputation.py`; 6 tests). Family **D EXECUTED** (risk-flagged)
→ `FLAG_RECOVERS_BLINDSPOT` (deterministic position-novelty flag catches 60% of the HIV EFV catalog blind spot
the learned rescue couldn't; `scripts/hiv_blindspot_position_novelty.py`; 7 tests). **ALL FOUR FAMILIES
EXECUTED.** Synthesis: additive suffices quantitatively (A), massive real joint linkage + imputation (C/C-deep),
the quantitative decoder is honestly calibrated (B), and a deterministic flag partially recovers the blind
spot the learned approach can't (D) — five honest world-model findings, all from data in hand, no new
labels/GPU/embeddings; frozen surface byte-unchanged throughout.

---

## 0 · The bigger picture (why the obvious moves are the wrong ones)

Every closed track in this project failed the SAME way: a **learned sequence model on a scarce label**, or a
**supervised decoder blocked by label circularity / censoring / sampling** (the 8 gates,
`wiki/negative_results_map_2026-06-13.md`). The verified thesis is *labels, not models*. Concretely:

- Foundation-model embeddings are **0-for-5** (learn lineage, not mechanism).
- ESM2 on AMR is **below chance** (antagonistic selection inverts likelihood models); ESM2 peaks at 650M.
- MIC-continuous from BV-BRC is **G1-circular + G6-censored**; pathotype is **label-blocked**; the provdisjoint
  grid is **saturated**. New independent labels need **paid acquisition** (user authority).

So "use our data creatively to improve the world model" must NOT re-run any of those. The trap is to reach for
a bigger GPU or another embedding — both are foreclosed. **The improvable frontier is the structure the
current additive catalog THROWS AWAY**, which we can model with data *already in hand* and no new labels:

1. **Interactions (epistasis).** The deployed rule is additive (count/OR). Reality has higher-order structure.
2. **Quantitative gradients.** The rule emits binary R/S. We hold *free, independent, continuous* labels (HIV
   fold-change, CRyPTIC MIC) that the binary output discards.
3. **Joint / linkage structure.** We run a *multi-cell* decoder (AMR × plasmid × virulence × serotype × MLST)
   over 744+ genomes — its own outputs encode a joint distribution we've never modelled.
4. **Bidirectionality (the "vice-versa").** Genotype→phenotype is inverted from phenotype→genotype only if we
   model the co-occurrence structure of determinants.

**Reframed north star:** the genome world model is a *mechanism-grounded model of the joint structure of
(determinants → interactions → quantitative effect → co-occurrence), bidirectional and self-aware of its
blind spots* — improved from the free continuous labels + our own decoder outputs we already possess, and
**explicitly not** from sequence embeddings (closed). This widens "world model" past the one definition that
kept failing (sequence→phenotype) to the definition the data can actually support.

---

## 1 · Decomposition — four families (screened against the 8 gates)

| Fam | Name | Substrate (IN HAND) | Improves the world model by | Gate screen |
|---|---|---|---|---|
| **A** ⭐ | **HIV epistasis + quantitative-effect model** | `data/raw/hiv/*` (2272/1867/2171/861 isolates, 94% multi-mutant, continuous fold-change) | modelling pairwise interactions the additive catalog misses + a calibrated continuous effect | HIV cleared all 8; **not** the embedding bet (mechanism-feature interaction model) |
| **B** | **Quantitative calibration layer** (R/S → calibrated continuous) | HIV fold + CoV-RDB fold + **CRyPTIC MIC** | honest continuous output + coverage-valid CI where a free continuous label exists | HIV/CRyPTIC free-continuous (bacterial BV-BRC MIC stays G1/G6-blocked → excluded) |
| **C** | **Joint co-occurrence / linkage world model + "vice-versa" inverter** | our OWN multi-cell decoder outputs (genome-map) over the 744+ committed accessions — self-distillation | modelling what-goes-with-what (co-resistance linkage, plasmid↔resistance) + P→G inversion | no new labels; associational (disclosed); novel |
| **D** | **Blind-spot / out-of-catalog self-awareness layer** | catalog-negative subsets + decoder residuals + the ΔΔG pocket finding (`hiv_blindspot_pocket_localization.py`) | flagging where the deployed rule is likely wrong | **RISK-FLAGGED:** the *likelihood/LM* blind-spot is a CLOSED negative on HIV — D must use POSITION-NOVELTY / selection, never a sequence LM |

### Flow-down + critical path
```
A (epistasis on HIV)  ──►  B (continuous calibration; needs A's interaction terms)
        │
        └─ residuals ──►  D (blind-spot; risk-flagged, partially closed — LAST)
C (joint co-occurrence)  ── independent ── runs in PARALLEL (different data/axis)
```
- **Critical path:** A → B. **Parallel:** C. **Optional/last (risk-flagged):** D.
- **Depth × breadth pairing (the recommended portfolio):** **A** = depth (quantitative + interaction on our
  best label) · **C** = breadth (joint structure across the whole multi-cell decoder). A proves the
  interaction-modelling thesis where data is richest; C makes it genome-wide.

### Honesty rails carried into every family
- **Frozen surface untouched.** All four are NEW non-frozen artifacts (scripts + wiki + tests). The frozen AMR
  decoder surface (`amr_rules.py` + `calibrated_amr_rules.json` + `shipped_decoder_surface.py`) stays
  byte-identical; prospective-lock `verify_lock` must stay green.
- **Beat the DOMAIN baseline, not a strawman.** A must beat the *strong additive OLS* already built
  (`hiv_nnrti_baseline.py`, the `hiv_*_mutant_catalog.py` deconfounded catalogs) — never a null.
- **Paired, held-out, nested-CV comparison** (per `feedback_paired_comparison_not_difference_of_medians`):
  interaction-vs-additive is a per-isolate paired delta on held-out folds with a bootstrap CI, never a
  difference of medians over different splits.
- **A negative is a shippable result.** "Additive suffices for HIV DR" is itself a world-model finding.

---

## 2 · `--plan` — Family A (flagship), ordered + optimized to a checkable bar

**Terminal claim:** a pairwise-interaction (epistasis) model of HIV drug-resistance fold-change **beats the
deployed additive mutant-specific catalog out-of-sample**, on a pre-registered majority of adequately-powered
drug-cells, WITHOUT overfitting — and emits a calibrated continuous fold-change with an honest CI.

**Pre-registered, DERIVED bar (frozen at run start, per R2 — thresholds derived, not asserted):**
- Nested 5×5 CV (outer = evaluation, inner = elastic-net λ selection). Metric = Spearman ρ (rank, matches
  Stanford `DRMcv.R`) AND R² on log10 fold-change.
- Interaction feature space is **constrained to control combinatorics**: pairwise terms only between positions
  that are (i) known DRM positions OR (ii) co-occur in ≥ K isolates (K derived from the per-class multi-mutant
  count); elastic-net (L1+L2) so unsupported pairs zero out. No triples in v0.
- **"Beat" = paired bootstrap (B=1000) on the outer-fold held-out predictions**, per drug: ΔSpearman
  (interaction − additive) 95%-CI lower bound > 0. Powered drug-cell = ≥ N_min isolates AND ≥ N_min in each
  fold class (N_min derived from the smallest reliably-scored existing HIV cell).
- **PASS** = CI-positive interaction gain on **≥ 50% of powered drug-cells** across NNRTI(5)+NRTI(6)+PI(8)+INI(5)
  (exact fraction frozen once the powered-cell count is computed in Step 2). **FAIL** = below that → bank the
  "additive suffices" negative with the per-cell table.

**Optimized step path (topo-ordered; independent steps batched):**

| # | Step | Class | Depends on |
|---|---|---|---|
| 1 | Load + normalize the 4 HIV DataSets into a shared `(isolate, {position:AA}, {drug:log10fold})` frame; reuse `hiv_rt_caller`/catalog parsers; unit-test the parse (censored-fold handling, `-`/`.` = WT) | edit-local-code | — |
| 2 | Compute the powered-cell census (isolates/drug, multi-mutant %, per-fold class balance) → **freeze N_min, K, the ≥50% fraction** into the pre-registration block; write it to the plan | run-tests | 1 |
| 3 | Re-fit the **additive baseline** in the SAME nested-CV harness (reuse the existing OLS/elastic-net catalog logic) → per-drug held-out ρ/R² = the number to beat | run-tests | 1,2 |
| 4 | Build the **constrained pairwise-interaction feature matrix** (DRM×DRM + co-occurring pairs ≥ K) + elastic-net fit in the SAME harness | edit-local-code | 2,3 |
| 5 | **Paired bootstrap** interaction-vs-additive on outer folds, per drug → CI table; **verify-in-batch:** inspect the top interaction terms — are they biologically real (known accessory/primary DRM pairs, e.g. NRTI TAMs, PI 82+54) or overfit noise? | run-tests | 3,4 |
| 6 | **Calibration head:** map the interaction model's output to a calibrated log10-fold with a held-out coverage check (does the 90% CI cover 90%?) | edit-local-code | 4,5 |
| 7 | Apply the frozen verdict (`compute` PASS/FAIL from the CI table) → write `wiki/hiv_epistasis_result_<date>.{md,json}` with the per-cell table + the top-interaction biology + the honest negative-or-positive headline | run-tests | 5,6 |
| 8 | Regression tests (parser, harness determinism, verdict function) + **assert frozen AMR surface + `verify_lock` unchanged** | run-tests | 7 |

**Gate/tool posture:** every step is `auto` (reads / in-cwd writes / local CPU tests — the HIV data is committed
and small; no Docker, no GPU, no network, no money). Commits to `main` = reversible-outward (push). No
genuinely-irreversible step. Estimated ~8–14 steps incl. recovery; well within the 100-step budget.

**What A explicitly does NOT do:** no sequence embedding, no ESM, no GPU, no new label. It is a mechanism-feature
statistical model on a label we already own — the exact niche that survives the regime boundary
(`feedback_g2p_decoder_regime_boundary`: molecular-property + *fitness-aligned continuous* label → learnable).

---

## 3 · `--plan` sketch — Family C (parallel breadth; lighter, more infra)

**Terminal claim:** a joint model over our multi-cell decoder outputs predicts a held-out functional feature
from the others better than its marginal prevalence, AND inverts to rank determinant sets for a target phenotype.

- **C1** Run the genome-map (Bakta + AMRFinder + VF + plasmid/serotype/MLST callers) over the 744+ committed
  accessions → a tidy `(genome × functional-feature)` matrix. *(Infra cost: Docker + hours; the one heavy part.)*
- **C2** Fit a joint model (e.g. regularized logistic / tree over the feature matrix); held-out predict each
  feature from the rest; compare to marginal-prevalence baseline with a paired CI. Falsifier = features
  conditionally independent (joint = marginal).
- **C3** Invert: given a target resistance, rank observed determinant sets by co-occurrence weight → the
  "vice-versa" P→G output; validate that the top-ranked set is the one actually seen in held-out R genomes.
- **Gate note:** self-distillation from our own decoder ⇒ associational, DISCLOSED (not a causal claim); no new
  label ⇒ no gate trip. The heavy step is C1 (Docker), so C lags A on wall-clock.

---

## 4 · Recommendation to the user (ratify one)

1. **Flagship-first (recommended):** ratify **Family A** and I execute its `--plan` (CPU-only, data in hand,
   ~1 session). Highest VOI, sharpest falsifier, directly extends shipped HIV cells, zero infra risk.
2. **Depth + breadth:** ratify **A + C** as the portfolio; A runs now, C's genome-map sweep runs in parallel
   (Docker).
3. **Full decomposition:** ratify all four; I `/project-init` a ledger per family (self-init cap OK) and
   sequence A→B, C parallel, D last (risk-flagged).

**On ratification I will:** freeze Family A's pre-registration (Step 2 derives the exact thresholds), then
execute — never before. No ledgers seeded, no build started, frozen surface untouched until you pick.
