# Reproducibility freeze — 2026-07-04 (refresh of 2026-06-13)

Consolidation snapshot of the validated decoder suite as of the 2026-07-04 consolidation pass. Supersedes
`wiki/reproducibility_freeze_2026-06-13.md` (which predates ~40 commits of viral-de-confounding + visible-trait
+ hybrid-imputation work). The **frozen AMR surface is byte-unchanged** since the 2026-06-13 freeze
(`dna_decode/eval/amr_rules.py` + `dna_decode/data/calibrated_amr_rules.json` last touched at `026488d`; the
`tests/test_tb_leak_guard.py` 9/9 guard has held across every commit since). This refresh does NOT re-freeze
the deployed rule surface — it records the ADDITIONAL validated artifacts layered on top of it.

## The frozen, deployed core (unchanged since 2026-06-13)
- `dna_decode/eval/amr_rules.py` + `dna_decode/data/calibrated_amr_rules.json` — the deployed bacterial AMR
  rule surface (6-drug DRUG_RULE + EXPRESSION_FLOOR abstainers). Byte-identical to 2026-06-13.
- `dna_decode/data/shipped_decoder_surface.py` — the deployed-claim surface.
- `dna_decode/data/mic_tiers.py` — per-drug breakpoints + tier classifier.
- The standing report cards: `wiki/decoder_validation_report_card.{md,json}` (bacterial provdisjoint,
  10 SCORED / clonality-disclosed) + `wiki/hiv_decoder_report_card.{md,json}` (25 cells, in-distribution).
- The 10 `wiki/provenance_disjoint_validation_*.json` + `wiki/external_validation_*` + lineage sidecars.

## Validated additions since 2026-06-13 (this refresh records them)
- **HIV within-subtype de-confounding — ALL 4 classes** (2026-07-03). NNRTI/PI/INSTI extended (NRTI held
  earlier): median within-B AUC 0.795 / 0.921 / 0.898; pooled−within-B ≈ 0 → the catalogs decode MECHANISM,
  not subtype structure. `scripts/hiv_within_subtype.py`; on the HIV report card (subtype-transfer column).
- **HIV v0.2 absolute-cutoff calibration** (2026-07-03). PI 8/8 + NNRTI 4/5 calibrated at DRMcv.R-SOURCED
  per-drug cutoffs; INSTI cutoff-walled (reported, not guessed). `scripts/hiv_absolute_cutoff_validate.py`.
- **Visible/organismal-trait breadth** (2026-07-03). 3 new deterministic single-locus openSNP cells —
  lactase (LCT), earwax (ABCC11), cilantro (OR6A2) — a calibrated spread (PILOT WIN / LABEL-LIMITED /
  CORRECTLY-FAILS-by-design). `dna_decode/data/single_snp_traits.py`.
- **The hybrid learned+deterministic arc — CLOSED** (`plans/Hybrid_Learned_Deterministic_Decoder_Plan.md`):
  V1 learned scoring (zero-shot + supervised head) beats zero-shot but LOSES to the curated catalog on
  held-out positions → fallback-only; V2 masked-genotype IMPUTATION PASSES (98.9% abstain-reduction) → the
  deployable win. Verdict: learned phenotype-prediction does not beat curated knowledge; learned
  input-completion does.
- **NEW deployable capability — the imputation pre-processor** (`dna_decode/imputation.py`, 2026-07-04):
  fail-closed (impute an uncallable determinant from a committed LD map only at tag-purity ≥ 0.90, else
  ABSTAIN), provenance-tagged (direct / imputed / abstain), NEVER touches the frozen surface. Committed real
  map `data/imputation/abo_rs8176719_from_rs657152.json` (LOO acc 0.985). Ready to plug in wherever a
  determinant with an LD map is added to a decode path.

## The standing conclusions (unchanged, now fully evidenced)
- **The deterministic decoder is the product.** Every learned-model expansion is closed with a recorded
  negative: embeddings 0-for-5 de-confounded; learned scoring fallback-only; molecular conservation competes
  on FUNCTION but resistance ⊥ conservation for pocket-evasion mechanisms.
- **The binding constraint is LABELS, not models** — exhaustively established.
- **Two non-foreclosed forward paths, both USER-gated:** (1) acquire a non-public/wet-lab/clinical label
  source (clears the wall by construction); (2) prospective-lock accrual (`scripts/prospective_lock_validate.py`,
  passive, time-gated). See `wiki/negative_results_map_2026-06-13.md` (the 8 reusable rejection gates — still
  current).

## Rerunnability
- Frozen-surface guard: `uv run pytest tests/test_tb_leak_guard.py -q` → 9/9 (byte-identity of the deployed
  rule surface).
- Full offline suite: **1939 passed, 7 skipped** as of this refresh (excluding
  `tests/test_models_foundation.py`, a host torch-paging limit). The only non-green entries are 10 ERRORs in
  `test_pipeline_predict_e2e.py` + `test_pipeline_predict_genome_input_e2e.py` — a pre-existing ENVIRONMENT
  gap (`xgboost` not installed in this venv; the fixtures call `train_xgboost_classifier`), NOT a logic
  failure. `uv sync` (or `uv pip install xgboost`) clears them.
- One-command HIV/bacterial report-card rebuild: `uv run python scripts/build_hiv_report_card.py` +
  `scripts/build_validation_report_card.py` (read-only roll-ups).
