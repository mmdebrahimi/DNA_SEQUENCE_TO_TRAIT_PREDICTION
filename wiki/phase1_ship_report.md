# Phase 1 Ship Report (GREEN)

**Date:** 2026-05-12
**Repo:** `mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`
**Tag:** `phase-1-shipped` (on the closeout commit chain)
**Verification gate:** GREEN — 287 passed, 0 failed, 1 skipped, pytest exit 0

## Waves shipped

**From `plans/Ecoli_G2P_Platform_Technical_Plan.md` (the 18-step technical plan):**

| Wave | Scope | Status |
|---|---|---|
| 0 | Project bootstrap (Step 1) | shipped |
| 1 | Data ingestion + foundation models + eval harness (Steps 0.5, 2, 3, 4, 5, 7, 10) | shipped |
| 1.5 | Hardening C1 + C2 + C3 + Tier 1-5 framework + Phase 2 backlog | shipped |
| 2 | Cohort catalog + HDF5 embedding cache (Steps 6, 8) | shipped |
| 2.5 | Hardening C4 + C5 + C6 (cohort + cache contract fixes) | shipped |
| 3 | Baseline classifiers + ISM + classical baselines (Steps 9, 11, 18) | shipped |
| 3.5 | Hardening C7 + C8 + M4 + M5 + use_label_encoder cleanup | shipped |

**From `plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md` (the 5+2-step ship-path contraction):**

| Step | Scope | Status |
|---|---|---|
| 15 | Smoke pipeline + synthetic fixtures | shipped |
| 14 | `scripts/pipeline.py` CLI with 4 subcommands | shipped |
| 13 | matplotlib viz + TSV export (pygenometracks → Phase 2) | shipped |
| 17 | Shell-loop leaderboard | shipped |
| 11.5 | Quantization-fidelity check | shipped |
| 16 | README + `docs/ARCHITECTURE.md` | shipped |
| 3.5 | Wave 3.5 hardening (also covered above) | shipped |

## Baseline test count + environment

```json
{
  "passed": 287,
  "failed": 0,
  "errors": 0,
  "skipped": 1,
  "pytest_exit_code": 0,
  "env": {"python": "3.11.5", "uv": "0.11.8", "platform": "win32"},
  "command": "uv run pytest tests/ -v"
}
```

Full record at `test_baseline.json` (committed).

**Closeout-time fix logged:** initial full-suite run returned 285 passed / 2 failed (smoke tests). Root cause: Wave 3.5's `MIN_TRAINING_SAMPLES = 10` production guard tripped on the smoke fixture's 4-strain cohort (LOSO trained each fold on N-1=3 samples). Fixed at the fixture layer (bumped to 12 strains in `scripts/smoke_pipeline.py`) — the production guard was not loosened. Also aligned CLI smoke threshold from 0.85 → 0.6 to match the test's "LOSO on small N is noisy" rationale, leaving 0.85 as the Phase 1 *real-data* gate documented in README.

## Doc updates

The prior session had already reconciled doc content but never committed. Closeout commit captured the reconciled state:

- `CLAUDE.md` — project-specific Claude guidance (architecture / common commands / gotchas / project workflow); reflects 18/18 steps shipped, no remaining "pending" markers.
- `TODOS.md` — archived Phase 1 done sections; preserved post-Phase-1 polish + Phase 2 entry criteria + Phase 2 real-data validation; added Phase 2 starter-genome trio (MG1655 / O157:H7 Sakai / ST131 EC958) with explicit clade-confound caveat.
- `README.md` — `Status: Phase 1 — code-complete (18/18 implementation steps + 3 hardening waves)`; CLI surface map; setup + quickstart + success criteria + optional cache routing.
- `LESSONS_LEARNED.md` — 10 retrospective bullets capturing the Phase-1 build pattern (inter-wave brainstorm gates, contract-gap categories, `/probe` calibration value, `/save-plan` after `/review`, annotation gene_id semantics, `.gitignore` anchor lesson, etc.).
- `FUTURE_FEATURES.md` — Phase 2/3/4+ capability backlog with `/probe` + `/brainstorm` synthesis (Attribution Refinement Engine; differentiable MLP + Captum IG; MIC regression head; pan-genome graph layer; multimodal long-term; AlphaFold-inspired arch with compute caveat).
- `docs/ARCHITECTURE.md` — one-page module map + data flow + key invariants.

## Archive status

- In-file status header added to `plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md` (commit `ae086b9`).
- Git tag `phase-1-shipped` attached to the ship-report commit + pushed to origin.
- `.claude/execute-plan-state/` added to `.gitignore`; the tracked `Ecoli_G2P_Platform_Technical_Plan.json` deleted in `aaed60e`.

## What Phase 1 ship gate validated

- All 18 implementation steps + 3 hardening waves compile + import + run from a clean `uv sync`.
- 287 unit tests pass against real Python deps (xgboost, sklearn, BioPython, h5py, pandas, numpy, scipy).
- `scripts/smoke_pipeline.py` runs end-to-end on synthetic 12-strain fixtures: cohort build → mean-pool embeddings → XGBoost LOSO-CV → ISM attribution → Tier 1-5 classification → markdown report. Top-1 attribution recovers the seeded signal gene (`g1` / `TAG_001`). AUROC well above the 0.6 noise floor.
- The 7 CLI entry points (`pilot_gate`, `pipeline {ingest, train, predict, attribute}`, `smoke_pipeline`, `leaderboard`, `quantize_fidelity_check`) have argparse + exit-code contracts unit-tested.
- The 2-phase production training guard (`MIN_TRAINING_SAMPLES = 10` + minority-class CV folds) survived the smoke regression without being loosened — the fix went to the fixture, not the guard.

## What Phase 1 ship gate did NOT validate

- **Prediction accuracy on real E. coli genomes.** No real data was ingested. The fluoroquinolone / β-lactam / tetracycline classifiers have never seen a real CARD / AMRFinder catalog entry or a real BV-BRC AST label. The Phase 1 AUROC targets (≥0.80 SLO / ≥0.85 stretch) and Tier 1-3 fractions (cipro ≥40%, ceftriaxone ≥25%, tet ≥30%) remain hypotheses.
- **Region-attribution biological correctness.** ISM Tier 1-3 hits against real resistance loci (gyrA, parC, blaCTX-M, blaSHV, tetA, etc.) have not been measured. Only synthetic-fixture top-1 attribution has been validated.
- **Multi-strain leaderboard performance.** `scripts/leaderboard.py` runs the shell-loop fan-out, but no real Evo / DNABERT-2 / Nucleotide Transformer / GENA-LM weights have been downloaded or compared.
- **Live BV-BRC AST + NCBI Datasets API integration.** `pilot.fetch_bvbrc_drug_counts` + `fetch_ncbi_assembly_quality` raise `NotImplementedError` without a local TSV / env var / config entry. Live REST resolution is deferred.
- **Mash CLI phylogeny clustering.** External binary dependency; Linux/WSL needs `apt install mash`, Windows needs manual binary install. Untested on this machine.
- **4-bit Evo quantization.** `bitsandbytes` is Linux/CUDA only; Windows users need WSL2. Untested.
- **`motif_recovery` is a Phase 2 placeholder** — currently returns the same high-impact position list for every motif name and warns on call.

## Phase 2 entry criteria

Per `TODOS.md`:

- **Real-data smoke:** one E. coli genome end-to-end with prediction + attribution captured (no ground-truth comparison required at entry).
- **Full test suite passes on real deps:** `uv run pytest tests/ -v` — 287 tests as of 2026-05-12. *Already passing at Phase 1 ship.*

Starter genome trio for the real-data smoke (infrastructure-only — not a model-quality test):
- K-12 MG1655 (`GCF_000005845.2`) — reference baseline
- O157:H7 Sakai — pathogenic comparison
- ST131 EC958 — multidrug-resistant clinical

Real model-quality validation requires within-lineage R/S pairs (or clade-stratified cohort), not the cross-lineage trio above — the trio is for ingestion / annotation / embedding-cache plumbing only. Cross-lineage R/S contrast confounds clade signature with resistance signal, which is exactly what the clade-only-baseline gate (`dna_decode/eval/clade_baseline.py:validation_gate`) is designed to detect.

## Open follow-ups (logged in TODOS / FUTURE_FEATURES)

- Pre-existing CLI/test smoke threshold drift was resolved this closeout; no remaining drift.
- `three_drug_intersection_target` is still enforced as a minimum, not a cap. Single-drug smoke special-case noted.
- Plan-text drift: technical plan + ship-path plan describe smoke as "5-strain corpus" while implementation uses 12. Historical text preserved; the implementation-of-record is the committed code + tests.
