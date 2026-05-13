# E. coli G2P — Phase 1 Ship-Path Plan

> **Status:** archived 2026-05-12 — Phase 1 shipped (all ship-path steps committed; full test suite green 287/0/1; `phase-1-shipped` tag created on the closeout commit).

> Contracted path to ship Phase 1 of `Ecoli_G2P_Platform_Technical_Plan.md`. Captures the `/review` synthesis verdict (HOLD scope + selective contraction within remaining steps) plus the deferred Wave 3.5 hardening fixes from the post-Wave-3 `/brainstorm`. Estimated remaining work: ~700 LOC across 5 implementation steps + 4 hardening edits.

---

## Problem Statement

13 of 18 implementation steps in `Ecoli_G2P_Platform_Technical_Plan.md` are complete (Waves 0 + 1 + 1.5 + 2 + 2.5 + 3). Five steps remain across Waves 4-7: Step 13 (viz), Step 14 (CLI), Step 15 (smoke + fixtures), Step 17 (leaderboard), Step 16 (docs).

Two parallel concerns:

1. **Wave 3 introduced 4 contract gaps** (post-Wave-3 `/brainstorm`): C7 calibration CV uses majority class; C8 attribution-report uses `locus_tag` first which never matches catalog gene symbols; M4 `motif_recovery` is a semantically non-discriminating placeholder; M5 `train_kmer_baseline` takes one concatenated string per strain (loses plasmid β-lactamase signal).

2. **The plan as written over-engineers remaining steps relative to Phase 1 verification gates** (`/review` synthesis): Step 17 dedicated module duplicates Step 14 `train` × N models; Step 13's pygenometracks dependency is heavy for what is a "nice-to-have" verification gate (#12); Step 14's 4-script decomposition could collapse to one `scripts/pipeline.py` with subcommands; Step 16's `HOW_TO_ADD_ORGANISM.md` is a Phase 2 concern.

This plan does NOT change the destination — Phase 1 success criteria stay identical (Tier 1-3 ≥40% cipro, ≥25% ceftriaxone, ≥30% tet, ≤20% Fail; clade-only baseline gap ≥0.10 on ≥75% of held-out clades; ≥3pp gap vs best classical baseline on ≥2 of 3 drugs). It only contracts the *path* and resolves the deferred contract gaps.

## Design Decisions

### D1: HOLD scope, do not expand

**Decision:** No new features added to Phase 1. The plan's "if classical baselines win on ≥2 drugs, Phase 2 must redesign" honesty stance is preserved.

**Rationale:** Both `/review` lenses agreed scope is dialed in. Multiple hardening passes have already removed cruft. Adding more before shipping invites further delay.

**Trade-off:** Selective addition rejected: tier rubric "Tier 1/2" targets may be unreachable in Phase 1 by construction (codon-level matching is Phase 2 work). Acknowledged but not gated on — Tier 3 is the realistic Phase 1 ceiling.

### D2: Reorder — Step 15 (smoke + fixtures) BEFORE Step 14 (CLI)

**Decision:** Write the smoke pipeline + synthetic 5-strain fixture FIRST, then build the CLI on top.

**Rationale:** The smoke pipeline exercises the full `cohort → cache → classifier → cv → mutagenesis → tier_classify` chain. Discovering integration bugs cheaply via mock data is faster than discovering them via real-data CLI runs. `MockFoundationModel` substrate is already in place (Step 7).

**Trade-off:** Plan's original order (14 → 15) would have shipped the CLI sooner but pushed integration discovery later.

### D3: Step 14 collapses to one `scripts/pipeline.py` with subcommands

**Decision:** Replace the 4 separate scripts (`ingest.py` / `train.py` / `predict.py` / `attribute.py`) with one `scripts/pipeline.py` exposing subcommands `ingest | train | predict | attribute`.

**Rationale:** `ingest` + `train` share ~80% setup (load config, load cohort, load cache). Argparse + exit-code pattern from `scripts/pilot_gate.py` is the template. ~250 LOC vs ~600 LOC for the original 4-script layout.

**Trade-off:** Individual scripts are simpler to reason about in isolation; subcommands add minimal dispatch overhead but reduce duplication.

### D4: Step 13 visualization uses matplotlib + TSV export, NOT pygenometracks

**Decision:** Ship a matplotlib line plot of `|prediction_delta|` per position + TSV export of `GeneEffectTable` / `PositionEffectTable`. Defer pygenometracks adapter to Phase 2.

**Rationale:** pygenometracks = external CLI binary + `.ini` config plumbing + subprocess fragility. Verification #12 (visual render) is "nice-to-have", not a Phase 1 ship gate. matplotlib path = ~40 LOC, no binary dep, works on Windows out of the box.

**Trade-off:** pygenometracks would produce publication-grade genome-browser figures; matplotlib gives engineer-grade plots only. Sufficient for Phase 1 introspection, insufficient for SME presentations.

### D5: Step 17 leaderboard collapses to a shell loop over `pipeline.py train`

**Decision:** Replace `dna_decode/eval/leaderboard.py` + `scripts/leaderboard.py` + dedicated tests with a ~30-line shell loop calling `pipeline.py train --model $model` for each model + a markdown table writer.

**Rationale:** The plan itself already softens at the original Step 17 ("Phase 1 ships with leaderboard initially run for Evo + DNABERT-2 only"). Dedicated module is over-engineering for what is fundamentally a fan-out + report.

**Trade-off:** Dedicated module would be cleaner for Phase 2+ when 4 models × 3 drugs × multiple cohorts compound. Phase 1 doesn't need that scale.

### D6: Step 16 docs trimmed to README + ARCHITECTURE.md only

**Decision:** Drop `HOW_TO_ADD_ORGANISM.md` from Phase 1 documentation. Keep README quickstart finalization + 1-page `ARCHITECTURE.md`.

**Rationale:** Multi-organism extension is Phase 2 work; documenting how to add organisms before Phase 2 exists is premature. README + ARCHITECTURE cover Phase 1 ship needs.

**Trade-off:** Future contributors will need to read code + plan history; not a blocker for a solo Phase 1 ship.

### D7: Apply Wave 3.5 hardening BEFORE Step 14 wiring fires

**Decision:** Resolve C7 / C8 / M4 / M5 + `use_label_encoder` cleanup before the CLI integrates Wave 1-3 modules.

**Rationale:** Step 14 `pipeline.py train` calls `train_xgboost_classifier` (C7 bug bites on rare-resistance drugs); `pipeline.py attribute` calls `build_attribution_report` (C8 bug undercounts Tier 1-3 hits on Bakta-annotated genomes). Wave 5 building on broken contracts repeats the Wave 1 / Wave 2 pattern.

**Trade-off:** ~115 LOC of fix-up work delays Step 14 by ~20 minutes; not fixing pushes the bugs into the integrated pipeline where they're harder to isolate.

### D8: Add quantization-fidelity micro-step (selective addition)

**Decision:** Add a 5-10-strain 4-bit vs full-precision ISM concordance check before declaring Phase 1 attribution-precision numbers final.

**Rationale:** Plan's compute default is 4-bit Evo via `bitsandbytes` on RTX 4090. If quantization distorts attribution maps materially, every Phase 1 attribution-precision number is quantization-conditional and the headline results become unreliable.

**Trade-off:** Adds one rented-A100 hour to Phase 1 budget. Cheap insurance vs publishing quantization-conditional results.

## Implementation Plan

Ordered by dependency + revised priority. Estimates are remaining-work LOC.

### Step 3.5: Wave 3.5 hardening pass (BEFORE Step 14)
Files: `dna_decode/models/classifiers.py`, `dna_decode/models/classical_baselines.py`, `dna_decode/interp/mutagenesis.py`, `tests/test_models_classifiers.py`, `tests/test_models_classical_baselines.py`, `tests/test_interp_mutagenesis.py`
Estimated: ~115 LOC

- **C7:** calibration CV uses minority-class count. `cv_folds = max(2, min(3, int(min((y==1).sum(), (y==0).sum()))))`. Skip calibration + warn when minority < 2.
- **C8:** `build_attribution_report` walks BOTH `gene_id` and `locus_tag` via `_best_tier_across_candidates(catalog, drug, *candidates)`; returns the lowest positive tier.
- **M4:** `motif_recovery` docstring marked PLACEHOLDER; emits `UserWarning` at first call so silent use is loud.
- **M5:** `train_kmer_baseline` accepts `dict[str, str | list[str]]`; concatenates list internally before k-mer extraction.
- **Cleanup:** remove `use_label_encoder=False` kwarg from `XGBClassifier` construction (xgboost 2.0+ ignores it but it's stale).

### Step 15: Smoke pipeline + fixtures (FIRST in remaining-work order)
Files: `scripts/smoke_pipeline.py`, `dna_decode/tests/fixtures/ecoli_mini/{genome.fna,annotations.gff3}`, `tests/test_smoke.py`, optional inline-string fixtures for `card_mini.json` / `amrfinder_mini.tsv` / `bvbrc_ast_mini.tsv` (reuse existing inline-fixture pattern from `tests/test_data_*.py`)
Estimated: ~250 LOC + 1 synthetic 5-strain corpus

- Synthetic 5-strain `ecoli_mini/` with 10 genes/strain; one gene carries a 50-bp seeded resistance signal; build labels accordingly.
- `scripts/smoke_pipeline.py`: end-to-end via `MockFoundationModel`. Steps: ingest fixture → annotate → build cohort → populate cache → train XGBoost → leave-one-out CV → ISM → tier classification → write `data/processed/smoke_report.md`. <60s on CPU.
- `tests/test_smoke.py`: invokes `smoke_pipeline.py`, asserts exit 0 + report exists + AUROC ≥0.85 on the seeded signal + attribution top-1 = the seeded gene.

### Step 14: Single `scripts/pipeline.py` with subcommands
Files: `scripts/pipeline.py`, `tests/test_pipeline_cli.py`
Estimated: ~250 LOC

- Subcommands: `ingest | train | predict | attribute`. Argparse subparsers; `pilot_gate.py` style.
- Shared setup helpers: `_load_config`, `_load_cohort`, `_load_cache`, `_resolve_paths`.
- `train` subcommand also runs Step 10's clade-only baseline + `validation_gate` per Verification #6.
- Exit codes: 0 = success, 1 = run-failure (e.g., AUROC below threshold flagged), 2 = config / IO error, 3 = missing dependency (e.g., Mash binary).
- Tests: mocked deps; exit-code coverage per `pilot_gate` pattern; verify `train --include-clade-baseline` runs the gate.

### Step 13: matplotlib visualization + TSV export
Files: `dna_decode/viz/browser.py`, `tests/test_viz_browser.py`
Estimated: ~80 LOC (down from ~200 with pygenometracks)

- `export_attribution_tsv(gene_effects, position_effects, output_path)`: dump `GeneEffectTable` + `PositionEffectTable` as side-by-side TSVs for downstream tooling.
- `render_attribution_plot(position_effects, output_path)`: matplotlib line plot of `|prediction_delta|` vs position; one subplot per gene in the top-K. Saves PNG.
- Tests: TSV column-count + matplotlib `Figure` smoke (no pygenometracks subprocess).
- pygenometracks adapter deferred to Phase 2 (add as separate `viz/pygenometracks_adapter.py` module when SME presentation work begins).

### Step 17: Shell-loop leaderboard
Files: `scripts/leaderboard.py`, `tests/test_leaderboard.py`
Estimated: ~80 LOC

- `scripts/leaderboard.py`: loops over a configurable model list (default `["evo", "dnabert2"]`; NT + GENA-LM added incrementally per the plan's softened position), calls `pipeline.py train` for each, aggregates the per-model metrics into a single markdown report at `data/processed/leaderboard.md`.
- Includes per-drug AUROC + AUPRC + clade-only-baseline gap + best-classical-baseline gap + attribution top-K precision.
- Tests: mocked `pipeline.py train` subprocesses; verify report structure + columns.

### Step 11.5: Quantization-fidelity check (selective addition per D8)
Files: `scripts/quantize_fidelity_check.py`
Estimated: ~50 LOC

- Run ISM on 5-10 strains at full precision (rented A100) and 4-bit (RTX 4090).
- Compute concordance: top-K=20 set intersection per (strain, drug) + Spearman rank correlation on prediction-delta values.
- Output `data/processed/quantize_fidelity_report.md` with concordance metrics + GO/NO-GO for using 4-bit attribution as the headline Phase 1 number.
- One-time cost (~1 A100 hour), runs after Step 14 training but before final Phase 1 results.

### Step 16: Trimmed documentation
Files: `README.md` (UPDATE), `docs/ARCHITECTURE.md`
Estimated: ~100 LOC of docs

- README: finalize quickstart with all subcommands of `scripts/pipeline.py`; add results-summary section reading from `data/processed/leaderboard.md`.
- `docs/ARCHITECTURE.md`: 1-page diagram of the data → embedding → classifier → ISM → tier → viz flow; module-by-module pointer with file paths.
- Drop `HOW_TO_ADD_ORGANISM.md` from Phase 1; mark as Phase 2 doc.

## Verification

Phase 1 success criteria UNCHANGED from `Ecoli_G2P_Platform_Technical_Plan.md`:

1. `uv run pytest tests/ -v` — all unit tests pass.
2. `uv run python -m scripts.pilot_gate --ast-tsv <bvbrc_ast_tsv>` — GO verdict.
3. `uv run python scripts/smoke_pipeline.py` — synthetic-fixture end-to-end completes in <60s with AUROC ≥0.85 on seeded signal.
4. `uv run python -m scripts.pipeline ingest --drugs ciprofloxacin,ceftriaxone,tetracycline` — real ingestion completes in <4 hours on home internet.
5. `uv run python -m scripts.pipeline train --drug ciprofloxacin --model evo` — LOMO-clade CV AUROC ≥0.80 SLO / ≥0.85 target.
6. `uv run python -m scripts.pipeline train --drug ciprofloxacin --model evo --include-clade-baseline` — embedding model AUROC ≥0.10 above clade-only baseline on ≥75% of held-out clades.
7. `uv run python -m scripts.pipeline attribute --strain-id <accession> --drug ciprofloxacin` — top-K=20 attribution hits with ≥40% Tier 1-3 (cipro), ≥25% (ceftriaxone), ≥30% (tet); ≤20% Fail-tier across all drugs.
8. `uv run python scripts/leaderboard.py --drugs ciprofloxacin,ceftriaxone,tetracycline --models evo,dnabert2` — best foundation model beats best classical baseline by ≥3pp AUROC on ≥2 of 3 drugs.
9. `uv run python scripts/quantize_fidelity_check.py` — 4-bit vs full-precision ISM concordance ≥0.7 Spearman on top-K=20 attribution-delta rank.

**Phase 1 ships when** verifications 1, 2, 3, 4, 5, 6, 7, 8 pass. Verification 9 may flag the headline numbers as "quantization-conditional" if concordance is below threshold; that's a documentation caveat, not a ship blocker.

## Open Clarifying Questions

1. **If Phase 1 ships with classical baselines beating embeddings on ≥2 drugs**, what's the project — pivot, or is "interpretability-via-ISM-on-foundation-models" still the deliverable independent of the AUROC gap?
2. **Is publication-grade visualization a Phase 1 requirement** (→ pygenometracks adapter back in scope), or is "engineer reads TSV / glances at matplotlib plot" sufficient (→ D4 stands)?
3. **Are NT + GENA-LM Phase 1 ship-blockers, or "Phase 1.5 / when GPU time permits"?** D5 assumes the latter; confirm.

Open questions do not block plan adoption; they refine the verification thresholds + Step 13 / Step 17 scope if answers change.
