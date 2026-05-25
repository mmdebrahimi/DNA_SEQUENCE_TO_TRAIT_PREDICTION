# DNA_SEQUENCE_TO_TRAIT_PREDICTION

Genotype-to-phenotype (G2P) inference platform — predicts phenotypic traits from genomic DNA sequences AND identifies which genomic regions are most strongly associated with those predictions. Biologically interpretable, not causal-claim-making.

## Status: Phase 1 — CLOSED 2026-05-17 (infrastructure + cross-drug architectural finding)

Phase 1 evidence collection closed 2026-05-17. Cross-drug architectural finding synthesis at `wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md`:

> At 12-strain smoke fidelity, frozen-NT-whole-genome-pooling PASSES on concentrated-signal AMR mechanisms (cipro QRDR point mutations: AUROC 0.750; cef plasmid acquired-gene β-lactamases: AUROC 0.833) AND FAILS on distributed mobile-element mechanisms (tet tet-family efflux + ribosomal protection: AUROC 0.400, anti-predictive). The architecture's failure mode appears mechanism-class-bounded, largely independent of drug identity at smoke fidelity.

EP1 cipro closed internally (`wiki/cipro_ep1_closeout_2026-05-17.md`) with a 4-tier adversarial audit infrastructure (mechanism × MIC × opacity merge with structurally-enforced SUSPEND gate). EP2 cef + tet smoke fired (cef PASS, tet FAIL, H17 falsified). No Databricks burst spent. External publication deferred per PC1=`internal_closeout`.

Phase 1 code: all 18 implementation steps shipped Wave 0-7 (2026-05-11 → 2026-05-12) + 3 hardening waves; cross-drug Evidence Packet evidence collection completed 2026-05-17 per the Evidence Packets framing reset 2026-05-15. **Phase 2 entry fired 2026-05-18**: BV-BRC strict-MIC 4-drug feasibility census ran (`scripts/bvbrc_strict_mic_4drug_census.py` + `wiki/bvbrc_strict_mic_4drug_census_2026-05-18.{md,json}`) — NO drug clears N=150 per-class at either strict-MIC or relaxed-MIC bars; structural bottleneck is `assembly_accession`. North star clarified: AI DNA decoder tool, not papers. v0 UX + success criteria LOCKED at `wiki/decoder_v0_ux_and_success_criterion.md` (CLI via `pipeline.py predict`, LOSO AUROC ≥ 0.70, cipro v0 / cef v0.1, JSON + markdown sidecar). 3 of 5 v0 criteria green via 24 new tests; 2 gated on Databricks N=147 cipro cache landing.

**Phase 2 in-flight (2026-05-22 → 2026-05-24)**: cipro interpretability audit completed on Precision 7780 (RTX 3500 Ada) by parallel Codex CLI session. Bounded-falsifier coordination plan + post-falsifier ship-path technical plan covered all 4 verdict branches × 3 gate states pre-committed. Codex on Precision 7780 ran the falsifier 2026-05-23 — **verdict = FAIL** (ranking-only rescue did not improve the ELX-family failure case on 12-strain Bucket B). Per the FAIL branch + north star, **v0 shipped 2026-05-24** as a cached-strain cipro predictor (`scripts/pipeline.py predict --strain-id ...`) with a documented scope-limit. v0 spec RELOCKED at `wiki/decoder_v0_ux_and_success_criterion.md` to match the implemented cached-strain surface (not the original genome-input decoder concept). Leakage-safe retrain on `leave_one_accession_out` CV yielded **AUROC 0.8697**. v0 closeout handoff: `wiki/dna_decoder_v0_closeout_handoff_2026-05-24.md`.

**v0.1 planning in flight (2026-05-24)**: ingestion contract plan covers BOTH candidate paths — real-genome-input cipro decode (`plans/v0.1_Ingestion_Contract_Plan.md` Path G) AND cef-cached expansion (Path C). Drug-agnostic mechanism audit (`scripts/drug_mechanism_audit.py`) enables Path C; cross-machine sync diagnostic (`scripts/cross_machine_sync_check.py`) prevents the divergence pattern that emerged during v0 closeout. 685 tests green (+316 vs Phase 2 entry).

See `plans/EP1_EP2_Cross_Drug_Synthesis_Plan.md` for the synthesis plan; `plans/Cipro_Decision_Bundle_Plan.md` + `plans/Cipro_Decision_Bundle_Technical_Plan.md` for the EP1 closeout planning chain; `plans/EP2_Cef_Tet_Smoke_Design_Plan.md` for the EP2 design. See `plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md` for the original Phase 1 contracted ship-path. See `plans/Ecoli_G2P_Platform_Technical_Plan.md` for the full Phase 1 plan with Tier 1-5 attribution-success framework. See `docs/ARCHITECTURE.md` for the module map.

What runs end-to-end today:

| Surface | Entry point | Notes |
|---|---|---|
| **Pilot gate (HARD)** | `python -m scripts.pilot_gate --ast-tsv <path>` | Validates per-drug strain counts before ingestion fires. Exit 0=GO, 1=NO-GO, 2=PilotGateError, 3=no source. |
| **Full pipeline** | `python -m scripts.pipeline {ingest, train, predict, attribute}` | Single CLI with 4 subcommands; shared config-driven path resolution. |
| **Smoke regression** | `python scripts/smoke_pipeline.py` | <60s synthetic-fixture end-to-end via MockFoundationModel; asserts AUROC ≥0.85 + top-1 attribution = seeded gene. |
| **Leaderboard fan-out** | `python scripts/leaderboard.py --drugs ... --models evo,dnabert2` | Loops pipeline.py train per (model × drug); writes `data/processed/leaderboard.md`. |
| **Quant-fidelity check** | `python scripts/quantize_fidelity_check.py --full-precision-attributions <manifest.json> --quantized-attributions <manifest.json>` | One-time 4-bit vs full-precision ISM concordance check; gates whether Phase 1 attribution numbers are quantization-conditional. |
| **Viz** | `dna_decode.viz.browser.render_attribution_plot` + `export_attribution_tsv` | matplotlib PNG + TSV export; pygenometracks deferred to Phase 2. |
| **BV-BRC strict-MIC 4-drug feasibility census** | `python -m scripts.bvbrc_strict_mic_4drug_census` | Phase 2 entry (2026-05-18). Per-drug feasibility at strict + relaxed bars for cipro/cef/tet/gent. Writes `wiki/bvbrc_strict_mic_4drug_census_<date>.{md,json}`. Imports from `dna_decode/data/mic_tiers.py` (shared per-drug catalogs). |
| **v0 decoder predict** | `python -m scripts.pipeline predict --strain-id X --model-path M.pkl --cache C.h5 --annotations G.gff3 --audit-merge-json A.json --output Y.json` | v0 schema per `wiki/decoder_v0_ux_and_success_criterion.md` (2026-05-18 LOCKED). Emits JSON + markdown sidecar with prediction + calibrated_probability + confidence_tier + top_k_attribution + audit_verdict (SUSPEND propagation) + provenance. |

Module map: `dna_decode/data/` (ingestion) + `dna_decode/models/` (foundation wrappers + cache + classifiers + classical baselines; `cache.verify_complete` integrity gate added 2026-05-15) + `dna_decode/interp/` (ISM + Tier 1-5 attribution) + `dna_decode/eval/` (CV + metrics + batched-call phylogeny + clade-only baseline) + `dna_decode/viz/` (browser) + `tools/` (Stage 2 bioinformatics-tool runner via Docker Desktop — Mash + AMRFinderPlus + Bakta).

## Phase 1 scope

| Aspect | Value |
|---|---|
| Organism | E. coli |
| Phenotypes | Ciprofloxacin, ceftriaxone, tetracycline binary resistance |
| Foundation models | Evo (primary), DNABERT-2, Nucleotide Transformer, GENA-LM (leaderboard) |
| Classical baselines | AMRFinder gene calls, k-mer logreg + XGBoost, gene-presence XGBoost (Step 18) |
| Baseline ML | Frozen foundation-model embeddings + XGBoost per drug |
| Attribution | In-silico mutagenesis (gene-level + nucleotide-level saturation) |
| CV | Leave-one-Mash-clade-out + clade-only baseline + per-clade reporting |
| Target | AUROC ≥0.80 SLO / ≥0.85 stretch; clade-baseline-gap ≥0.10 on ≥75% of held-out clades; ≥3pp gap vs best classical baseline on ≥2 of 3 drugs |
| Horizon | 3 months Phase 1; 12 months Phase 1+2+3 |
| Compute | Local GTX 860M (4 GiB Maxwell, NT v2 only — verified 2026-05-13) + Databricks burst for larger cohorts. 4-bit Evo unavailable (bitsandbytes requires CC ≥ 7.0). Original target was RTX 4090 + 4-bit Evo; never materialized. |

## Long-term vision

Multimodal genotype-phenotype platform — start with bacterial AMR (Phase 1), expand toward eukaryotes + image-paired phenotype data in later phases. NOT a direct stepping stone to "DNA → animal image" prediction; that would require a parallel multimodal track.

## Setup

```bash
# 1. Install uv (if not already on PATH)
#    Windows PowerShell:  irm https://astral.sh/uv/install.ps1 | iex
#    Linux/macOS:         curl -LsSf https://astral.sh/uv/install.sh | sh

# 2. Sync deps. pytest is in default deps (Wave 1.5 hardening fix):
uv sync

# 3. Run the test suite
uv run pytest tests/ -v

# 4. Optional: install dev tooling (ruff + pytest-cov)
uv sync --extra dev

# 5. (Advanced, gated on hardware) install bitsandbytes for 4-bit Evo quantization
#    Requires CC ≥ 7.0 GPU (Ampere / Ada / Hopper). NOT compatible with the project's
#    actual GTX 860M (CC=5.0). Skip unless running on A100+ or similar.
uv sync --extra quantize
```

## Phase 1 quickstart

End-to-end Phase 1 run (assumes `uv sync` + BV-BRC AST TSV downloaded + Mash CLI installed):

```bash
# 1. HARD gate: confirm you have enough labeled strains per drug
uv run python -m scripts.pilot_gate \
  --drugs ciprofloxacin,ceftriaxone,tetracycline \
  --target-per-drug 150 \
  --ast-tsv path/to/bvbrc_ast.tsv

# 2. Smoke test: <60s end-to-end on synthetic fixtures (sanity check before real run)
uv run python scripts/smoke_pipeline.py

# 3. Ingest: build cohort + download cohort genomes
#    For real-data runs, pass --assembly-metadata-csv pointing at the BV-BRC
#    Genomes-tab export (CSV). The adapter at dna_decode/data/bvbrc_genome.py
#    feeds contig_count + N50 + MLST + assembly_accession into
#    candidates_from_bvbrc_ast. Both --assembly-metadata (legacy YAML) and
#    --assembly-metadata-csv are mutually exclusive.
uv run python -m scripts.pipeline ingest \
  --drugs ciprofloxacin,ceftriaxone,tetracycline \
  --ast-tsv path/to/BVBRC_genome_amr.csv \
  --assembly-metadata-csv path/to/BVBRC_genome.csv \
  --download-genomes

# 4. Populate the embedding cache (deferred — see ARCHITECTURE.md for the wiring;
#    embedding cache populate is invoked from a Phase 2 helper script that hasn't
#    shipped yet; Phase 1 callers populate the cache externally via cache.populate()).

# 5. Train per-drug classifier + run CV + emit clade-only baseline + validation gate
uv run python -m scripts.pipeline train \
  --drug ciprofloxacin --model evo --include-clade-baseline

# 6. Run ISM attribution + Tier 1-5 classification for one strain
uv run python -m scripts.pipeline attribute \
  --strain-id <bvbrc-strain-id> \
  --drug ciprofloxacin \
  --card-path path/to/card.json \
  --amrfinder-path path/to/amrfinder.tsv \
  --output data/processed/attribution_report.json

# 7. Build leaderboard across foundation models + classical baselines
uv run python scripts/leaderboard.py \
  --drugs ciprofloxacin,ceftriaxone,tetracycline \
  --models evo,dnabert2

# 8. (Optional, gated on CC ≥ 7.0 GPU) Validate that 4-bit Evo attribution matches full-precision
uv run python scripts/quantize_fidelity_check.py \
  --full-precision-attributions full_manifest.json \
  --quantized-attributions quantized_manifest.json \
  --drug ciprofloxacin
```

## Decoder v0 quickstart (Phase 2 in-flight)

The v0 AI DNA decoder operates on **cached strains** — a strain whose NT embeddings already live in the HDF5 cache (built by `pipeline ingest` + the Databricks N=147 cipro populate). UX + success criteria locked in `wiki/decoder_v0_ux_and_success_criterion.md`.

```bash
# Predict cipro R/S for a cached strain, with top-K attribution + audit-verdict propagation.
uv run python -m scripts.pipeline predict \
  --model-path data/processed/models/ciprofloxacin_nucleotide_transformer.pkl \
  --strain-id 562.12345 \
  --cache D:/dna_decode_cache/embeddings/nt_n147_cipro.h5 \
  --annotations D:/dna_decode_cache/refseq/GCF_xxx.x/annotations.gff3 \
  --audit-merge-json wiki/cipro_mechanism_phenotype_merge_2026-05-17.json \
  --output result.json
```

Writes `result.json` + `result.md` (markdown sidecar) per the v0 schema:

- `prediction` (R/S) + `calibrated_probability` + `confidence_tier` (HIGH/MEDIUM/LOW)
- `top_k_attribution` — gene-level ISM hits with resistance-catalog tier labels (Tier 1–5)
- `audit_verdict` — propagated from the merge gate; explicit `suspend_gate_fired` flag + verdict explanation when training cohort had `SUSPEND_CONDITION_4`
- `provenance` — model, training cohort, LOSO AUROC, trained-on date

**Not a clinical decision support tool.** Audit verdict + provenance must accompany any downstream interpretation. See `wiki/decoder_v0_ux_and_success_criterion.md` for full v0 schema + success criteria.

## Phase 1 success criteria

Phase 1 ships when:

- Smoke pipeline passes (`scripts/smoke_pipeline.py` returns exit 0)
- LOMO-clade-out CV AUROC ≥0.80 SLO / ≥0.85 target per drug
- Embedding model AUROC ≥0.10 above clade-only baseline on ≥75% of held-out clades
- Top-K=20 attribution-tier distribution: cipro ≥40% Tier 1-3 hits; ceftriaxone ≥25%; tet ≥30%; all ≤20% Fail
- Best foundation model beats best classical baseline by ≥3pp AUROC on ≥2 of 3 drugs
- Quantization-fidelity check returns GO (mean Spearman ≥0.7, intersection ≥0.6)

Phase 2 redesign trigger: classical baselines win on ≥2 drugs. The Step 18 classical-baselines control wires this empirically — see `plans/Ecoli_G2P_Platform_Technical_Plan.md` validation-gate section.

## Pilot gate alternate inputs

- `BVBRC_AST_TSV=path/to/ast.tsv` env var
- `bvbrc_ast.local_tsv_path: path/to/ast.tsv` in `config/datasources.yaml`

## Optional: route caches to a USB drive

Phase 1 runtime needs ~25GB (foundation models + strain genomes + embeddings). If your C: drive is tight, route caches to external storage:

```bash
# Replace E: with your drive letter (Windows) or /mnt/usb (Linux)
export HF_HOME=E:/hf_cache               # HuggingFace tokenizer + model cache
export DNA_DECODE_CACHE_ROOT=E:/dna_decode_cache
```

Then edit `config/datasources.yaml` to point `cache_dir` fields at the USB-backed path.

## Project workflow

Built using a personal Claude Code skill ladder for project planning:
- `/idea-anchor` → `/project-init` → `/brainstorm` ×3 → `/technical-plan` → `/probe` → `/execute-plan`
- Project ledger maintained via `/project-state` skill
- Execution state tracked in `.claude/execute-plan-state/Ecoli_G2P_Platform_Technical_Plan.json`
- All planning artifacts captured as audit trail

See `project_state/dna-decode-2026-05-11.md` for full decision history (17 hypotheses, 12+ decisions made, 54+ action-log entries as of 2026-05-17; "Phase 1 / 2 / 3" labels retrospective-only — new work tracked as Evidence Packets per the 2026-05-15 framing reset; Phase 1 evidence collection CLOSED 2026-05-17 with the cross-drug architectural finding synthesis).
