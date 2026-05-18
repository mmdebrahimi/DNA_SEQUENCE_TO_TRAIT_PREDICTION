# DNA_SEQUENCE_TO_TRAIT_PREDICTION

Genotype-to-phenotype (G2P) inference platform — predicts phenotypic traits from genomic DNA sequences AND identifies which genomic regions are most strongly associated with those predictions. Biologically interpretable, not causal-claim-making.

## Status: Phase 1 — CLOSED 2026-05-17 (infrastructure + cross-drug architectural finding)

Phase 1 evidence collection closed 2026-05-17. Cross-drug architectural finding synthesis at `wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md`:

> At 12-strain smoke fidelity, frozen-NT-whole-genome-pooling PASSES on concentrated-signal AMR mechanisms (cipro QRDR point mutations: AUROC 0.750; cef plasmid acquired-gene β-lactamases: AUROC 0.833) AND FAILS on distributed mobile-element mechanisms (tet tet-family efflux + ribosomal protection: AUROC 0.400, anti-predictive). The architecture's failure mode appears mechanism-class-bounded, largely independent of drug identity at smoke fidelity.

EP1 cipro closed internally (`wiki/cipro_ep1_closeout_2026-05-17.md`) with a 4-tier adversarial audit infrastructure (mechanism × MIC × opacity merge with structurally-enforced SUSPEND gate). EP2 cef + tet smoke fired (cef PASS, tet FAIL, H17 falsified). No Databricks burst spent. External publication deferred per PC1=`internal_closeout`.

Phase 1 code: all 18 implementation steps shipped Wave 0-7 (2026-05-11 → 2026-05-12) + 3 hardening waves; cross-drug Evidence Packet evidence collection completed 2026-05-17 per the Evidence Packets framing reset 2026-05-15. **Phase 2 entry:** BV-BRC strict-MIC 3-drug feasibility census (per `project_state/dna-decode-2026-05-11.md` Candidate-next-actions row 1) — deferred to a fresh session; open via `/idea-anchor` + `/project-init` first per the synthesis's narrow reopen rule.

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
