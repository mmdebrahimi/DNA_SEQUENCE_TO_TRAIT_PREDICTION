# CLAUDE.md

Project-specific guidance for Claude Code (claude.ai/code) when working in this repo.

## What this is

Personal solo Phase 1 of an E. coli genotype-to-phenotype prediction platform. Predicts antibiotic resistance (ciprofloxacin / ceftriaxone / tetracycline) from genomic DNA + identifies which genomic regions drive predictions. Biologically interpretable; **not** causal-claim-making.

Long-term vision: multimodal genotype-phenotype platform expanding to eukaryotes + image-paired phenotype data. Phase 1 is foundation infrastructure, NOT a stepping stone to "DNA → animal image" prediction.

## Architecture (one level of depth that matters)

Layered Python package under `dna_decode/`:

- `data/` — ingestion + preprocessing layer
  - `pilot.py` — Step 0.5 HARD-gate (counts per-drug labeled strains before downstream fires)
  - `refseq.py` — NCBI Datasets ZIP downloader + unpacker (`genome.fna` / `annotations.gff3` / `annotations.gbk`)
  - `annotations.py` — GFF3 + GenBank parsers; CDS + intergenic extractors
  - `resistance_db.py` — CARD + AMRFinder loaders → unified `ResistanceCatalog`
  - `ast_data.py` — BV-BRC AST TSV loader with broth-microdilution + organism filters
  - `cohort.py` — drug-first cohort builder + MLST round-robin balancing + assembly-quality threshold + save/load parquet + `candidates_from_bvbrc_ast` adapter
- `models/` — foundation models + classifiers
  - `foundation.py` — Evo + DNABERT-2 + Nucleotide Transformer + GENA-LM wrappers (lazy-load); `MockFoundationModel` for tests
  - `cache.py` — HDF5 embedding cache; opens once per `populate(model, strain_genomes, annotations)` batch; version-mismatch refusal
  - `classifiers.py` — XGBoost + sigmoid calibration on mean-pooled strain embeddings; NaN-aware aggregation
  - `classical_baselines.py` — AMRFinder + k-mer + Bakta-gene-presence baselines (the "would k-mer beat embeddings?" control)
- `interp/` — interpretability
  - `mutagenesis.py` — gene-level ISM + saturation mutagenesis; Tier 1-5 attribution-success framework; motif-recovery placeholder
- `eval/` — evaluation harness
  - `cv.py` — LOSO + LOMO + leave-one-Mash-clade-out CV
  - `metrics.py` — AUROC / AUPRC / Brier / ECE + attribution-precision + per-clade aggregation
  - `phylogeny.py` — Mash distance + ANI clustering (subprocess wrapper)
  - `clade_baseline.py` — clade-only baseline classifier + `validation_gate` (Phase 1 ship gate)
- `viz/` — visualization (`browser.py` matplotlib + TSV export; pygenometracks deferred to Phase 2)

Plus `scripts/` (CLI entry points) + `tests/` (pytest) + `config/datasources.yaml` (declarative data-source registry) + `plans/` (technical plan + ship-path plan) + `project_state/` (Bellman-inspired project ledger).

## Common commands

```bash
# Install + sync (pytest is in default deps; ruff/pytest-cov in dev extras)
uv sync
uv sync --extra dev          # optional: ruff + pytest-cov
uv sync --extra quantize     # optional: bitsandbytes (Linux/WSL only)

# Test suite
uv run pytest tests/ -v

# Pilot gate (Step 0.5 HARD gate). Download BV-BRC AST TSV first from ftp.bvbrc.org.
uv run python -m scripts.pilot_gate \
  --drugs ciprofloxacin,ceftriaxone,tetracycline \
  --target-per-drug 150 \
  --ast-tsv path/to/bvbrc_ast.tsv
# Exit: 0=GO, 1=NO-GO, 2=PilotGateError, 3=no AST source provided

# Phase 1 full pipeline
uv run python -m scripts.pipeline ingest --drugs ciprofloxacin,ceftriaxone,tetracycline
uv run python -m scripts.pipeline train --drug ciprofloxacin --model evo --include-clade-baseline
uv run python -m scripts.pipeline attribute --strain-id <accession> --drug ciprofloxacin
```

## Phase 1 success criteria (per technical plan + Tier 1-5 rubric)

Phase 1 ships when all of these pass:

- LOMO-clade-out CV AUROC ≥0.80 SLO / ≥0.85 target per drug
- Embedding model AUROC ≥0.10 above clade-only baseline on ≥75% of held-out clades (mechanistic-signal validation gate)
- Top-K=20 attribution tier distribution: cipro ≥40% Tier 1-3 hits; ceftriaxone ≥25%; tet ≥30%; all drugs ≤20% Fail-tier
- Best foundation model beats best classical baseline by ≥3pp AUROC on ≥2 of 3 drugs

Phase 2 redesign trigger: classical baselines win on ≥2 drugs. The plan's Step 18 wires this control explicitly; **don't soften it**.

## Gotchas

- **BV-BRC strain_id ≠ NCBI assembly_accession.** `CandidateStrain` carries both. `download_cohort_genomes` resolves via `assembly_accession` (NCBI accession like `GCF_000005845.2`); falling back to `strain_id` would 404 NCBI Datasets.
- **HDF5 cache opens ONCE per `populate()` batch.** Don't refactor `cache.populate` to call `cache.put()` in a loop — that's the 2.25M-file-open footgun that Wave 2.5 fixed.
- **`build_attribution_report` tries BOTH `gene_id` and `locus_tag`** against the resistance catalog via `_best_tier_across_candidates` helper. Locus tags (`b0001`, `TAG_001`) never match gene symbols (`gyrA`); `gene_id` may or may not be a gene symbol depending on annotation source.
- **Calibration CV folds use minority-class count, NOT majority.** Rare-resistance drugs (e.g., low ceftriaxone resistance prevalence) break the original `max(positives, negatives)` formulation. Fixed in Wave 3.5 hardening.
- **Mash CLI is an external binary dependency** (Step 10 phylogeny clustering). Linux/WSL2: `apt install mash`. Windows: download binary from Mash GitHub releases. `pyani` fallback exists in spec but is unwired — confirm before first real-data run.
- **4-bit Evo via bitsandbytes is Linux/CUDA only.** Windows users need WSL2 for the quantization extras. `[platform_system != 'Windows']` marker in `pyproject.toml` keeps the dep optional.
- **`pilot.fetch_bvbrc_drug_counts` raises NotImplementedError without a local TSV / env var / config entry.** Live API integration is deferred until first real-data run. Workaround: download BV-BRC AST TSV from `ftp.bvbrc.org` + pass via `--ast-tsv`.
- **`pilot.fetch_ncbi_assembly_quality` stays scaffolded INTENTIONALLY.** Phase 2 introduced a separate CSV adapter (`dna_decode/data/bvbrc_genome.py`) that bypasses it via `pipeline ingest --assembly-metadata-csv path/to/BVBRC_genome.csv`. The scaffold's 2-key return (contig_count + n50) is wrong-shaped for the real CSV's 7+ richer fields; do not "implement the NotImplementedError." Live NCBI Datasets REST integration is Phase 3 work and will replace `fetch_ncbi_assembly_quality` then.
- **`motif_recovery` is a placeholder** — returns the same high-impact position list for every motif name. Phase 2 work; not a Phase 1 ship gate.
- **GFF3 annotation source variance.** `parse_gff3` collapses `ID=` / `Name=` / `gene=` attributes into one `gene_id` column. Bakta emits `ID=g1` + `locus_tag=TAG_001`; RefSeq emits `gene=gyrA`. Downstream code matching against a gene-symbol catalog needs to try both `gene_id` + `locus_tag`.
- **`.gitignore`'s `/data/` is anchored** (leading slash = repo root only). Bare `data/` would also match `dna_decode/data/` subpackage source.

## Project workflow

Built using a personal Claude Code skill ladder:
- `/idea-anchor` → `/project-init` → `/brainstorm` ×N → `/technical-plan` → `/probe` → `/save-plan` → `/execute-plan`
- Project ledger maintained via `/project-state` skill at `project_state/dna-decode-2026-05-11.md`
- Execution state tracked at `.claude/execute-plan-state/Ecoli_G2P_Platform_Technical_Plan.json`
- Plan index at `wiki/plans-index.md`

Run `/brainstorm` between every wave; don't skip "to save time." Three brainstorm rounds caught 3 grounded contract gaps each — that's the pattern.
