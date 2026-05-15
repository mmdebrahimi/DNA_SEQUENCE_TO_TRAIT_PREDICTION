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
  - `cache.py` — HDF5 embedding cache; opens once per `populate(model, strain_genomes, annotations)` batch; version-mismatch refusal. `verify_complete(expected_genes_by_strain) → CompletenessReport` is the consumer-side integrity gate (added 2026-05-15): every mean-pool consumer MUST call it and bail unless `report.all_complete` is True (half-flushed strain becomes valid-looking mean otherwise).
  - `classifiers.py` — XGBoost + sigmoid calibration on mean-pooled strain embeddings; NaN-aware aggregation
  - `classical_baselines.py` — AMRFinder + k-mer + Bakta-gene-presence baselines (the "would k-mer beat embeddings?" control)
- `interp/` — interpretability
  - `mutagenesis.py` — gene-level ISM + saturation mutagenesis; Tier 1-5 attribution-success framework; motif-recovery placeholder
- `eval/` — evaluation harness
  - `cv.py` — LOSO + LOMO + leave-one-Mash-clade-out CV
  - `metrics.py` — AUROC / AUPRC / Brier / ECE + attribution-precision + per-clade aggregation
  - `phylogeny.py` — Mash distance + ANI clustering. **Batched 2-call refactor 2026-05-15**: `compute_mash_distances` now issues 1 `mash sketch -o sketch.msh <all-fastas>` + 1 `mash dist sketch.msh sketch.msh` (vs prior N*(N-1)/2 calls — 10,731 at N=147 → 2). New `use_docker=True` kwarg routes through `tools/docker_runner.run` for Windows hosts without a native mash binary.
  - `clade_baseline.py` — clade-only baseline classifier + `validation_gate` (Phase 1 ship gate)
- `viz/` — visualization (`browser.py` matplotlib + TSV export; pygenometracks deferred to Phase 2)

Plus `tools/` (Stage 2 bioinformatics-tool runner via Docker Desktop — `docker_runner.py` exposes `run(image, args, mounts, env, capture_output, check, timeout)` with FileNotFoundError + TimeoutExpired wrapping as DockerRunnerError; routes Mash + AMRFinderPlus + Bakta on Windows hosts that lack native binaries) + `scripts/` (CLI entry points) + `tests/` (pytest) + `config/datasources.yaml` (declarative data-source registry) + `plans/` (technical plan + ship-path plan + EP design plans) + `project_state/` (Bellman-inspired project ledger).

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

# Phase 2 smoke gate (12-strain cipro mini-cohort; 3 variants: NT-XGBoost + k-mer + gene-presence)
# Uses data/processed/mini_cipro_nt_cache.h5 + gate_b_mini_cohort.parquet.
# AMRFinder deferred to Stage 1 prep (cohort's persisted resistance-gene fields are empty for the mini).
# IMPORTANT: NT path uses calibrate=False — CalibratedClassifierCV overcorrects at N≤20 and inverts predictions.
uv run python scripts/smoke_gate_12strain_cipro.py
# Writes wiki/smoke_gate_12strain_cipro_<date>.md result packet. Exit 0=PASS, 1=FAIL.

# Phase 2 Stage 1 engineering screen (N=40 cipro; 4 variants: NT-XGBoost gate + NT-logreg sanity + k-mer-XGB classical + NT+k-mer-fusion-logreg diagnostic)
# Uses D:/dna_decode_cache/embeddings/nt_n40_cipro.h5 (populate first) + gate_b_n40_cipro_cohort.parquet.
# Refactored 2026-05-14: thin orchestration over `dna_decode/eval/cv.py` (NT variants via `leave_one_strain_out_cv`)
# + `dna_decode/eval/loso_kmer.py` (k-mer + fusion runners shared with the smoke gate).
# Gate formula: max(NT-XGBoost, NT-logreg) AUROC − k-mer-XGB AUROC ≥ 3 pp under LOSO.
# Verdict bucket (frozen pure function of point gap): ≥5 pp CLEAN PASS / 3-5 pp NOISY PASS / <3 pp FAIL.
# Stage 2 action (decision-layer, separate field): BURST_STAGE_2 / HOLD_STAGE_2_CI_DEGENERATE /
# ALTERNATIVE_POOLING_RERUN / PIVOT_TO_BAKTA — computed from (verdict, ci_lo, fusion behavior) per the
# plan's Verdict-Time Pre-Commitments table.
# All variants run with calibrate=False (pinned by test at all 3 gate-bearing call sites);
# paired bootstrap CI (B=1000) on the gap surfaces clean-vs-noisy.
# `compute_gate_outcome` validates NT-vs-k-mer strain_ids alignment (raises on mismatch).
# Fusion variant is DIAGNOSTIC ONLY — does NOT count for gate.
HF_HOME=D:/hf_cache uv run python scripts/stage1_n40_cipro.py
# Writes wiki/stage1_n40_cipro_<date>.md result packet. Exit 0 ONLY when stage2_action == BURST_STAGE_2; exit 1 otherwise.

# Populate NT cache for the N=40 cohort (long-running on GTX 860M; ~5-7 hr wallclock)
# WARNING: external Seagate Portable D: drive USB hiccups have crashed populates mid-run
# with HDF5 EOA-truncation (errno=13, unrecoverable). For N=150 use Databricks burst.
HF_HOME=D:/hf_cache uv run python scripts/populate_cache.py \
  --cohort data/processed/gate_b_n40_cipro_cohort.parquet \
  --model nucleotide_transformer \
  --cache D:/dna_decode_cache/embeddings/nt_n40_cipro.h5 \
  --refseq-cache D:/dna_decode_cache/refseq \
  --device cuda

# Stage 2 cohort builder (N=150 target; effective ~147 after BV-BRC ceilings).
# Label-stratified MLST-balanced selection (per-class R/S separately); raises on
# missing MLST or imbalance outside slack. Hardcoded for cipro; pass --drug to change.
uv run python scripts/build_stage2_n150_cohort.py \
  --ast-tsv "C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv" \
  --assembly-metadata-csv "C:/Users/Farshad/Downloads/BVBRC_genome (1).csv" \
  --drug ciprofloxacin --target-total 150 --per-class 75 --balance-slack 15

# BV-BRC MLST-gap diagnostic. Layer-by-layer inspection of why cohort R-counts
# come up short. Surfaces whether the bottleneck is MLST, assembly_accession,
# AST coverage, or assembly-quality filters.
uv run python scripts/diagnose_bvbrc_mlst_gaps.py

# NT embedding-cache integrity probe (run BEFORE Stage 1 if the populate may have
# been interrupted — half-flushed strain at crash time becomes valid-looking mean
# embedding via Stage 1's gene-list-then-mean aggregation without this gate).
uv run python scripts/probe_nt_cache.py \
  --cohort data/processed/gate_b_n40_cipro_cohort.parquet \
  --cache D:/dna_decode_cache/embeddings/nt_n40_cipro.h5 \
  --refseq-cache D:/dna_decode_cache/refseq
# Exit: 0 = ALL_COMPLETE, 1 = INCOMPLETE (partial/absent/corrupt strains), 2 = cache fault.

# Stage 2 bioinformatics toolchain (Mash + AMRFinderPlus + Bakta) via Docker.
# DBs installed at C:/Users/Farshad/dna_decode_stage2/{amrfinder,bakta}_db/.
# Install + smoke artifact: wiki/stage2_install_artifact_2026-05-15.md.
# Pinned images: quay.io/biocontainers/mash:2.3--hb105d93_10 / ncbi/amr:4.2.7-2026-03-24.1 / oschwengers/bakta:v1.11.4.
# Direct docker invocation from Git Bash REQUIRES MSYS_NO_PATHCONV=1 (see Gotchas).
# Bakta bakta_db REQUIRES bash-c entrypoint wrapper for conda-init (see Gotchas).
MSYS_NO_PATHCONV=1 docker run --rm \
  -v C:/Users/Farshad/AppData/Local/Temp/mash_smoke:/data \
  quay.io/biocontainers/mash:2.3--hb105d93_10 \
  mash sketch -o /data/sketch /data/*.fna
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
- **Mash CLI is an external binary dependency** (Step 10 phylogeny clustering). **As of 2026-05-15: Docker route is the canonical install path** via `tools/docker_runner.py` + pinned image `quay.io/biocontainers/mash:2.3--hb105d93_10`. `compute_mash_distances(..., use_docker=True)` routes through it. Native binary on PATH still works (Linux/WSL2: `apt install mash`; Windows: GitHub releases binary); the Docker path eliminates host-environment variance. `pyani` fallback in spec remains unwired.
- **4-bit Evo via bitsandbytes is unavailable on this project's hardware.** Actual GPU is GTX 860M (Maxwell, CC=5.0, 4 GiB VRAM); bitsandbytes requires CC ≥ 7.0. Linux/CUDA caveat is moot here. The `[platform_system != 'Windows']` marker in `pyproject.toml` keeps the dep optional. Real-data runs use NT v2 100M only; DNABERT-2 + Evo are not callable on this machine.
- **`pilot.fetch_bvbrc_drug_counts` raises NotImplementedError without a local TSV / env var / config entry.** Live API integration is deferred until first real-data run. Workaround: download BV-BRC AST TSV from `ftp.bvbrc.org` + pass via `--ast-tsv`.
- **`pilot.fetch_ncbi_assembly_quality` stays scaffolded INTENTIONALLY.** Phase 2 introduced a separate CSV adapter (`dna_decode/data/bvbrc_genome.py`) that bypasses it via `pipeline ingest --assembly-metadata-csv path/to/BVBRC_genome.csv`. The scaffold's 2-key return (contig_count + n50) is wrong-shaped for the real CSV's 7+ richer fields; do not "implement the NotImplementedError." Live NCBI Datasets REST integration is Phase 3 work and will replace `fetch_ncbi_assembly_quality` then.
- **`motif_recovery` is a placeholder** — returns the same high-impact position list for every motif name. Phase 2 work; not a Phase 1 ship gate.
- **GFF3 annotation source variance — `gene_id` vs `gene_symbol`.** As of 2026-05-14, `parse_gff3` populates THREE separate columns: `gene_id` (from GFF3 `ID=` / `Name=` — strain-unique by construction, e.g., `gene-b0001`; this is the embedding-cache key in `extract_cds_sequences`), `gene_symbol` (from GFF3 `gene=` — cross-strain, e.g., `gyrA`; populated for only ~11% of CDSs in typical RefSeq GFF3), and `locus_tag` (also strain-unique). For gene-family / presence-absence features that need to generalize across strains, use `gene_symbol` first. Do NOT use `gene_id` as a cross-strain identifier — it WILL cause 0% LOSO vocab overlap and AUROC=0.000 (mechanism: held-out rows all-zero → XGBoost predicts training class prior → prior inverts against held-out label on balanced LOSO). See `plans/Gene_Presence_AUROC_Bug_Fix_Plan.md`. `parse_genbank` populates `gene_id` + `gene_symbol` to the same value (the qualifier `gene`) — known asymmetry; downstream consumers should still prefer `gene_symbol`. The smoke runner gates on median per-fold vocab overlap (`min_median_vocab_overlap=0.20`) and reports `INDETERMINATE_IDENTIFIER_OOV` rather than a misleading AUROC when degenerate.
- **`parse_gff3` has TWO-PASS parent→CDS gene_symbol propagation** (added 2026-05-14 PM for Bakta compatibility). Bakta-style GFF3 puts `gene=gyrA` on the parent `gene` row; CDS rows link via `Parent=` attribute. Pass 1 builds `parent_id → gene_symbol` map from gene-type rows; Pass 2 inherits onto CDS rows lacking their own `gene=`. Same-row `gene=` always wins (RefSeq-style unaffected). If you find yourself adding a Roary / Bakta / orthogroup pipeline, this is the parser layer that already knows how to handle parent-linked annotations.
- **`leave_one_*_out_cv` raises on unassigned strain_ids by default** (`dna_decode/eval/cv.py`, 2026-05-14 PM defensive fix). The old silent `__unassigned__` bucketing was a Stage 2 risk — silent bucketing warps fold structure (one giant fold for all unassigned). Pass `allow_unassigned=True` for the legacy silent-bucket behavior; only do this if you've explicitly audited the unassigned set.
- **BV-BRC's real cohort bottleneck is `assembly_accession`, NOT MLST.** Raw BVBRC_genome.csv has 96% MLST coverage. The parser drops 35,790 of 85,114 rows that lack `assembly_accession` (NCBI Datasets API can't fetch them). For cipro: 568 R total in AST, but only 77 R have downloadable accessions, and 72 of those pass quality+MLST filters. If a cohort comes up short on R strains, the lever is NOT to chase MLST — it's to find R strains with NCBI-downloadable assemblies (or change source databases). See `scripts/diagnose_bvbrc_mlst_gaps.py` for the layer-by-layer diagnostic.
- **`build_cohort`'s default MLST-balanced selection does NOT stratify by label.** Stage 2 cohort builder (`scripts/build_stage2_n150_cohort.py`) uses per-class `_mlst_balanced_selection` for R + S separately. The naive single-pool call left ~28 available R strains on the table at small cohorts (49R/101S vs achievable 72R/75S). Always per-class for binary-classification cohorts.
- **`.gitignore`'s `/data/` is anchored** (leading slash = repo root only). Bare `data/` would also match `dna_decode/data/` subpackage source.
- **`EmbeddingCache.populate(skip_existing=True)` skips at GENE-DATASET path level, NOT strain level** (`cache.py:296`). Stage 1's loader admits a strain on ≥1 cached gene + mean-pools whatever's there. A mid-flush crash leaves a partial gene set that LOOKS like a complete strain to Stage 1. f9ed79f flush patch bounds loss to ~1 strain on crash but doesn't eliminate. **Consumer-side defense:** every mean-pool consumer MUST call `cache.verify_complete(expected_genes_by_strain) → CompletenessReport` and bail unless `report.all_complete` is True. Standalone audit CLI at `scripts/probe_nt_cache.py`. 8 unit tests at `tests/test_models_cache.py` pin the 4 status buckets + the `all_complete=False` rule on partial.
- **Git Bash silently breaks Docker `-v <host>:/<container>` volume mounts via MSYS path conversion** (`/db` → `C:/Program Files/Git/db`, silent unbound ephemeral mount). For direct `docker run` from Git Bash: prefix every invocation with `MSYS_NO_PATHCONV=1`. Symptom on failure: command exits 0, downloads succeed inside the container, host directory is 0 bytes after `--rm`. Python `subprocess.run` (e.g., `tools/docker_runner.run`) is unaffected — no shell to munge the path.
- **Bakta `bakta_db download` AMRFinder dep-check fails when `--entrypoint /opt/conda/bin/bakta_db` skips the image's bash shell init**. Conda env stays un-activated; `bakta_db`'s deeper-than-`which-amrfinder` dep check then fails despite the binary being on PATH (v1.11.4 + v1.12.0 both affected). **Workaround:** invoke via `--entrypoint /bin/bash -c "bakta_db download --output /db --type light"`. DB v6.0 light at `C:/Users/Farshad/dna_decode_stage2/bakta_db/db-light/` (4.0 GB; Zenodo DOI 10.5281/zenodo.14916843). See `wiki/stage2_install_artifact_2026-05-15.md` §B1.
- **Project framing: Phase 1 / 2 / 3 labels are now retrospective-only** (2026-05-15 framing reset per `project_state/dna-decode-2026-05-11.md` Pending Decisions row 8). New work is tracked as Evidence Packets (EP1 cipro Stage 2 / EP2 cef + tet smoke / EP3 attribution audit / EP4 clade-shift / EP5 cef + tet Stage 2) in the Mid-term table. Parked horizon tracks: 2nd-organism portability, MIC continuous head, pan-genome graph, multimodal DNA+image. The "Phase 1 shipped" label in historical artifacts means infrastructure-complete only — scientific ship criteria gate at Stage 2 N=150 evaluation per `wiki/phase1_ship_report.md:77-80`.

## Project workflow

Built using a personal Claude Code skill ladder:
- `/idea-anchor` → `/project-init` → `/brainstorm` ×N → `/technical-plan` → `/probe` → `/save-plan` → `/execute-plan`
- Project ledger maintained via `/project-state` skill at `project_state/dna-decode-2026-05-11.md`
- Execution state tracked at `.claude/execute-plan-state/Ecoli_G2P_Platform_Technical_Plan.json`
- Plan index at `wiki/plans-index.md`

Run `/brainstorm` between every wave; don't skip "to save time." Three brainstorm rounds caught 3 grounded contract gaps each — that's the pattern.
