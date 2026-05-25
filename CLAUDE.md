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

# EP1 cipro audit infrastructure (shipped 2026-05-17 as part of Phase 1 closeout).
# All 4 scripts emit both .md (narrative) + .json (machine-readable) sidecars.
# Cohort-provenance caveat: cohort uses the N=38 cipro cohort; for cef/tet adapt
# the cohort path (mini cohorts at data/processed/gate_b_mini_{cef,tet}_cohort.parquet).

# Cipro raw BV-BRC MIC rejoin audit (label-noise diagnostic; ~5 min).
# Tiers each strain under CLSI + EUCAST breakpoints; classifies HIGH_R / HIGH_S /
# DECISIVE / BORDERLINE / AMBIGUOUS / CONFLICT / NO_MIC.
uv run python scripts/cipro_mic_audit.py \
  --cohort data/processed/gate_b_n40_cipro_cohort.parquet \
  --ast-csv "C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv"

# Cipro AMRFinderPlus mechanism audit (per-strain QRDR/plasmid/efflux/regulatory hit
# detection across all R + S strains). ~95s/strain × 38 ≈ ~1 hour first run; cached.
# Runs AMRFinder via tools.docker_runner.run; AMRFINDER_DB at dna_decode_stage2.
# Filter: CIPRO_RELEVANT_AMR_CLASSES = {QUINOLONE, FLUOROQUINOLONE, MULTIDRUG}.
# Synonymous-SNP filter + main.tsv/mutations.tsv dedup pinned.
uv run python scripts/cipro_mechanism_audit.py
# (Detached batch wrapper at run_mechanism_audit_detached.bat for ~1 hr CPU.)

# Cipro mechanism × MIC merge (joins the audit outputs above into a single
# noise_class table per strain). Emits structurally-enforced SUSPEND_CONDITION_4
# gate verdict + recommended_next_step.
uv run python scripts/cipro_mechanism_phenotype_merge.py
# Verdict thresholds: SIGNAL_DOMINATES (>=0.70 clean) / MIXED (0.40-0.70) /
# NOISE_DOMINATES (<0.40 clean). Today's verdict 2026-05-17: NOISE_DOMINATES
# (signal quality 0.17, opacity_count=0 — AMRFinder is NOT the bottleneck).

# Cipro curated AMR baseline (2-layer verdict; load-bearing for PIVOT TRIGGER
# condition 4). LR + XGB over (AMRFinder POINT + acquired + k-mer + MLST) LOSO.
# Emits original_condition_4 (frozen) + amended_condition_4 (no_POINT >= 0.773
# OR mechanism_only >= 0.80) + given_suspended_gate=INFORMATIONAL_ONLY field.
# Refuses to fire when the merge gate's clean_count < 10 (SUSPEND_CONDITION_4).
uv run python scripts/cipro_curated_baseline.py
# Today (2026-05-17): gate fires SUSPEND_CONDITION_4; script honored gate + did
# not produce a misleading AUROC verdict on uninterpretable labels.

# BV-BRC strict-MIC 4-drug feasibility census (shipped 2026-05-18; Phase 2 entry).
# Counts per-drug feasibility of building N=150 cohorts at TWO label-quality bars:
#   strict-MIC = HIGH_R + HIGH_S only (4x safety margin)
#   relaxed-MIC = strict + DECISIVE_R + DECISIVE_S (audit framework as downstream gate)
# 6-stage pipeline: AST rows -> distinct genomes -> with MIC -> classification pass
#   -> with assembly_accession -> passing assembly QC.
# Headline 2026-05-18: NO drug clears N=150 per-class at either bar; assembly_accession
# is the structural bottleneck (96% drop on cipro strict-MIC).
uv run python -m scripts.bvbrc_strict_mic_4drug_census
# Writes wiki/bvbrc_strict_mic_4drug_census_<date>.{md,json}.

# v0 decoder predict (shipped 2026-05-18; UX + criteria LOCKED at wiki/decoder_v0_ux_and_success_criterion.md).
# Emits v0 JSON + markdown sidecar with prediction + calibrated_probability +
# confidence_tier + top_k_attribution (gene-level ISM + Tier 1-5 catalog tier)
# + audit_verdict (SUSPEND propagation) + provenance.
uv run python -m scripts.pipeline predict \
  --model-path data/processed/models/ciprofloxacin_nucleotide_transformer.pkl \
  --strain-id <BV-BRC-strain-id> \
  --cache D:/dna_decode_cache/embeddings/nt_n147_cipro.h5 \
  --annotations D:/dna_decode_cache/refseq/GCF_xxx.x/annotations.gff3 \
  --audit-merge-json wiki/cipro_mechanism_phenotype_merge_2026-05-17.json \
  --top-k 10 --output result.json
# Writes result.json + result.md sidecar.
```

```bash
# EP2 cef + tet smoke gate (fired 2026-05-17). Uses cipro mini-cohort smoke
# infrastructure with --drug arg + new mini cohorts. Today's verdicts:
#   cef: NT-XGBoost 0.833 = k-mer 0.833 (PASS); H17 partially preserved on cef.
#   tet: NT-XGBoost 0.400 anti-predictive vs k-mer 0.722 (FAIL); H17 falsified.
# Re-fire the smoke (e.g., after smoke-runner bug fix) with:
HF_HOME=D:/hf_cache uv run python scripts/smoke_gate_12strain_cipro.py \
  --cohort data/processed/gate_b_mini_cef_cohort.parquet \
  --nt-cache D:/dna_decode_cache/embeddings/nt_n40_cipro.h5 \
  --refseq-cache D:/dna_decode_cache/refseq \
  --drug ceftriaxone
# Replace --drug ceftriaxone with --drug tetracycline + mini_tet_cohort.parquet
# for tet. Output strings now templated on --drug (output filename auto-named).
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
- **Project framing: Phase 1 / 2 / 3 labels are now retrospective-only** (2026-05-15 framing reset per `project_state/dna-decode-2026-05-11.md` Pending Decisions row 8). New work is tracked as Evidence Packets (EP1 cipro / EP2 cef + tet smoke / EP3 attribution audit / EP4 clade-shift / EP5 cef + tet Stage 2) in the Mid-term table. Parked horizon tracks: 2nd-organism portability, MIC continuous head, pan-genome graph, multimodal DNA+image. **Phase 1 evidence collection CLOSED 2026-05-17** with cross-drug architectural finding synthesis at `wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md`: NT-frozen-whole-genome-pooling PASSES on concentrated-signal mechanisms (cipro QRDR 0.750 + cef plasmid β-lactamases 0.833) and FAILS on distributed mobile-element mechanisms (tet 0.400 anti-predictive). EP1 cipro closed internally (`wiki/cipro_ep1_closeout_2026-05-17.md`; no Databricks burst). EP2 H17 falsified. External publication deferred per PC1=`internal_closeout`. Phase 2 entry: BV-BRC strict-MIC 3-drug feasibility census (`project_state/dna-decode-2026-05-11.md` Candidate-next-actions row 1) — deferred to a fresh session per the synthesis's narrow reopen rule (reopen ONLY for internal contradiction or factual mismatch).
- **Anti-predictive AUROC at N=12 with `calibrate=False` is data-shape divergence, NOT a plumbing bug** (extends 2026-05-14 LESSON on calibration overcorrection). The 2026-05-14 LESSON applies specifically to `CalibratedClassifierCV` with isotonic regression at N=11 training (symmetric two-value output around 0.5). NT-XGBoost on the tet 12-strain smoke 2026-05-17 returned 0.400 vs k-mer 0.722 with `calibrate=False` already set — that's genuine architectural mismatch on distributed mobile-element resistance, not a calibration bug. The two failure-mode patterns are distinct: (a) calibrate=True at N≤20 → AUROC ≈ 0 with symmetric two-value scores (calibration bug); (b) calibrate=False + distributed-mechanism resistance → AUROC 0.3-0.45 with non-degenerate score distribution (data-shape divergence). Diagnose by checking the calibrate flag + score-distribution shape BEFORE attributing to either.
- **Cef + tet mini cohorts filtered from the cipro N=38 cohort, not built BV-BRC-wide.** `data/processed/gate_b_mini_cef_cohort.parquet` + `gate_b_mini_tet_cohort.parquet` were built 2026-05-17 via inline pandas filter on the existing N=38 cipro cohort (6R/6S each, 12 unique MLSTs, full assembly availability). This reuses the populated NT cache at `D:/dna_decode_cache/embeddings/nt_n40_cipro.h5` without re-populate. **Provenance caveat:** mechanism distribution may reflect cipro-cohort selection artifacts, not cef/tet biology — any conclusions are scoped to "within this reused mini-cohort" unless a BV-BRC-wide cef/tet cohort is built (deferred). For Stage 2 + publication-facing claims, build the cohort from raw BV-BRC AST directly.
- **Smoke runner had a silent-variant-drop bug** (fixed in `plans/Cef_Mechanism_Audit_Plan.md` Step 1; reduced plan post-/brainstorm round 3). `scripts/smoke_gate_12strain_cipro.py:374-382` had `try/except Exception` that swallowed any variant runner's exception silently. Today's cef + tet 2026-05-17 smoke reports both show only 2 variants (NT-XGBoost + k-mer) — gene-presence raised `FileNotFoundError` (missing GFF3 on one strain) and got silently dropped. Fix is queued in the cef plan: re-raise on unexpected errors; render `INDETERMINATE_<reason>` rows when known exception types fire. Same script may have hidden other variant failures in earlier runs — review existing smoke artifacts if you find a 2-variant table where 3+ were expected.
- **AMRFinder cipro-relevant Class filter:** keep MULTIDRUG, not just QUINOLONE / FLUOROQUINOLONE. `scripts/cipro_mechanism_audit.py` filters mutations.tsv to `CIPRO_RELEVANT_AMR_CLASSES = QUINOLONE_CLASSES | {"MULTIDRUG"}` — regulatory mutations like `marR_V84WfsTer` + `acrR_S30HfsTer` come through with AMRFinder Class=MULTIDRUG (not QUINOLONE). Filtering on QUINOLONE-only would drop real cipro-affecting regulatory frameshifts. The same discipline applies to any future cef / tet mechanism audit (cef-relevant: BETA-LACTAM + CARBAPENEM + CEPHALOSPORIN + MULTIDRUG; tet-relevant: TETRACYCLINE + MULTIDRUG). **Per-drug AMRFinder Class filter is now centralized in `dna_decode/data/mic_tiers.py::DRUG_AMRFINDER_CLASSES`** (added 2026-05-18) — future drug audits should `amrfinder_classes_for(drug)` rather than hardcoding.
- **Shared MIC-tier classifier + per-drug catalogs at `dna_decode/data/mic_tiers.py`** (added 2026-05-18). Single source of truth for: `breakpoints_for(drug)` (CLSI 2024 + EUCAST 14.0 E. coli — cipro/cef/tet/gent), `classify_tier(mics, distinct_calls, breakpoints)` (HIGH_R/HIGH_S/DECISIVE_R/DECISIVE_S/BORDERLINE/AMBIGUOUS/CONFLICT/NO_MIC), `amrfinder_classes_for(drug)`, `loci_by_mechanism_for(drug)` (mechanism → loci catalog), `primary_mechanisms_for(drug)`, `classify_gene_symbol(drug, symbol)` (tolerant prefix-match like `qnrB19` → `qnrB`). Constants: `DRUG_BREAKPOINTS`, `DRUG_LOCI_BY_MECHANISM`, `DRUG_PRIMARY_MECHANISMS`, `DRUG_AMRFINDER_CLASSES`, `CO_RESISTANCE_MECHANISMS` (efflux/regulatory/porin_loss — shared cross-drug modifiers), `STRICT_MIC_TIERS`, `RELAXED_MIC_TIERS`. 76 unit tests at `tests/test_mic_tiers.py`. `scripts/bvbrc_strict_mic_4drug_census.py` uses it; `scripts/cipro_mic_audit.py` + `scripts/cipro_mechanism_audit.py` + `scripts/cipro_mechanism_phenotype_merge.py` + `scripts/cipro_curated_baseline.py` left as cipro-specific by design (have 80+ dedicated tests; migration is drift-prevention only).
- **BV-BRC strict-MIC 4-drug feasibility verdict 2026-05-18 — NO drug clears N=150 per-class at either strict-MIC OR relaxed-MIC.** Cipro 17R/4S, cef 66R/2S, tet 1R/0S, gent 2R/132S. Structural bottleneck is `assembly_accession` (96% drop on cipro strict-MIC); relaxed-MIC barely helps because the few DECISIVE strains also lack downloadable accessions. **The categorical N=147 cipro cohort is THE viable v0 substrate** (already built at `data/processed/stage2_n150_cipro_cohort.parquet`; Databricks burst populating NT cache 2026-05-18). Cef may be v0.1 substrate if cef-S labels are backfilled from PATRIC/NARMS/EuSCAPE (deferred). Tet + gent dropped from Phase 2 candidate list. 4th-mechanism-class falsifier (gent) substrate also infeasible — needs a new candidate. See `wiki/bvbrc_strict_mic_4drug_census_2026-05-18.{md,json}` + `scripts/bvbrc_strict_mic_4drug_census.py`.
- **Audit framework relaxes the cohort-construction bar** (architectural insight emerged from the BV-BRC infeasibility finding). Pre-audit-framework, cohort construction needed strict-MIC labels because nothing downstream caught dirty ones. Post-audit-framework (mechanism × MIC × opacity merge with SUSPEND gate, shipped Phase 1 closeout), permissive cohort + downstream audit gate is the right shape. Strict-MIC at 4× safety margin is paper-grade label-quality; pre-filtering on it duplicates the audit framework's job AND discards usable training data. v0 decoder UX (`wiki/decoder_v0_ux_and_success_criterion.md`) reflects this: train on categorical labels, propagate audit verdict (SUSPEND framing) into every prediction's JSON output.
- **AMD Polaris (RX 4xx/5xx) is NOT viable for modern PyTorch** (verified empirically 2026-05-18). ROCm dropped Polaris (gfx803) support after v4.5 (2021); modern PyTorch (2.x) requires ROCm 5.x+. User had a 12-GPU RX 570/580 mining rig that couldn't be used — visual confirmation via the AORUS RX 580 8GB box label gated the decision in ~5 minutes vs hours of failed ROCm install. **Verify GPU architecture tier BEFORE committing to setup work**: AMD Polaris/Vega/RDNA1 = avoid; RDNA2 = unofficial-but-workable via `HSA_OVERRIDE_GFX_VERSION=10.3.0`; RDNA3 = officially supported. NVIDIA: Pascal (GTX 10xx, CC 6.1) works but no Tensor Cores / no bitsandbytes; Turing+ (CC≥7.5, RTX 20xx+) = full support. Cross-project memory at `~/.claude/projects/.../memory/feedback_gpu_arch_verify_before_setup.md`.
- **Databricks DBFS FUSE mount `mkdir` fails with `OSError [Errno 5] Input/output error`** on modern runtimes (especially Unity-Catalog-enabled). Error trace lands at `pathlib.py:1313 os.mkdir(self, mode)`. The `/dbfs/` FUSE mount is read-mostly for native Python `mkdir`; reads from `/dbfs/...` paths still work. **NEVER use `Path("/dbfs/...").mkdir(...)` or `os.makedirs("/dbfs/...")` from notebook code.** Workaround: use `/local_disk0/` for working storage during the run + push final outputs to `dbfs:/...` via `dbutils.fs.cp` at end (REST API path, bypasses FUSE). Set cluster auto-terminate ≥240 min for ML inference workloads. Cross-project memory at `~/.claude/projects/.../memory/feedback_databricks_dbfs_fuse_mkdir_bug.md`. Working notebook template: `C:/Users/Farshad/Downloads/Stress_Load_decoder_fixed_2026-5-18.py`.
- **Project north star clarified 2026-05-18**: "AI DNA decoder tool, not papers. Failure-tolerant iteration." Supersedes the Phase 2 framing brainstorm's 3-candidate research/audit-tier slate. Frame architectural findings (mechanism-class-bounded NT-pooling, SUSPEND gate) as engineering inputs to the decoder, not standalone research outputs. Don't propose paper/blog/arXiv deliverables as Phase 2+ anchors. v0 decoder spec at `wiki/decoder_v0_ux_and_success_criterion.md` (LOCKED). See project memory at `~/.claude/projects/.../memory/project_dna_decode_north_star.md`.
- **`build_cohort` asserts assembly_accession uniqueness by default** (added 2026-05-22). Raises `CohortConstructionError` if any two candidate strains share an `assembly_accession` (excluding empty strings) — prevents same-genome LOSO leakage by construction. Override via `build_cohort(..., allow_duplicate_accessions=True)` for the rare case of intentional same-genome registration. Triggered by `GCA_025200635.1` duplicated as both `562.109860` AND `562.111036` in the N=147 cipro cohort. See `tests/test_cohort_build_dedup.py` (6 regression tests). Lesson at `LESSONS_LEARNED.md` 2026-05-22.
- **`pipeline.py predict` emits `attribution_scope_confidence` field** (added 2026-05-22; pre-falsifier defaults to `INDETERMINATE`). Helper `_classify_attribution_scope(prefix, saturated, all_negative_delta, falsifier_verdict, falsifier_pass_passes_high) -> {HIGH, PARTIAL, INDETERMINATE}`. Locus-tag-prefix proxy (ERS = HIGH; ELX/ELY/ELV/ELU/ELT = PARTIAL; saturated/all-negative-delta = INDETERMINATE) derived from the 2026-05-21 cipro audit JSON. Lower fidelity than the future Mash-clade-based field (PASS-path artifact). 12 unit tests at `tests/test_pipeline_predict_v0.py`.
- **`pipeline.py predict` provenance fields RELOCKED 2026-05-23** (per Codex's v0 closeout on Precision 7780). Canonical fields: `cv_strategy` (e.g., `leave_one_accession_out`) + `cv_auroc` (primary CV AUROC for that strategy). Legacy `loso_auroc` preserved for older bundles; markdown sidecar renders BOTH when both are present. `reporting_mode = canonical_audit_aware` when `--audit-merge-json` is supplied, else `debug_internal`. Schema additive, backward-compatible — existing tests + older training pickles continue to work. v0 spec at `wiki/decoder_v0_ux_and_success_criterion.md` (RELOCKED 2026-05-23, narrower than the 2026-05-18 LOCKED version — surface is now cached-strain only, NOT genome-input).
- **v0 closed 2026-05-24** per the post-falsifier plan's FAIL branch. Bounded falsifier on Precision 7780 returned FAIL (ranking-only rescue didn't improve ELX-family failure cases). v0 shipped anyway with `wiki/dna_decoder_v0_closeout_handoff_2026-05-24.md` + scope-limit doc (lives on Precision 7780 as `reports/cipro_v0_scope_limit_decision_2026-05-23.md` — not yet on origin). Leakage-safe retrain on `leave_one_accession_out` CV: AUROC 0.8697. v0.1 first question: real-genome-input cipro decode (Path G) vs cef-cached expansion (Path C) — see `plans/v0.1_Ingestion_Contract_Plan.md`.
- **Cross-machine sync diagnostic at `scripts/cross_machine_sync_check.py`** (shipped 2026-05-24). Detects drift between origin and local across 5 axes (commit gap, working-tree dirtiness, Downloads/ recent artifacts, spec-divergence spot-check, pytest --collect-only). Run after any cross-machine handoff. 7 unit tests at `tests/test_cross_machine_sync_check.py`. `KNOWN_DIVERGENCE_TARGETS` pinned for the 2026-05-23 Codex relock (decoder spec RELOCKED, pipeline.py cv_strategy + leave_one_accession_out markers).
- **Drug-agnostic mechanism audit at `scripts/drug_mechanism_audit.py`** (shipped 2026-05-24). Generalized version of `scripts/cipro_mechanism_audit.py` using `mic_tiers.py` per-drug catalogs. Takes `--drug` arg. Enables cef + tet + gent audits without refactoring the cipro-hardcoded script (which has 80+ dedicated tests). Verdicts: PRIMARY_DOMINANT / MIXED_MECHANISMS / MOSTLY_UNKNOWN / EMPTY_R_SET. 14 unit tests at `tests/test_drug_mechanism_audit.py` (cipro QRDR / synonymous filter / class filter / drug-switch / POINTX dedup / cef bla CTX-M).
- **Bounded-falsifier coordination contract lives at `wiki/cipro_bounded_falsifier_coordination_plan_2026-05-22.md` + `wiki/cipro_bounded_falsifier_subset_2026-05-22.json`** (shipped 2026-05-22). Codex on Precision 7780 owns runner mechanics; Claude on the GTX 860M laptop owns subset selection + leakage check + diagnostic export schema. Runner draft at `scripts/cipro_bounded_falsifier.py` with 15 contract tests pinning the verdict matrix (PASS / FAIL / RUNNER_REGRESSION / INDETERMINATE_BUCKET_C). Leakage check at `scripts/leakage_check_dup_accession.py` (< 5s gate; blocks falsifier interpretation if `loso_leakage_present=True`). Post-falsifier ship-path technical plan at `plans/Cipro_Post_Falsifier_Ship_Path_Technical_Plan.md` (449 LOC, 4 verdict branches × 3 gate states pre-committed before results land per the verdict-vs-budget LESSON 2026-05-14).
- **Mash-cluster orchestration at `scripts/mash_cluster_n147.py`** (RECLASSIFIED 2026-05-24 from PASS-path artifact → v0.1 infrastructure). The post-falsifier plan tagged it as PASS-only; falsifier shipped FAIL → v0 closed WITHOUT firing Mash. But Mash-cluster N=147 is now load-bearing for BOTH v0.1 paths: (a) real-genome-input decode needs phylogenetic context for new strains; (b) cef-cached expansion needs clade-balanced cohort selection. Threshold sweep over `(0.02, 0.03, 0.04, 0.05, 0.07, 0.10)` picks lowest intra/inter variance ratio among thresholds satisfying min-clades ≥ 3 AND max-clade-fraction < 0.60. Fallback 0.05 if no candidate qualifies. Pure logic (`score_threshold` / `pick_best_threshold` / `per_clade_label_balance`) tested locally with 8 tests at `tests/test_mash_cluster_n147.py`; Docker-dependent Mash invocation reuses `dna_decode/eval/phylogeny.py::compute_mash_distances(..., use_docker=True)` on a Docker-equipped host (Precision 7780).
- **`scripts/pipeline.py predict` is v0-compliant** (rewrite 2026-05-18). Emits the v0 JSON + markdown sidecar schema per `wiki/decoder_v0_ux_and_success_criterion.md`: `prediction` + `calibrated_probability` + `confidence_tier` (HIGH/MEDIUM/LOW from `_confidence_tier(proba)`) + `top_k_attribution` (gene-level ISM + Tier 1-5 catalog labels) + `audit_verdict` (SUSPEND propagation from merge-gate JSON sidecar via `--audit-merge-json`) + `provenance` (model, training_cohort, loso_auroc, trained_on — populated from train pickle's enriched provenance fields). New CLI args: `--annotations`, `--card-path`, `--amrfinder-path`, `--audit-merge-json`, `--top-k`, `--output-md`, `--no-attribution`. 16 unit tests at `tests/test_pipeline_predict_v0.py` + 6 E2E integration tests at `tests/test_pipeline_predict_e2e.py` (uses synthetic fixtures — no Databricks cache needed for regression).

## Project workflow

Built using a personal Claude Code skill ladder:
- `/idea-anchor` → `/project-init` → `/brainstorm` ×N → `/technical-plan` → `/probe` → `/save-plan` → `/execute-plan`
- Project ledger maintained via `/project-state` skill at `project_state/dna-decode-2026-05-11.md`
- Execution state tracked at `.claude/execute-plan-state/Ecoli_G2P_Platform_Technical_Plan.json`
- Plan index at `wiki/plans-index.md`

Run `/brainstorm` between every wave; don't skip "to save time." Three brainstorm rounds caught 3 grounded contract gaps each — that's the pattern.
