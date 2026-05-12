# DNA_SEQUENCE_TO_TRAIT_PREDICTION

Genotype-to-phenotype (G2P) inference platform — predicts phenotypic traits from genomic DNA sequences AND identifies which genomic regions are most strongly associated with those predictions. Biologically interpretable, not causal-claim-making.

## Status: Phase 1 — Wave 1.5 hardened

Phase 1 scope ships in 17 implementation steps across 8 waves. **Currently shipped: Wave 0 + Wave 1 + Wave 1.5 hardening pass** (7/17 implementation steps complete; ~3700 LOC + tests; ~14 commits).

What's actually runnable today:

| Module | Status | Notes |
|---|---|---|
| `dna_decode/data/pilot.py` | **runnable** | `python -m scripts.pilot_gate --ast-tsv <bvbrc-ast-tsv-path>` runs the HARD-gate against a downloaded BV-BRC AST TSV |
| `dna_decode/data/refseq.py` | **runnable** | Downloads NCBI Datasets ZIP + unpacks into `genome.fna` / `annotations.gff3` / `annotations.gbk` |
| `dna_decode/data/annotations.py` | **runnable** | Parses GFF3 + GenBank; extracts CDS + intergenic sequences |
| `dna_decode/data/resistance_db.py` | **runnable** | Loads CARD JSON + AMRFinder TSV |
| `dna_decode/data/ast_data.py` | **runnable** | Loads BV-BRC AST TSV with broth-microdilution + organism filters |
| `dna_decode/models/foundation.py` | **scaffolded** | Wrappers for Evo + DNABERT-2 + Nucleotide Transformer + GENA-LM (lazy-load; needs HuggingFace download + GPU at first `embed()` call). MockFoundationModel works without weights/GPU. |
| `dna_decode/eval/` | **runnable** | LOSO + LOMO + leave-one-Mash-clade-out CV, metrics (AUROC/AUPRC/Brier/ECE/attribution-precision), Mash phylogeny clustering, clade-only baseline classifier |

What's NOT yet implemented (waves 2-7):

- Step 6 Strain/AST cohort catalog (drug-first; depends on Step 0.5 + 2 + 4 + 5)
- Step 8 HDF5 embedding cache
- Step 9 Baseline classifiers (XGBoost on frozen foundation embeddings)
- Step 11 In-silico mutagenesis (gene-level + saturation)
- Step 13 Genome-browser visualization (pygenometracks)
- Step 14 CLI entry points (ingest / train / predict / attribute)
- Step 15 Smoke pipeline + fixtures (the synthetic-data end-to-end smoke test referenced earlier was Wave 6 work; not yet shipped)
- Step 17 Leaderboard (foundation models + classical baselines on same cohort + CV protocol)
- Step 18 Classical baseline benchmark (AMRFinder + k-mer + gene-presence)
- Step 16 Documentation finalization

See `plans/Ecoli_G2P_Platform_Technical_Plan.md` for the full plan with Files + Depends on metadata.

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
| Compute | Single RTX 4090 with 4-bit Evo (Linux/WSL); A100 fallback |

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

# 5. Optional: install bitsandbytes for 4-bit Evo quantization (Linux/WSL only)
uv sync --extra quantize
```

## Running the pilot gate (Step 0.5 HARD gate)

The pilot gate validates that you have enough labeled isolates per drug BEFORE downstream ingestion fires.

```bash
# Download a BV-BRC AST TSV from ftp.bvbrc.org first, then:
uv run python -m scripts.pilot_gate \
  --drugs ciprofloxacin,ceftriaxone,tetracycline \
  --target-per-drug 150 \
  --ast-tsv path/to/bvbrc_ast.tsv

# Exit codes:
#   0 = GO (all per-drug counts >= target AND 3-drug intersection >= target)
#   1 = NO-GO (some count below target; see failure_reasons in the report)
#   2 = PilotGateError (config or filesystem issue)
#   3 = NotImplementedError (no AST source provided + live API path is deferred)
```

Alternative ways to provide the BV-BRC TSV:
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

See `project_state/dna-decode-2026-05-11.md` for full decision history (11 hypotheses, 4 decisions made, 14 action-log entries).
