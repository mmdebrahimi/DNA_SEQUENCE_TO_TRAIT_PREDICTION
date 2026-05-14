# Phase 2 Preflight Status

**Captured:** 2026-05-12 (post Phase 1 ship)
**Host:** Windows 10, `C:\Users\Farshad\PythonProjects\dna_decode`
**Purpose:** machine-state snapshot before Path 1 (`populate_cache.py`) and Gate A (one-genome smoke).

## Disk

| Mount | Total | Used | Free | % Used |
|---|---|---|---|---|
| `C:\` | 237 GB | 202 GB | **35 GB** | 86% |

**Concern:** Phase 1 README projects ~25 GB for full real-data run (foundation-model weights + strain genomes + embedding cache). 35 GB free is uncomfortably close. **Recommendation:** route caches to external storage before any HuggingFace download.

```powershell
# Recommended cache routing (Windows PowerShell). Replace E: with your drive.
[Environment]::SetEnvironmentVariable("HF_HOME", "E:\hf_cache", "User")
[Environment]::SetEnvironmentVariable("DNA_DECODE_CACHE_ROOT", "E:\dna_decode_cache", "User")
```

## Environment variables

| Var | Status | Purpose |
|---|---|---|
| `HF_HOME` | **unset** | HuggingFace tokenizer + model cache location |
| `DNA_DECODE_CACHE_ROOT` | **unset** | RefSeq + cohort + embedding cache root |
| `BVBRC_AST_TSV` | **unset** | Path to BV-BRC AMR TSV (alternative to `--ast-tsv` flag) |

Action: set `HF_HOME` and `DNA_DECODE_CACHE_ROOT` to non-`C:` paths before downloading any foundation-model weights.

## External data

| Asset | Status | Action |
|---|---|---|
| **BV-BRC genome metadata CSV** | ✅ present at `C:\Users\Farshad\Downloads\BVBRC_genome.csv` (2157 rows; 152 E. coli) | Use for the scaffolded `pilot.fetch_ncbi_assembly_quality` path (CheckM completeness / N50 / contig count are present) |
| **BV-BRC AMR Phenotypes TSV** | ❌ NOT downloaded — required for pilot gate + cohort + train | Go to https://www.bv-brc.org/ → search "Escherichia coli" → **AMR Phenotypes tab** (not Genomes) → export. Expect ~30K-60K strain×drug rows. |
| **CARD catalog** | ❌ not cached at `data/cache/card/` | Auto-fetched by `resistance_db.py` on first use (URL in `config/datasources.yaml`) |
| **AMRFinder catalog** | ❌ not cached at `data/cache/amrfinder/` | Same — auto-fetched on first use |
| **NCBI Datasets reachability** | ✅ via Python REST in `refseq.py` (no CLI needed) | `datasets` CLI not on PATH but unused — REST API is the canonical path |

**CSV/TSV gap:** `dna_decode/data/ast_data.py:96` hardcodes `sep="\t"` (TSV). BV-BRC defaults to CSV export. **Path 1 will add CSV/TSV auto-detect** — one-line fix using `sep=None, engine="python"`.

## Compute

| Asset | Status | Notes |
|---|---|---|
| Python | ✅ 3.11.5 in `.venv/` | matches `.python-version` pin |
| `uv` | ✅ 0.11.8 at `~/.local/bin/uv.exe` | not on bash `$PATH` by default; PowerShell finds it |
| `torch` | ✅ 2.11.0 (CPU build) installed | usable for DNABERT-2 / Nucleotide Transformer / GENA-LM (CPU-feasible) |
| `transformers` | ✅ 5.8.0 | HuggingFace model loading works |
| `huggingface-hub` | ✅ 1.14.0 | download path works |
| `xgboost` | ✅ 3.2.0 | classifier head ready |
| `biopython` | ✅ 1.87 | annotation parsing ready |
| `h5py` | ✅ 3.16.0 | embedding cache backend ready |
| GPU / CUDA | ✅ GTX 860M, CUDA 11.8 (4 GiB VRAM, Maxwell, CC=5.0) | NT v2 100M works; 4-bit Evo via bitsandbytes requires CC ≥ 7.0 → unavailable on this GPU |
| WSL2 | ✅ installed (Docker Desktop distros visible) | available if Evo path becomes a Phase 2 priority |
| `bitsandbytes` | n/a Windows | Linux/CUDA only; install inside WSL2 if needed |
| Mash CLI | ❌ not on PATH | needed for `phylogeny.py` clustering (Step 10); install via `apt install mash` in WSL2 or manual binary on Windows |

## Test gate

| Item | Status |
|---|---|
| Phase 1 full test suite | ✅ 287 passed / 0 failed / 1 skipped (`test_baseline.json` recorded) |
| Synthetic smoke pipeline | ✅ runs end-to-end on 12-strain fixture |
| Phase 1 ship tag | ✅ `phase-1-shipped` pushed to origin |

## Foundation model decision tree

For Gate A and Gate B (Phase 2 entry smoke), pick a model based on what's reachable today:

| Model | Params | Setup needed | Feasible today? | Use for |
|---|---|---|---|---|
| **MockFoundationModel** | 0 (deterministic hash) | none | ✅ yes | Plumbing-only smoke — confirms ingest/cache/train/attribute wiring works on real metadata; biology is meaningless |
| **DNABERT-2** | 117M | HuggingFace download (~500 MB) | ✅ yes on Windows CPU | First real-embedding smoke; CPU-feasible without WSL2 |
| **Nucleotide Transformer** | 100M (v2) | HuggingFace download (~400 MB) | ✅ yes on Windows CPU | Alternative real-embedding smoke |
| **GENA-LM** | ~110M | HuggingFace download (~450 MB) | ✅ yes on Windows CPU | Alternative real-embedding smoke |
| **Evo 4-bit** | 7B | WSL2 + CUDA GPU + bitsandbytes install | ⚠️ needs WSL2 setup | Highest-capacity option; deferred until earlier models prove pipeline |

**Recommendation:** mock for Gate A (plumbing only), DNABERT-2 for Gate B (first real-embedding cohort run). Evo is Phase 2 / Phase 3.

## Path 0 blockers / unblockers summary

**No blockers for Path 1 (cache populate driver).** Can start immediately.

**Two blockers for Gate A (one-genome smoke):**
1. `pipeline ingest` works without BV-BRC AMR TSV if invoked with just `--accessions GCF_000005845.2` (need to verify CLI accepts this path).
2. `populate_cache.py --allow-mock` must be built first (Path 1).

**Three blockers for Gate B (mini-cohort dry-run):**
1. BV-BRC AMR TSV download (user action).
2. CSV/TSV auto-detect in `ast_data.py` (Path 1 scope).
3. `populate_cache.py` (Path 1 scope).

**Suggested env-var setup before any download:**
```powershell
mkdir E:\hf_cache E:\dna_decode_cache  # or another non-C: drive
[Environment]::SetEnvironmentVariable("HF_HOME", "E:\hf_cache", "User")
[Environment]::SetEnvironmentVariable("DNA_DECODE_CACHE_ROOT", "E:\dna_decode_cache", "User")
```

## Recommended next-action order

1. **User action:** download BV-BRC AMR Phenotypes for E. coli (parallel with Path 1 work).
2. **User action:** set `HF_HOME` + `DNA_DECODE_CACHE_ROOT` env vars to non-`C:` paths.
3. **Code work:** Path 1 — `scripts/populate_cache.py` standalone driver + CSV/TSV auto-detect in `ast_data.py`.
4. **Smoke run:** Gate A — MG1655 ingest + mock embed.
5. **Smoke run (gated on user BV-BRC download):** Gate B — mini-cohort + DNABERT-2.
