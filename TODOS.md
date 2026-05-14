# TODOs

Deferred work and post-Phase-1 polish.

See also: `plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md` for the contracted ship-path plan (archived 2026-05-12); `FUTURE_FEATURES.md` for Phase 2+ ideas.

## Phase 2 — Entry criteria (gate before any Phase 2 work)

- [ ] Real-data smoke: one E. coli genome end-to-end with prediction + attribution captured (no ground-truth comparison required at entry)
- [ ] Full test suite passes on real deps (`uv run pytest tests/ -v` — 287 tests as of 2026-05-12)

## Phase 2 — Real-data validation (deferred from Phase 1)

- [ ] Real-data pilot gate against a downloaded BV-BRC AST TSV — confirm ≥150 strains per drug for cipro / ceftriaxone / tet after broth-microdilution + assembly-quality filters.
- [ ] End-to-end `pipeline.py ingest → train → attribute` on the pilot-validated cohort; record CV AUROC + Tier 1-3 fractions per drug.
- [ ] Mash CLI install on the target machine (Linux/WSL2 via `apt install mash`, OR `pyani` Python fallback if Mash unavailable on Windows).
- [x] GPU constraint observed 2026-05-14: actual hardware is GTX 860M (4 GiB Maxwell, CC=5.0, 2014). 4-bit Evo via bitsandbytes requires CC ≥ 7.0 — structurally unavailable. DNABERT-2 blocked by Triton 3.x state_dict incompat. NT v2 100M works locally (slow, ~9 min/strain). Larger-cohort runs (≥150 strains) gated on Databricks burst.

## Phase 2 — Starter genome set (infrastructure smoke; NOT a model-quality test)

For the one-genome / few-genome end-to-end smoke that gates Phase 2 entry. Useful for proving the ingestion → annotation → embedding → cache → train → attribute pipeline runs against real NCBI data. **Not** suitable as a model-quality test — cross-lineage R/S contrast confounds clade signature with resistance signal, which is exactly the trap the clade-only-baseline gate is supposed to detect.

| Strain | Assembly accession | Role |
|---|---|---|
| K-12 MG1655 | `GCF_000005845.2` | Reference baseline (non-pathogenic, well-annotated) |
| O157:H7 Sakai | `GCF_000008865.2` (verify via NCBI Datasets) | Pathogenic comparison |
| ST131 EC958 | `GCF_000285655.3` (verify via NCBI Datasets) | Multidrug-resistant clinical lineage |

Download via:
```bash
datasets download genome accession GCF_000005845.2
```

- [ ] Phase 2 real model-quality test needs **within-lineage** susceptible+resistant pairs (or clade-stratified cohort), NOT the cross-lineage trio above. Use the trio for infrastructure smoke only.

## Phase 2.5 perf hardening (deferred from Gate B prep)

- [ ] **`load_bvbrc_ast` is slow on real BV-BRC AMR exports** — Python pandas engine (sep=None) + iterrows() over 50K+ rows takes minutes. Audit cohort generator skips --ast on real data because of this. Fix: peek at first line to detect separator → use C engine; replace iterrows with vectorized ops. Surfaced 2026-05-13 during audit_cohort.py run.
- [x] **Embedding cache populate is per-sequence on GPU** — RESOLVED 2026-05-13. Added `_embed_window_batch` to `FoundationModel` base + override in `NucleotideTransformerModel` with mask-aware mean pooling. `cache.populate()` now chunks pending pairs into `EMBED_BATCH_SIZE=4` groups (tuned for 4 GiB VRAM on GTX 860M; raise to 32+ on larger GPUs). Regression test `test_nt_embed_window_batch_matches_per_sequence` confirms numerical equivalence with per-sequence path within rtol=1e-4. Empirical speedup on GTX 860M Maxwell was modest (~20%) — kernel-launch overhead is small fraction of compute-bound forward pass on this GPU; speedup would be 5-25× on Ada/Ampere hardware.
- [ ] **[BLOCKED] NT AutoModel swap** (`foundation.py:239`) — `AutoModel.from_pretrained` fails at load on the NT v2 100M `trust_remote_code` checkpoint with state_dict shape mismatch (`Linear[4096, 512]` vs `Linear[512, 2048]`). InstaDeep's `trust_remote_code` defines architecturally distinct AutoModel vs AutoModelForMaskedLM variants. Deferred indefinitely per `plans/Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md` Commit 2 gating rule.
  - Environment at failure: `transformers==5.8.0`, `torch==2.11.0`, `python==3.11.5`, os=`Windows 10`, cuda=`available`, NT revision=`f34324c6fde36a4f635f0f1f06cac5d25acd6798` (HF cache at `F:\hf_cache\hub\` via `HF_HOME` env var), loader=`AutoModel.from_pretrained(trust_remote_code=True)`.
  - Revisit when: InstaDeep ships an AutoModel-compatible checkpoint, OR project drops `trust_remote_code` in favor of a manual reimplementation.
- [ ] **[OPEN] NT hidden-state hardening** (`foundation.py:259`) — current code calls `model(**inputs, output_hidden_states=True)` then takes only `hidden_states[-1]`. Drop intermediate-layer materialization for a low-risk speedup. AutoModel swap is blocked (see above); diagnostic spike needed to check whether the loaded `AutoModelForMaskedLM` exposes a base-encoder accessor (`.base_model` / `.encoder` / `.nucleotide_transformer`) returning final hidden state directly. Defer until NT becomes critical-path during Phase 2 smoke.
- [ ] **[OPEN] Cache-integrity validation before `.complete` sentinel** (`refseq.py:_unpack_ncbi_datasets_zip` + `download_genome`) — current code writes extracted files non-atomically then drops `.complete` sentinel. `is_cache_complete()` checks only sentinel presence, not file content. Surfaced 2026-05-13: `GCA_003204155.1/genome.fna` had a 4,100-byte block of binary garbage starting at byte 1,159,168 (likely partial overwrite from a prior interrupted run); `.complete` was still present so subsequent `download_cohort_genomes` calls skipped it for weeks. Suggested fix BEFORE `.complete` is written: (a) `SeqIO.parse(fasta)` smoke-parse succeeds, (b) sequence alphabet ⊆ `{A,C,G,T,N}`, (c) file size within plausible range for the organism, (d) cross-check extracted `genome.fna` bytes against the genomic `.fna` member inside the cached `package.zip`. Recovered manually for the 12-strain mini-cohort run; permanent fix deferred.

## Phase 2 — Decision gate prerequisites and deferred classifier variants

Per `plans/Phase2_Decision_Gate_Plan.md` (saved 2026-05-13). Tiered N=50→N=150 staged path with Option-C threshold (≥5 pp AUROC delta AND top-K attribution includes gyrA / parC / parE).

- [ ] **[OPEN] Real clade-only baseline (Stage 1 prerequisite, NOT smoke gate)** (`scripts/pipeline.py:285`, `dna_decode/eval/clade_baseline.py`, `dna_decode/data/cohort.py`) — originally framed as smoke-gate prerequisite; **rescoped 2026-05-14** to Stage 1 N=50+ per B-B decision lock (smoke gate drops clade-only entirely; 12 unique MLST → singleton clades → AUROC ≈ 0.5 regardless of fix). Stage 1 scope:
  - **Replace `hash(s.mlst) % 10` placeholder** with `mlst_to_clade_id` helper in `dna_decode/data/cohort.py` (deterministic via `zlib.crc32`, scheme-aware; see `plans/Sidework_Sequence_Ship_Path_Plan.md` D7 pseudocode).
  - **Auto-disable with warning on sparse cohorts.** Empirically verified 2026-05-14 on `gate_b_cohort.parquet`: **61 unique MLSTs across 67 strains** (unique_mlst_fraction=0.910). At Stage 1 N=50 sampled from this, sparsity will remain ≥0.85 — clade-only stays degenerate even after the parsing fix. Pattern: `if unique_mlst_fraction > 0.8: logger.warning("Clade-only baseline disabled: insufficient lineage reuse (fraction={f:.2f} > 0.8)"); skip_clade_baseline = True`.
  - **Coarser binning escape hatch.** ST-level grouping is too fine for the cohort's MLST diversity; consider binning to (a) scheme name only (`ecoli_achtman_4` vs `Escherichia_coli_1` — only 2-3 groups), OR (b) Mash clusters at ANI≥99% (true lineage grouping). Without coarser binning, clade-only stays degenerate at any cohort size up to gate_b's 67-strain limit.
  - **Strain-keying semantics fix** in `predict_clade_only` / `per_clade_baseline` — current code assigns same clade AUROC to every strain in a clade group regardless of LOSO-fold structure. Separable concern flagged by /review 2026-05-13.
- [ ] **[OPEN] Add Random Forest wrapper to `classifiers.py`** — for Stage 2 decision gate (N=150). Field-standard downstream classifier on frozen DNA-FM embeddings per the 2025 Nature Communications benchmark (research_outputs/sota-bacterial-amr-prediction-small-cohorts-2026-05-13.md). Straightforward sklearn-style wrapper alongside existing XGBoost path.
- [ ] **[OPEN] Add TabPFN wrapper to `classifiers.py` with PCA-to-≤2000-feature step; pin package version** — for Stage 2 decision gate. NT outputs 512-dim → fits TabPFN-2.6's 2000-feature envelope directly, but PCA still recommended for stability + reproducibility (and to fit older TabPFN-2's 500-feature limit if package falls back). McElfresh NeurIPS 2024 shows TabPFN beats XGBoost on ≤1250-sample regime (rank 4.88 vs 8.30). Pin exact TabPFN version for reproducibility.
- [ ] **[OPEN] Add SNP-table feature variant — parse AMRFinderPlus output (`--organism Escherichia --mutation_all`) for POINT* rows** — for Stage 2 decision gate. AMRFinderPlus already emits point mutations as `gyrA_S83I` / `parC_S80I`-style rows (Method column = `POINT`, `POINTP`, `POINTX`, or `POINTN`) when run with the right flags. Treat as binary features alongside existing gene-presence. Cipro-specific classical baseline: gyrA / parC / parE rows are the textbook resistance signal. Scope: hours (parsing existing AMRFinder output), not days (no variant-calling pipeline needed).
- [ ] **[OPEN] XGBoost calibration overcorrection at N≤20 training** (`dna_decode/models/classifiers.py:train_xgboost_classifier(..., calibrate=True)`) — surfaced 2026-05-14 during 12-strain cipro smoke gate. At N=11 training (LOSO at N=12), `CalibratedClassifierCV` produces perfectly anti-predictive output (NT AUROC=0.000 vs 0.750 without calibration; output collapses to symmetric values around 0.5). Smoke gate runner now uses `calibrate=False`. Decision needed before Stage 1 N=50: (a) re-enable at N=50 where calibration CV has enough samples per fold, OR (b) require N≥30 before calibration kicks in (add a guard in `train_xgboost_classifier`). Related Wave 3.5 C7 fix already handles minority-class folds; this is a separate "calibrator's own training data is too small at LOSO N=12" issue.
- [ ] **[OPEN] Gene-presence + XGBoost returns AUROC=0.000 at N=12 even without calibration** (`scripts/smoke_gate_12strain_cipro.py:run_gene_presence_xgboost`) — surfaced 2026-05-14. Symmetric anti-predictive pattern (similar to the calibration bug) but `calibrate=False` is already set. Likely XGBoost overfitting at N=11 train × high-dim sparse binary features (~thousands of gene IDs from full Bakta-style annotation). NT and k-mer don't show this — NT is dense continuous, k-mer is dense count. Suggestions: (a) regularize XGBoost more aggressively (lower max_depth, higher min_child_weight, L2 reg) when feature_dim > n_train, (b) feature pre-select top-K most-variant genes before training, (c) use simpler classifier (logistic regression with L2) for binary-sparse features. Investigate before Stage 1 N=50 where gene-presence will run again.

## Pre-existing known limitations (not bugs)

- **Live BV-BRC API integration**: `pilot.fetch_bvbrc_drug_counts` raises `NotImplementedError` when no `--ast-tsv` flag / env var / config entry is provided. Live REST endpoint resolution deferred until first real-data run. Workaround: download an AST TSV/CSV from BV-BRC.
- **`fetch_ncbi_assembly_quality`** in `pilot.py` stays scaffolded intentionally. Phase 2 ships a separate CSV adapter (`dna_decode/data/bvbrc_genome.py`) that bypasses it via `pipeline ingest --assembly-metadata-csv`. Live NCBI Datasets REST integration is Phase 3 work.
- **Annotation source variance**: `parse_gff3` collapses `ID=` / `Name=` / `gene=` into one `gene_id` column. Annotation-source-aware extraction (gene_symbol as separate column) is Phase 2 cleanup work.
