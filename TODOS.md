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

## Pre-existing known limitations (not bugs)

- **Live BV-BRC API integration**: `pilot.fetch_bvbrc_drug_counts` raises `NotImplementedError` when no `--ast-tsv` flag / env var / config entry is provided. Live REST endpoint resolution deferred until first real-data run. Workaround: download an AST TSV/CSV from BV-BRC.
- **`fetch_ncbi_assembly_quality`** in `pilot.py` stays scaffolded intentionally. Phase 2 ships a separate CSV adapter (`dna_decode/data/bvbrc_genome.py`) that bypasses it via `pipeline ingest --assembly-metadata-csv`. Live NCBI Datasets REST integration is Phase 3 work.
- **Annotation source variance**: `parse_gff3` collapses `ID=` / `Name=` / `gene=` into one `gene_id` column. Annotation-source-aware extraction (gene_symbol as separate column) is Phase 2 cleanup work.
