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
- [ ] GPU sanity check: confirm RTX 4090 has ≥24GB VRAM for 4-bit Evo + bitsandbytes installs cleanly (Linux/WSL2 only — Windows skips quantize extras).

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

## Pre-existing known limitations (not bugs)

- **Live BV-BRC API integration**: `pilot.fetch_bvbrc_drug_counts` raises `NotImplementedError` when no `--ast-tsv` flag / env var / config entry is provided. Live REST endpoint resolution deferred until first real-data run. Workaround: download an AST TSV/CSV from BV-BRC.
- **`fetch_ncbi_assembly_quality`** in `pilot.py` stays scaffolded intentionally. Phase 2 ships a separate CSV adapter (`dna_decode/data/bvbrc_genome.py`) that bypasses it via `pipeline ingest --assembly-metadata-csv`. Live NCBI Datasets REST integration is Phase 3 work.
- **Annotation source variance**: `parse_gff3` collapses `ID=` / `Name=` / `gene=` into one `gene_id` column. Annotation-source-aware extraction (gene_symbol as separate column) is Phase 2 cleanup work.
