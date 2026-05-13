# Gate A: One-Genome Plumbing Smoke (GREEN)

**Date:** 2026-05-12
**Outcome:** SUCCESS — ingest + parse + cache populate paths all work end-to-end on real NCBI data
**Runtime:** 5 seconds (post-download; download itself ~3 seconds for 5 MB ZIP)
**Strain:** E. coli K-12 substr. MG1655 (`GCF_000005845.2`)

## What was validated

| Layer | Path | Result |
|---|---|---|
| Network | `refseq.download_genome('GCF_000005845.2', ...)` → NCBI Datasets v2 API | ✅ 5.2 MB ZIP downloaded |
| Unzip | `_unpack_ncbi_datasets_zip` → `genome.fna` + `annotations.gff3` + `annotations.gbk` | ✅ all three files extracted |
| Manual cohort | `CandidateStrain` + `StrainCohort` + `save_cohort` → parquet | ✅ 1-strain cohort written |
| Annotation parse | `parse_gff3` on real RefSeq GFF3 (2.5 MB) | ✅ 4318 CDS rows parsed |
| CDS extract | `extract_cds_sequences(genome.fna, annotations)` | ✅ all 4318 CDS sequences extracted |
| Mock embedding | `MockFoundationModel.embed_batch()` × 4318 sequences | ✅ deterministic 128-d float32 vectors |
| HDF5 write | `EmbeddingCache.populate()` opens cache file ONCE, writes 4318 datasets | ✅ `gate_a_cache.h5` populated |
| Sample read | `f['strains/MG1655/cds-NP_414542.1']` (the first CDS — `thrL`) | ✅ shape=(128,), float32 |

## Reproducibility

```bash
cd C:/Users/Farshad/PythonProjects/dna_decode
$env:PATH = "$HOME\.local\bin;$env:PATH"   # if uv not on PATH

mkdir -p data/cache/refseq data/processed

# 1. Download MG1655 (idempotent — re-running short-circuits via .complete sentinel)
uv run python -c "from dna_decode.data.refseq import download_genome; download_genome('GCF_000005845.2', 'data/cache/refseq')"

# 2. Build a 1-strain cohort (assembly_accession is the load-bearing field)
uv run python -c "
from dna_decode.data.cohort import CandidateStrain, StrainCohort, save_cohort
s = CandidateStrain(strain_id='MG1655', assembly_accession='GCF_000005845.2', mlst='ST10', contig_count=1, n50=4641652, ast_labels={'ciprofloxacin': 0})
c = StrainCohort(strains=[s], per_drug_strain_ids={'ciprofloxacin': ['MG1655']}, three_drug_intersection=['MG1655'])
save_cohort(c, 'data/processed/gate_a_cohort.parquet')
"

# 3. Populate the embedding cache (mock model, no biology)
uv run python scripts/populate_cache.py \
  --cohort data/processed/gate_a_cohort.parquet \
  --cache data/processed/gate_a_cache.h5 \
  --refseq-cache data/cache/refseq \
  --model mock --allow-mock

# Expected output:
#   [populate_cache] loaded cohort: 1 strains
#   [populate_cache] resolved 1 strain(s); skipped 0
#   [populate_cache] model=mock (embedding_dim=128, device=cpu)
#   [populate_cache]   MG1655: wrote 4318 / 4318 embeddings
#   [populate_cache] DONE: 4318 new embeddings across 1 strains -> data\processed\gate_a_cache.h5
```

## What this does NOT validate

- **Foundation model biology** — mock embeddings are deterministic hash output; no biological signal. The plumbing layer is real; the embedding semantics are not.
- **Train / predict / attribute paths** — 1 genome cannot reach `pipeline train` (needs ≥10 samples + both R/S classes). Gate B exercises those paths with a 12-20-strain BV-BRC-derived cohort.
- **Phylogeny clustering** — Mash CLI isn't installed; the LOMO-clade-out CV path is untested on real data.
- **Quantization** — Evo 4-bit needs Linux/CUDA via WSL2; not exercised here.

## Bug fixed during Gate A

- `scripts/populate_cache.py:196` — unicode arrow `→` in the final success print broke Windows cp1252 console output (`UnicodeEncodeError`). Replaced with ASCII `->`. The 4318 embeddings had already been written successfully when the print crashed; the exit code was non-zero only because of the crash. **Fix included in the Gate A commit.**

## Cache footprint

| File | Size |
|---|---|
| `data/cache/refseq/GCF_000005845.2/genome.fna` | 4.7 MB |
| `data/cache/refseq/GCF_000005845.2/annotations.gff3` | 2.5 MB |
| `data/cache/refseq/GCF_000005845.2/annotations.gbk` | 11.9 MB |
| `data/cache/refseq/GCF_000005845.2/package.zip` (original) | 5.2 MB |
| `data/processed/gate_a_cohort.parquet` | ~3 KB |
| `data/processed/gate_a_cache.h5` (4318 × 128 × 4 bytes + HDF5 overhead) | ~3 MB |
| **Total** | **~28 MB** |

Trio (MG1655 + Sakai + ST131) projected: ~85-100 MB. Well within the 35 GB free on C:.

## Gate A status: PASS

The plumbing for real E. coli genome ingest → annotation parse → CDS extract → cache populate works. Phase 2 entry's first criterion — "real-data smoke: one E. coli genome end-to-end with prediction + attribution captured" — is partially satisfied (ingest + embed; train + predict + attribute are Gate B scope).

## Next gate

**Gate B (labeled mini-cohort dry-run)** — requires BV-BRC AMR Phenotypes TSV download. Pending user action per `wiki/PHASE2_PREFLIGHT.md`.
