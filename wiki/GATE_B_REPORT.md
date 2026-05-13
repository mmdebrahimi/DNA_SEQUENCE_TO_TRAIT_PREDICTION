# Gate B Report — Mini-Cohort Real-Data Dry-Run

**Status:** ⏸ pending populate restart (external storage disconnected mid-run on full 67-strain attempt 2026-05-13 ~10:30 UTC; pivoting to 12-strain mini-cohort per /brainstorm Option C)

**Cohort:** 12 strains (6 cipro-R + 6 cipro-S), all 12 distinct MLST clades, all assemblies pass Phase 1 QC (contig_count ≤ 5, N50 ≥ 324 K for one strain — others ≥ 4.6 Mbp).
**Embedding model:** Nucleotide Transformer v2 100M multi-species (`InstaDeepAI/nucleotide-transformer-v2-100m-multi-species`, 512-dim).
**GPU:** NVIDIA GeForce GTX 860M (4 GB VRAM, compute capability 5.0, Maxwell), CUDA 11.8 via torch 2.7.1+cu118.
**Storage:** `F:\dna_decode_cache\` (external 4 TB drive — see Reliability Notes below).
**Gate B scope:** infrastructure dry-run — can the real-data pipeline run ingest → embed → train → predict → attribute end-to-end without crashing? Not a model-quality test.

---

## Selected mini-cohort

| Strain ID | MLST | contig_count | N50 | cipro |
|---|---|---|---|---|
| 562.28805 | ST410 | 2 | 4,695,904 | R |
| 562.30362 | ST4 | 4 | 5,213,052 | R |
| 562.28563 | ST156 | 4 | 4,877,957 | R |
| 562.17621 | ST167 | 4 | 4,864,149 | R |
| 562.17721 | ST1284 | 4 | 4,746,578 | R |
| 1328433.3 | **ST131** | 5 | 5,076,304 | R |
| 562.28389 | ST1276 | 2 | 4,987,160 | S |
| 562.7627 | ST4554 | 2 | 661,558 | S |
| 562.16326 | ST1809 | 3 | 4,658,583 | S |
| 562.16325 | ST1408 | 3 | 4,595,238 | S |
| 562.7641 | ST5543 | 3 | 324,659 | S |
| 562.52722 | ST29 | 4 | 5,570,476 | S |

Notable: ST131 (resistant) is the globally dominant fluoroquinolone-resistant E. coli clade — appropriate positive control.

---

## Execution plan (post-F:-reconnect)

```bash
# Pre-requisite: F: drive remounted; F:\dna_decode_cache\refseq\ may need re-download
# of the 12 mini-cohort genomes.

cd C:/Users/Farshad/PythonProjects/dna_decode
export HF_HOME=F:/hf_cache
export PATH=~/.local/bin:$PATH

# 1. Re-ingest (no AST processing — just genome downloads for the 12 strains).
#    Cohort already exists; we only need the genome files.
uv run python -m scripts.pipeline ingest \
  --drugs ciprofloxacin \
  --ast-tsv "C:/Users/Farshad/Downloads/BVBRC_genome_amr (1).csv" \
  --assembly-metadata-csv "C:/Users/Farshad/Downloads/BVBRC_genome (1).csv" \
  --target-per-drug 6 \
  --intersection-target 0 \
  --cohort-out "data/processed/gate_b_mini_cohort_regen.parquet" \
  --download-genomes

# Or shortcut: download_cohort_genomes called directly with the 12-strain parquet.

# 2. Populate cache with NT on GPU (~30 min expected).
uv run python scripts/populate_cache.py \
  --cohort data/processed/gate_b_mini_cohort.parquet \
  --cache F:/dna_decode_cache/embeddings/nt_mini.h5 \
  --refseq-cache F:/dna_decode_cache/refseq \
  --model nucleotide_transformer \
  --device cuda

# 3. Train cipro classifier.
uv run python -m scripts.pipeline train \
  --drug ciprofloxacin \
  --model nucleotide_transformer \
  --cohort data/processed/gate_b_mini_cohort.parquet \
  --cache F:/dna_decode_cache/embeddings/nt_mini.h5 \
  --include-clade-baseline \
  --min-auroc 0.5

# 4. Predict on a held-out strain.
uv run python -m scripts.pipeline predict \
  --model-path data/processed/models/ciprofloxacin_nucleotide_transformer.pkl \
  --strain-id 1328433.3 \
  --cache F:/dna_decode_cache/embeddings/nt_mini.h5

# 5. Attribute (ISM + Tier 1-5) on the same strain.
uv run python -m scripts.pipeline attribute \
  --model-path data/processed/models/ciprofloxacin_nucleotide_transformer.pkl \
  --strain-id 1328433.3 \
  --cache F:/dna_decode_cache/embeddings/nt_mini.h5 \
  --annotations F:/dna_decode_cache/refseq/<accession>/annotations.gff3 \
  --output data/processed/gate_b_attribution.json
```

---

## Results (to fill in after run completes)

### Step 1 — Genome download
- Wallclock: TBD
- Strains downloaded: TBD / 12
- Disk: TBD MB on F:

### Step 2 — Embedding cache populate
- Wallclock: TBD
- Embeddings written: TBD
- HDF5 size: TBD MB
- GPU utilization peak: TBD
- Issues encountered: TBD

### Step 3 — Train
- LOSO CV AUROC: TBD
- Clade-only baseline AUROC: TBD
- Validation gate: TBD (PASS / FAIL — note Gate B isn't gated on this)
- Issues: TBD

### Step 4 — Predict
- Held-out strain: 1328433.3 (ST131, R)
- Predicted P(R): TBD
- Issues: TBD

### Step 5 — Attribute
- Top-K (K=20) genes with highest |delta|: TBD
- Tier 1-3 hits: TBD / 20
- Known cipro resistance loci present in top-K: TBD (gyrA / parC / qnrA / qnrB / qnrS / aac(6')-Ib-cr)
- Issues: TBD

---

## What this report does and does NOT claim

**Does claim:**
- The Phase 2 pipeline can run end-to-end on a real (small) labeled E. coli cohort
- The data ingestion + cohort assembly + genome download + NT embedding + classification + attribution path is wired correctly

**Does NOT claim:**
- Real biological insight — 12-strain cohort is too small for statistical claims about cipro resistance mechanisms
- Reproducible AUROC numbers — the cohort is hand-picked top-assembly-quality strains, not a random sample
- Production-quality attribution — Tier 1-3 overlap at K=20 is informational only on a 12-strain cohort
- Anything about ceftriaxone or tetracycline — those drugs are out of mini-cohort scope

---

## Reliability notes

**External storage F: disconnected during first attempt (2026-05-13):**
- Full 67-strain NT populate ran for ~2 hr wallclock before F: drive vanished
- HDF5 write failures cascaded; cache file lost
- ~430 MB of downloaded genomes also on F: — also lost when drive disconnected
- HF model weights (1.5 GB) on F: — also lost (will redownload)

**Mitigations applied:**
- `scripts/populate_cache.py` now emits stderr-flushed progress per strain (visible recovery checkpoint even on stdout buffering / disconnect)
- 12-strain mini-cohort scope reduces single-run blast radius from 2.5 hr → 30 min
- Cohort manifest stays on C: (`data/processed/gate_b_mini_cohort.parquet`) so a disconnect doesn't lose the strain-selection state

**Open follow-ups (post Gate B):**
- Cache versioning with `pooling_strategy` field (from /brainstorm D-2) — prevents hybrid-cache bugs if NT pooling is ever refactored to batched mask-aware
- True batched inference in `embed_batch` + `cache.populate` collection — Phase 2.5 hardening (~2-3 hr refactor with proper equivalence tests)
- Optional `output_hidden_states=True` removal (low-risk speedup; see /brainstorm M2)

---

## /brainstorm decision log

Two adversarial rounds during Gate B execution:
1. After full 67-strain cohort built: chose Option C (mini-cohort) over Option B (refactor batched inference) — Option B has hidden pooling-semantics + hybrid-cache risk that Mock-only tests don't catch
2. Pre-existing M2 (drop `output_hidden_states=True`) deferred to Phase 2.5 alongside the batching refactor — same review surface, same correctness considerations

Reference: in-conversation /brainstorm Round 1 + Round 2 against `dna_decode/models/foundation.py` + `dna_decode/models/cache.py` paths.
