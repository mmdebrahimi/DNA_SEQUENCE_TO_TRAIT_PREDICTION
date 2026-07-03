# AlphaMissense vs the deterministic floor — the measured deterministic→strong-predictor gap (molecular phenotype)

**Date:** 2026-07-02
**Script:** `scripts/dms_alphamissense_benchmark.py` (`wiki/dms_alphamissense_benchmark_scores.json`)
**Data (free):** AlphaMissense precomputed scores (Cheng et al. 2023, *Science*; Zenodo 8208688, 1.15 GB, all
human missense) joined to the 24-assay MAVEDB DMS benchmark. 13 assays / 11 genes had a UniProt id + offset →
**offset-corrected join match_rate = 1.0 on every assay.** On D:.

## What this measures

The multi-gene DMS benchmark established the **deterministic** floor (BLOSUM62 ≈ 0.22 median) and *cited*
published learned models at ~0.4–0.5 — not run (GPU). This **turns that cited contrast into a measured one**
using AlphaMissense: a free, **precomputed** variant-effect predictor (no GPU) covering all human missense
variants. AM is a *learned* (AlphaFold-distilled) predictor — **not** the project's deterministic ethos; it is
benchmarked here as the "strong usable predictor" to size the gap. Both are put on the same
"positive = predicts functionality" axis (AM is a pathogenicity score → sign-aligned via each assay's nonsense
polarity anchor + UniProt-offset position mapping).

## Result — AlphaMissense ≈ 2.1× the deterministic BLOSUM floor

| | polarity-corrected Spearman (median, 12 clean assays) |
|---|---|
| **AlphaMissense (precomputed, free)** | **0.515** |
| BLOSUM62 deterministic floor (paired) | 0.240 |
| **measured gain (AM − BLOSUM)** | **+0.265** |

Per-assay (AM vs BLOSUM, both polarity-corrected; join match_rate = 1.0):

| gene | modality | n | **AM** | BLOSUM |
|---|---|---|---|---|
| CYP2C19 | abundance | 7,830 | **+0.676** | +0.339 |
| p53 | function | 4,069 | **+0.671** | +0.222 |
| NUDT15 | function | 2,934 | **+0.602** | +0.155 |
| CYP2C9 | abundance | 6,370 | **+0.598** | +0.333 |
| TPMT | abundance | 3,689 | **+0.551** | +0.240 |
| PTEN | function | 7,260 | **+0.539** | +0.182 |
| MLH1 | abundance | 4,802 | **+0.515** | +0.273 |
| ASPA | abundance | 5,843 | **+0.486** | +0.274 |
| PTEN | abundance | 5,083 | **+0.477** | +0.217 |
| STK11 | function | 6,021 | **+0.372** | +0.085 |
| KRAS | binding | 1,190 | +0.340 | +0.274 |
| KRAS | abundance | 1,188 | +0.118 | +0.200 |

(MSH2 excluded — polarity-unknown, no nonsense anchor → sign ambiguous.)

## What it means for the project

- **The molecular-phenotype modality has a usable, free, no-GPU strong predictor** — AlphaMissense recovers
  measured protein-variant effect at ~0.52 median, ~2× the deterministic substitution floor. Where the earlier
  benchmark could only *cite* the learned ~0.4–0.5, this **measures** it (0.515) and identifies the deployable
  tool. AM is strongest exactly where BLOSUM is weakest — **function** assays (p53 +0.67, NUDT15 +0.60, PTEN
  function +0.54): the position-specific signal BLOSUM can't see is precisely what the distilled model adds.
- **The boundary map is now fully quantified with hard numbers:**

  | regime | winner | measured |
  |---|---|---|
  | curated high-effect catalog | **deterministic** | AMR/PGx/ben-1 (sens/spec) |
  | organism-level polygenic | *neither* (learned 0-for-5 de-confounded) | yeast…DGRP |
  | **molecular property, no catalog** | **learned/distilled** | **AM 0.52 ≫ BLOSUM 0.24 (this benchmark)** |

- **PGx bridge, sharpened:** CYP2C19/CYP2C9/TPMT/NUDT15 — the very genes the project's PGx star-allele decoder
  calls at the diplotype level — are here scored per-variant at AM ~0.55–0.68. A future PGx v0.2 could use AM as
  a per-variant functional prior for the non-core alleles the deterministic star-allele proxy currently withholds.

## Honest scope

- **AlphaMissense is a learned/AlphaFold-distilled predictor**, benchmarked here as the strong usable tool — it
  is explicitly NOT claimed as deterministic. The finding is "the strong precomputed predictor ~doubles the
  deterministic floor," measured on the same assays.
- The join is **offset-corrected** (UniProt pos = DMS pos + MAVEDB offset); match_rate = 1.0 on all 13 assays
  validates the mapping. 11 of 24 assays lacked a UniProt id in MAVEDB metadata → not joined (not a failure,
  just missing metadata).
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9).

## Reproduce

```bash
# AlphaMissense_aa_substitutions.tsv.gz (Zenodo 8208688) on D:; filter to the benchmark UniProts, then:
uv run python scripts/dms_alphamissense_benchmark.py
uv run pytest tests/test_dms_alphamissense_benchmark.py -q   # 3 offline synthetic tests
```
