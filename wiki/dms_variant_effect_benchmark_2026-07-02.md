# Multi-gene DMS variant-effect benchmark — the deterministic floor across the human molecular-phenotype landscape

**Date:** 2026-07-02
**Script:** `scripts/dms_variant_effect_benchmark.py` (`wiki/dms_variant_effect_benchmark_scores.json`)
**Data (free, unit-level):** 24 human Deep Mutational Scanning assays from **MAVEDB**, 17 genes, 3 phenotype
modalities (abundance / function / binding), ~180k scored missense variants total. On D:.

## What this establishes

Extends the single PTEN cell to a **systematic benchmark**: does a deterministic substitution-severity score
(BLOSUM62) capture measured protein-variant effect, across many human proteins, and does it depend on the
phenotype modality? This maps the deterministic floor of the molecular-phenotype regime (the tractable
"more human G2P" frontier).

**Methodology fix (load-bearing):** DMS assays have opposite polarities — high score = functional in some,
= damaging in others. Each assay is polarity-anchored by its **nonsense variants** (maximally damaging), and
the reported Spearman is signed so **positive always means "conservative substitution → preserved function"** —
making all 24 assays comparable + aggregatable. (2 assays lacked a nonsense anchor → polarity=unknown, excluded
from the aggregate.)

## Result

| | polarity-corrected BLOSUM62 Spearman |
|---|---|
| **overall median (22 assays)** | **0.219** |
| abundance (10 assays) | **0.243** |
| binding (3 assays) | 0.274 |
| function (9 assays) | **0.163** |

**Top-line finding — the deterministic substitution-severity rule captures a *consistent, modest* fraction of
protein-variant effect (~0.22 median), and it is MODALITY-DEPENDENT:** it predicts **abundance** (structural
destabilization) better than **function** (0.24 vs 0.16). That is mechanistically sensible — a chemically
disruptive substitution destabilizes a fold (abundance) more reliably than it perturbs a position-specific
catalytic/interface role (function).

Per-assay highlights (all positive except a few near-zero/negative function assays):

| assay | modality | n | Spearman |
|---|---|---|---|
| CYP2C19 | abundance | 7,830 | +0.339 |
| CYP2C9 | abundance | 6,370 | +0.333 |
| CYP2C9 | function | 6,142 | +0.320 |
| MLH1 / ASPA / TPMT / SOD1 | abundance | 3–6k | +0.24–0.27 |
| PTEN | abundance | 5,083 | +0.217 |
| LDLR | abundance | 15,011 | +0.073 |
| TP53 | function | 7,446 | −0.067 |

The PGx-relevant genes (CYP2C19, CYP2C9, TPMT, NUDT15, GCK) are a **bridge to the project's PGx lane** — the
same genes the deterministic star-allele decoder calls at the *diplotype* level are here scored at the
*per-variant molecular* level.

## The boundary, now quantified across 17 genes

This is the **deterministic floor** of the molecular-phenotype modality. The published ProteinGym leaderboard
puts zero-shot learned models (ESM-1v / EVE / TranceptEVE) at **~0.4–0.5 Spearman** on such assays (cited,
**not** run here — GPU) — roughly **2× the BLOSUM floor**. So the boundary the project has been mapping now has
a systematic, multi-gene anchor:

| regime | winner | evidence |
|---|---|---|
| curated high-effect catalog | **deterministic** | AMR, TB, fungal, HIV, SARS, PGx, ben-1 |
| organism-level polygenic | *neither* (learned 0-for-5 de-confounded) | yeast, Arabidopsis, cipro, pathotype, DGRP |
| **molecular property, no catalog** | **learned** (~0.45) ≫ deterministic (~0.22) | **this 24-assay benchmark** |

## Honest scope

- **BLOSUM62 is a deliberately simple deterministic baseline** (BLOSUM45/80 give near-identical numbers — in
  the JSON). Position-specific conservation would do better but needs per-protein alignments; a distilled
  predictor (AlphaMissense) would approach the learned range — both are follow-ups.
- **The learned ~0.45 is cited from ProteinGym, not reproduced here** (GPU). No fabricated figure.
- 2 assays (MSH2, NUDT15) lacked a clean nonsense polarity anchor → excluded from the aggregate (flagged in
  the per-assay table).
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9).

## Reproduce

```bash
uv run python scripts/dms_variant_effect_benchmark.py   # reads D:/dna_decode_cache/proteingym/benchmark_manifest.json
uv run pytest tests/test_dms_variant_effect_benchmark.py -q   # 2 offline synthetic tests
# assays fetched from https://api.mavedb.org/api/v1/score-sets/<urn>/scores
```
