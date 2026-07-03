# DMS variant-effect benchmark — definitive N (63 assays / 40 genes) + the modality×predictor interaction

**Date:** 2026-07-03
**Scripts:** `scripts/dms_variant_effect_benchmark.py` + `scripts/dms_alphamissense_benchmark.py`
(`wiki/dms_benchmark_big_scores.json`, `wiki/dms_am_big_scores.json`)
**Data (free):** 63 human MAVEDB DMS assays / 40 genes (BLOSUM); 48 assays / 34 genes with a UniProt+offset
join to AlphaMissense (Cheng 2023). Scales the 2026-07-02 N=11/17-gene benchmark to a definitive N.

## Headline (confirmed at 3–4× the genes)

| | median polarity-corrected Spearman |
|---|---|
| deterministic BLOSUM62 floor (63 assays / 40 genes) | **0.210** |
| AlphaMissense (48 assays / 34 genes, join match-gated) | **0.423** |
| measured gain (AM − BLOSUM, paired) | **+0.210 (~2×)** |

The last night's suggestive result (BLOSUM 0.22 / AM 0.52 / ~2× at N=11) holds **definitively** at 34–40 genes:
the deterministic substitution-severity floor is ~0.21, a free precomputed predictor (AlphaMissense, no GPU)
roughly doubles it.

## The new conclusion — AM's gain is MODALITY-DEPENDENT (needs the large N)

Per-modality, on the 44 paired assays with a known polarity:

| modality | n | BLOSUM floor | AM | **AM gain** |
|---|---|---|---|---|
| function | 3 | 0.182 | 0.539 | **+0.357** |
| abundance | 12 | 0.231 | 0.479 | +0.244 |
| binding | 12 | **0.268** | 0.396 | **+0.149** |
| other | 17 | 0.188 | 0.367 | +0.231 |

**An inverse relationship: AlphaMissense adds LEAST exactly where the deterministic floor is HIGHEST.**
Binding assays have the highest BLOSUM floor (0.268) and the smallest AM gain (+0.149); function assays have
the lowest floor and the largest gain. Mechanistically interpretable: a chemically-disruptive substitution at
a binding interface is partly captured by substitution-severity (so BLOSUM already works, AM adds little),
whereas position-specific catalytic-function effects are invisible to substitution type (so BLOSUM fails, and
the learned/evolutionary signal AM carries is what closes the gap).

- **Robust part:** abundance (n=12) vs binding (n=12) — both well-powered; binding floor higher, binding gain
  lower. Clean.
- **Under-powered:** the **function** row is n=3 (few pure-function assays cleared the polarity+join gates) —
  its +0.357 gain is directional, flagged, not definitive.

## What it means for the project

- The molecular-phenotype boundary is now mapped with a **definitive** N: deterministic substitution-severity
  is a real but ~0.21 floor; the deployable strong tool (AlphaMissense, free, no GPU) ~doubles it; and **the
  deterministic approach is *least* deficient for binding-interface effects, *most* deficient for
  position-specific function.** That tells you exactly where a deterministic decoder is "good enough" (binding /
  structural destabilization) vs where you must reach for the learned predictor (catalytic function).
- This directly motivates the **deterministic-conservation breakthrough candidate** (dedicated session): a
  position-specific conservation score (per-protein MSA) targets exactly the *function* deficit BLOSUM can't
  touch — the one place it could close the gap to AM without a learned model.

## Honest scope

- AlphaMissense is a learned/AF-distilled predictor (benchmarked as the strong usable tool; not the
  deterministic ethos). Gains are measured on the same assays, join match-gated (≥0.5), offset-corrected.
- 15 of 63 assays lacked a UniProt id/clean join → BLOSUM-only (not AM-joined); polarity-unknown assays
  excluded from medians.
- Supersedes (extends) `wiki/dms_variant_effect_benchmark_2026-07-02.md` +
  `wiki/dms_alphamissense_benchmark_2026-07-02.md` (N=11 versions). Frozen AMR surface byte-unchanged
  (leak guard 9/9).

## Reproduce

```bash
uv run python scripts/dms_variant_effect_benchmark.py --manifest D:/dna_decode_cache/proteingym/big_manifest.json --out wiki/dms_benchmark_big_scores.json
uv run python scripts/dms_alphamissense_benchmark.py --manifest D:/dna_decode_cache/proteingym/big_manifest.json --am D:/dna_decode_cache/proteingym/am_filtered_big.tsv --blosum wiki/dms_benchmark_big_scores.json --out wiki/dms_am_big_scores.json
```
