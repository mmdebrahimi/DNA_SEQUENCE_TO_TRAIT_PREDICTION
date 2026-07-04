# Phase 2 — supervised ESM head vs catalog vs zero-shot (2026-07-03)

**The decisive "does learning beat the curated catalog" test.** A supervised head
(StandardScaler→PCA(30)→balanced LogReg) on `[zero-shot LLR, BLOSUM62, ESM WT per-position embedding]`,
**leave-one-POSITION-out** CV (generalize to UNSEEN resistance sites, not memorize), scored against BOTH
baselines on the SAME held-out HIV single-mutant variants (n=293, 102 R). Script:
`scripts/esm_supervised_resistance_head.py`.

## Result — `PARTIAL`: beats zero-shot, does NOT beat the curated catalog

| set | n (R) | head AUC | zero-shot AUC | head balacc | **catalog balacc** |
|---|---|---|---|---|---|
| **pooled** | 293 (102) | **0.691** | 0.587 | 0.590 | **0.783** |
| NNRTI (RT/EFV) | 57 (33) | 0.348 | 0.244 | 0.388 | 0.790 |
| INSTI (IN/RAL) | 175 (54) | 0.758 | 0.724 | 0.596 | 0.769 |
| NRTI (RT/3TC) | 28 (14) | 0.707 | 0.821 | 0.750 | 0.750 |
| PI (PR/NFV) | 33 (1) | 0.031* | 1.0* | — | 0.938 |

*PI n_R=1 → degenerate, ignore.

## What it means (pre-committed verdict branch: PARTIAL)

- **The learned head BEATS zero-shot ESM** (pooled AUC 0.691 vs 0.587): training on fold labels extracts
  MORE resistance signal from the embeddings than the fitness-only zero-shot score, and it generalizes across
  unseen positions (leave-one-position-out). It even nudges the anti-predictive NNRTI zero-shot (0.244) up
  toward chance (0.348).
- **But it does NOT beat the deterministic catalog** (pooled balacc 0.590 vs 0.783). On KNOWN resistance
  positions, **curated determinant knowledge beats learning-from-embeddings** — the catalog already encodes
  exactly what the head must struggle to infer, and it wins on every powered class. The head does not recover
  the NNRTI pocket signal to a deployable level (balacc 0.388 vs catalog 0.79) — pocket-evasion resistance is
  not a cross-position-generalizable embedding pattern.
- **Per the plan's pre-committed PARTIAL branch:** the learned scorer is a **novel-variant FALLBACK only** —
  it may add value scoring variants the catalog is SILENT on (uncatalogued positions), but it is NOT a
  catalog replacement. Integrate cautiously in Phase 4 (fallback path, fail-closed to the catalog); do NOT
  route known-position calls through it.

## Consequence for the hybrid plan

- **V1 (learned scorer) = bounded fallback, not a replacement.** The deterministic catalog stays primary.
- **The grand "learned beats deterministic" vision does NOT hold** on independent-label, held-out-position
  tests — consistent with the whole project's finding that the binding constraint is curated knowledge +
  labels, not model capacity. This is the honest ceiling of the protein-level learned branch.
- **Phase 4 refinement (named):** the head's true value is on the catalog-SILENT subset — evaluate it there
  specifically (does it rank catalog-missed true-R above catalog-correct true-S?) before wiring any fallback.
- **Phase 5 (genotype world-model) stays deferred** — a protein-level supervised head already tops out below
  the catalog; a genotype-level model faces the (stronger) LD confound.

## Honest scope
- Leave-one-position-out is the HARD honest CV; a random split would inflate the head via position
  memorization. N=293 single-mutants, in-distribution HIVDB fold.
- The catalog is a strong baseline (known-DRM lookup); the fair comparison it wins is "curated vs learned on
  known positions" — it does not, by construction, cover novel positions (the head's intended niche).
- **verify-in-batch caught a real bug pre-verdict:** the catalog baseline was passing the fold COLUMN CODE
  (`EFV`/`3TC`) instead of the drug NAME to `call_hiv_observed` → INDETERMINATE → false all-S (balacc 0.5),
  which would have falsely inflated the head's win. Fixed (`_DRUG_NAME` map) before the reported numbers.
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9); ESM weights on D:.

## Reproduce
```bash
HF_HOME=D:/hf_cache TORCH_HOME=D:/dna_decode_cache/torch \
  uv run python scripts/esm_supervised_resistance_head.py --model esm2_t33_650M_UR50D
uv run pytest tests/test_esm_supervised_head.py -q   # 3 offline helper tests (no model)
```
