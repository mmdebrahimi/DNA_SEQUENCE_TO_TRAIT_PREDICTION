# Full-manifest ESM2-650M — the definitive large-N R2 prospective number (2026-07-22)

**Status:** ✅ the definitive, large-N prospective result for the R2 molecular cell. ESM2-650M over **2383**
leakage-free held-out MaveDB DMS assays (842 genes, 978 human) → **median |Spearman| = 0.478** (shuffled
control 0.019). Ran overnight as 2 concurrent Kaggle T4 shards; merged 2026-07-22. Frozen surface byte-unchanged.

## Result
| metric | value |
|---|---|
| assays scored | **2383** (of 2569 held-out; 80 transient fetch-fails over the ~6h run, 22 too_long >1022aa, 84 too_few <20 missense) |
| distinct genes | **842** |
| **median \|Spearman\|** | **0.478** |
| shuffled negative control | **0.019** (near-zero → the signal is real) |
| **human-only median** | **0.492** (978 human assays) |
| top per-gene | MAGI3 0.786 · RSP5 0.752 · WWTR1 0.748 · TP63 0.743 · NEDD4 0.732 · NCF1 0.733 … |

## Reading it honestly
- **This SUPERSEDES the 86-assay 0.503 as the headline.** The 86-assay single page was a slightly favorable
  small sample; the tight large-N estimate over 2383 assays is **0.478** — and it still matches the published
  ESM2-650M ProteinGym field number (~0.48). Small-sample optimism, corrected by scale.
- **Leakage-free by construction** (target genes NOT in the ProteinGym benchmark ESM2/the hybrid was tuned on)
  and **confound-free** (R2 = designed mutant libraries, no population-structure axis). So 0.478 is real
  generalization at scale, not benchmark overfit.
- **Human at scale:** 0.492 median over 978 human-protein assays — the north-star "decode human" crossing, now
  measured on ~1000 human proteins, for free, in the regime where learned models work.
- The prior narrower results stand as consistent slices: pharmacogene subset 0.547 (16 assays), the 84-assay
  90%-paired-BLOSUM-beat (p=5e-15). The regime claim (learned beats the cheap baseline on fitness-aligned
  molecular properties) holds from n=16 to n=2383.

## Provenance
- Overnight run: 2 Kaggle T4 kernels `mavedb-full-esm2-shard{a,b}` (1285 + 1284 assays), masked-marginal core
  byte-faithful to `scripts/esm_zeroshot_dms.py`. Merged by `scripts/mavedb_full_esm2_collect.py`.
- Manifest: `wiki/mavedb_prospective_holdout_full_2026-07-21.json` (2569 held-out). Merged JSON:
  `wiki/mavedb_full_esm2_2026-07-22.json`. Regime lens: `plans/Trait_Decoding_Roadmap.md`.
