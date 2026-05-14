# Smoke Gate — 12-strain cipro cohort (2026-05-14)

> **This is an engineering smoke / falsification gate, NOT a powered classifier comparison.**
> At N=12 with balanced LOSO the per-strain noise floor is ±8.3% and 95% AUROC CI width is ~±0.19.
> The 15-percentage-point gap acceptance bar is an engineering heuristic to catch "NT obviously broken" —
> NOT statistically powered ranking. Multiple-comparison statistical power at N=12 forbids classifier ranking.
> Real decision gate scheduled at Stage 1 N=50 (local engineering screen) → Stage 2 N=150 (Databricks).

**Cohort:** `data\processed\gate_b_mini_cohort.parquet` (12 strains, 6R/6S cipro)
**Drug:** ciprofloxacin
**Gap threshold:** ≥15 pp (NT-XGBoost AUROC ≥ best-classical AUROC − 15 pp)
**Verdict:** PASS (NT not obviously worse)

## Per-variant LOSO results

| Variant | AUROC | AUPRC | N | Label balance |
|---|---:|---:|---:|---|
| NT-XGBoost (nucleotide_transformer) | 0.750 | 0.692 | 12 | 6R / 6S |
| k-mer (k=8) + XGBoost | 0.694 | 0.663 | 12 | 6R / 6S |
| Gene-presence + XGBoost | 0.000 | 0.394 | 12 | 6R / 6S |

## Gap analysis

- NT-XGBoost AUROC: **0.750**
- Best classical baseline (k-mer (k=8) + XGBoost): **0.694**
- Gap: **-5.6 pp** (best_classical − NT)
- Acceptance bar: gap < 15 pp → PASS

## Notes

- 3 variants run (NT-XGBoost + k-mer + gene-presence). AMRFinder deferred to Stage 1 prep — cohort's persisted `plasmid_resistance_genes` / `chromosome_resistance_genes` fields are empty for the mini cohort, and no per-strain AMRFinderPlus CLI infrastructure exists yet.
- k-mer and gene-presence both use within-fold vocabulary rebuild (training-set-only) to prevent held-out leakage.
- Top-K attribution genes (Tier 1-5 classification; gyrA/parC/parE presence check) NOT included in this smoke. Run as a separate `pipeline.py attribute` step if smoke passes.
- Locked decisions reflected: B-B (smoke = 4→3 variants; clade-only dropped), `--per-class 20` for Stage 1 cohort (N=40), deterministic `zlib.crc32` for any future MLST hashing.

## Next action

- **PASS** → proceed to Stage 1 N=50 local engineering screen (~4 hours of GTX 860M time).
- **FAIL** → NT is broken on real data. Demote NT track. Classical baselines become project spine.
- **INDETERMINATE** → fix missing variant or rerun smoke.