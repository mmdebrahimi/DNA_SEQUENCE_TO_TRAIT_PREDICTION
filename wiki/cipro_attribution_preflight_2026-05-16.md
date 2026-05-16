# Cipro attribution preflight v2 — N=38 cohort, mean-pool NT-XGBoost (2026-05-16)

**Purpose:** test whether NT embeddings carry cipro-resistance signal at all (independently of Stage 1's verdict gate).
**Method:** train final NT-XGBoost on all 38 strains (no LOSO); per cipro-R strain, run gene_level_mutagenesis; dual aggregations (sum positive delta + frequency-in-top-K) across 19 cipro-R strains.
**Pooling:** mean-pool (512-dim) — matches Stage 1 (NOT Stage 1b mean+max).
**Signed-delta filter:** only positive deltas (knockout LOWERED R-prob = gene supporting R) count toward sum aggregation.
**Loci tracked:** expanded set across 5 mechanism classes (QRDR target-alteration + plasmid protect/modify + efflux + porin + regulatory). v1 narrow {gyrA, parC, parE} preserved as TEXTBOOK_QRDR_LOCI subset.

**Verdict:** INCONCLUSIVE_MISS
- Textbook QRDR found: []
- Expanded cipro-loci found: []
- Mechanism classes hit: []

## Cohort-wide top-20 by SUM POSITIVE delta

| rank | gene_symbol | sum(+delta) | freq | mechanism |
|---:|---|---:|---:|---|
| 1 | wcaA | 0.2707 | 2 |  |
| 2 | prpD | 0.1587 | 1 |  |
| 3 | htpG | 0.1520 | 1 |  |
| 4 | miaB | 0.1505 | 1 |  |
| 5 | rmuC | 0.1407 | 2 |  |
| 6 | surA | 0.1405 | 2 |  |
| 7 | holD | 0.1404 | 1 |  |
| 8 | yibH | 0.1404 | 1 |  |
| 9 | selA | 0.1404 | 1 |  |
| 10 | mltF | 0.1398 | 1 |  |
| 11 | plsB | 0.1398 | 1 |  |
| 12 | prfC | 0.1398 | 1 |  |
| 13 | ykfG_1 | 0.1369 | 1 |  |
| 14 | cydA | 0.1307 | 1 |  |
| 15 | potA_1 | 0.1299 | 1 |  |
| 16 | hisC | 0.1283 | 1 |  |
| 17 | lacI | 0.1283 | 1 |  |
| 18 | hemA | 0.1283 | 1 |  |
| 19 | gsiA_4 | 0.1283 | 1 |  |
| 20 | fdhE | 0.1283 | 1 |  |

## Cohort-wide top-20 by FREQUENCY in per-strain top-K

| rank | gene_symbol | freq | sum(+delta) | mechanism |
|---:|---|---:|---:|---|
| 1 | wcaA | 2/19 | 0.2707 |  |
| 2 | rmuC | 2/19 | 0.1407 |  |
| 3 | surA | 2/19 | 0.1405 |  |
| 4 | yncL | 2/19 | 0.0972 |  |
| 5 | yqiK | 2/19 | 0.0397 |  |
| 6 | hisB | 2/19 | 0.0276 |  |
| 7 | prpD | 1/19 | 0.1587 |  |
| 8 | htpG | 1/19 | 0.1520 |  |
| 9 | miaB | 1/19 | 0.1505 |  |
| 10 | holD | 1/19 | 0.1404 |  |
| 11 | yibH | 1/19 | 0.1404 |  |
| 12 | selA | 1/19 | 0.1404 |  |
| 13 | mltF | 1/19 | 0.1398 |  |
| 14 | plsB | 1/19 | 0.1398 |  |
| 15 | prfC | 1/19 | 0.1398 |  |
| 16 | ykfG_1 | 1/19 | 0.1369 |  |
| 17 | cydA | 1/19 | 0.1307 |  |
| 18 | potA_1 | 1/19 | 0.1299 |  |
| 19 | hisC | 1/19 | 0.1283 |  |
| 20 | lacI | 1/19 | 0.1283 |  |

## Mechanism-class breakdown

| mechanism | loci tracked | found in any top-K |
|---|---|---|
| QRDR_target_alteration | 4 | (none) |
| plasmid_protect_modify | 8 | (none) |
| efflux | 7 | (none) |
| porin | 2 | (none) |
| regulatory | 5 | (none) |

## Verdict interpretation rubric

- **STRONG_POSITIVE:** ≥1 QRDR locus (gyrA/B, parC/E) AND ≥1 mechanism class hit; positive signed delta in multiple R strains. Architecture has real biology.
- **WEAK_POSITIVE:** any expanded-set locus appears in top-K, but no QRDR or only single-strain. Suggestive but not load-bearing.
- **INCONCLUSIVE_MISS:** no cipro-loci recovered. DO NOT treat as architecture mismatch — mean-pool dilution (1/N_genes per knockout) + symbol-coverage limits (RefSeq ~11% CDSs carry `gene=`) are confounders. Next escalation: refactor `gene_level_mutagenesis` to mean+max + retest matching Stage 1b's pooling.
- **DAMNING_MISS:** v2 + mean+max preflight both return no cipro-loci with positive signed delta. THEN architecture mismatch is defensible.

## Notes

- N=38 cohort; 19R/19S balance per Stage 1 verdict.
- Mean-pool 512-dim features (matches Stage 1, NOT Stage 1b's 1024-dim mean+max). `gene_level_mutagenesis` hardcodes mean re-aggregation after knockout; mean+max attribution requires refactor.
- Empty `gene_symbol` entries are dropped from cohort aggregation (only ~11% of CDSs in RefSeq GFF3 carry `gene=` per CLAUDE.md gotcha). Genes-without-symbol may still rank high per-strain via gene_id but won't aggregate cross-strain. Per-strain JSON sidecar preserves raw gene_id ranking for offline inspection.
- Classifier exit-code-agnostic: the preflight uses the model regardless of Stage 1 verdict. If NT-XGBoost AUROC was FAIL at LOSO but attribution recovers QRDR, the architecture has signal that small-N LOSO can't surface.
- Per-strain top-K table persisted at sidecar JSON: `cipro_attribution_preflight_2026-05-16.per_strain.json`.