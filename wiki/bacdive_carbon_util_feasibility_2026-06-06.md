# BacDive carbon-utilization substrate feasibility census — 2026-06-06

> EP-6 entry gate. Which E. coli carbon source has a de-confoundable
> utilizer/non-utilizer cohort big enough to test NT embeddings vs gene-content
> on a sampling-INDEPENDENT label with NO curated catalog?

- Organism: `Escherichia coli`
- Carbon sources in export: 1
- Min-strains floor: 100 · min-minority-frac: 0.15
- MLST sidecar provided (de-confound gate runnable): False
- **Feasible carbon sources: 0**

## Feasible (ranked: balance, then downloadable cohort size)

| carbon source | N | +/- | minority frac | with-accession | de-confound | verdict |
|---|---:|---|---:|---:|---|---|
| (none) | | | | | | |

## All assessments

| carbon source | N | +/- | with-acc | verdict |
|---|---:|---|---:|---|
| tryptophan | 27 | 25/2 | 0 | INFEASIBLE_TOO_FEW |

## Interpretation

- `FEASIBLE` = clears count + balance + accession floors AND de-confound gate
  returned DE_CONFOUNDED (within-MLST utilizer/non-utilizer contrast exists).
- `FEASIBLE_PENDING_DECONFOUND` = count floors clear; MLST not supplied so the
  de-confound gate could not run. Supply `--mlst` to resolve.
- `BLOCKED_CONFOUNDED` = no within-lineage contrast → predicting utilization would
  predict lineage, not metabolism (the Li et al. 2023 phylogeny-dominance trap).
- Ranking prioritizes HIGHER minority_fraction (balanced/harder sources) because
  easy carbon sources are already near-ceiling for gene-content RF out-of-clade —
  the embedding niche is the hard ones.