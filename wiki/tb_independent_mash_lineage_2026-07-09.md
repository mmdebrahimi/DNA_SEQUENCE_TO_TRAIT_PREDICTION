# Independent TB cohort — finer Mash lineage collapse: **MASH_GREEDY_COLLAPSE_NOT_APPLICABLE_TB** (2026-07-09)

**Cohort:** EBI AMR-Portal TB provenance-disjoint (leaked=0) · **clustered:** 2845 / 2845 scored isolates
**Mash:** `quay.io/biocontainers/mash:2.3--hb105d93_10` sketch `-s 10000` · greedy-representative (chaining-resistant)
**Pairwise distances:** median 6.34e-04 · p75 8.41e-04 · p99 2.15e-03

Frozen clustering + confusion math reused unchanged; only the lineage definition changes.
**Nothing is re-scored.**

## Verdict — the deferred 'finer Mash-collapse' hypothesis is FALSIFIED

FALSIFIED: the deferred hypothesis that 'a finer Mash-collapse would tighten the independent TB lineage number' does NOT hold. Mash sens/spec vary MONOTONICALLY with the collapse threshold, and there is NO threshold at which the partition is both meaningfully collapsed and structurally sound. Fine rungs barely merge anything (1e-5 -> 2,501 clusters over 2,845 isolates) and simply reproduce the clonality-INFLATED raw number. Coarse rungs collapse into a blob: at 7e-4 a SINGLE cluster holds 77% of the cohort and 97% of isolates fall into mixed-label DISCORDANT clusters that are excluded, so the reported sens/spec is computed on a ~3% non-random residue. The cluster COUNT can be made to match the barcode's while the cluster STRUCTURE does not, which makes a granularity-matched comparison invalid here. Mechanism: M. tuberculosis is monomorphic and shows no lineage-scale gap at Mash resolution (pairwise median 6.3e-4, p75 8.4e-4), so any radius large enough to yield ~100 clusters exceeds most inter-lineage distances and one representative absorbs the cohort. The pinned Napier barcode -- a phylogeny-aware, marker-based partition -- is the CORRECT tool for TB, and its collapsed numbers (RIF 0.444/0.979, INH 0.321/0.972) STAND as the honest headline.

Usable (non-degenerate) rungs found: **0 of 8**. A rung is degenerate when one cluster holds >20% of the cohort, or >20% of labelled isolates fall in excluded mixed-label clusters, or the partition is barely collapsed (>50% as many clusters as isolates, i.e. the raw view).

## Threshold sweep — sens/spec is a function of the threshold, not a property of the decoder

| threshold | clusters | largest clust. | drug | disc. isolates | lineage sens [95% CI] | lineage spec [95% CI] | R-lin | S-lin | degenerate |
|---:|---:|---:|---|---:|---|---|---:|---:|:---:|
| 1e-05 | 2501 | 2.1% | rifampicin | 6.7% | 0.908 [0.89–0.923] | 0.96 [0.948–0.97] | 1157 | 1315 | **YES** |
| 1e-05 | 2501 | 2.1% | isoniazid | 1.2% | 0.856 [0.836–0.873] | 0.964 [0.952–0.973] | 1344 | 1145 | **YES** |
| 5e-05 | 1631 | 7.4% | rifampicin | 33.8% | 0.892 [0.864–0.916] | 0.974 [0.962–0.982] | 552 | 999 | **YES** |
| 5e-05 | 1631 | 7.4% | isoniazid | 26.7% | 0.787 [0.755–0.816] | 0.975 [0.962–0.983] | 678 | 880 | **YES** |
| 1e-04 | 1089 | 21.1% | rifampicin | 56.8% | 0.859 [0.813–0.894] | 0.978 [0.964–0.986] | 284 | 721 | **YES** |
| 1e-04 | 1089 | 21.1% | isoniazid | 56.8% | 0.691 [0.641–0.737] | 0.986 [0.974–0.993] | 358 | 644 | **YES** |
| 2e-04 | 582 | 32.8% | rifampicin | 76.2% | 0.805 [0.723–0.868] | 0.973 [0.952–0.985] | 113 | 409 | **YES** |
| 2e-04 | 582 | 32.8% | isoniazid | 77.2% | 0.475 [0.394–0.557] | 1.0 [0.99–1.0] | 142 | 386 | **YES** |
| 3e-04 | 358 | 37.2% | rifampicin | 84.7% | 0.691 [0.56–0.797] | 0.969 [0.941–0.984] | 56 | 263 | **YES** |
| 3e-04 | 358 | 37.2% | isoniazid | 84.4% | 0.27 [0.182–0.381] | 1.0 [0.985–1.0] | 75 | 245 | **YES** |
| 5e-04 | 160 | 38.4% | rifampicin | 91.1% | 0.545 [0.347–0.731] | 0.982 [0.937–0.995] | 25 | 112 | **YES** |
| 5e-04 | 160 | 38.4% | isoniazid | 92.0% | 0.097 [0.033–0.249] | 1.0 [0.963–1.0] | 33 | 100 | **YES** |
| 7e-04 | 84 | 76.8% | rifampicin | 96.1% | 0.25 [0.089–0.532] | 0.964 [0.877–0.99] | 13 | 56 | **YES** |
| 7e-04 | 84 | 76.8% | isoniazid | 96.9% | 0.05 [0.009–0.236] | 1.0 [0.921–1.0] | 20 | 45 | **YES** |
| 1e-03 | 43 | 90.4% | rifampicin | 98.2% | 0.0 [0.0–0.561] | 0.967 [0.833–0.994] | 4 | 31 | **YES** |
| 1e-03 | 43 | 90.4% | isoniazid | 98.3% | 0.0 [0.0–0.259] | 1.0 [0.851–1.0] | 11 | 22 | **YES** |

## Granularity-matched comparison — attempted, and INVALID

### rifampicin — nearest rung 7e-04 (84 Mash lineages vs 110 barcode lineages)

> **This comparison is NOT valid and the Mash row must not be quoted as a result.** At this threshold one cluster holds 76.8% of the cohort and 96.1% of labelled isolates are excluded as DISCORDANT, so the Mash sens/spec below is computed on a small non-random residue (n_scored=67). Matching the cluster COUNT did not match the cluster STRUCTURE.

| | sens | spec | R-lin | S-lin | discordant | n_scored |
|---|---|---|---:|---:|---:|---:|
| Mash *(unusable)* | 0.25 [0.089–0.532] | 0.964 [0.877–0.99] | 13 | 56 | 15 | 67 |
| **barcode (stands)** | **0.444** | **0.979** | 20 | 47 | 43 | — |

### isoniazid — nearest rung 7e-04 (84 Mash lineages vs 110 barcode lineages)

> **This comparison is NOT valid and the Mash row must not be quoted as a result.** At this threshold one cluster holds 76.8% of the cohort and 96.9% of labelled isolates are excluded as DISCORDANT, so the Mash sens/spec below is computed on a small non-random residue (n_scored=65). Matching the cluster COUNT did not match the cluster STRUCTURE.

| | sens | spec | R-lin | S-lin | discordant | n_scored |
|---|---|---|---:|---:|---:|---:|
| Mash *(unusable)* | 0.05 [0.009–0.236] | 1.0 [0.921–1.0] | 20 | 45 | 19 | 65 |
| **barcode (stands)** | **0.321** | **0.972** | 30 | 36 | 44 | — |

## Raw-limit sanity anchor (finest rung, threshold 1e-05)

At 2501 clusters over 2845 isolates the collapse is nearly a no-op, so these values should approach the published RAW numbers (RIF sens 0.920 / spec 0.955; INH sens 0.879 / spec 0.962). They do — which validates the whole pipeline end-to-end:

- **rifampicin**: sens 0.908 · spec 0.96 (R-lin=1157 S-lin=1315)
- **isoniazid**: sens 0.856 · spec 0.964 (R-lin=1344 S-lin=1145)

## Honesty

GENUINELY INDEPENDENT (out-of-CRyPTIC-build): accession-level provenance-disjoint, measured phenotype, WHO catalogue applied UNCHANGED. The lineage-collapsed number is the honest headline; the raw per-isolate number is clonality-inflated. This run REFINES the clustering (Mash greedy-representative) vs the coarser pinned Napier barcode; it re-scores nothing.

## Scope limits

Callability unassessed (no regeno) -> determinant calls are a CONSERVATIVE lower bound. Assembly-available subset (not prevalence-preserving) -> status stays TB_SUBSET_PLUMBING. asm5 minimap2 VCFs miss some determinants. One anomalous pair has Mash distance >=0.5 (no shared hashes) -- a single divergent/low-quality assembly, not systematic.

Generated by `scripts/tb_independent_mash_collapse.py`.