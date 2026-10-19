**The reuse seam, read in source (`dna_decode/eval/genomic_prediction.py`).**

- `cv_ridge_gp(X, y, *, trait, n_splits=5, alphas=(1.0,10.0,100.0,1000.0), seed=0, n_perm=100) -> GPResult`
  and `cv_model_gp(X, y, *, model="gbm", ..., n_perm=0)` are the two public entry points. Both build folds
  via `_kfold_indices(n, k, seed)` and then consume a `folds` list-of-index-arrays.
- **The permutation null already inherits the fold structure.** `cv_model_gp`'s null is
  `_cv_r_with(X, rng.permutation(y), folds, fp)` — the *same* `folds` object as the measured run. So a
  grouped `folds` makes the null grouped **by construction**, with no second change. This is the single
  most important fact for the design: it means the null never silently bounds a different quantity than
  the metric it is compared against.
- `GPResult` already carries `null_r_mean`, `null_r_p95`, `beats_null` — the per-substrate, instrument-native
  threshold this plan uses in place of any asserted number.
- **A latent leak in the ridge path.** `_cv_predictive_r` picks `alpha` on an inner split taken as
  `train_idx[:n_in]` / `train_idx[n_in:]` — a positional slice with no group awareness. Under grouped
  outer folds, inner train and inner validation can hold members of the same group, so the alpha choice
  leaks within the fold. Small (alpha selection only, on a 4-value grid) but it is exactly the class of
  defect this plan exists to remove, so Step 2 fixes it rather than documenting it.

**The grouping primitive, read in source (`dna_decode/eval/clonality.py`).**

- `greedy_representative_clusters_from_matrix(distance_matrix, threshold) -> dict[str, int]` — pick a
  representative, drop everything within `threshold`, repeat. **Chaining-resistant**, deterministic via
  `sorted(distance_matrix.strain_ids)`. This is reused, not reimplemented: CLAUDE.md records that
  `phylogeny.cluster_by_ani`'s single-linkage **chains** (57.8% largest cluster vs greedy-rep's 31.3% on
  the same data).

**The worked precedent (`scripts/crossaxis_lineage_deconfound.py`).** Already does
Mash-cluster → `GroupKFold` vs naive `StratifiedKFold` with a collapse-toward-chance reading, for AMR
cohorts. It establishes the *technique*; it does not generalise (bacterial Mash distances, binary R/S
labels, no permutation null). Its reading convention — compare against what chance looks like under the
*same* fold structure — carries over.

**Consumers that unlock (no edits in this plan).** `scripts/yeast_bloom_gp_arm.py`
(`D:/dna_decode_cache/bloom/geno_v2.txt` + `BYxRM_PhenoData.txt`), `scripts/arabmagic_gp_arm.py`
(`D:/dna_decode_cache/arabmagic/arabmagic_{geno,pheno}.csv`), `scripts/bxd_gp_arm.py`,
`scripts/dgrp_learned_decoder.py`.

**Substrate reality, measured not assumed (this is a correction to the stated motivation for F3).**
`data/arabidopsis/` contains **phenotypes only** — `FT10_pheno_261.csv`, `FT16_pheno_262.csv`,
`flowering_1001g/pheno_{261,262,29_FLC}.csv`, each carrying accession lat/long/country. There is **no
cached 1001-Genomes genotype matrix.** The only cached Arabidopsis genotypes are `arabmagic`, a
**constructed** multi-parent panel. So the natural-population positive control I had assumed was on disk
is not, and the plan must *select* its control pair from what is actually intact rather than name it in
advance.

### Reusable-Code Survey

Searched: `dna_decode/eval/`, `dna_decode/glm/`, `scripts/*_arm.py`, `scripts/*deconfound*`,
`scripts/*lineage*`, `graphify-out/GRAPH_REPORT.md` (absent), `src/`|`lib/`|`utils/` (absent — this repo
is a flat package under `dna_decode/`).

| Candidate | Verdict |
|---|---|
| `eval/genomic_prediction.py` → `_kfold_indices`, `_cv_r_with`, `_cv_predictive_r`, `GPResult` | **EXTEND.** Thread `groups`; the `folds` seam carries the null for free. |
| `eval/clonality.py` → `greedy_representative_clusters_from_matrix`, `DistanceMatrix` | **REUSE AS-IS.** Chaining-resistant; the alternative is recorded as defective. |
| `eval/clonality.py` → `cluster_weighted_confusion` | **REJECT.** Collapses a *classification* confusion matrix; this arm is continuous. |
| `eval/phylogeny.py` → `cluster_by_ani` | **REJECT.** Single-linkage, chains (recorded, measured). |
| `scripts/crossaxis_lineage_deconfound.py` | **PATTERN ONLY.** Technique reused; code is cohort/binary-specific. |
| `scripts/source_diverse_validate.py` | **PATTERN ONLY** — the refuse-and-emit-no-metrics shape for a substrate that fails its own bar. |
| `sklearn.model_selection.GroupKFold` | **REJECT as the fold builder.** It balances group *sizes*, not fold *counts*, and silently reduces effective splits when one group dominates; a local builder that can **refuse** is what the over-split and dominant-group guards need. Used only as a cross-check in tests. |
