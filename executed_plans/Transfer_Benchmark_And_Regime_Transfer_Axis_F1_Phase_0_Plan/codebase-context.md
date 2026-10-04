### The package v1 missed, which this plan now builds on

`dna_decode/deconfound/` is canonical, exported and tested (`tests/test_deconfound_package.py:17::test_public_api_imports` pins the public surface). Signatures read from source:

| Function | Signature / behaviour | Role here |
|---|---|---|
| `cv_r2(X, y, groups=None, alpha=10.0, seed=0) -> float` | `deconfound.py:24`. With `groups`, leave-one-group-out, pooled out-of-fold r². **SILENTLY SKIPS a group when `tr.sum() < 5 or te.sum() < 1` (`:30-31`, `continue`)**, and returns `nan` when `ok.sum() <= 10`. | the zero-shot held-out-group score |
| `within_group_r2(X, y, groups, min_n=30, alpha=10.0) -> (float, int)` | `:43`. Pooled out-of-fold r² computed INSIDE each group of ≥`min_n`, scored on group-centered residuals with test truth centered by the **TRAIN** mean. Returns `(r2, n_groups_used)`. | **the honest de-confounded metric** |
| `group_centered_association(y, x, groups) -> (global_rho, within_group_rho, n_used)` | `:90`. NaN-aware; returns the confounded AND de-confounded association in one call. | the pooled-vs-within contrast |
| `permutation_null(y, x_resid, groups, n=200) -> np.ndarray` | `:138`. Shuffles y **INSIDE each group**, re-centers, correlates against fixed `x_resid`. Deterministic (per-draw seed = loop index); **no seed parameter**. | the correct null for the within-group estimand |
| `cluster_from_distance(dist, k) -> np.ndarray` | `:70`. **AVERAGE**-linkage (`:75`, `linkage(..., method="average")`), K clusters via `fcluster(..., criterion="maxclust")`. | clade folds |
| `r2(y_true, y_pred)` | `:18`. Plain r², can go negative. | shared scorer |

Also `scorecard.py` (`Candidate`, `GATE_KEYS`, `score`) — dataset-candidate gating, **not used by this plan**.

### Import posture — a constraint, measured

`dna_decode/eval/__init__.py` is **docstring-only**, no eager submodule imports. `eval/regime.py` imports in **0.005 s** and `eval/genomic_prediction.py` in **0.002 s** because the latter imports sklearn **lazily inside functions** (`:49`, `:97`, `:109`). `dna_decode/deconfound/deconfound.py` imports `scipy.cluster.hierarchy`, `scipy.spatial.distance`, `scipy.stats` and `sklearn.linear_model` at **module scope**.

Therefore every `dna_decode.deconfound` import in this plan is **lazy, inside the function that uses it**, matching `genomic_prediction.py`'s existing convention. A module-scope import would give `eval/transfer.py` an eager scipy+sklearn load and slow every test collection that touches `eval/`.

### Why the gauntlet lives in `eval/`, not in `deconfound/` — measured, not preference

`dna_decode/deconfound/__init__.py` imports its **whole export surface eagerly** (`from dna_decode.deconfound.deconfound import cluster_from_distance, cv_r2, …`). So placing the gauntlet inside that package would put transfer-benchmark code in the import path of all **8** existing consumers — `dataset_candidate_scorecard`, `depmap_decoder`, `depmap_fusion_methylation`, `depmap_multimodal`, `dgrp_learned_decoder`, `gdsc_fusion_decoder`, `yeast_cnv_attribution`, `yeast_growth_decoder` — none of which has anything to do with transfer. The de-confounding **mathematics** still has exactly one home; what lives in `eval/` is orchestration with zero statistics of its own, and the Step 5 AST guard is what enforces that rather than a promise.

### The binding group-size floor is 30, not 10 — measured

`within_group_r2(X, y, groups, min_n=30, …)` skips any group with fewer than `min_n` members (`deconfound.py:52`, `continue`) and returns `(nan, used)`; it also returns `nan` when fewer than 11 points were scored. Empirically, on 6 groups × 20 members:

| call | result |
|---|---|
| `within_group_r2(X, y, g)` | **`nan`, `n_groups_used=0`** |
| `cv_r2(X, y, groups=g)` | **`−0.4230`** — a confident-looking number on the same data |
| `within_group_r2` at 40/group | `0.5939`, `n_groups_used=6` |

Since PASS is read off the **within-group** cell, the usable k ceiling is governed by `max(MIN_SCORED, 30)`, and a `nan` / `n_groups_used == 0` must be surfaced as `unscorable` rather than silently falling back to the pooled cell.

### Substrate reachability — verified on disk

| Arm | Substrate | State |
|---|---|---|
| A1 scoring-core control | Bloom BYxRM | `D:/dna_decode_cache/bloom/BYxRM_GenoData.txt` (1,048,576 B) + `BYxRM_PhenoData.txt` (797,659 B). One cross. |
| B synthetic negative | generated in-process | none needed; deterministic |
| C verdict-replay regression | Arabidopsis FT10 G2 result | `wiki/phase2_arabidopsis_result_2026-06-12.md` — **markdown table only, no JSON sidecar exists.** The per-seed numbers are not machine-readable today. |
| (not an arm) | BXD mouse, Arabidopsis MAGIC | `D:/dna_decode_cache/bxd/`, `D:/dna_decode_cache/arabmagic/` (677 lines × 1260 markers × 8 traits). Single panels; disjoint marker spaces. |

`wiki/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.json` carries the three-metric *shape* (`forced_readout.r2_only_improvement_present: true`, `within_group_clean_pass_present: false`) but **not** the per-seed values.

### Reusable-Code Survey

**Adopted:** `dna_decode/deconfound/` — `cv_r2`, `within_group_r2`, `group_centered_association`, `permutation_null`, `cluster_from_distance`, `r2` (all six; no equivalent is re-written) · `eval/genomic_prediction.py` lazy-import convention · `tests/test_deconfound_package.py:23::test_within_group_r2_deconfounds` as the **fixture template for the synthetic negative arm** (group offsets plus a within-group signal on one feature — the negative arm is that fixture with the signal term set to zero) · the repo's `ContractNumberBinding` verified-not-trusted pattern for Step 2.

**Found and deliberately NOT called:** `eval/clade_baseline.py` + `eval/cv.py::leave_one_clade_out_cv` + `fba/nulls.py::curveball_shuffle` — classification/binary-matrix shaped (AUROC dicts, `train_fn`/`predict_fn`, 0/1 matrices) while both live substrates are continuous · `eval/clonality.py::greedy_representative_clusters_from_matrix` — v1 preferred it over `cluster_from_distance` on the grounds that "single-linkage chains"; that is true of `eval/phylogeny.py::cluster_by_ani:232` (union-find) but **`cluster_from_distance` is average-linkage**, so the argument never applied to the incumbent being replaced. Retracted; `cluster_from_distance` is used.

**Searched:** `dna_decode/eval/` (23 modules) · `dna_decode/deconfound/` (3) · `dna_decode/fba/nulls.py` · `scripts/{yeast_bloom_gp_arm,bxd_gp_arm,arabmagic_gp_arm,yeast_growth_decoder,dgrp_learned_decoder}.py` · `pyproject.toml` · `graphify-out/GRAPH_REPORT.md` (**absent**) · `src/utils`, `src/lib`, `src/common`, `utils`, `lib`, `common` (**none exist**).

### Prior-work replay (HIGH-salience, overlapping modules)

- *2026-08-31* — "DERIVE the blocker from your own committed artifacts." This plan exists because v1 asserted a substrate it never verified.
- *2026-08-26* — "pin a non-vacuity count so it cannot go green on an empty set"; "scope a guard to the cases the rule actually applies to." Both drive Step 5's guard design.
- *2026-06-10* — shared-key silent-overwrite trap → `wiki/transfer_benchmark_*` namespace, verified free.
- *2026-05-24* — pre-commit the verdict-conditional response BEFORE results land → Step 6's frozen bar above the data load.
- *2026-09-11* — "parameterise the blocked value; build the seam" → the replay fixture (Step 2) is data, not a hard-coded constant.
