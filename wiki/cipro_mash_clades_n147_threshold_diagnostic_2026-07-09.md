# N=147 cipro clade partition — threshold + selection-rule diagnostic (2026-07-09)

**Companion to** `wiki/cipro_mash_clades_n147_2026-07-09.{md,json}` (the first real run of
`scripts/mash_cluster_n147.py`, all 147 strains). **Nothing here changes that artifact.** This memo records
two grounded defects in how its threshold was chosen, surfaced by the same structural diagnostics that
falsified the TB Mash collapse the same day (`wiki/tb_independent_mash_lineage_2026-07-09.md`).

**Status: FINDING + RECOMMENDATION. The selection rule is NOT changed here** — it feeds downstream
contracts (clade-balanced cohort selection; leave-one-clade-out CV folds), so the fix is surfaced for
ratification rather than applied unilaterally.

## The shipped result

The run chose **threshold 0.020 → 6 clades, largest clade 57.8%** of the cohort. It was the **only**
qualifying rung in the pinned sweep `(0.02, 0.03, 0.04, 0.05, 0.07, 0.10)`.

## Defect 1 — the chosen threshold is a grid-boundary artifact

`0.020` is the **floor** of the pinned candidate grid. Every coarser rung fails (`max_clade_fraction ≥ 0.60`),
so "the sweep chose 0.020" really means "0.020 is the smallest number we tried." Finer thresholds were never
evaluated, and they are strictly better-structured. Measured on the real 147-genome Mash matrix
(E. coli pairwise distance: median 0.0262, p25 0.0162, p75 0.0320):

| threshold | clades | largest clade | singletons | variance_ratio | qualifies? |
|---:|---:|---:|---:|---:|:---:|
| 0.002 | 101 | 17.7% | 56.5% | 0.0812 | yes |
| 0.005 | 69 | 18.4% | 27.2% | 0.0941 | yes |
| **0.008** | **45** | **19.7%** | **14.3%** | 0.1735 | yes |
| **0.010** | **36** | **19.7%** | **12.2%** | 0.2351 | yes |
| 0.015 | 11 | 57.8% | 3.4% | 0.3988 | yes |
| **0.020 (chosen)** | **6** | **57.8%** | **0.7%** | 0.3990 | yes |
| 0.030–0.100 | 2 | 99.3% | 0.7% | 0.1185 | no (max-fraction) |

A 57.8% clade means **leave-one-clade-out CV is effectively leave-58%-out** on one fold — the very fold
structure the clade partition exists to provide. At 0.008–0.010 the largest clade is 19.7%.

Note also that `0.015` (variance_ratio 0.3988) would **beat the chosen 0.020** (0.3990) on the rule's own
criterion if it were merely added to the grid. The current answer is decided at the fourth decimal place by
where the grid happens to stop.

## Defect 2 — the selection rule is monotonically biased toward over-splitting

`score_threshold` computes `variance_ratio = mean(intra-clade distance) / mean(inter-clade distance)`, and
`pick_best_threshold` takes the **minimum** among qualifying rungs. As the threshold falls, clusters shrink,
intra-clade distances → 0, and so `variance_ratio → 0` **by construction**. It is not a well-behaved
objective: it prefers finer partitions without bound.

So the grid floor is **load-bearing by accident**. Extending the grid downward without fixing the rule
picks `0.002` (lowest variance_ratio, 0.0812) → **101 clades, 56.5% of strains are singletons** — a
degenerate over-split. The two defects have been masking each other: a rule that always wants finer, and a
grid that refuses to offer finer.

The qualifying criteria (`n_clades ≥ 3`, `max_clade_fraction < 0.60`) guard only against a **dominant
clade**. Nothing guards against **over-splitting**. This is the exact structural blind spot that made the
TB granularity-matched comparison invalid: *cluster count and a scalar ratio do not describe partition
shape.*

## Defect 3 (contributing) — `cluster_by_ani` is single-linkage, and it chains

`mash_cluster_n147.py` clusters with `phylogeny.cluster_by_ani` (single-linkage union-find), whereas
`clonality.py` deliberately uses **greedy-representative** clustering precisely because single-linkage
chains `A~B~C` into one cluster. Same matrix, same thresholds:

| threshold | single-linkage clades / largest | greedy-rep clades / largest |
|---:|---|---|
| 0.008 | 45 / 19.7% | 49 / 19.7% |
| 0.010 | 36 / 19.7% | 41 / 19.7% |
| **0.015** | **11 / 57.8%** | **15 / 31.3%** |
| 0.020 | 6 / 57.8% | 9 / 55.8% |
| 0.030 | 2 / 99.3% | 3 / 70.7% |

At 0.015 the 57.8% blob is **substantially a chaining artifact** (greedy-rep: 31.3%). At 0.020 the cohort
genuinely has one large group at that radius, so both methods blob.

## Also fixed (safe, unambiguous)

The module docstring claimed the sweep picks the threshold with "intra-clade-vs-inter-clade variance ratio
**maximized**". The implementation and its tests **minimize** it (lower = tighter, better-separated clades).
Docstring corrected; no behavior change.

## Recommendation (for ratification — not applied)

1. Add an **over-split guard** to the qualifying criteria before extending the grid, e.g. a maximum
   singleton fraction (~≤20%) and/or a minimum clade size. Without it, a finer grid selects `0.002`.
2. Then extend `CANDIDATE_THRESHOLDS` downward (`0.005, 0.008, 0.010, 0.015, …`). With the guard, the
   best-structured rung on this cohort is **0.008–0.010** (36–45 clades, largest 19.7%, singletons 12–14%).
3. Prefer `clonality.greedy_representative_clusters_from_matrix` over `cluster_by_ani` for clade definition
   (chaining-resistant; already the method the frozen lineage layer trusts).
4. Report **partition structure** (largest-clade fraction, singleton fraction) alongside clade count in the
   artifact, so a degenerate partition can never be reported as a clean one.

Until ratified, the shipped `cipro_mash_clades_n147_2026-07-09` artifact stands **with this caveat
attached**: its 6-clade / 57.8%-largest partition is coarse and boundary-selected, and should not be used
for leave-one-clade-out CV folds without re-running under the guard above.

## Reproduce

```bash
uv run python scripts/mash_cluster_n147.py \
  --cohort data/processed/stage2_n150_cipro_cohort.parquet \
  --refseq-cache D:/dna_decode_cache/refseq --drug ciprofloxacin
```
Requires all 146 assembly accessions present in the refseq cache — the script **WARN-and-skips** a missing
FASTA, so an incomplete cache silently clusters a subset (43/146 before this run's fetch).
