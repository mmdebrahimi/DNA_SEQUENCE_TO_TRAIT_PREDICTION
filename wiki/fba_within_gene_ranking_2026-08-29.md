# Where the model varies at all, it varies in the right direction — the switch failure is silence, not error

The last untried cheap lever on the open conditional-switch cell was "score the continuous knockout growth
ratio as a ranking instead of thresholding at 1%." It is now run, and the answer is a **bounded PASS**
that sharpens the diagnosis rather than reopening the cell.

## Why the existing number could not answer this

`continuous_readout` scores the ratio as a ranking **pooled over every gene × condition cell** and reports
AUROC ~0.59. That reads as "the cutoff is discarding signal" — but a pooled ranking is dominated by the
**gene main effect**. A gene essential in all four media has a low ratio in all four and contributes four
correctly-ranked positives without the model ever having switched. **A pooled 0.59 is reachable with
exactly zero within-gene signal.**

The switch question is strictly *within* gene: for one gene, is the ratio lower in the conditions where
that gene is actually essential? Conditioning on the gene removes its main effect by construction — the
same de-confounding idiom already used here for lineage, clonality and ancestry, applied to the one axis
it had never been applied to.

## Pre-registered, then run

| | |
|---|---|
| primary | mean within-gene AUROC over **non-flat** conditionally-essential genes |
| PASS | > 0.60 and permutation p < 0.05 |
| FAIL | ≈ 0.5 — variation is uninformative even where it exists; lever closes |
| must-hold | flat fraction reproduces the 2026-08-12 artifact's ~64% |

## Result

| metric | value | n |
|---|---:|---:|
| pooled (existing, gene effect **not** removed) | ~0.59 | 268 cells |
| **within-gene, non-flat — pre-registered primary** | **0.7308** | 26 genes |
| within-gene, all genes | 0.5896 | 67 genes |
| permutation p (within-gene label shuffle) | **0.001** | 2000 perms |
| flat fraction (must-hold) | 0.6119 ✓ | 41/67 |

**Deterministic:** 0.7308 on every repeat, spread 0.0 against a margin of 0.131 over the bar.

Per-gene distribution over the 26 non-flat genes: **11 rank perfectly (1.0)**, 2 at 0.833, 9 at 0.667,
1 at 0.333, **3 fully inverted (0.0)**. The signal is not carried by one or two genes.

## What it means

**The model is not wrong about direction — it is silent.** Where its growth ratio varies across media at
all, it points the right way 73% of the time. But **61% of these genes emit one identical number for all
four media**, contributing exactly 0.5 each and dragging the all-genes mean to 0.5896 — which is why the
pooled readout looked like weak-but-real signal.

This independently corroborates the existing `MIS_CONDITIONED = 0` finding on a completely different
metric: the model **fires roughly right or stays silent, never in the wrong place**. Two measurements,
two axes, same shape.

## The ceiling, and why this does not reopen the cell

An **oracle** relative rule — handed each gene's true number of essential conditions *k*, and calling the
*k* lowest-ratio conditions essential — gets **11 of 67** genes' exact pattern right, against the deployed
threshold's 3. Better, and still small.

Three honest bounds on that 11:

1. **It is handed *k*.** A deployed rule must infer how many conditions a gene is essential in, and
   nothing here estimates that. So this **ranks, it does not call** — the same shape as the `inverse`
   cell, which ranks edits and refuses to dose them.
2. **11 is exactly the count of AUROC-1.0 genes**, as it must be: top-*k* selection is right precisely
   when every essential condition ranks below every dispensable one. The ceiling is not independent
   evidence; it is the same fact in the project's own metric.
3. **The flat 61% are unreachable by any readout change whatsoever.** The model emits one number; no
   thresholding, ranking, or calibration recovers a distinction that was never computed.

So the 2026-08-12 verdict — *"the readout costs real signal; most of the deficit is still the model"* —
**stands, and is now quantified**: the readout's share is bounded at 11/67 with an oracle, and 61% of the
deficit is structural flatness.

## Four defects, all in my own metric, all flattering

Every one made the result look better, and all four are the same root: **this data is full of exact ties,
and at every tie an arbitrary choice can masquerade as a result.**

| # | defect | effect |
|---|---|---|
| 1 | flatness used a 1e-9 tolerance; the tie test used exact float equality | 36 of 41 flat genes scored LP noise as signal |
| 2 | flat genes counted as oracle hits via stable-sort index order | ceiling inflated 11 → **23** |
| 3 | a tie straddling the top-*k* boundary counted as a hit | ceiling inflated 11 → 13 |
| 4 | determinism tested by strict equality | apparent ±0.013 run-to-run variance — which was **defect 1**, not LP degeneracy |

Defect 1 was caught by arithmetic, not by inspection: 41 genes at exactly 0.5 plus 26 at 0.718 cannot
average 0.6045. Fixing it made the whole run bit-deterministic across repeats.

This is the same trap recorded for BLOSUM62 mid-ranks in the resistance-conservativeness probe, where
`sorted()`-order tie-breaking silently shifted a median and moved a p-value from 0.682 to 0.614. **Second
independent instance in this repo. Use tolerance-aware ties, and refuse rather than break them.**

All four are pinned by `tests/test_fba_within_gene_ranking.py` (9 tests, offline).

Reproduce: `uv run python scripts/fba_within_gene_ranking.py --repeat 2` (cobrapy, ~2 min/repeat,
single-process by design).
