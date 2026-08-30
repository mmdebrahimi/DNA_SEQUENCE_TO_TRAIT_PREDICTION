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

## Does the same confound bite anywhere else? Audited: no

The failure mode generalises — *a pooled ranking cannot answer a within-group question* — so every other
AUROC in the package was checked for a repeated-measure grouping where within-group is the real question.

**No second instance.** The others are either single-measure per unit (`essentiality`: one call per gene,
so there is no within-gene structure to confound) or already group-aware by construction (clade baseline,
per-clade metrics, HIV within-subtype transfer, the lineage-collapse layer). Gene × condition is the only
repeated-measure design in the repo, which is exactly why this cell was the gap.

An honest negative, recorded so the check does not get repeated.

---

# Replication on the 25-source carbon axis — the direction result strengthens, the LEVER nearly vanishes

The result above rests on 26 varying genes across **4** conditions, where a within-gene AUROC can take
only a handful of discrete values. The repo already has an independent, far better-powered axis: the
Fitness Browser **Keio carbon panel — 25 mappable sources, 217 conditionally-essential genes**, with a
different substrate and a different label source (transposon fitness, not Orth's curated E/N calls).

Same metric, same code path (`--axis carbon`).

| | media4 (Orth) | carbon (Keio) |
|---|---:|---:|
| conditions | 4 | **25** |
| conditionally-essential genes | 67 | **217** |
| flat (one ratio for every condition) | 61.2% | **68.2%** |
| **within-gene AUROC, non-flat** | 0.7308 (n=26) | **0.8133 (n=69)** |
| permutation p | 0.001 | **0.0005** |
| within-gene AUROC, all genes | 0.5896 | 0.5996 |
| deployed exact-set (that axis's own) | 3/67 = 4.5% | 23/217 = 10.6% |
| oracle relative-rule ceiling | 11/67 = 16.4% | 27/217 = 12.4% |
| **headroom** | **+8 genes (+11.9 pp)** | **+4 genes (+1.8 pp)** |

Deterministic: 0.8133 on both repeats, spread 0.0.

## Two findings, pulling opposite ways

**1. The direction result replicates and strengthens.** Where the model's ratio varies across conditions,
it points the right way **81%** of the time on 69 genes — better than the 4-media 73% on 26, on different
data with a different label source. And the flatness finding replicates and *worsens*: **68% of these
genes emit one identical growth ratio across all 25 carbon sources.** "Silence, not error" is now a
two-substrate result.

**2. The practical lever nearly vanishes on the better-measured axis.** This corrects the impression the
4-media ceiling gave. There, ranking looked like a 3→11 win (~3.7×). On carbon the deployed absolute
threshold already reaches 23/217 and the oracle ceiling is 27/217 — **+4 genes, +1.8 pp, with an oracle
that is handed each gene's true essential-condition count.**

The 4-media "3.7×" was a **small-axis artifact**: with only four conditions the absolute threshold does
badly, so a relative rule looks like a large relative win over a tiny base. Given 25 conditions the
threshold does much better and the headroom collapses.

**So: build the relative rule? No.** A deployable version must also infer *k*, would recover at most 4 of
217 genes, and inferring *k* is the original problem restated. The honest conclusion is that the readout
is **not** where the deficit lives — measured on the axis best able to say so.

## What this does and does not change

- **Unchanged:** the switch cell stays open, and its bottleneck stays the one already measured — the
  conditioning signal is not measured in the conditions the phenotype data uses.
- **Sharpened:** "the readout costs real signal" is now bounded at **+1.8 pp** on the better axis, against
  +11.9 pp on the smaller one. Quote the carbon number.
- **Strengthened:** flatness is the dominant term on both axes and grows with condition count
  (61% → 68%). A model that emits one number for 25 different carbon sources is the thing to fix.

## One defect, same family as the other four

The first carbon run printed its oracle ceiling against **`deployed 3/67`** — the *4-media* baseline —
which made the lever look roughly four times more valuable than it is. `deployed_exact_set(axis)` now
reads each axis's own committed artifact. It also had to tolerate schema drift between two generations of
the same producer (`n_scored_exact_set` on carbon, only `n_conditionally_essential` on the older 4-media
artifact); reading one key alone reported a silent `unknown`. 3 tests added (12 total).

Reproduce: `uv run python scripts/fba_within_gene_ranking.py --axis carbon --repeat 2`
(needs `D:/dna_decode_cache/fitness_browser/feba.db`).
