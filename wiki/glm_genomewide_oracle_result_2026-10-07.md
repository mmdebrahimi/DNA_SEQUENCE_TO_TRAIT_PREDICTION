# The oracle's 0.59 is GRID-SPECIFIC — and the feature ranking INVERTS on real genomic sequence

**Date:** 2026-10-07 · **Artifact:** `wiki/glm_genomewide_oracle_2026-10-07.json`
**Reproduce:** `uv run python scripts/glm_genomewide_oracle.py --limit 30000` (CPU, ~4 min)
**Substrate:** GEO **GSE144621** — 321,123 sheared *E. coli* MG1655 genomic fragments in **LB and M9**, plus
46,713 peak tiles. Downloaded and measured, not read from the paper.

**This retracts the generality of today's earlier headline.** The oracle reached ρ 0.59 on the designed σ70
grid and captured 79% of achievable selection value there. Its own honest-limit section said a fully-crossed
designed grid is not a natural promoter distribution. **That limit has now been tested and it binds.**

## The result: the ranking does not merely weaken, it reverses

| feature | designed grid (leave-element-out) | **genomic peak tiles** (leave-peak-out) |
|---|---|---|
| GC content (1 feature) | 0.1414 | **+0.3000 — best** |
| 4-mer (256) | 0.4137 | +0.2769 |
| 6-mer (4,096) | 0.5872 | +0.2022 |
| positional one-hot (600) | **0.5914 — best** | +0.2691 |

On the grid, GC was worst by a factor of four and positional one-hot won. On real genomic promoter tiles
**GC alone beats every learned feature**, and the one-hot feature that won on the grid drops from 0.5914 to
0.2691. Verdict `DOES_NOT_GENERALISE` on the pre-registered bar (best learned must beat GC by ≥0.05).

**And the data is not the excuse.** The tiles' own replicates agree at 0.7945, so the ceiling on any
sequence model here is **0.9410** — a *cleaner* assay than the frag library. The best learned feature reaches
**29.4% of that ceiling.** A clean assay plus a weak model is a model finding, not a data finding.

**Mechanism, and it is the obvious one:** a fully-crossed grid varies a handful of defined elements at
**fixed offsets**, so a positional model can learn "this −10 variant at this offset → this strength". Real
genomic promoters sit at arbitrary offsets inside arbitrary flanking sequence, so position-specific weights
have nothing stable to attach to. What survives is bulk composition — and the GC signal is **negative
(ρ = −0.3253), i.e. AT-rich tiles are more active**, which is correct *E. coli* promoter biology (the −10
box TATAAT is AT-rich) rather than an artifact.

## The random-fragment arm is a different question, and that matters

| feature | random genomic fragments (position-blocked) |
|---|---|
| GC | +0.1023 |
| 4-mer | +0.1313 |
| 6-mer | +0.0919 |

Near-zero, and **not from lack of data** — a size ladder at n = 3,000 / 10,000 / 30,000 gives 0.099 / 0.154 /
0.131 for 4-mers, so it plateaus rather than climbing, and the negative is not power-limited.

The reason is that it asks a different question. **96.7% of random fragments sit below 2× the median**,
packed at a ~0.87 baseline with a thin tail to 400. A sheared-genome library asks *"is there any promoter in
this arbitrary piece of genome?"* — rare-event **detection**. The grid asked *"how strong is this promoter?"*
Quoting one against the other would be comparing two tasks. That is why the peak-tile arm exists and is the
like-for-like comparator.

## Growth condition carries real signal — the bacterial analogue of cell-state conditioning

The same fragments measured in two media, so this is directly measurable. **The comparison must be
noise-matched, and getting it wrong inverts the answer.**

| comparison | ρ |
|---|---|
| within LB, replicate 1 vs 2 | 0.4578 |
| within M9, replicate 1 vs 2 | 0.4774 |
| **cross LB–M9, single-vs-single (mean of 4 pairings)** | **0.4062** |
| cross LB–M9, average-vs-average (the *unmatched* form) | 0.5025 |

Naively, LB and M9 agree (0.5025) *better* than either medium agrees with itself (0.458/0.477) — which reads
as "condition does not matter". That is backwards: `RNA_exp_ave` is a 2-replicate mean, so comparing it
across media is average-vs-average while replicate agreement is single-vs-single, and averaging cuts noise by
~√2. Matched, **cross (0.4062) falls below within (0.4676), gap −0.0614 on 297,599 shared fragments** →
`CONDITION_CARRIES_SIGNAL`. Disattenuated, the **true cross-condition correlation is 0.7887**, so roughly
38% of the real variance is condition-specific.

**So condition-conditioning is worth building, measured rather than assumed.** This is the one item from the
AlphaGenome-critique triage that was both new and actionable, and it now has a number behind it.

## Six traps, and three of them were my own errors caught in flight

1. **The delimiter is MIXED within one GEO series.** The two frag files are space-delimited; the peak-tile
   file is tab-delimited. Either assumption applied globally parses the other to **zero rows**, and a
   zero-row parse that is then averaged reads as "no signal" rather than "no data". The loader sniffs per
   file and refuses on an implausible column count. **I had told a collaborating session these files were
   space-delimited, as fact, from the sibling dataset's behaviour. That was wrong for one of the three.**
2. **`RNA_exp_ave` is already DNA-normalised** — measured ρ with `DNA_ave` = **−0.039**. Dividing by
   `DNA_ave` again, the obvious MPRA move, would double-normalise and inject the DNA column's noise.
3. **91.5% of consecutive fragments overlap**, measured. A random split puts near-duplicates on both sides,
   so only position-blocked and leave-peak-out splits are offered — there is deliberately no `random_split`
   in the module to reach for.
4. **I assumed peak tiles were a fixed 150 bp. They are 48–150.** The module's own `positional_onehot` guard
   refused the ragged input instead of padding it, which is the only reason a silently-anchored, plausible,
   meaningless number was not produced. 97.4% are exactly 150, so filtering to those costs 2.6% and makes
   the positional arm genuinely valid.
5. **I normalised the tile arm by the FRAG file's ceiling**, and it flattered the model: 34.9% of 0.7928
   where the honest figure is 29.4% of 0.9410. The two subsets differ hugely in reproducibility (0.458 vs
   0.795 single-replicate), so each arm must carry its own ceiling.
6. **1,512 designed controls** (pos/neg/random) are chosen to be extreme and were dropped — leaving them in
   credits the model with spread it never had to earn.

## Honest limits

- **The grid-vs-genomic comparison crosses two different studies**, so part of the gap could be assay or lab
  difference rather than grid-versus-genomic. The *within-substrate* finding — GC beating every learned
  feature on tiles, under one split, one assay — does not depend on that comparison and is the stronger half.
- **GSE108535 has no derivable ceiling** (its two published tables correlate at exactly 1.0 — one monotone
  transform of the other, not replicates), so its 0.59 cannot be ceiling-normalised and the two
  fraction-of-ceiling figures are **not** comparable. Only the raw ρ values are.
- Frag model arms run on a 30,000 subsample (6-mers are 4,096 features); the ceiling and the condition
  effect use the full 321k/318k.
- Ridge on fixed features only. **No learned representation was trained**, which is exactly the gap the
  next step addresses.
- A fragment straddling a held-block boundary is assigned by midpoint, so a little sequence crosses the
  seam: 2 seams over 4.64 Mb.
- Leave-peak-out holds out genomic locations, which is positional, but it is not a phylogenetic split.

## What this changes

1. **Today's 79%-selection-value result is grid-specific and must be quoted that way.** The scoring head is
   validated on a designed grid and *not* on real genomic sequence.
2. **It raises the value of a learned representation rather than lowering it.** Fixed features demonstrably
   fail to transfer; whether a CNN transfers is now an open and well-posed question, with a derived ceiling
   (0.941) to measure it against and a crude baseline (0.300) it must beat.
3. **The sharpest consequence for the generative loop:** if the generator produces grid-like sequence the
   oracle can score it; if it produces genome-like sequence, the oracle is at 29% of ceiling and loses to GC
   content. **Generator and oracle must be validated on the same distribution**, and right now they are not.
