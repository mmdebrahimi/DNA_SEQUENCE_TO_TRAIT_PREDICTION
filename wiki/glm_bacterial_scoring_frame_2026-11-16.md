# The bacterial promoter number and the human number measure different tasks

**The root cause of the 20-25x gap is, at least in large part, the SCORING FRAME -- not bacterial
biology.** Measured 2026-11-16 on `data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz`,
locally, CPU, using the same `leave_peak_out` splitter and the same 10 seeds as the published arms.

## The headline: the comparison flips sign, unanimously

| frame | n | GC | 4-mer ridge | margin | learned wins |
|---|---|---|---|---|---|
| **pooled** (the published frame) | 44,106 | **+0.3112** | +0.1438 | **-0.1793** | **0 of 10 seeds** |
| **active only** (the strength question) | 14,913 | **-0.0464** | +0.1320 | **+0.1820** | **10 of 10 seeds** |

On the frame the published arms were scored on, a learned sequence representation is *worse than one
GC number on every single seed*. Restricted to the active tiles, the same representation *beats GC on
every single seed* by +0.18.

## Why: GC is a detector, not a strength predictor

> **RETRACTION, same day.** This section originally argued that *"a single binary split of the label
> recovers 0.8194 of the 0.941 ceiling -- 87% of the achievable ranking signal"*. **That statistic is
> vacuous and the inference is withdrawn.** Spearman is Pearson on ranks, and a threshold indicator is 1
> on exactly the top *k* ranks, so the quantity is a closed form in the split proportion alone:
> `rho = sqrt(3k(n-k)/(n^2-1))`, which for `n=44,106`, `k=14,913` equals **0.8194 exactly**. Four
> unrelated distributions at that split -- uniform, Gaussian, heavy-tailed, and a deliberately bimodal
> mixture -- all return 0.8194 to four decimals, and the maximum is 0.866 at any median split of any
> continuous label. It measured the split, never the label. **The finding below does not depend on it**
> (see "what survives"). A valid replacement would be a quantile-binned reproducibility curve, which is
> distribution-sensitive and threshold-free; it has not been run.

The `active` flag is a **pure threshold on the label itself** -- `expn_med >= 0.5778`, with **zero
overlap** on either side. The label is split 29,193 tiles in a narrow band `[0, 0.578]` against 14,913 in
a wide one `[0.578, 158.3]`.

**What survives, and it is the part that was always doing the work.** GC content separates those two
groups and does nothing else. Pooled it scores +0.311; **within actives it
scores -0.046, i.e. nothing.** The mechanism is unremarkable: AT-rich DNA melts more easily, so GC
predicts whether a window can form an open complex -- a threshold property -- and says nothing about
how strong a real promoter is.

The learned representation has the mirror-image profile: a *weaker* detector (0.144 against GC's
0.311) and the *only* carrier of within-active strength signal (+0.132 against GC's -0.046).

**So the published +0.0122 conv-vs-GC margin is the near-cancellation of two opposite effects, not an
absence of learnable signal.** The pooled metric rewards precisely what GC is good at and penalises
precisely what the learned model is good at.

## Reproducibility is not uniform either

| subset | replicate Spearman | implied ceiling |
|---|---|---|
| all tiles | 0.7945 | 0.9410 |
| active | 0.8849 | 0.9690 |
| inactive | 0.3927 | 0.7510 |

Two-thirds of the pooled scoring set carries a label whose own two replicates agree at **0.393**.
(Both subsets are threshold-defined, so range restriction attenuates these; they are not directly
comparable to the pooled figure.)

## The objection to this analysis, and why the finding survives it

Restricting to actives **conditions on the outcome**, which attenuates correlations. But it attenuates
*both* predictors identically, and they move in opposite directions: GC collapses from +0.311 to
-0.046 while the 4-mer ridge holds at 0.144 -> 0.132, a drop of 8%. A shared range restriction cannot
produce a differential of that shape unless the two predictors carry different *kinds* of signal --
which is the conclusion, not an assumption.

## The project had already found this failure mode and the fix was incomplete

`dna_decode/glm/genomewide.py`'s `Tile` docstring records it: a random-fragment library asks a
rare-event **detection** question while the designed grid asked a **strength** question, and peak tiles
were adopted to ask the grid's question on real genomic sequence. The step was right and it did not go
far enough -- **nobody measured the active fraction within peak tiles, and it is 34%.**

## What this does and does not do to the 20-25x headline

It makes it **suspect for a specific, testable reason. It does not refute it.** The eukaryotic arm's
own first honest limit says it cannot say which of four differences -- assay, element length, kingdom,
library construction -- carries the effect. This says **library construction is a live candidate
carrying part of it**: the bacterial side is a genome-wide detection-heavy scan, the human side a
designed library of candidate elements with a continuous `log2(RNA/DNA)` readout and no binary active
flag anywhere in the pipeline. Weakly corroborating: the human GC baselines (0.2362, 0.2200) sit well
*below* the bacterial pooled GC (0.3112), which is what less detection-dominance looks like.

**No arm artifact, verdict field or frozen file was edited.** Re-interpreting a published verdict is an
authority call; recording a new measurement beside it is not.

## Honest limits

1. **The decisive measurement is missing.** Whether the human label is similarly bimodal could not be
   checked: the data sits on a detached volume and its loader needs `openpyxl`, deliberately not
   installed. Until then, the frame-artifact share of the 20-25x is **unknown**.
2. **The learned proxy is a 4-mer ridge, not the conv.** The conv scored 0.3234 pooled, *above* GC;
   the ridge scored 0.1438, well below. So the conv is a materially better detector than this proxy,
   and its active-only score is unmeasured.
3. The active frame conditions on the outcome. It is a legitimate way to ask the strength question,
   but it is not "the correct frame" and its absolutes are not comparable to pooled ones.
4. **Even on the active frame, +0.132 against a 0.969 ceiling is small.** Bacterial promoter strength
   is still hard. This corrects the *diagnosis*, not the difficulty.
5. One substrate, one organism, one assay series. Nothing here generalises to bacterial regulation.

## Cheapest next measurements, in order

1. **Re-run the published conv (6/20/75, 13,217 params) on the 14,913 active tiles**, same 10
   registered seeds. Architecture-matched, local data, CPU. This is the number that settles whether
   the bacterial substrate supports a human-sized learned gain.
2. Measure the human lentiMPRA activity distribution for bimodality.
3. Re-score the best bacterial arm (PromoterAtlas embeddings, pooled 0.5076, pooled margin over GC
   already +0.1964) on actives only -- the realistic upper bound on this substrate.

## Pipeline validation (run before any of the above was read)

`n=44,106`, `active=14,913`, `inactive=29,193`, single-replicate `0.7945`, ceiling `0.9410`, GC
`0.3110` against the published `0.3112` -- all matched. The script **refuses** the decomposition on any
mismatch, and that refusal fired once on a real bug: the `active` column holds the literal strings
`active`/`inactive`, not `0`/`1`, so on the first run every tile landed in one bucket.
