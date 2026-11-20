# The +0.0271 "fine-tuning made it worse" delta is NOT established; losing to a 3-mer Markov chain IS

**Status:** `DIAGNOSTIC_METHOD_VARIANT`. Pre-registered at
`wiki/glm_generator_paired_noise_bar_2026-11-20.json` (frozen 2026-11-20, before any paired number
existed). Grid: `wiki/glm_generator_paired_sweep_2026-11-20.json` (40/40 cells, committed verbatim
before this reading). Verdict recomputed reproducibly by
`scripts/glm_paired_sweep_verdict.py` -> `wiki/glm_generator_paired_sweep_verdict_2026-11-20.json`.

## The correction

The carried claim was that LoRA fine-tuning moved the generator's distinguishability the WRONG way
by **+0.0271** (zero-shot 0.7787 -> fine-tuned 0.8058, mode-both, lower is better). Both numbers were
a single draw at the falsifier's default seed 7.

Scored at **ten seeds**, with every model scored at the same seeds so each difference is paired:

| statistic | value |
|---|---|
| mean paired delta | **+0.0458** |
| spread of the paired deltas | **0.2072** |
| seeds against the committed sign | **3 of 10** (seeds 1, 4, 7) |
| spread / committed effect | **7.7x** |
| delta at seed 7, this run | **-0.0142** — the OPPOSITE sign to the committed +0.0271 |

Frozen rule: same sign on every seed AND mean greater than the deltas' own spread.
Neither holds -> **`DELTA_NOT_ESTABLISHED`**.

**This does NOT establish the opposite.** The mean delta is positive, and slightly larger than the
committed value, so the DIRECTION is weakly consistent with the original reading. What fails is the
claim that the effect is resolved: a third of the seeds reverse it, and the spread is nearly eight
times the effect.

## What SURVIVES, decisively

Against its **own seed-matched** 3-mer Markov null, the generator is worse at **10 of 10 seeds**,
mean gap **+0.1186** (range +0.0794 to +0.1660). It never once beats it.

So the embarrassing half of the original finding is robust and the other half is not. The two were
reported together; they should not have been, and the distinction is the measured one: a 0.119 gap
survives a 0.207 spread, a 0.027 gap does not.

## The trajectory: unresolved, with a pattern I am NOT calling a finding

**`TRAJECTORY_UNRESOLVED`.** No consecutive epoch step is sign-consistent across seeds, and the
delta spread (0.2072) exceeds the total epoch-1-to-epoch-3 range (0.1257).

The suggestive part, recorded because suppressing it would be its own distortion: step means rise
monotonically (+0.0016 -> +0.0187 -> +0.0255) and positive counts rise (6 -> 8 -> 9), which is what
progressive over-training would look like. The frozen rule says UNRESOLVED and it is **not
upgraded** — its own text prescribes more draws, not a conclusion. The first epoch does essentially
nothing (mean +0.0016, 6 positive / 4 negative).

## Two defects in my own instrumentation, both found by this run

**1. The in-kernel gate was wrong, and the artifact says so.** The grid's own `verdict` field reads
`INSTRUMENT_INVALID` with "the Markov null did not reproduce". That is a gate defect, not a fact
about the data. The gate required every cell's null within 1e-6 of `BAR_BOTH = 0.6679` — but 0.6679
is the null's value **at seed 7 only**. The falsifier takes a single `--seed`
(`scripts/glm_generator_falsifier.py:124`, default 7) and threads it into `build_natural(...,
a.seed)` at line 185, which at line 57 does `random.Random(seed).shuffle(wins)` and so selects the
natural windows and the fit/eval split, and into `MarkovGenerator(order=..., seed=a.seed).fit(...)`
at line 284. The null is therefore refit per seed, and in this grid takes **10 distinct values from
0.5674 to 0.7024**. The corrected gate anchors reproduction at seed 7 alone.

**2. The pre-registration's band was mislabelled.** It calls the across-seed spread
"generation-sampling variation". Because the seed also moves the corpus split and the null refit,
the spread is **measurement variation**: generation + corpus split + null refit. The PAIRED delta is
unaffected — within one seed both models see the same windows, the same split and the same null, so
the confound cancels inside each delta, which is a second and unanticipated justification for
pairing — but the band's name was narrower than what the seed actually moves. The bar is amended
with this correction rather than edited, since it was frozen before this was known. It also left the
precedence between `FLAT_NOT_EPOCHS` and `TRAJECTORY_UNRESOLVED` unstated when both hold; the more
specific one is reported and the omission is recorded.

## The variance that dominates is the one this design does not bound

The bar recorded, in advance, that seed pairing varies **generation only** and leaves TRAINING
variance unmeasured. That caveat turned out to name the largest term.

- Zero-shot at seed 7 reproduces the committed **0.7787 EXACTLY** (delta 0.0). The scoring path is
  deterministic and the instrument is validated.
- Epoch 3 at seed 7 is **0.7645** against the committed **0.8058** — a drift of **0.0413**, which is
  *larger than the 0.0271 effect under test*.
- The loss curves differ: **7.43638 / 7.16777 / 6.91858** here versus **7.4361 / 7.1655 / 6.9148**
  committed. Same method, same data, same seed — so training is nondeterministic.

A single training run therefore cannot support a delta claim of this size, whatever the seed count.
Bounding it needs repeated training at different torch seeds, which is explicitly out of scope here.

## What this changes, and what is a user decision

Changed: nothing in the frozen surface. `scripts/kaggle/glm_finetune_kernel.py` is unchanged
(sha256 `6bdbd02641abedc5`), `verdict()` and the five method constants are untouched, and the two
committed arm artifacts keep their measured numbers — what changes is what may be CONCLUDED from
their difference.

**Open and NOT taken here:** whether to propagate this correction to `CLAUDE.md`, the
`features/glm-generator-finetune/` record, and the project memory, all of which currently carry
"fine-tuning made it worse" as settled. `CLAUDE.md` edits are user authority.

## Honest limits

- One organism (E. coli K-12), one corpus, **one training trajectory**.
- An upstream window is a PROXY for a promoter, not a mapped transcription start site.
- A k-mer discriminator sees LOCAL composition while the published failures are LONG-RANGE, so a low
  score is necessary-not-sufficient — and this run produced no low scores.
- `n=120` generated sequences per cell, the committed geometry, kept for comparability rather than
  because it is the right choice for variance. Raising it would attack the spread at its source.
- The spread is measurement variation, not generation variation (see defect 2).
- 10 seeds bound the spread weakly; the mean delta's own uncertainty is not reported because the
  frozen rule is a sign-and-spread rule, not an interval.
