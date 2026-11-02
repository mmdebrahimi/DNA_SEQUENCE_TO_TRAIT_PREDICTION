# The positional PromoterAtlas arm: does position-resolved latent structure beat the mean-pool 0.5076

<!-- doc-status: active -->
<!-- measured-on: 2026-11-01 -->

Graded by `wiki/positional_arm_acceptance_bar.json`, frozen before any number below existed.

## DECISIONS NEEDED

**None from this arm's verdict.** `ceiling_verdict` reads `blocked:untried-arm`, so the criterion is
unchanged by this arm and no scope decision falls out of it.

**It stays open on 3 declared-but-unrun arm(s), which are named rather than left to
be inferred:**

- `asgiven-position-resolved-promoteratlas-150bp`
- `finetuned-promoteratlas-encoder-150bp`
- `shipped-expression-heads-zeroshot-150bp`

`best` is unmoved at `promoteratlas-base-embeddings-frozen-150bp-padded-to-200`. The criterion cannot read `met` while any of these is
untried, which is the whole point of declaring them: an arm declared only once its
cheap sibling succeeds is an arm that evaporates.

## CORRECTIONS

**Two scoping corrections, both of which NARROW this arm from how its ledger rationale
declared it.** They were fixed in the bar before any number existed, not written after.

1. This arm CANNOT test the -35/-10 spacing claim its own ledger rationale asserts. genomewide.Tile exposes seq, expression, peak, strand, category, active, rep1, rep2 -- there is no TSS, no promoter register and no peak-relative offset, so a bin index is a position in an ARBITRARY 150 bp window, not a promoter coordinate. What is testable is narrower: does a position-resolved reduction beat a position-averaged one.

2. Because of that, the readout matching the stated spacing mechanism is NOT bins. Spacing is translation-invariant, so the anchor-free autocorrelation is a declared CO-PRIMARY testing a DIFFERENT hypothesis -- not a second shot at the same one.

**The orientation convention is no longer unresolved, and it resolved AGAINST orienting.**
The frozen bar records it as a 1-bit inference *"from a performance difference between
two rungs"* and lists it under what this arm does not establish. This run supplied that
difference: as-given scores **0.5086** against oriented **0.4893**, so reversing the bin axis
on minus-strand tiles **destroys about 0.0193 of signal rather than aligning the cohorts**.
The bar's own text is left UNEDITED — it was correct when frozen — and this resolution is
recorded beside it rather than written back into it.

**The consequence is load-bearing and uncomfortable: the arm's declared co-primary was the
ORIENTED rung, so the eligible primary was handicapped by the transform this same run
shows to be the wrong one.** The verdict stands unedited anyway, because as-given does not
clear the margin either — gain -0.0042, CI [-0.0079, +0.0017], which CONTAINS zero, so it TIES
the incumbent rather than beating it. A variant with as-given as the primary would be a
NEW ARM under R6, never a re-read of this one.

**And this is the multiplicity rule doing real work rather than decorating the bar.**
`bins10_mean_asgiven` at 0.5086 is the highest rung on the ladder and sits ABOVE the
incumbent 0.5076. It is a declared SECONDARY and INELIGIBLE as this arm's result, so it
is not the metric and does not set the criterion's `best`. Reporting the best of six
rungs would have claimed a gain — and even that rung's own interval contains zero.

**A further correction, to an inherited surface rather than to this arm.** INHERITED from the incumbent arm's IMPLEMENTATION, not from its prose. That bar's text reads '|median gain over GC| < 0.02' but its kernel computes no GC term at all. The implementation is the statistic that produced the fragment arm's +0.0517 and the atlas arm's 0.0203; a literal arm_shuffled-minus-GC would sit near -0.31 and could never clear a 0.02 bar. Discrepancy recorded rather than silently resolved.

## The result

**Verdict: `POSITIONAL_NO_GAIN`** — both co-primary intervals lie entirely BELOW the margin +0.0200

The incumbent is the mean-pool arm at median Spearman **0.5076**, read from its committed
per-seed table and never re-run, so this arm cannot move the number it is graded against.
The margin is **+0.02**, inherited verbatim from the tile-rung arms.

| rung | role | features | median rho | gain vs atlas | corrected CI | seed wins |
|---|---|---|---|---|---|---|
| `bins10_mean_oriented` | co-primary | 1280 | 0.4893 | -0.0187 | [-0.0219, -0.0160] | 0/10 |
| `autocorr` | co-primary | 2560 | 0.3975 | -0.1089 | [-0.1227, -0.1032] | 0/10 |
| `bins10_mean_asgiven` | secondary | 1280 | 0.5086 | -0.0042 | [-0.0079, +0.0017] | 4/10 |
| `bins2_mean_oriented` | secondary | 256 | 0.5000 | -0.0052 | [-0.0088, -0.0040] | 0/10 |
| `bins5_mean_oriented` | secondary | 640 | 0.4970 | -0.0108 | [-0.0142, -0.0080] | 0/10 |
| `bins10_max_oriented` | secondary | 1280 | 0.4954 | -0.0127 | [-0.0192, -0.0061] | 1/10 |

`gain` is the PAIRED difference against the committed atlas value on the same seed. The
interval is a paired bootstrap over the 10 seeds at the Bonferroni-corrected level: the two
co-primaries at alpha/2, the four secondaries at alpha/6.

**Only the two co-primaries are eligible as this arm's result, and that is enforced rather
than promised** — the verdict function cannot see the secondaries at all. A secondary rung
scoring well is a diagnostic, never the headline.

## The gates, which ran before any rung was scored

| gate | result | what it establishes |
|---|---|---|
| split equality | PASS | all 10 recomputed splits match the committed `n_test_per_seed`, so the arm is scored on the incumbent's held-out peaks |
| parity | PASS | the global mean DERIVED from the stored 10-bin means reproduces the committed embedding to max |dev| 1.058e-06 and its per-seed Spearmans (median 0.5076 vs 0.5076), so every rung is comparable to 0.5076 rather than merely adjacent to it |
| shuffled-label null | PASS | median +0.0053 over 10 seeds against the bar `|median| < 0.02` |

## What this does NOT establish

- No promoter anchor exists on this substrate, so a NO_GAIN says that position bins relative to an ARBITRARY window add nothing -- NOT that positional structure is absent. The negative must never be over-read as closing positional readouts generally.
- The orientation convention of `seq` is UNRESOLVED and is inferred from a performance difference between two rungs. That is a 1-bit inference, not established biology: only 576 of 22,285 minus sequences (2.6%) have their reverse complement among the plus sequences, and ZERO of those share a peak, so plus and minus tiles are distinct genomic windows rather than two orientations of one element.
- 689 duplicate sequences are known on this substrate and are UNCHANGED by this arm. A train/test-straddling duplicate affects this arm exactly as it affects the incumbent, so the PAIRED comparison is the right shape while the absolute levels inherit the issue.
- The encoder is FROZEN and its hyperparameters were chosen elsewhere; this is not a tuned positional model and a NO_GAIN does not bound what a fine-tuned one could do.
- Ten seeds is the inherited pairing. The fragment-pyramid arm needed 270 pooled seeds before its null cleared, so ten may be too few here for the same reason.
- The rung ceiling is 0.941 and replicate-derived. No tile-rung number may be read as a fraction of it without saying so.

**My forecast, recorded before any number existed:** `POSITIONAL_NO_GAIN` at p=0.45. 3 of 7. The fragment-null P1 was falsified twice, P2 mislocated its interval, and the atlas arm's own bar put p=0.55 on LOSES_ALONE and p=0.08 on what actually happened. DO NOT defer to these priors.

## Reproduce

```bash
uv run python -m scripts.positional_arm_probe
uv run python scripts/freeze_positional_arm_bar.py
uv run python -m scripts.positional_arm_extract
uv run python -m scripts.positional_arm_score
uv run python -m scripts.positional_arm_memo
```
