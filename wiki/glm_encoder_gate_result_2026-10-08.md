# G-A: a real conv encoder does NOT reliably beat one GC number on genomic promoter tiles

**Frozen verdict `WEAK`; paired analysis says that WEAK is not distinguishable from
`DOES_NOT_GENERALISE`.** `wiki/glm_encoder_gate_2026-10-08.json`, schema `glm-encoder-gate-v1`. 44,106
peak tiles (14,913 active / 29,193 inactive), leave-peak-out, registered 10-seed list, CPU.

| arm | median Spearman | AUROC (active) | params |
|---|---|---|---|
| `gc` — one GC number | **+0.3112** | 0.6672 | 2 |
| `encoder-w20` | +0.3163 | 0.6739 | 2,657 |
| `encoder-multi` (6/20/75) | **+0.3234** | 0.6774 | 13,217 |
| `encoder-shuffled` (null) | −0.0240 | 0.4892 | 13,217 |

Median margin **+0.0122** against a pre-registered 0.05 bar. `registered_protocol_used=true`,
`null_clean=true`, `reconcile_ok=true` — and the reconcile anchor landed on **0.3000 exactly**, so this is
the same split and data as the published baseline.

## The median is not the instrument. The paired test is.

Every arm sees the *identical* leave-peak-out split per seed, so the margins are genuinely paired:

| arm | wins | sign-test p | mean margin | 95% CI |
|---|---|---|---|---|
| `encoder-multi` | **6/10** | 0.754 | +0.0125 | **[−0.0099, +0.0349]** — includes 0 |
| `encoder-w20` | **5/10** | 1.000 | +0.0010 | [−0.0202, +0.0222] — includes 0 |
| `encoder-shuffled` | 0/10 | 0.002 | −0.3288 | [−0.3597, −0.2980] — excludes 0 |

`encoder-multi`'s per-seed margins are **+.024 +.038 +.022 +.032 −.032 −.022 −.020 −.020 +.018 +.085** —
positive on six seeds, negative on four, sign-flipping. **The positive median is carried by seed
variation, not by a reliable edge.** `encoder-w20` is a literal coin flip at 5/10.

**The null's CI is what makes this a finding rather than a blunt instrument**: the same paired test
detects the shuffled arm's deficit at p=0.002 with a CI nowhere near zero. So the test *can* resolve a
real difference on 10 seeds; it does not resolve one here.

**Reusable, and this project had already recorded it**: a difference of medians is not a paired
comparison. The frozen verdict rule compares medians and returned `WEAK`; the paired instrument on the
same numbers says *not distinguishable from zero*. The verdict field is deliberately left as the frozen
rule's mechanical output with the paired block beside it under its own key — re-sizing a pre-registered
bar after seeing the result is an authority call, not an analysis choice.

## So Step 3's stop gate was substantively right, and the build was still worth doing

The encoder was built on an explicit user decision **against** this gate's prediction
(`wiki/glm_tile_headroom_2026-10-08.json`, position gain −0.0068 against +0.02, verdict
`NO_POSITIONAL_HEADROOM`). That decision was correct to make, and the result is stronger than either
outcome alone would have been:

- **The gate's two honest limits were real and are now closed.** Its proxy was positional GC in 10 bins
  under a *linear* ridge, so it could not see a motif at a precise offset and could not represent a
  composition×position interaction. A conv encoder escapes both. That is why the gate recommended
  ratification instead of closing G-A itself.
- **Having escaped them, the encoder still does not win.** The negative is no longer a prediction from a
  weak proxy — it is a measurement from the model class the proxy stood in for.

A caution against over-reading the direction: the median margin *is* positive, and on the first two seeds
(+0.024, +0.038) it looked like a clean win. Ten seeds and a paired test are what resolved it. Two
favourable seeds are not a result.

## The encoder is not broken, and that is established separately

A negative is only interpretable if the encoder can learn position *at all*. Three independent checks say
it can:

1. **Planted motif.** A 6 bp motif at a fixed offset — exactly the signal the gate's proxy is blind to —
   is learned to Spearman **0.84** against a ceiling of 1.0 (`tests/test_glm_encoder.py`). The ceiling is
   computed, not assumed; see below.
2. **The null collapses.** Shuffled training labels give −0.0240 median and lose on 10/10 seeds, so the
   architecture is not scoring by some pipeline artifact.
3. **It tracks the incumbent closely rather than failing.** 0.3234 vs 0.3112 is a model that has learned
   roughly what GC content already carries — not one that failed to train.

**Both framings agree**, which matters because regression may be the wrong frame here: with 29,193
inactive tiles a Spearman over everything can be substantially ranking noise, so the gate also reports
AUROC on the `active` column. `encoder-multi` 0.6774 vs `gc` 0.6672 — the same small, unreliable edge. A
finding that appeared in only one framing would have been about the framing.

## Defects found while building, all mine

- **The AUROC second framing would have measured nothing.** `fit_encoder` did not return predictions, so
  the gate's `r.pop("auroc_active", None)` was always `None` — a column that renders as a framing while
  containing no data, and the plan's own risk flag requires that framing. Both the encoder arms and the
  incumbent now return predictions; scoring only the encoder would have left the second framing with
  nothing to compare against.
- **My planted-signal bar was mis-specified and the model took the blame first.** It failed at 0.306 and
  the first reading was "the encoder cannot learn position". Probing instead of patching: 0.306 was
  under-training (75 Adam steps), and more steps reached 0.71 — then plateaued. The *ceiling* was the real
  story: with within-class noise, truth has n distinct ranks while a binary detector emits two values, so
  a **perfect** detector is tie-capped at **0.8661**, and I had set the bar at 0.8 — asking for 92% of an
  unreachable maximum. The model's 0.71 was 82% of the achievable ceiling. Removing the noise makes the
  ceiling exactly 1.0 and the encoder reaches 0.84. **Third mis-specified bar in this plan**, same root
  cause each time: a threshold frozen before the mechanism analysis that says what a correct result can
  achieve.
- **The inner-split grouping was hardcoded to `peak`.** Now a caller-supplied `group_key`, and
  `tile_peak_key` REFUSES a record with no peak rather than defaulting to `""` — a silent default would
  collapse every row into one group, and a later convenience fallback to a random inner split would leak
  on a substrate that needs a blocked one.
- Two of my *test* expectations were wrong while the implementation was right: the reconcile message
  starts a sentence, and `auroc([.1,.2,.3,.4],[0,1,0,1])` is 0.75 not 0.5 (interleaved is not a coin
  flip — the positives sit at ranks 2 and 4).

## Honest limits

- **This is THIS encoder at THIS size** — 32 channels per branch, 13,217 parameters — on THIS substrate.
  Strong evidence against the conv-encoder lever on genomic promoter tiles; **not** a statement about
  learned representations in general, and not about G-B.
- The encoder is deliberately small because **capacity has already hurt here**: 6-mers (4,096 features)
  scored worst of every arm at 0.2022. A larger encoder is a different experiment. The shuffled null is
  what would detect it over-fitting, and it is clean.
- Early stopping fires at 8–16 of 40 epochs on every seed, so the result is not an under-training
  artifact — but it does mean the inner validation split, not the epoch cap, is choosing the model.
- **Everything is read against the derived ceiling of 0.941**, never against 1.0.
- The 0.05 generalises bar and the 0.05 null bar are **ASSERTED**, not derived.
- Peak-blocked splits hold out whole peaks, but a tile can still share a little sequence with a training
  tile across a peak boundary.

Reproduce: `uv run python scripts/glm_encoder_gate.py --seeds 10 --no-download` (~35 min CPU, needs
`D:/dna_decode_cache/mpra`). Tests: `uv run pytest tests/test_glm_encoder.py tests/test_glm_encoder_gate.py -q`
(55 tests; the gate's are offline and torch-free).
