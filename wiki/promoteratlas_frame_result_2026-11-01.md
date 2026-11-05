# Corpus-window POSITION takes a material share of PromoterAtlas's 0.5076 — and outside those windows the gain is indistinguishable from zero

<!-- doc-status: active -->
<!-- measured-on: 2026-11-05 -->

<!-- Execution record for plans/Frame_Conditioning_Does_Corpus_Window_Position_Explain_PromoterAtlas_0_5076_Plan/
     (candidate 7 of project_state/glm-learned-representation-2026-10-07.md), all 7 steps.
     Bar: wiki/promoteratlas_frame_acceptance_bar.json (frozen by code before any number existed).
     Result: wiki/promoteratlas_frame_conditioning_2026-11-01.json. Per-seed detail is gitignored at
     data/promoteratlas/frame_conditioning_detail_2026-11-01.json.
     THE FILENAME DATE IS THE PLAN'S, NOT THE MEASUREMENT'S: the plan was written 2026-11-01 and
     names these artifacts, and the bar plus the arm artifact's `memo` pointer already reference
     them, so the names are kept for referential integrity. The run is dated 2026-11-05 inside every
     artifact. -->

**Verdict: `FRAME_PARTIAL` — and this run does NOT resolve `FRAME_PARTIAL` versus `FRAME_EXPLAINS`.**
The distinction was decided by a point median landing **0.0004** above a 0.02 bar on a statistic
whose ten seeds split **5–5** across that bar with **four of them negative**. That is a coin flip,
not a finding. Read the bound, not the label.

**The bounded form, which is what may be quoted:** the corpus-window rule *as replicated* — 200 bp
upstream of an annotated gene start, 5 bp intergenic gate, sorted ordering, one organism's
annotation — absorbs a material part of the atlas margin, and **the advantage over GC outside those
windows is not distinguishable from zero.** It may be written neither as *"position excluded"* nor as
*"position explains it"*. `arms.native150 = 0.5076` keeps its measured value and the ceiling verdict
is untouched.

---

## DECISIONS NEEDED

**1. The bar grades a point median where this statistic needs an interval. Should it be changed?**
This is a change to a frozen bar, so it is yours, not mine. C3's weaker stratum margin is
`0.0230` by median of per-seed differences (`0.0204` as the difference of medians — the gap between
those two readings is itself a symptom), against the inherited `0.02`. Per-seed it is
`[0.042, 0.050, -0.020, 0.027, -0.035, -0.049, 0.043, -0.001, 0.019, 0.108]`. **My
recommendation: yes** — require an interval on the stratum margins, and until then treat
PARTIAL-vs-EXPLAINS as open. I did **not** re-size the bar after seeing the result; that is the one
move a frozen bar exists to prevent, so the verdict above is reported exactly as the frozen grader
computed it.

**2. Ratify or reject the eighth outcome.** The plan names seven outcomes; the grader has eight.
The addition is `INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE`, built because the plan's own
pre-execution review held that *shrinkage is not absorption unless it exceeds what a
matched-informativeness, non-positional feature produces* — and a control that cannot change a
verdict is not a control. With seven outcomes the grader would have been forced to report
`FRAME_EXPLAINS` over a measurement it could already see was confounded. **Recommendation: ratify**,
and amend the plan rather than the grader. The plan file is left unedited.

**3. The matched-strength control did not hold its own premise. Build the clean one?**
Its premise is *matched informativeness*, and the non-positional comparators turned out **weaker
standalone** than the frame: dinucleotide `0.2895` and k-mer4 `0.2811` against the frame's `0.3170`.
So `attributable_absorption = +0.0949` is **weakly controlled** — part of it may measure the frame
being a *better* feature rather than a *positional* one. **Recommendation: a shuffled-window
control** — identical feature construction with the window set randomly relocated along the genome.
That is matched by construction rather than by hoping two feature families happen to be equally
strong, it is the direct analogue of the subtractive arm's random-direction control, and it is
cheap. Not built here: it is a new control the frozen bar does not name, and adding one after seeing
the result needs your call.

**4. Candidate 8 is stale in the ledger and I did not touch it.** The G-A candidate table lists it
`pending`, while the arms ledger already records
`positional-non-mean-pooled-promoteratlas-150bp` as **measured**. That is another session's work,
and silently marking someone else's candidate done is the opposite of what the journal is for. Only
candidate 7 was retired.

---

## CORRECTIONS

**1. The plan's definition of `in_window_full` disagreed with the number it pinned, and the number
wins.** The plan's prose says tiles *"entirely inside one"* window; the 8,239 figure it pins was
measured by `scripts/frame_plan_feasibility_probe.py` as `c == TILE_LEN` against a single **union**
coverage mask. Those are different predicates — two overlapping windows can jointly cover a tile
neither contains. The measured definition is implemented, and `tests/test_frame_features.py`
constructs exactly that two-window case so the two cannot drift under one name.

**2. My own grader had two validity gates in the wrong order.** It checked frame informativeness
before the non-strengthened-baseline gate; the plan specifies the reverse (*"a non-strengthened
baseline beats both"*), and the plan is also right on the merits — a baseline that never received
the feature can itself make the frame look uninformative. Caught by the ordering test, not by
reading the code.

**3. The parent exposure memo's honest limit is AMENDED, not replaced.** Its bullet now carries a
positional clause with the original wording preserved above it. The linear-removal finding is **not**
wrong: removing a detection *direction* and conditioning on window *membership* are different
interventions, and the second is the one that bites.

**4. My forecast missed, and in the opposite direction from my recorded bias.** The frozen bar's
modal prediction was `FRAME_EXCLUDED_CONDITIONAL` at `p = 0.44`; the outcome was `FRAME_PARTIAL`
(`p = 0.25`). Standing record **2 of 7 → 2 of 8**. The recorded bias is *pessimism* — talking the
project out of runs that then succeed. This error runs the other way: over-confidence that the arm's
gain would hold up. A new data point, not another instance of the old pattern.

---

## What was measured

The harness was proved first. On the committed splits the cached embedding pass reproduced
`arms.native150` at **delta 0.0000**, so nothing downstream rests on a drifted representation. Every
frame-bearing number below is measured on the **located subset** (43,698 of 44,106 tiles) with its
own baseline re-measured there, so the atlas figure is `0.5108` rather than the full-cohort `0.5076`
and the two must not be compared.

**C1 — the encoder plainly knows where the corpus is.** Predicting window membership from the
128-dim embedding reaches AUROC **0.8346**, against **0.5342** for GC alone, which is barely above
chance. This is the positional restatement of the already-measured 0.9156 control-vs-tile
separability. It establishes that the information is *present*; it does not by itself say the
expression gain is *made of* it.

**C2 — handing the frame to the baseline absorbs a material share, and the gain still survives.**
GC alone scores `0.3159`; GC plus the five frame features scores `0.3846`. The atlas margin falls
from **0.1950** to **0.1262**, so the frame absorbs **0.0687** — well past the 0.02 bar, but with
0.1262 still standing. The frame alone scores **0.3170**, essentially matching GC, so it is a real
feature and the vacuity gate never fires.

The control series behaved in a way worth reporting rather than smoothing: adding non-positional
composition columns **did not shrink the margin at all — it widened it.** GC plus dinucleotide
(`0.2897`) and GC plus k-mer4 (`0.2814`) both score *below* plain GC on held-out data, so those
baselines lose ground under leave-peak-out while the frame gains it. Directionally that is
informative — the frame's absorption is not merely "more columns" — but it is also precisely why the
control fails its matched-informativeness premise, per DECISION 3.

**C3 — the cell the plan pre-committed as the stronger design, and the one that matters.** It
controls for window membership *by construction* rather than by regression, so unlike C2 it does not
require the frame to be linearly modellable. Within corpus-window tiles the atlas advantage over GC
is **0.2689**, tight across seeds at `[0.2361, 0.3213]`. Outside them it is **0.0230**, with four of
ten seeds **negative** and a seed range of `[-0.0492, +0.1081]`.

That collapse is not an artifact of stratifying. A **size-matched random stratification** — same
stratum sizes, membership shuffled — leaves a minimum margin of **0.1860** `[0.1609, 0.2350]`. So
splitting the cohort in two does not by itself destroy the margin; splitting it *on corpus-window
membership* does.

The two strata are not exchangeable and the spreads ship beside the margins for exactly that reason:
in-window expression has interquartile range `0.4735` against out-of-window `0.2710`, with the
standard deviations inverting (`3.380` versus `5.535`). Range restriction attenuates a rank
correlation for **both** channels, which is part of why the out-of-window margin is noisy — and why
reading it as a clean zero would be as wrong as reading it as a clean 0.02.

**The two cells disagree in strength, and the plan said in advance which to trust.** C2 says most of
the gain survives a frame-strengthened baseline; C3 says the gain is concentrated almost entirely in
corpus-window tiles. Both can hold at once, because C2's baseline is a *linear* model over five
crude features and cannot express what window membership does to expression, while C3 conditions on
membership directly. The plan called C3 the strongest cell before any number existed. It is also the
cell whose statistic this run could not resolve.

## Honest limits

- Scoped to the corpus rule **as specified** (200 bp, 5 bp gate, sorted ordering), never *as the
  authors ran it*: the upstream extractor threads its gate in **file order** over a GenBank file
  this repo does not have, and the two orderings differ by 26 windows and 146 tiles.
- One organism, one annotation build, one encoder, one rung.
- 325 ambiguous and 31 not-found sequences are excluded and counted, so 408 of the 44,106 scored
  tiles carry no frame. Emitted plus excluded reconciles exactly.
- 689 of the 44,106 rows share a sequence with another row, so a sequence duplicated across two
  peaks can land on both sides of a leave-peak-out split. A property of the substrate and of every
  arm on this rung, reported rather than fixed here.
- The frame features are crude by construction. Failing to absorb the whole margin bounds **these
  five features**, not every possible positional encoding.
- `attributable_absorption` is weakly controlled (DECISION 3).
- The verdict's deciding statistic is unresolved (DECISION 1).

## What did not change

`arms.native150 = 0.5076` keeps its measured value, per the bar's pre-commitment (d) — a measured
metric is not retro-edited to match a later finding. The ceiling verdict stays
`blocked:untried-arm`; this is a **limit on an existing measured arm**, not a re-read of it and not
a variant of it, and a positional or fine-tuned PromoterAtlas variant remains a NEW ARM under R6.
The frozen five and `uv.lock` are byte-unchanged, as are `genomewide.py` and `expression.py` (both
are embedded base64+gzip in the arm kernel with a byte-identity assert). The arm artifact's update is
purely additive: 53 insertions, zero deletions, and `arms{}` asserted byte-identical across the
write by the script that performed it. No float dump reached `wiki/`.
