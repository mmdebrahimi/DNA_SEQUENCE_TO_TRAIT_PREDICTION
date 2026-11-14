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

### MEASURED UPDATE 2026-11-14 -- two of the four decisions now rest on measurement, not judgement

Everything above is the original 2026-11-01 text and is unchanged. This block is additive. Nothing
frozen was edited to produce it: the five pinned files, `uv.lock`, the acceptance bar, the
conditioning artifact and the arm artifact are all byte-unchanged, the published `FRAME_PARTIAL`
still carries `resolved: false`, and the 20 tests that rebuild the grader's input from the
committed artifact and re-grade it still pass.

**DECISION 3 is answered, and BOTH halves of the answer are measured.** The recommendation above --
a shuffled-window control -- is **dead, and would have produced a confident wrong answer.**
`frame_features(located, mask, cds_rows, *, tile_length)` takes window information ONLY as the
union mask, so relocating windows changes three of its five columns and leaves the two
`cds_rows`-derived columns byte-identical between treatment and control. Measured on the same ten
committed folds (`wiki/promoteratlas_frame_column_split_2026-11-05.json`): those two CDS columns
**alone** score `0.3183` standalone against the full five-column table's `0.3170`, and the
comparison has to be paired within seed to mean anything because the seed spread exceeds the gap --
the CDS pair is at least the full table on **6 of 10** seeds, while the three window columns are
below it on **10 of 10**. Handed to the GC baseline the CDS columns absorb `+0.0705` against the
frame's `+0.0687`, so the control would have absorbed MORE than the treatment, driven
`attributable_absorption` negative, and fired `INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE` for a
reason that is an artifact of shared columns rather than a fact about position. Cost of finding
that out: about seven seconds of compute.

The replacement is a **circular rotation of the whole annotation**
(`dna_decode/glm/annotation_shift.py`), which regenerates BOTH halves of the table from the rotated
rows so nothing is invariant. It is matched by construction on gene count, the gene-length
multiset, strand composition, total CDS-covered bases and coding density -- exactly, not
approximately, and the runner refuses any offset that fails to preserve them. Independent
per-window relocation cannot do this: at `0.8626` coding density most relocated windows land inside
genes, so `in_window_any` would start tracking "inside a gene body" instead of "upstream of a gene
start", substituting a different confound for the one being removed.

**It was run, and the attribution survives it** (`wiki/promoteratlas_frame_shift_null_2026-11-14.json`,
five offsets, the same ten folds for every arm, true arm recomputed in the same run and reproducing
`0.3170` and `+0.0687` exactly):

| quantity | true annotation | rotated (5 offsets) |
|---|---|---|
| frame absorption | `+0.0687` | `-0.0009` (`-0.0050` to `+0.0015`) |
| frame features standalone | `0.3170` | `0.0114` |
| in-window margin | `0.2689` | `0.1774` |
| out-of-window margin | `0.0230` | `0.1938` |
| in/out collapse | `+0.2459` | `-0.0164` (`-0.0751` to `+0.0003`) |

So the absorption is **specific to true annotation positions**, not a property of handing a ridge
five annotation-shaped columns: `attributable_absorption = +0.0696` against a control whose premise
holds by construction instead of by hope. The rotated table being predictively dead on its own
(`0.0114`) is the sanity check that the rotation does what it claims. And the in/out collapse does
**not** reproduce under rotation -- both rotated strata sit near `0.19`, which independently agrees
with the committed size-matched random null's `0.1860` and with the located subset's own margin
over GC of `0.1950`. Two differently-built nulls landing in the same place is worth more than
either alone, and the rotation is the stronger of the two because it matches the spatial clustering
of the strata and not merely their size.

**This is REPORTED, not graded, and it does not resolve the verdict.** The frozen bar does not name
this control; letting a control added after the fact select a verdict is the move the bar exists to
prevent. PARTIAL-versus-EXPLAINS still turns on the out-of-window stratum margin sitting `0.0004`
from its bar, which this control does not address. What it settles is the ATTRIBUTION half.

**DECISION 1: my recommendation above is now withdrawn as stated, because it is unachievable.**
"Require an interval, and until then treat PARTIAL-vs-EXPLAINS as open" means *permanently* open on
this design. On the ten per-seed values the 95 percent t-interval is `[-0.0149, 0.0517]`, a
half-width of `0.0333` -- **wider than the `0.02` threshold it must discriminate**, and containing
both zero and the bar. Reaching a half-width below the bar needs about **21** seeds; separating the
point estimate from the bar needs about **3,253**. Worse, the arithmetic mean is `0.0184`, **below**
the bar, while both estimators reported above (`0.0230`, `0.0204`) are above it -- so the spread of
defensible point estimates already straddles the threshold, which is a stronger statement than
"0.0004 away". A sign test is `5` of `10`, two-sided `p = 1.000`. More seeds is therefore the wrong
lever: the design is the limit, and "until then" would be never.

There is also a statistic mismatch that makes "an interval" under-specified until it is resolved:
the grader reads a **difference of medians** (`min_stratum_margin = 0.0204`) while the stability
block reports the spread of **per-seed paired margins** (median `0.0230`). An interval built from
the second does not bound the first.

**Revised recommendation: edit no frozen artifact.** Bars here are per-run pre-registrations --
there are 23 `*_acceptance_bar.json` files in `wiki/`, one per arm -- and this one has already
discharged its function. Freezing a NEW bar for the successor run, which names ONE statistic and
grades it, carries a bootstrap interval on that statistic, names `UNRESOLVED` as a first-class
outcome rather than forcing a choice between two substantive ones, and adds the missing gate below,
requires no outcome-neutrality proof and no revert pre-commitment, because nothing is edited. The
frame bar would take at most a non-grading `superseded_by` pointer, which touches no predicate. That
dissolves the problem instead of satisfying it. **Drafted, not frozen** -- the acceptance criteria
are yours.

**DECISION 2: still ratify, but not for the reason given above, and the record should say so.** The
eighth outcome was **INERT on this run's data**: re-grading from the committed artifact gives
`absorbs_unattributably = False` (C2 never fired -- `margin_over_gc_plus_frame 0.1262` is not below
`0.02` -- and C3's random null cleared the bar at `0.186`). With seven outcomes the grader would
still have reported `FRAME_PARTIAL`, **not** `FRAME_EXPLAINS`. So it should be ratified on
prospective grounds -- it sits before every substantive outcome in `CHECK_ORDER`, so it can only
ever demote, which satisfies the narrow checkable rule that a post-hoc outcome may only be an
abstention inserted ahead of all substantive outcomes -- and it should not be credited with having
caught this run's confound, because it did not.

**The gate that WOULD have caught it does not exist.** The artifact records
`matched_control_premise.control_is_at_least_as_strong_as_the_frame = False`, and a grep for that
field across `dna_decode/glm/frame_verdict.py` and `scripts/freeze_frame_bar.py` returns **zero**
hits: it is computed and discarded. The eighth outcome guards *"absorption did not exceed the
matched control"*; nothing guards *"the control was not matched"*. The missing gate is an
`INDETERMINATE_CONTROL_PREMISE_UNMET` reading a boolean the pipeline already produces. That is this
repo's own inert-control pattern in a new place, and it belongs in the successor bar rather than in
the frozen one.

**DECISION 4: leaving the row alone is right; the gap is that nothing DETECTS the divergence.**
Seven test files reference `arms_ledger` and none reconciles a candidate row's `pending` against the
arms ledger's `measured`. Meanwhile `advance_ranker`'s `next_action` is `table[0]` unconditionally,
so a stale `pending` row actively misdirects the next `--advance` run -- which already happened once
in this arc. The ownership argument against marking another session's candidate done is unchanged
and correct. A reconciliation check that reports divergence without editing either side is the
addition, and it belongs to no single session.


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
