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

---

# Third axis, and a mechanism: flatness is predicted by the axis's own dynamic range

Nitrogen was run to test a **pre-registered prediction**, not to repeat a result. Flatness had risen
61.2% → 68.2% from 4 media to 25 carbon sources; if a gene's ratio is flat because the *axis* carries
little dynamic range, nitrogen should be flattest of all — six of its thirteen conditions were already
known to give identical wildtype growth.

**It is.** And the direction result holds a third time.

| axis | conditions | genes | **flatness** | within-gene AUROC (non-flat) | p |
|---|---:|---:|---:|---:|---:|
| media4 (Orth) | 4 | 67 | 61.2% | 0.7308 (n=26) | 0.001 |
| carbon (Keio) | 25 | 217 | 68.2% | 0.8133 (n=69) | 0.0005 |
| **nitrogen (Keio)** | 13 | 155 | **75.5%** | 0.7088 (n=38) | 0.0005 |

## The predictor, measured on all three with one yardstick

The prediction was fired with dynamic range measured on nitrogen *only*, which makes it an ordering with
one measured predictor rather than a relationship. `scripts/fba_axis_dynamic_range.py` measures it on all
three — cheaply, since it needs **wildtype** growth per condition only (42 LP solves, no deletion panel):

| axis | distinct growths | distinct fraction | CV | flatness |
|---|---:|---:|---:|---:|
| media4 | 4 / 4 | **1.00** | 0.574 | 61.2% |
| carbon | 21 / 25 | 0.84 | 0.473 | 68.2% |
| nitrogen | 8 / 13 | **0.615** | **0.308** | 75.5% |

**Monotonic on both summaries.** Two are reported because they answer different questions —
`distinct_fraction` asks whether the model can tell the conditions apart *at all*, CV asks *by how much*;
an axis can score 1.0 on the first while barely spreading.

## What this buys

Flatness is the dominant term in the switch failure — it is 61–76% of the genes on every axis, and it is
unreachable by any change of readout. It now has a **predictor that can be measured before committing to
an axis, in seconds, without running a single deletion.**

That is the honest, quantified version of "axis choice is a free lever": **to look for conditional signal,
pick an axis whose wildtype growth actually spreads.** Nitrogen was the worst available choice and is now
measured as such rather than suspected.

**Honest limits.** n=3 axes cannot establish a relationship — this says only that the pre-registered
direction survives a common yardstick. The three axes also differ in substrate, label source and gene set,
any of which could drive flatness instead; the dynamic-range summaries are simply the one thing now
measured identically across all three. And the causal story is unverified: a plausible alternative is that
axes the model resolves poorly are also axes whose *genes* are peripheral to its stoichiometry.

Reproduce: `uv run python scripts/fba_axis_dynamic_range.py` (seconds) and
`uv run python scripts/fba_within_gene_ranking.py --axis nitrogen` (needs `feba.db` on D:).

---

# The denominator was wrong: the model is SILENT on ~90% of genes and ~70% right when it commits

Everything above quotes exact-set match as *"23 of 217"* on carbon, 10.6%. That denominator scores the
model against a target it **structurally cannot hit** for most of those genes.

Every conditionally-essential gene is, by definition, essential in some conditions and not others. So a
**constant** prediction can never match one — not on a good day, not with a better threshold. And the
model's call is constant for any gene whose ratio is flat, plus any gene whose ratio moves but never
crosses the 1% cutoff.

Verified before relying on it, in the committed carbon artifact: `commit_strata.predicted_constant` =
**184 genes, 0 exact-set matches**. My independent recomputation reproduces it exactly (184 constant-call
genes → 0 hits; 33 committing genes → 23 hits).

## Three nested strata, all three axes

| axis | flat (one ratio for every condition) | varies, never crosses cutoff | **call varies — model commits** | exact when committing | commit rate |
|---|---:|---:|---:|---:|---:|
| media4 | 41 | 22 | 4 | 3/4 = 75% | 6% |
| **carbon** | 148 | 36 | **33** | **23/33 = 70%** | **15%** |
| nitrogen | 117 | 22 | 16 | 8/16 = 50% | 10% |

**media4's 75% is on four genes and means nothing on its own.** Carbon (n=33) is the powered one;
nitrogen (n=16) sits between.

## Is 70% impressive? Anchored, not asserted

An exact-set "hit" means naming *precisely* the right subset of 25 conditions. The null: keep the model's
own count of essential conditions but place them at random. A gene whose predicted count differs from the
truth's cannot match at any placement and contributes zero.

| axis | observed exact hits | chance expectation |
|---|---:|---:|
| media4 | 3 | 0.75 |
| **carbon** | **23** | **0.78** |
| nitrogen | 8 | 0.62 |

Carbon is **~30× chance**. The placement is real.

## What this does and does not say

**It does not say the model is good.** It says the model is a **high-precision, very-low-coverage**
conditional predictor: it declines to answer for 85–94% of genes, and is right about 70% of the time on
the rest. Coverage is the deficit, and coverage is capped by flatness — which is structural, unreachable
by any readout change, and predicted by the axis's own dynamic range (previous section).

**It does say the published headline mis-states which thing is broken.** *"10.6% exact-set"* reads as *the
model is bad at conditional essentiality*. The measured decomposition is *the model is accurate when it
speaks and almost always silent*. Those imply different next moves: the first says fix the predictor, the
second says the predictor's reachable set is small and the question is whether it can be widened at all.

It is also the **fourth** independent corroboration of `MIS_CONDITIONED = 0` — the model does not fire in
the wrong place. Four measurements, four different metrics, same shape.

## Honest limits

- The strata are defined by the model's own output, not by the truth, so they are not a post-hoc
  cherry-pick — but "when it commits" is still a **conditional** accuracy and must never be quoted as the
  model's accuracy.
- Committing-stratum sizes are small on two of three axes (4 and 16). Only carbon's 33 is reasonably
  powered.
- The chance null holds the model's own predicted count fixed. A null that also let the count vary would
  be even lower, so this is the conservative choice.

Reproduce: `uv run python scripts/fba_within_gene_ranking.py --axis carbon` — the anatomy block prints
with every run.

---

# Adversarial review, 2026-08-30: the mechanism was unnamed, one null was too easy, and one axis fails

An adversarial pass over the committed conclusions produced three findings. Two are corrections to
what I published; one concern I raised turned out to be refuted by measurement.

## 1. ACCEPTED, and the most valuable finding: the mechanism is INFEASIBILITY, and I did not name it

Verified in `wiki/fba_conditional_carbon_2026-08-13.json` — an independently produced, previously
committed artifact:

| | |
|---|---:|
| committing genes | 33 |
| committing genes touching a nonoptimal/infeasible cell | **32** |
| exact-set matches | 23 |
| exact matches whose essential calls are **ALL** suspect | **23 / 23** |

**Every single exact match is driven entirely by infeasibility events.** Under this repo's prior
finding (an infeasible deletion LP is genuine essentiality — hard ATPM floor, 38 of 39 such cells
experimentally confirmed) that is not an artifact. It is the **mechanism**, and calling the result
"a high-precision conditional predictor" invited a general reading it does not support.

**Corrected phrasing:** *the model's conditional signal is essentially all hard feasibility breaks —
deleting the gene makes growth infeasible on some media and not others. Where that happens it names
the exact set ~70% of the time. Everywhere else it is silent.*

This also explains the rest better than my framing did: "flat" means no feasibility break anywhere on
the axis, and the signal is **binary, not graded** — which independently predicts that a ranking rule
cannot help. It **strengthens** "do not build the relative rule" rather than weakening it.

## 2. ACCEPTED: the chance null treated conditions as interchangeable

The original null asks only "how many placements of *k* essential conditions exist" — `1/C(n,k)`. If
true essentiality concentrates in a few substrates *and* the model tends to break on those same
substrates, a model could beat that null by matching the marginal shape rather than the per-gene
placement.

Added a **second, strictly harder null**: shuffle the TRUTH matrix preserving **both** margins (every
gene keeps its essential-condition count, every condition keeps its essential-gene count) with the
model's predictions held fixed. Reuses the repo's tested `nulls.curveball_shuffle`, which raises if a
margin ever breaks rather than returning an invalid null. 200 draws.

| axis | observed | interchangeable null | **both-margins null** (mean / p95 / max) | verdict |
|---|---:|---:|---|---|
| media4 | 3 | 0.75 | 1.275 / 2.0 / **3.0** | **NOT distinguishable — the null's max equals the observed** |
| carbon | **23** | 0.78 | 1.035 / 3.0 / 4.0 | clears it decisively |
| nitrogen | **8** | 0.62 | 1.615 / 3.0 / 4.0 | clears it |

The concern was directionally right — the harder null is higher on all three axes — and **immaterial
on carbon** (23 against a 200-draw maximum of 4). But it **changes the verdict on media4**: its 3
observed hits sit exactly at the stricter null's own maximum. I had already flagged media4's "75%" as
meaningless on n=4; it is now quantitatively confirmed as not significant, and should not be quoted.

## 3. REFUTED BY MEASUREMENT: the tiny-spread concern runs the other way

The worry: `FLAT_EPS = 1e-9` may count numerical noise as signal, inflating the primary. Sweeping the
spread floor (mean within-gene AUROC over non-flat genes):

| axis | ≥1e-9 | ≥1e-6 | ≥1e-4 | ≥1e-2 |
|---|---|---|---|---|
| media4 | 0.7308 (26) | 0.7308 (26) | 0.7308 (26) | **0.7500** (24) |
| carbon | 0.8144 (68) | 0.8144 (68) | 0.8163 (67) | **0.8345** (62) |
| nitrogen | 0.7088 (38) | 0.7088 (38) | 0.7280 (37) | **0.8196** (29) |

**Tightening the floor RAISES the primary on all three axes** — the near-noise genes were diluting the
result, not inflating it, so the published numbers are conservative. Carbon's non-flat spreads have
median 0.178 and max 1.0; only 2 of 69 fall below 1e-4.

(The `≥1e-9` column shows carbon n=68 against the published n=69 because the artifact's `spread` field
is rounded to 6dp and one gene's true spread lies between 1e-9 and 5e-7. A reporting-precision
artifact, not a second number.)

## 4. PARTIALLY ACCEPTED: demote distinct-growth fraction, lead with CV

**Distinct-growth fraction is mechanically confounded with condition count** — with 4 conditions,
all-distinct is nearly automatic; with 25 it is not. That half of the dynamic-range evidence is weak
and is demoted. **CV is scale-free and not mechanically tied to condition count**, and is
independently monotonic (0.574 / 0.473 / 0.308), so it carries the claim.

The prose is also softened: **"pick an axis whose wildtype growth spreads" is a pre-run TRIAGE SCREEN,
not a predictor.** Knockout flatness can also be driven by gene-set composition, substrate biology and
label source, none of which n=3 axes can separate. The script's own `honest_limit` already said n=3
cannot establish a relationship; the prose now stays at that level.

## Still open

- Of the 10 carbon committing MISSES, are they also infeasibility-driven but on the wrong substrate,
  or are they the genuinely graded threshold cases?
- Are condition marginals highly uneven across the 25 carbon sources, especially for singleton true
  essentiality? The margin-preserving null bounds the answer's consequence but does not describe it.
- Does ratio-spread filtering hide absolute-growth issues, where a low wildtype growth makes a small
  absolute change look like a large ratio change?
