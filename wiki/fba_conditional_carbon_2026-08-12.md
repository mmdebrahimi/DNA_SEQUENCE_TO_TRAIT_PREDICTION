# FBA commits rarely, and is accurate when it commits — 25 carbon sources (2026-08-12)

> **This memo was rewritten the same day it was written.** Its first version claimed the 4-media
> "the model is not switching" finding had **reversed** — headline "0.0% constant across 20 shapes". That
> figure was **an artifact of a bug**, and the reversal claim is **withdrawn**. See *The withdrawn claim*
> below. What survives is real but narrower, and the corrected framing is more useful than either
> previous headline.

## The question

Does iML1515 predict *conditional* gene essentiality — a gene dispensable on one carbon source and
required on another? Measured two days ago on the Orth 2011 4-media screen (268 gene × condition cells),
the answer looked like a flat no. This re-asks it on the Fitness Browser RB-TnSeq compendium:
**25 carbon sources, 217 conditionally-essential genes, 5,425 cells**.

## Result — the honest framing

| | 4 media (Orth) | **25 carbon sources** |
|---|---|---|
| conditionally-essential genes | 67 | **217** |
| exact-set match | 3/67 = 4.5% | **23/217 = 10.6%** |
| per-cell agreement | 0.5709 | **0.7368** |
| best constant null | 0.5588 | 0.6623 |
| lift over null | +0.0121 | **+0.0745** |
| mean per-condition MCC | 0.0611 | **0.3864** |
| **model predicts a CONSTANT pattern** | **94.0%** | **84.8%** |

**The model is constant for the large majority of switching genes on BOTH substrates.** It is *less*
collapsed on the wider panel (84.8% vs 94.0%) and *more* accurate on every metric, but this is a
difference of degree, not a reversal.

### The sharper reading — commit rate

A constant prediction can never exactly match a two-sided gene, so **all 23 exact matches must come from
the 33 genes where the model predicts a varying pattern**:

- the model **commits** to a varying pattern for **33/217 = 15.2%** of genes
- where it commits, it is **exactly right 23/33 = 70%** of the time

That is the most informative statement available from this data: **FBA commits rarely, but is accurate
when it commits.** Neither "it is not switching" (the old claim) nor "it is switching" (the withdrawn one)
captures it.

**With the honest caveat that its commitments are on the easy end.** 25 of the 33 varying predictions call
essentiality in exactly 1 of 25 conditions, and the true labels are sparse in the same way — 54 genes are
truly essential in exactly 1/25, 22 in 2/25, so 76/217 (35%) are essential in ≤2 conditions. The model
commits where the answer is most concentrated.

## The withdrawn claim, and why it survived a self-check

The first version reported **0.0% constant across 20 shapes** and concluded "the model is switching."

`pattern_distribution` tested for constancy with `if p in ("....", "EEEE")` — **hardcoded to four
characters**. On a 25-condition panel, all-dispensable is 25 dots and all-essential is 25 E's; neither
matched, so the count came out zero. The true figure is **184/217 = 84.8%** (145 all-dispensable + 39
all-essential).

This is painful for a specific reason: **earlier the same day I fixed a different hardcoded-4 bug in this
same function** — `pattern_distribution` was iterating the four media names — and that fix created false
confidence that the function was clean. Generalising the condition *keys* was not enough; the constant
*test* was hardcoded too. One fix in a function is not a clearance of that function.

The earlier bug was caught by an impossible self-contradiction in the run's own output. This one was not,
because the artifact recorded only the *predicted* pattern distribution — the true distribution, which
would have shown the mismatch immediately, was never persisted. It is now.

**Blast radius: contained.** At 4 conditions the old literal test and the correct one agree, so the
published 4-media figure (62/68 = 94.0%) needs no revision. Pinned by
`test_the_four_media_constant_numbers_are_UNAFFECTED_by_the_fix`.

## What this substrate is — and is not

It is **not a re-measurement** of the 4-media result. Several things differ at once:

| | 4 media | 25 carbon sources |
|---|---|---|
| label modality | binary growth/no-growth of individual deletion strains | pooled RB-TnSeq competitive fitness, thresholded at `fit < -2` |
| oxygen | aerobic + anaerobic | **aerobic only** |
| strain | MG1655-based model vs Orth screen | Keio / BW25113 |
| reproduction gate | the paper's own iJO1366 FBA columns | **none** |

So the accuracy improvement **cannot be attributed to condition breadth alone** — label modality,
aerobicity, condition distribution and strain context all moved together. The defensible claim is
substrate-bounded: *on a larger aerobic carbon-source thresholded-fitness substrate, the model is
measurably more accurate and somewhat less collapsed.*

## Honest limits

- **Label sensitivity: RESOLVED 2026-08-13 — the headline is not a lucky cutoff.**
  *(`wiki/fba_label_threshold_sweep_2026-08-13.md`; supersedes the "no threshold sweep was run" gap.)*
  Sweeping the fitness bar over `{-1.0, -1.5, -2.0, -2.5, -3.0}` × a per-cell confidence bar over
  `{none, |t|>=2, 3, 4}` gives **positive lift over the constant null at all 20 sound settings**, rising
  monotonically as the label bar tightens (+0.0396 at `fit<-1.0` -> **+0.0745 shipped** -> +0.1413 at
  `fit<-3.0`). The constant fraction stays 85-87% throughout. The shipped `< -2` sits mid-range, not at a
  maximum.
  A 21st-through-35th setting used a first-draft confidence filter requiring `|t| >= bar` in EVERY
  condition; **all 15 collapsed to zero commitments** because that rule is ANTI-SELECTIVE for switchers
  (a switcher is confidently essential in one condition and confidently NEUTRAL in the rest). That is the
  instrument failing, not the claim; `per_cell` is now the default mode.
- **Exact-set is a weak headline metric here** given how concentrated the true patterns are. The commit-rate
  decomposition above is the more honest summary.
- **The null moved too** (0.5588 → 0.6623), so "lift" is being compared across different base rates.
- **The 39 non-optimal solves ARE the commitments — and they are correct.** *(resolved 2026-08-13; see
  `wiki/fba_infeasibility_finding_2026-08-13.md`)* All five deletion scripts now audit solver status per
  cell. The audit showed a total concentration: **23 of 23 exact-set matches, and 32 of 33 committed
  genes, have essential calls that are entirely non-optimal solves** (25 genes × 1 cell + 7 × 2 = 39, the
  whole suspect set), while the 184 predicted-constant genes touch none. That looked damning. It is not:
  all 39 re-solve `infeasible` deterministically, each is the canonical catabolic gene for exactly that
  carbon source (galT/galE/galK on galactose, malEFGK/malQ on maltose, xylA on xylose, sdhABCD on
  succinate, kgtP on α-ketoglutarate …), and **38 of 39 are experimentally essential in that condition**.
  With `ATPM lower_bound = 6.86`, deleting the sole catabolic route leaves no feasible flux distribution
  at all — infeasibility is the model's signature for *this carbon source cannot be used*. The
  NaN-to-essential coding is correct.
  **So the commit claim now has a mechanism, not just a count** — and a sharper caveat: those
  commitments are structurally easy, because sole-route catabolism is the one place a stoichiometric
  model with a maintenance floor cannot help but be right. It says nothing about the 184 constant genes,
  where the conditional deficit still lives.
- **The 184 constant genes are now diagnosed — the deficit is a MODEL problem, not a readout one.**
  *(2026-08-13, `wiki/fba_constant_gene_diagnostic_2026-08-13.md`.)* Of 1,083 missed essential cells,
  **833 (76.9%) are `flat`** — deleting the gene changed predicted growth by *nothing* (median missed
  ratio **1.0**); only **10.7% are readout-recoverable** (23% even counting the generous `slight`
  bucket). 110 of the 145 all-dispensable genes are flat in EVERY condition where truth says essential,
  and only **29% of those are isozyme-redundant by GPR** — so the obvious mechanism explains less than a
  third. The 39 over-callers behave oppositely: 315 wrong cells, **0 from infeasibility**. Infeasibility
  only ever produces CORRECT calls here; both error directions live in the finite-ratio regime.
  **Consequence:** threshold retuning (`continuous_readout`'s oracle, `deployable_threshold`) addresses
  at most ~11% of the deficit and should stop being proposed as the next lever; the redundancy that
  pFBA-restriction attacks is the real target.
- **True patterns are far more diverse than predicted ones**: 141 distinct true shapes among 217 genes,
  against 20 predicted shapes. That asymmetry, not the aggregate agreement, is the clearest statement of
  the model's conditional resolution.

## Reproduce

```bash
uv run python scripts/fba_conditional_carbon_validate.py
```

Needs `feba.db` (7.4 GB, not committed — figshare `10.6084/m9.figshare.25236931`, CC BY 4.0).
Sidecar: `wiki/fba_conditional_carbon_2026-08-12.json` (now carries both predicted **and** experimental
pattern distributions). Tests: `tests/test_fba_conditional_essentiality.py`.
