# F3 — phylogeny de-confounding: the instrument is BUILT and NOT CALIBRATED

**Verdict: `blocked:positive_control_not_fold_compatible`.** A blocked verdict here is a
pre-registered outcome of the plan, not a partial result. **No calibration metrics were
emitted**, and none should be quoted from this run.

**Scope of the claim.** A grouped-CV path now exists and is tested
(`eval/genomic_prediction.py`, `groups=`). What does **not** exist is evidence that it
*responds differentially to real population structure*. Until that is measured, a grouped
number from it is not evidence of de-confounding, and `eval/genomic_prediction`'s module
docstring and `eval/regime.py` both say so.

## What was built

| unit | status |
|---|---|
| `eval/genomic_prediction.py` grouped-fold seam (`groups=`) | shipped, 21 tests |
| `eval/panel_structure.py` IBS + guard-selected threshold + `pc1_variance_fraction` | shipped, 33 tests |
| `scripts/f3_substrate_census.py` integrity + structure census | shipped, 22 tests, artifact |
| `scripts/f3_deconfound_calibration.py` double control + 8 frozen verdicts | shipped, 23 tests, artifact |

`groups=None` is a **bit-for-bit no-op**, proven against the pre-change module extracted from
`git HEAD` rather than against hardcoded goldens: ridge over 3 seeds, gbm, and the NaN-filter
path all return identical values (anchor `predictive_r = 0.4723343055050292`).

The permutation null consumes the **same `folds` object** as the measured run, so a grouped run
has a grouped null **by construction** rather than via a second code path. That is asserted
directly by capturing `id(folds)` across the measured run and every null draw, for both the
ridge and the model path — a null that shares the metric's fold structure is the only kind
that bounds the metric it is compared against.

## Why it is blocked

All four cached panels are **intact** (bloom 1008×11623, arabmagic 703×1260, bxd 198×7320,
dgrp 205×12365; full parse, not a presence check). The blocker is not missing data.

**The one natural-population panel cannot support grouped CV at all.** A grouping can clear
the panel-level 60% dominant-cluster bound and still be unusable for *k*-fold CV, where a group
above *n/k* makes one fold exceed its share. Across the whole percentile ladder, dgrp has no
rung that is **both** fold-compatible (largest ≤ 1/5) **and** not mostly singletons:

| rung | largest cluster | singletons | usable at k=5 |
|---|---|---|---|
| p1 – p30 | 0.024 – 0.195 | **0.658 – 0.836** | no — mostly singletons |
| p40 – p50 | 0.302 – 0.342 | 0.568 – 0.600 | no — both bounds fail |
| p60 | **0.522** | 0.481 | no — above the 1/5 fold share |

bloom, arabmagic and bxd each have such a rung. **Unblocking is a user decision**, because it
means acquiring a substrate: a natural-population panel whose relatedness is not concentrated
in one clump — e.g. the 1001-Genomes genotype matrix, which is *not* on disk (see below).

## Two premises the plan was written on, both corrected by measurement

**1. There is no cached natural-population Arabidopsis genotype matrix.** `data/arabidopsis/`
holds phenotypes only (`FT10_pheno_261.csv`, `FT16_pheno_262.csv`, `flowering_1001g/pheno_*`);
the only cached Arabidopsis genotypes are `arabmagic`, a **constructed** multi-parent panel.
The plan's motivating premise came from reading a directory listing, which is the failure its
own Step 1 exists to prevent.

**2. `pc1_variance_fraction` cannot select the positive control.** The plan pre-registered
"the intact panel with the highest measured `pc1_variance_fraction`". Measured, that rule is
**self-contradictory on this panel set**:

| ordering | result |
|---|---|
| raw `pc1` | bxd 0.0354 > bloom 0.0261 > dgrp 0.0197 > arabmagic 0.0173 |
| excess over its own within-marker permutation null, matched marker count | **bloom 7.91** > bxd 3.84 > arabmagic 3.52 > **dgrp 2.00** |

Under **both** normalisations the one natural population ranks last or third, and under excess
the top panel is **bloom — the bi-parental cross that must be the negative control**.

The mechanism is **structural, not statistical**: `pc1` measures total genotypic correlation,
which is linkage disequilibrium **plus** population structure, and a bi-parental cross
*maximises LD by construction* while carrying no population structure at all.

**So the control pair is selected BY CONSTRUCTION**, from each panel's `design` field, declared
in the census `PANELS` table before any number was read and independent of every measured
statistic. `negative = bloom` (bi-parental), `positive = dgrp` (one natural *Drosophila*
population). A test pins that the selector does not read `pc1`.

**Cluster shape is independent corroboration and agrees.** dgrp sits at 66–84% singletons from
p1 to p50 and passes only at p60; all three crosses pass from p1–p2 with low singletons. That
is the cryptic-relatedness signature — a few tight clumps among unrelated lines — reached by a
different route than `design`, so the agreement is corroboration rather than circularity.

## The criterion, and why it is differential

A natural panel may carry structure **and** real functional signal, and a *correct*
de-confounder removes the former while leaving the latter. An absolute "the positive control
falls to its own null" bar would therefore **condemn a working instrument**. What a
structure-de-confounder must do is respond **more** where structure exists than where it is
randomised by construction — so the criterion is the **contrast** between the two controls,
which is what makes them a pair rather than two runs. It asserts no magnitude.

Both adverse verdicts are genuinely reachable, which is what separates this from an experiment
with a foregone outcome: `INVALID_OVER_AGGRESSIVE` (the grouping destroys signal on the cross,
where the partition is near-arbitrary) and `INVALID_INSENSITIVE` (IBS clustering fails to
recover structure that is measured to exist). All 8 verdict states have a test that reaches
them.

## A measured property of the threshold rule, worth carrying

A **percentile threshold always finds something, and the guards cannot see it.** On an i.i.d.
binary panel with no structure, a 10th-percentile cut still produced 19 well-shaped clusters
over 50 individuals (largest 14%), clearing both guards — ~p% of pairs fall below any quantile
by construction. **A clean cluster shape is not evidence of structure.** Two consequences: read
`pc1` (no free parameter) for existence, and note that this is a *feature* for the negative
control, since an arbitrary grouping on a cross is exactly what must be exercised there.

A single fixed percentile is also wrong in the other direction: at p10 a genuinely structured
two-block panel (`pc1` 0.894) was **refused** for 84% singletons. So the threshold is swept
over a ladder of the panel's own distance percentiles and selected **by the guards** — the fix
CLAUDE.md prescribes for the Mash defect, *"add an over-split guard, THEN extend the grid"*, in
that order. Measured: the structured panel then selects p50 and recovers **both planted blocks
exactly** (2 clusters, 30/30, zero contamination) where the fixed rung shattered them.

`select_threshold` takes the **smallest** passing rung and its bias is named rather than hidden:
smallest under-merges, sitting closer to random CV, which biases **against** detecting a
de-confounding response. A response surviving that choice is not an artifact of it; the
opposite rule would flatter the instrument.

## Honest limits

- A blocked run bounds **nothing** about the instrument's quality. It is neither validated nor
  invalidated.
- `design` is recorded **by construction**, not measured. `pc1_variance_fraction` is the
  measured statistic, and it is the one shown above to be unfit for selection.
- `pc1_excess` uses 3 permutations at 1200 subsampled markers — enough to establish an
  ordering, not a precise effect size.
- The two controls would have had **different traits**, so the between-control comparison is
  two-sample, not paired; only the within-control sign tests are paired.
- The positive control's structure is **measured** while the negative control's randomisation
  is **by construction** — a deliberate asymmetry in evidence strength.
- IBS-greedy-representative is **one** grouping. An `INVALID_INSENSITIVE` verdict would not
  distinguish the grouping's blindness from the instrument's insensitivity.
- `MIN_TRAITS_FOR_SIGN_TEST`, `MAX_TRAITS`, `IBS_PERCENTILE_LADDER`,
  `MAX_LARGEST_CLUSTER_FRACTION` and `MAX_SINGLETON_FRACTION` are **asserted, not derived**.
- The frozen AMR surface is byte-unchanged; the prospective lock still verifies.

## Reproduce

```
uv run python scripts/f3_substrate_census.py          # read-only, exit 0 always
uv run python scripts/f3_deconfound_calibration.py
```

Artifacts: `wiki/f3_substrate_census_2026-10-09.json`,
`wiki/f3_deconfound_calibration_2026-10-09.json`.
