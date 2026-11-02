# The eukaryotic arm: the tile-chain ceiling was the SUBSTRATE's, not the method's

<!-- doc-status: active -->
<!-- measured-on: 2026-11-01 -->

Graded by `wiki/eukaryotic_arm_acceptance_bar.json`, frozen before any number below existed.

## DECISIONS NEEDED

**One, and it is structural rather than about this result.** This arm's metric is a HUMAN
lentiMPRA Spearman sitting in a criterion whose other ten metrics are BACTERIAL, and
`ceiling_verdict` picks `best` as a max over metrics. It is harmless today — 0.4853 is below the
bacterial best 0.5076, so `best` is unmoved — but a future human arm scoring above that would
silently become `best` on a bacterial criterion. The arm's real contribution is the CONTRAST,
not the level. Moving it to its own criterion is a scope change to a criterion declared
2026-10-09, so it is yours; a guard
(`test_the_human_arms_metric_is_not_read_as_best`) now fires if it ever stops being harmless.

**The ceiling criterion does NOT close, and that is correct rather than an oversight.** This
was the last `planned` arm, but under R6 the positional (non-mean-pooled) PromoterAtlas arm
that G-A candidate 8 names was declared at the moment this one ran — because an arm declared
only once its cheap sibling succeeds is an arm that evaporates. So `ceiling_verdict` stays
`blocked:untried-arm` on a genuine untried candidate: 11 measured, 1 planned, `best` unmoved.

## CORRECTIONS

**The standing handoff called this arm a funding decision. It was neither a funding nor an
access problem.** The substrate was on an attached drive, MD5-verified, with a working loader,
and the arm ran locally in 38 minutes of wall clock. A prior session the same day had already
established it was *not* drive-blocked and judged it cost-gated via a Kaggle fetch; the drive
being attached makes even that unnecessary.

## The result

| | GC-only ridge | conv encoder | gain |
|---|---|---|---|
| human HepG2 (median over 10 folds) | 0.2362 | 0.4853 | **+0.2491** |
| bacterial tiles (committed arms) | 0.3112 | 0.3234 | **+0.0122** |

Paired bootstrap over the 10 published folds: mean +0.2471, 95% CI [+0.2399, +0.2536], against a
bacterial anchor of +0.0122.

**Verdict: `SUBSTRATE_LIMITED`** — the CI for the human gain lies entirely ABOVE the bacterial gain 0.0122

## What the gain is actually made of

A gain-over-baseline contrast only crosses a substrate boundary if the baseline means the same
thing on both sides, and here it does not — GC is weaker on human regulatory elements than on
bacterial tiles. Separating the two sources:

- the conv encoder moves **+0.1619** (0.3234 → 0.4853)
- the GC baseline moves **+0.0750** in the direction that widens the gap (0.3112 → 0.2362)
- together they account for the **+0.2369** difference in gain

So 68% of the widened gain is the encoder doing better and 32% is GC doing worse. **Both point the
same way**, which is what makes the reading hold: on this substrate a learned encoder finds
signal that a single composition feature cannot, and on bacterial tiles it essentially could
not. Had the gain widened *only* because GC collapsed, the arm would say nothing about the
encoder.

## What this does NOT establish

**The published incumbents are NOT a comparator, and that was pre-committed.** Table S6 gives
per-fold values for six methods on these exact folds and never states the metric; the sheet was
checked directly, including the rows the parser skips (`unparsed_rows()` returns 0 for
`Benchmarking results`, verified rather than trusted). What can be said is weaker and still
useful: whatever those values measure, six independent methods cluster at roughly 0.72–0.87 and
order the way a performance metric would. That is an inference about an unnamed quantity.

- THE BIGGEST ONE: this is ONE substrate change, not a substrate sweep. A human lentiMPRA panel differs from bacterial tiles in assay, in element length (230 bp vs 150 bp), in kingdom and in how the library was constructed. The arm cannot say WHICH of those differences carries any effect it finds.
- The bacterial anchor 0.0122 is treated as a FIXED POINT; its own uncertainty is not propagated into the CI, so a verdict that lands near the boundary is weaker than the predicate makes it look.
- The conv encoder's hyperparameters were chosen on the BACTERIAL rung. Reusing them on human sequence is the honest apples-to-apples choice for a method-class comparison, but it is NOT a tuned human model and must never be reported as one.
- A constructed panel is the regime this project already knows learned methods work in (constructed -> molecular, Spearman 0.761 on TEM-1). So a SUBSTRATE_LIMITED result CORROBORATES an existing finding rather than discovering something; its value is that it bounds the tile-chain ceiling, not that it is news about constructed panels.
- Nothing here licenses a claim about eukaryotic g->p generally. The project's organism-level eukaryotic embedding arm is a CLOSED NEGATIVE on natural populations (Arabidopsis, de-confounded) and this arm does not reopen it: constructed variation is a different regime, which is the whole reason this arm is scoped to a constructed panel.

## Reproduce

```bash
uv run python -m scripts.eukaryotic_arm_probe          # validity first: join, folds, degeneracy
uv run python scripts/freeze_eukaryotic_arm_bar.py     # refuses if the anchors drift
uv run python -m scripts.eukaryotic_arm_run --cell HepG2
uv run python -m scripts.eukaryotic_arm_memo --cell HepG2
```
