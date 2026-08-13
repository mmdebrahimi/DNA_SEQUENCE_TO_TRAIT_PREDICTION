# The regulatory restriction is not a lever — it makes the conditional metric WORSE at scale (2026-08-13)

> The last surviving candidate direction for the FBA conditional deficit, tested on the wide panel and
> **rejected**. Also the run where my own verdict function fired a false `CONFIRMED` and was caught before
> publication.

## What was tested

pFBA restriction — in each condition solve for the cheapest way to grow, then force off every
gene-associated reaction carrying no flux — was the only lever still standing:

- **gap-filling** — ruled out (154 flips of 5,425, exact-set −1)
- **threshold retuning** — ruled out quantitatively (≤11% of the deficit)
- **route-constraining** — the direction the constant-gene diagnostic points at (76.9% of missed cells
  are `flat`, i.e. redundancy)

It had scored well on 4 media (per-cell 0.5709 → 0.6157, p = 0.0 vs the rate-matched null), then dropped
to **p = 0.06** once a margin-preserving null was built. This re-runs the identical intervention on
**25 carbon sources × 217 genes = 5,425 cells**, ~20× the data.

## Result — it degrades every axis

| | per-cell | exact-set | mean MCC | TP | FP | precision |
|---|---|---|---|---|---|---|
| baseline (all routes) | **0.7368** | 23/217 | **0.3864** | 721 | 317 | **0.6946** |
| pFBA-restricted | **0.6839** | 23/217 | 0.3081 | 988 | **871** | 0.5315 |
| **delta** | **−0.0529** | **0** | −0.0783 | +267 | **+554** | −0.163 |

**Verdict: `REGULATORY_RESTRICTION_MAKES_IT_WORSE`.**

The mechanism is visible in the columns. Forcing a unique route makes far more genes look essential
(+267 TP but **+554 FP**). On the 4-media substrate the baseline was badly *under*-calling — TP 10, FP 6
on 268 cells — so any push toward more calls helped. On the wide panel the baseline already makes 1,038
calls at precision 0.69, and the restriction just buries it in false positives. **The 4-media "lift" was
the intervention compensating for an under-calling baseline, not adding condition-specificity.**

Exact-set is unchanged at 23/217 — the restriction does not reproduce a single additional switch pattern.

## The false positive my verdict function produced

The run first printed **`REGULATORY_LIFT_CONFIRMED_ON_WIDE_PANEL`**, because `verdict_for` compared the
observed value only against the **nulls** and never against the **baseline arm it was meant to improve
on**. Both null comparisons genuinely pass:

| null | mean | max | p vs observed 0.6839 |
|---|---|---|---|
| rate-matched (weak) | 0.5510 | 0.5666 | 0.0 |
| margin-preserving (strong) | 0.6750 | 0.6794 | **0.0** (0/200 reach it) |

Those are not wrong — they are answering a different question. A margin-preserving null is built from the
**restricted arm's own margins**, so beating it means *"given how many calls this arm made, they are well
placed."* That is entirely compatible with **making that many calls being a bad idea in the first place**.

On 4 media the restriction also improved the baseline, so the missing comparator was invisible. The
baseline check now precedes the null check and is pinned by three tests.

This is the same defect class as everything else caught in this cluster: **a metric that doesn't compare
against the right reference.** It is worth naming that the error kept appearing in *new* places after
each fix — condition keys, then the constant test, then the confidence filter, now the comparator.

## Against the pre-registered expectation

Registered before the run: *"it will NOT clear p<0.05."*

**Half right, and wrong about the mechanism.** It cleared the strong null decisively (p = 0.0) — my
reasoning about a small-sample ceiling was wrong. But the intervention still fails, for a reason the
forecast never considered: it is worse than doing nothing. Recording that the prediction was right for
the wrong reason, which is not the same as being right.

## Where this leaves the FBA conditional cell

All three levers are now measured and closed:

| lever | verdict | evidence |
|---|---|---|
| add reactions (gap-fill) | ruled out | 154/5,425 flips, exact-set −1 |
| retune the threshold | ruled out | ≤11% of the deficit is readout-recoverable |
| constrain routes (pFBA) | **ruled out as a METHOD** | −0.0529 per-cell, +554 FP, exact-set +0 |

**The diagnosis stands; the treatment does not exist yet.** 76.9% of missed cells are deletions that
changed nothing, and that number owes nothing to pFBA — the redundancy is real. But the one available way
of removing redundancy makes things worse, because it removes it *indiscriminately*: ~1,868 reactions
forced off per condition, chosen by cheapest-flux rather than by biology.

What that implies for a next attempt, stated as a hypothesis rather than a finding: a useful restriction
would have to be **selective** — targeting the specific redundant routes that make a gene flat — rather
than a blanket parsimony sweep. The ~78 flat genes whose redundancy is *not* explained by isozyme GPR
structure are where such a thing would have to start.

## Honest limits

- pFBA picks **one** optimal-flux solution; alternate optima of equal cost would force off different
  routes. Untested, and it could change which genes look essential.
- Forcing off ~69% of gene-associated reactions is a crude proxy for regulation. This result rejects
  **that proxy**, not "regulation matters" as a hypothesis.
- All 25 conditions are aerobic carbon sources. No oxygen axis — and the 4-media substrate that produced
  the original positive *did* carry one.

## Reproduce

```bash
uv run python scripts/fba_regulatory_carbon_test.py
```

Needs `feba.db` (7.4 GB, not committed — figshare `10.6084/m9.figshare.25236931`, CC BY 4.0).
Sidecar: `wiki/fba_regulatory_carbon_test_2026-08-13.json` (carries both the corrected verdict and the
wrong first one, under `verdict_first_version_WRONG`).
Tests: `tests/test_fba_regulatory_carbon.py` (8).
