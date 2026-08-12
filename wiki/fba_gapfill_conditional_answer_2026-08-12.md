# Does label-blind gap-filling move the conditional metric? No. (2026-08-12)

The open question left by `wiki/fba_conditional_essentiality_2026-08-12.md`. Answered by experiment, with
the prediction registered before the run.

> **Prediction, written first:** gap-filling should make the conditional metric **worse**. Conditional
> essentiality requires the *absence* of an alternative route in one medium; gap-filling *adds* routes, so
> it should push genes toward "dispensable everywhere" — the constant pattern the model already
> over-produces.

The prediction was wrong in an interesting way. Gap-filling doesn't make it worse. It does **nothing at
all.**

## Result

Three label-blind arms against iML1515, donor = Salmonella pan-reactome iYS1720 (enterobacterial, 1,125
reactions absent from iML1515). No arm consults the essentiality labels.

| arm | reactions added | exact-set | per-cell | deployed MCC | **binary call flips** |
|---|---|---|---|---|---|
| baseline | 0 | 3/67 | 0.5709 | 0.0918 | — |
| random, k=25 | 25 | 3/67 | 0.5709 | 0.0918 | **0** |
| random, k=100 | 100 | 3/67 | 0.5709 | 0.0918 | **0** |
| random, k=400 | 400 | 3/67 | 0.5709 | 0.0918 | **0** |
| targeted (dead-end-closing) | 28 | 3/67 | 0.5709 | 0.0918 | **0** |
| **maximal (every donor reaction)** | **1,125** | 3/67 | 0.5709 | 0.0918 | **0** |

Identical to four decimals, at every dose, across seeds, including the maximal arm.

## The model really is changing — that's what makes this interesting

This is not a no-op that failed to apply. The augmented model is substantially different:

- reactions **2,712 → 3,837**; genes **1,516 → 2,111**
- anaerobic wild-type growth **0.1575 → 1.3529** (8.6×)
- knockout ratios change for **12–17 of 67 genes** (78 of 268 gene × condition cells)

**And not one of those changes crosses the essentiality threshold.** The ratios move — mostly *upward*,
toward 1.0000 (`b0115` 0.9774 → 1.0000, `b2277` 0.8860 → 1.0000), i.e. toward *more dispensable*, which is
the direction the prediction anticipated. They just never move far enough to flip a call.

The effect also **saturates immediately**: 78 cells move at k=400 and exactly the same 78 at k=1,125. A
small handful of donor reactions accounts for the entire effect, and adding a further 700 changes nothing.

## Threshold-free, so the cutoff isn't hiding a win

| | AUROC | deployed MCC | oracle MCC |
|---|---|---|---|
| baseline | 0.645 | 0.0918 | 0.271 |
| +1,125 donor reactions | 0.654 | 0.0918 | 0.254 *(worse)* |

The AUROC difference (+0.009) is **inside the noise** — this metric carries several points of run-to-run
and invocation-shape variation from degenerate LP optima, so a delta that small is not evidence of
anything. The oracle threshold gets slightly *worse*. There is no hidden win.

## Why — and it is the same mechanism as the original finding

The conditional deficit was diagnosed as **64% of these genes having a perfectly flat knockout ratio**
across all four media. A flat ratio means the model already has alternative routes that work identically
regardless of carbon source. **Adding more routes cannot create a medium-specific dependency** — it can
only deepen the redundancy that produced the flatness.

Gap-filling answers "what biochemistry is missing?". The conditional deficit is not missing biochemistry.

## Two supporting findings from the same run

**The yeast premise does not transfer to E. coli.** Gap-adjacency among conditionally-essential genes is
**27–40%** in iML1515 versus **66%** in yeast/iMM904, and iML1515 has 7.4% dead-end metabolites against
yeast's 18.4%. iML1515 is simply a much better-curated model, so there was less to fill — which is part of
why the maximal arm changed so little.

**The reproduction-gate numbers are partly in-sample, which makes the headline stronger.** Orth et al.
state plainly that this screen was used to *fix* the model: *"By comparing model predicted growth
phenotypes to the measurements, errors in the reconstruction were found and several updates were made."*
So iJO1366 was corrected using this very dataset, and iML1515 inherits those corrections. The ~5%
conditional switch score is therefore an **optimistic ceiling measured on tuned-against data**, not a
floor.

## What this rules in and out

**Ruled out:** adding reactions. Any form of it — random, dead-end-targeted, or exhaustive.

**Not tested, and now the leading candidate:** regulatory and uptake constraints. The flat ratios say the
model keeps the same routes available in every medium. Real cells don't — they repress and induce. A
method that switches *which reactions are available* per condition (regulatory FBA, or condition-specific
expression constraints) attacks the actual mechanism, whereas gap-filling attacks a mechanism that isn't
the problem here.

Also untested: a different donor. The maximal arm exhausts *this* donor, but another reconstruction could
in principle carry a reaction that matters.

## Reproduce

```bash
uv run python scripts/fba_gapfill_conditional_test.py
```

Sidecar: `wiki/fba_gapfill_conditional_test_2026-08-12.json`. Tests: `tests/test_fba_gapfill_conditional.py`.
