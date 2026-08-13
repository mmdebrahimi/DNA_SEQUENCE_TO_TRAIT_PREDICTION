# Does label-blind gap-filling move the conditional metric? No. (2026-08-12)

> ## ⚠ PARTIALLY SUPERSEDED THE SAME DAY — read this first
>
> Everything below was measured on the Orth **4-media** substrate. A re-test on **25 carbon sources**
> (`scripts/fba_gapfill_carbon_recheck.py`) splits this document's conclusion in two.
>
> **Second correction (same day):** an earlier version of this banner said the 4 media had "no room to
> move". That was based on a "0% constant" figure which was itself a bug — a constant-pattern test
> hardcoded to four characters. Corrected, the model is constant for **84.8%** of switching genes on the
> 25 sources against **94.0%** on the 4 media. The wider panel gives *somewhat* more room, not a
> transformation. See `wiki/fba_conditional_carbon_2026-08-12.md`.
>
> | claim below | status on the wider panel |
> |---|---|
> | **Mechanism:** "adding reactions cannot change a single call" (0 flips at every dose) | **FALSIFIED** — the maximal arm changes **154 of 5,425** calls |
> | **Practical:** "gap-filling does not help" | **CONFIRMED, on much better evidence** — exact-set goes **23 → 22** (down 1) and per-cell moves **+0.0003** |
>
> **Third note (2026-08-13): the negative is now solver-audited and survives cleanly.** Of the 154
> flips on the 25-source panel, only **4** involve a cell whose solve was non-optimal in either arm — so
> the "changed calls are noise, not signal" reading is not itself a solver artifact. Separately, the
> NaN-as-essential coding these scripts use was suspected of manufacturing results and was tested and
> **refuted**: it is deterministic ATPM-maintenance infeasibility, i.e. genuine essentiality. See
> `wiki/fba_infeasibility_finding_2026-08-13.md`.
>
> So the *headline answer is unchanged and now better supported* — gap-filling still does not improve
> conditional essentiality — but the **"it cannot move anything at all" framing was substrate-bounded**:
> on a wider panel it moves 154 calls. The changed calls are noise, not signal.
>
> The mechanism paragraph below ("adding more routes cannot create a medium-specific dependency") is the
> part to distrust: it is a clean argument that the 4-media substrate could not put under strain. Note the
> wider panel differs in label modality, oxygen coverage and strain as well as breadth, so the difference
> is **not** attributable to condition count alone.

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

**Ruled IN — tested, not asserted:** regulatory constraints. See below.

Also untested: a different donor. The maximal arm exhausts *this* donor, but another reconstruction could
in principle carry a reaction that matters.

## What *does* move it — the regulatory constraint

Rather than leave "regulation is the leading candidate" as an assertion, it was tested the same day
(`scripts/fba_regulatory_conditional_test.py`). The intervention is a **parsimonious (pFBA) restriction**:
in each medium, solve for the cheapest way to grow, then force off every gene-associated reaction carrying
no flux in that solution. The surviving route differs *by medium* — the condition-specificity the flat
ratios lack. pFBA never sees the labels, so it is label-blind in the same sense as the gap-fill arms.

| | exact-set | per-cell | mean MCC | AUROC | TP | FP | precision |
|---|---|---|---|---|---|---|---|
| baseline | 3/67 | 0.5709 | 0.0611 | 0.612 | 10 | 6 | 0.625 |
| **pFBA-restricted** | **5/67** | **0.6157** | **0.217** | 0.576 | **56** | 40 | 0.583 |

**Lift over the constant null goes from +0.012 to +0.057 — nearly 5×**, and mean per-condition MCC more
than triples.

**The control is what makes that meaningful.** Forcing a unique route makes far more genes look essential,
so the gain could be nothing but a better-matched base rate. Scoring a **rate-matched random** predictor —
calling the same 96 cells essential, at random, 200 times — gives mean **0.5172** (sd 0.028, max 0.5933).
**Zero of 200 draws reach the observed 0.6157** (empirical p < 0.005). The gain is genuine
condition-specific signal.

**But read it as a direction, not a method:**

- It is a **crude** proxy for regulation — it forces off ~**1,865 of 2,712** gene-associated reactions
  (~69% of the model).
- The gain is **recall only**. Precision does not improve (0.625 → 0.583), and the threshold-free **AUROC
  gets worse** (0.612 → 0.576) — the continuous ranking degrades even as the binary calls improve.
- It still reproduces only **5 of 67** switches. Most of the deficit remains.
- pFBA picks *one* optimal-flux solution; alternate optima of equal cost would force off a different route.

So the honest summary of the pair: **adding reactions cannot move this metric at all; constraining which
reactions are available moves it about 5× over null.** That is a clear signpost for where the FBA cell's
condition-awareness has to come from — and it points at regulation, not at biochemistry coverage.

## Reproduce

```bash
uv run python scripts/fba_gapfill_conditional_test.py
```

Sidecar: `wiki/fba_gapfill_conditional_test_2026-08-12.json`. Tests: `tests/test_fba_gapfill_conditional.py`.
