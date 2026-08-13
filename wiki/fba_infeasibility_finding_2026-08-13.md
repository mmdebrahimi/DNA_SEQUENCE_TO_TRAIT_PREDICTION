# An `infeasible` FBA solve is not a bug — it is how iML1515 says "essential here" (2026-08-13)

> **This memo records a remediation that reversed its own premise.** A plan was written, reviewed and
> executed on the hypothesis that a published result was a solver artifact. The hypothesis was wrong,
> and proving it wrong produced a better-grounded result than the one it set out to retract.

## The suspicion

Every FBA deletion script in this repo codes a NaN growth as *essential*:

```python
scripts/fba_regulatory_conditional_test.py:157   d[gid] = 0.0 if g != g else g / wt   # 0.0 <= 0.01 → essential
scripts/fba_gapfill_carbon_recheck.py:82         d[gid] = (g != g) or (g < FRAC * wt)
```

cobrapy returns NaN exactly when a solve is non-optimal. So the scripts appeared to code **"the solver
failed" as "the gene is required"** — silently, since none of them read the `status` column.

Two published claims sat directly on top of that:

- the **pFBA regulatory lift** (per-cell 0.5709 → 0.6157, ~5× the constant-null lift, 0/200 rate-matched
  draws reaching it), measured in an arm that forces off ~69% of gene-associated reactions before
  deleting anything — the condition most likely in this codebase to drive LPs infeasible;
- the **carbon-panel commit claim** ("FBA commits rarely, and is accurate when it commits").

The result's shape fit the suspicion uncomfortably well: the regulatory arm's gain was **recall-only**,
with precision falling 0.625 → 0.583 and threshold-free AUROC falling 0.612 → 0.576. That is what you
would see if infeasible solves were being converted into "essential".

## What the audit found — and why it looked damning

A per-cell solver audit (`dna_decode/fba/solver_audit.py`) was added to all five deletion scripts.

| arm | suspect cells | of | fraction |
|---|---|---|---|
| pFBA-restricted, 4 media | **74** | 268 | **27.6%** |
| carbon panel, 25 sources | 39 | 5,425 | 0.7% |

All 74 restricted-arm cells were `infeasible`. Under abstention the regulatory arm's TP collapsed
**56 → 12**, and its per-cell score (0.6237) no longer beat a rate-matched null recomputed on the
abstained denominator (max 0.6546, empirical p = 0.125). The pre-committed rule fired:
`REGULATORY_LIFT_IS_A_SOLVER_ARTIFACT`.

On the carbon panel the concentration was even starker — and the arithmetic is exact:

- **23 of 23** exact-set matches have essential calls that are *all* suspect cells
- **32 of 33** committed genes likewise
- 184 predicted-constant genes touch **zero** suspect cells
- 25 genes × 1 cell + 7 genes × 2 cells = **39** — the entire suspect set, with nothing left over

Read naively: every varying prediction the model makes is the LP failing.

## The check that reversed it

Before writing that up, the 39 cells were re-solved one gene at a time at `processes=1`
(`scripts/fba_infeasibility_probe.py`). Two things came back.

**1. They are deterministic.** 39 of 39 re-solve `infeasible`. Zero returned real growth. A numerical or
threading artifact would not reproduce identically.

**2. They are the right genes.** Every single one is the canonical catabolic gene for exactly the carbon
source it fails on:

| carbon source | genes | pathway |
|---|---|---|
| D-Galactose | galT, galE, galK | Leloir |
| D-Maltose | malE, malF, malG, malK, malQ | ABC transport + amylomaltase |
| D-Mannose | manX, manY, manZ, manA | PTS + P-isomerase |
| D-Glucosamine | manXYZ, nagB | PTS + deaminase |
| N-Acetyl-D-Glucosamine | nagA, nagB | |
| D-Galacturonate / D-Glucuronate | uxaC, uxaA, uxuA, kdgK, eda | hexuronate / Entner–Doudoroff |
| D-Xylose | xylA | isomerase |
| D-Ribose | rbsK | ribokinase |
| D-Sorbitol | srlD | sorbitol-6-P dehydrogenase |
| D-Mannitol | mtlD | mannitol-1-P dehydrogenase |
| L-Fucose | fucK | fuculokinase |
| Glycolate | glcD, glcE, glcF | glycolate oxidase |
| Succinate | sdhA, sdhB, sdhC, sdhD | succinate dehydrogenase |
| α-Ketoglutarate | kgtP | α-KG transporter |

**38 of the 39 are experimentally essential in that very condition** under the RB-TnSeq labels (97%).

## The mechanism

iML1515 carries a hard maintenance floor: `ATPM lower_bound = 6.86`. Delete the only route to catabolise
the **sole** carbon source and there is no flux distribution that meets maintenance at all — so the LP is
genuinely **infeasible**, not feasible-with-zero-growth.

Ordinary lethality looks different. On glucose, 29 of 150 gene knockouts returned `growth ≈ 0` with
`status == optimal`, and **zero** returned infeasible. Infeasibility is the specific signature of *this
carbon source cannot be used at all*, which is exactly the condition-specific essentiality being measured.

So the NaN-to-essential coding is **correct**. It was never a bug.

## What this does to the claims

| claim | status |
|---|---|
| pFBA regulatory lift (0.6157, p < 0.005) | **STANDS.** Abstention removes true positives; its 0.6237 is a biased lower bound, not a cleaner number. |
| Carbon commit claim ("accurate when it commits") | **STANDS, and is now mechanistically explained** rather than merely observed. |
| Gap-filling does not help | **STANDS, cleanly.** 154 flips of 5,425, only **4** involving a suspect solve. |
| Carbon headline numbers | **Reproduced exactly**: exact-set 23/217, per-cell 0.7368, constant 184/217 = 84.8%. |

The commit claim is materially stronger than when it was written. "The model commits for 15.2% of genes
and is right 70% of the time where it commits" was an observation about counts. It is now a statement
with a named mechanism and a checkable gene list: **the model commits precisely where a catabolic route
is the sole route, and it names the right enzyme.**

The caveat sharpens rather than disappears. Those commitments are *structurally easy* — sole-route
catabolism on a single carbon source is the one place a stoichiometric model with a maintenance floor
cannot help but be right. It says nothing about the 184 genes where the model stays constant, which is
where the real conditional deficit still lives.

## The pre-committed rule was right to fire, and its label was wrong

The verdict rule was authored before the run, per this project's verdict-vs-budget discipline, and the
discipline worked as intended: it fired without me getting to choose the answer after seeing the number.

But **pre-committing a rule does not pre-commit its premise.** The rule encoded an unexamined assumption
— that abstaining on non-optimal solves is the conservative move. That assumption was false, so a
correctly-fired rule produced a false label. The artifact now records all three: what the rule said, what
the probe found, and that the premise was falsified.

The generalisable lesson: **a pre-registered decision rule protects against motivated reasoning, not
against a wrong model of the mechanism.** When a rule fires against a published result, the cheap
confirmation — 39 re-solves, a gene-name lookup — comes *before* the retraction, not after.

## Reproduce

```bash
uv run python scripts/fba_conditional_carbon_validate.py    # writes solver_audit + commit_strata
uv run python scripts/fba_infeasibility_probe.py            # the decisive re-solve
uv run python scripts/fba_regulatory_conditional_test.py    # consumes the probe's verdict
```

Sidecars: `wiki/fba_infeasibility_probe_2026-08-13.json`,
`wiki/fba_conditional_carbon_2026-08-13.json`,
`wiki/fba_regulatory_conditional_recheck_2026-08-13.json`.
Tests: `tests/test_fba_solver_audit.py` (23).
