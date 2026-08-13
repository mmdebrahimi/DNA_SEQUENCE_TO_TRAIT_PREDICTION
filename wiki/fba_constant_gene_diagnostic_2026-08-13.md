# Why iML1515 won't switch: 77% of its misses are deletions that changed *nothing* (2026-08-13)

> The conditional deficit is a **model** problem, not a readout problem. A perfect threshold retune
> recovers at most ~11% of it. Measured, not argued.

## The question

The model commits to a varying essentiality pattern for only 33 of 217 conditionally-essential genes, and
those commitments are now fully explained — sole-route catabolism detected through ATPM-maintenance
infeasibility (`wiki/fba_infeasibility_finding_2026-08-13.md`). That left **184 constant genes** where the
whole conditional deficit lives, and nothing had ever asked why.

The question separates two diagnoses with **opposite fixes**:

- **MODEL problem** — the knockout growth ratio is flat at 1.0 in the conditions where the gene is truly
  essential. The model has an alternative route the real cell does not use. No readout change recovers
  this; it needs regulation or gene-content curation.
- **READOUT problem** — the ratio is materially depressed but sits above the 1% cutoff. The model *sees*
  the defect and the binary threshold discards it. A graded metric recovers it for free.

## Result

| stratum | n |
|---|---|
| predicted **all-dispensable** | **145** |
| predicted **all-essential** | 39 |
| commits to a varying pattern | 33 |

**1,083 missed essential cells**, split by what the model actually saw:

| cause | cells | share | fixable by a better readout? |
|---|---|---|---|
| **`flat`** — ratio ≈ 1.0, the deletion changed **nothing** | **833** | **76.9%** | **No** |
| `slight` — ratio 0.90–1.0 | 134 | 12.4% | Barely |
| `material` — ratio 0.10–0.90, a real defect discarded by the cutoff | 116 | 10.7% | **Yes** |
| `near_threshold` — ratio just above 1% | 0 | 0% | Yes |

**Median missed ratio: 1.0.** More than three-quarters of the time the model misses a conditionally
essential gene, deleting that gene changes predicted growth by *literally nothing*.

**110 of the 145** all-dispensable genes are flat in **every** condition where the gene is truly
essential. The model cannot call those essential under any threshold.

**Verdict: `DEFICIT_IS_A_MODEL_PROBLEM_NOT_A_READOUT_ONE`** — readout-recoverable share **10.7%**.

**Robustness of that verdict:** even counting the entire `slight` bucket as recoverable (a deliberately
generous reading — a 2–10% growth dip is not really a missed essentiality call) gives 250/1,083 = **23%**,
still far below the 50% bar for "mostly a readout problem". The conclusion does not turn on where the
`slight`/`material` boundary sits.

## The mechanistic cross-check, and the part that surprised me

The obvious explanation for a flat ratio is **isozyme redundancy** — another gene covers the same
reaction via an `or` in the gene-reaction rule, so a single deletion cannot block anything.

Only **32 of the 110 flat genes (29%)** are isozyme-redundant by GPR structure.

So the obvious mechanism accounts for **less than a third**. The remaining ~78 genes are flat for reasons
that live above the GPR: alternative *pathways* rather than isozymes on the same reaction, transport
redundancy, or a biomass objective that never demands what the gene produces. That is a harder and more
interesting class than "the model has duplicate enzymes", and it is not addressable by GPR curation alone.

## The over-calls behave oppositely

The 39 predicted-essential-everywhere genes over-call **315 cells** where truth says dispensable — and the
cause split is **`{infeasible: 0, sub_threshold: 315}`**. Not one over-call comes from infeasibility.

That is a clean asymmetry against the commitments, where infeasibility carried 54–66% of the calls
(`wiki/fba_label_threshold_sweep_2026-08-13.md`). **Infeasibility only ever produces correct calls in this
data; the errors come entirely from the finite-ratio regime.** The boolean signal is trustworthy; the
graded one is where the model goes wrong in both directions.

## What this rules in and out

**Ruled out — threshold retuning as a route to conditional accuracy.** `continuous_readout`'s oracle
threshold and `deployable_threshold`'s retune address at most ~11% of the deficit. That line is
quantitatively dead for *this* metric, and should stop being proposed as the next lever.

**Ruled out (already, independently) — gap-filling.** Adding reactions cannot help a model whose problem
is that deletions do nothing; more routes make flatness worse, which is exactly what the gap-fill arms
measured (154 flips of 5,425, exact-set −1).

**Pointed at — constraining which routes are available.** 77% flatness *is* the redundancy the pFBA
restriction attacks. Two independently-produced results point the same way: **the deficit is too many
routes, not too few, and not a badly-chosen cutoff.**

**But the pFBA result itself got weaker the same day, and the direction should not borrow its old
strength.** Its named weakness — a null that sampled cells independently while real gene patterns are
correlated — was closed by building the margin-preserving null (`dna_decode/fba/nulls.py`, Curveball swap
randomization preserving every gene's and every condition's essential-call count):

| null | mean | max | p vs observed 0.6157 |
|---|---|---|---|
| rate-matched (only the grand total fixed) | 0.5172 | 0.5933 | **0.0** — 0/200 draws reach it |
| **margin-preserving (both margins fixed)** | **0.5946** | **0.6157** | **0.06** — 12/200 reach it |

So the published "~5× lift over null, p < 0.005" was **substantially an artifact of a null that was too
easy to beat**. Against the strong null the pFBA arm sits at roughly the 94th percentile and **does not
clear p < 0.05**.

The honest position: **this diagnostic establishes the direction on its own evidence** (77% of misses are
deletions that changed nothing — that number owes nothing to pFBA). What it does **not** do is validate
pFBA-restriction as the *method*; that intervention is now suggestive-but-not-significant, on 4 media,
with a crude 69%-of-reactions-off proxy. Direction: supported. Method: open.

## Honest limits

- `flat` uses an exact-1.0 tolerance of 1e-6. Degenerate LP optima shift mid-range ratios between runs,
  but a ratio *at* 1.0 is stable — the deletion genuinely changed nothing.
- Isozyme redundancy is read from GPR `or` structure. That is a **model** property: it says the model
  cannot call the gene essential, not that the real cell is redundant.
- Calling `material` cells "readout-recoverable" is an **upper bound**. Recovering them needs a threshold
  that does not simultaneously destroy precision elsewhere — untested here, and the 4-media oracle retune
  bought MCC at no precision cost only on that much smaller substrate.
- All 25 conditions are aerobic carbon sources. No oxygen axis.

## Reproduce

```bash
uv run python scripts/fba_constant_gene_diagnostic.py
```

Needs `feba.db` (7.4 GB, not committed — figshare `10.6084/m9.figshare.25236931`, CC BY 4.0).
Sidecar: `wiki/fba_constant_gene_diagnostic_2026-08-13.json`.
Tests: `tests/test_fba_constant_gene_diagnostic.py` (8).
