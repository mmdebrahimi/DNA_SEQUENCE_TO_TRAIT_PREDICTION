# The "FBA can't switch" finding REVERSES at 25 carbon sources (2026-08-12)

Two days ago, on the Orth 2011 4-media screen, iML1515 reproduced the medium-dependent essentiality switch
for **3 of 67 genes (4.5%)** and predicted a **constant** pattern for **94%** of the genes whose
essentiality actually depends on the medium — a lift of **+0.012** over a constant-predictor null. The
conclusion recorded was that the model "is not switching at all."

Re-asked against the Fitness Browser RB-TnSeq compendium — **25 carbon sources, 217 conditionally-essential
genes**, ~125× the cell count — that conclusion does not survive.

| | 4 media (Orth 2011) | **25 carbon sources (Fitness Browser)** |
|---|---|---|
| conditionally-essential genes | 67 | **217** |
| exact-set match | 3/67 = 4.5% | **23/217 = 10.6%** |
| per-cell agreement | 0.5709 | **0.7368** |
| best constant null | 0.5588 | 0.6623 |
| **lift over null** | **+0.0121** | **+0.0745** (6×) |
| mean per-condition MCC | 0.0611 | **0.3864** (6×) |
| **genes predicted with a CONSTANT pattern** | **94.0%** | **0.0%** |
| distinct predicted shapes | 3 (truth had 12) | **20** |

**The headline correction: the model is switching.** Zero of 217 genes get a constant prediction, across
20 distinct patterns. The earlier "94% constant / not switching at all" was a real measurement of a real
model — but on a substrate too small and too narrow to see the behaviour.

## Why the 4-media result was pessimistic

Not because it was wrong, but because of what those four media *are*. Two of the four are glucose
(differing only in oxygen), and the other two are gluconeogenic (lactate, succinate). A model whose
conditional resolution runs through **carbon-catabolic** pathways has very little room to express itself
across that set — three of the four conditions enter metabolism at nearly the same point.

The 25 carbon sources span hexoses, pentoses, sugar acids, amino sugars, a disaccharide, polyols and
organic acids. Different entry points into central metabolism, so different genes become load-bearing —
and the model tracks that.

Wild-type growth alone shows the range is real: maltose **1.78 /h** (a disaccharide, ~2× glucose, exactly
as expected), glucose 0.877, glycerol 0.495, acetate **0.210**, glycolate **0.153**.

## The bug this run caught — and why the contradiction was the tell

The first run printed `constant pattern for 217/217 (1.0) across 1 shape` **alongside** 23 exact-set
matches. Those cannot both be true: an exact-set match requires the predicted set to equal the true set,
and a two-sided gene's true set is by definition neither empty nor full, so a constant prediction can
never match one.

Cause: `pattern_distribution` still iterated the four hardcoded media names while the prediction was keyed
by carbon source. Every lookup missed, every gene read as `....`, and the function reported a fabricated
"100% constant". `switch_accuracy` and `constant_baselines` had been generalised; this one had not.

Had the contradiction not been checked, this run would have published **"the model is 100% constant"** —
a dramatic and completely false claim, and the exact opposite of the truth. Pinned by
`test_pattern_distribution_MUST_use_the_callers_conditions`.

> **Reusable:** when two metrics from a single run disagree in a way that is logically impossible, that is
> a bug signal, not a finding. Check the harness before writing the interpretation.

## Honest limits

- **No reproduction gate.** The Orth substrate ships the paper's own iJO1366 FBA columns, so the pipeline
  could be checked against a published prediction before any new number was trusted. The Fitness Browser
  has no such column. This number has no equivalent prior check.
- **Aerobicity is not varied.** All 25 are aerobic carbon-source assays, so the oxygen axis the 4-media set
  carried is absent here. This measures **carbon-source specificity**, not oxygen response — the two
  results are not strictly the same question, and the difference cuts both ways.
- **RB-TnSeq fitness is a pooled competition readout**, not a growth/no-growth call. The `fit < -2` cutoff
  is inherited from the shipped Keio validation (Bernstein 2023) for comparability, not re-derived here.
- **10.6% exact-set is still low in absolute terms.** Six times the lift is a real improvement over a
  near-null result; it is not a model that reproduces conditional essentiality well.
- Replicates per carbon source are averaged (62 experiments over 25 sources).

## What this changes

The 4-media conclusion — "adding reactions can't help, regulation is the lever" — was drawn from a
substrate where the model had almost no room to demonstrate conditional resolution. Both prior findings
should be re-read in that light:

- The **gap-filling negative** (0 binary flips at any dose) was measured on the same 4 media. Worth
  re-running here before treating it as settled.
- The **pFBA regulatory positive** (~5× lift over null) was also 4-media. If the model already switches on
  a carbon panel, the headroom that intervention was filling may be smaller than it appeared.

## Reproduce

```bash
uv run python scripts/fba_conditional_carbon_validate.py
```

Needs `feba.db` (7.4 GB, not committed — figshare `10.6084/m9.figshare.25236931`, CC BY 4.0).
Sidecar: `wiki/fba_conditional_carbon_2026-08-12.json`. Tests: `tests/test_fba_conditional_essentiality.py`.
