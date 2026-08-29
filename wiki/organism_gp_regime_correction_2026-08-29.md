# The organism-level g→p line is neither closed nor unexplored — it is the most-worked line here, with a measured bottleneck

I told the user twice this session that organism-level polygenic genotype→phenotype was a **"closed
negative."** Both times were wrong, and the second was wrong in a new way. This records the correction,
because a wrong strategic summary in a file that loads every session misdirects every future run — which is
the exact failure mode this session spent days measuring.

## Error 1 — over-compression (third occurrence)

"❌ closed negative" dropped the scope. The accurate statement:

- **Zero-shot foundation-model embeddings** on **natural populations** fail (0-for-5, de-confounded) —
  they learn population structure, not causal signal.
- **Constructed-variation designs succeed.** The yeast segregant cross (Bloom 2013) decoded **12/12
  quantitative traits at r 0.46–0.80** — recorded here as *the project's first clean g→p positive*. Blaming
  organism complexity was the wrong reading; the discriminating variable is **population design**.

My own memory already carried this warning and noted I had compressed it twice. This was the third.

## Error 2 — proposing as novel the line that is already deepest

I then proposed "an E. coli knockout panel × measured growth, with FBA as the mechanistic baseline" as the
innovation direction. **That is already built and heavily worked** — ~25 `wiki/fba_*` artifacts dated
2026-08-03 → 08-22, with pre-registrations, results, an adversarial review, and a retraction.

What exists:

- **Substrate, committed locally:** Orth 2011 Table S1 — **1,075 E. coli K-12 genes × 4 minimal media**,
  experimental E/N per cell, **68 conditionally essential**. Two-sided by construction (a gene is its own
  control; only the medium moves), and it ships the paper's own iJO1366 calls as a reproduction gate.
- **Fitness Browser `Keio`** is also already wired: carbon (28 conditions), nitrogen (32), stress (55).

So the substrate I was about to recommend has been in use for a month.

## What the line actually found

**Per-condition essentiality is solved well.** iML1515 MCC **0.70–0.74** across the four media, ~+0.04 over
the published iJO1366 gate.

**The conditional SWITCH — the property strain design depends on — is not.**

| | exact-set match | per-cell |
|---|---:|---:|
| null (always dispensable) | 0/68 | 0.5588 |
| iJO1366 (paper's own) | 4/68 = 5.9% | 0.5735 |
| **iML1515 (ours)** | **3/67 = 4.5%** | **0.5709 (+0.012 over null)** |

And the mechanism is diagnosed, not guessed: **the model is not switching at all.** Turning each gene into
a 4-character pattern over conditions, experimental truth has **12 distinct shapes**; iML1515 produces
**3**, predicting a constant pattern for **94%** of genes. On this subset it very nearly *is* the constant
predictor.

**`MIS_CONDITIONED = 0`, replicated on an independent axis** (nitrogen: 155 genes / 13 conditions, vs
carbon: 217 / 25; class shares within ~2 points). The model never fires in the *wrong* place — **it fires
roughly right or stays silent.** The failure is silence, not error, which is a materially different problem
from "the model is confused."

## The bottleneck, measured rather than asserted

Interventions tried: expression-gated GPR (recall 0.338 → 0.648), composed selector, biomass completion,
demand completion, e-flux bridge (**retracted**), cross-organism *P. putida*.

The composed selector recovered **1 gene of 8**, and the reason is now quantified:

| | |
|---|---:|
| Fitness Browser `Keio` carbon conditions | 28 (25 map to an iML1515 exchange) |
| PRECISE-1K distinct carbon values | 18 |
| **intersection** | **11** |
| unmatched conditions present in PRECISE-1K under any name | **0** |

Of 1,035 expression samples, **621 are glucose**; **six of the eleven scored conditions rest on ≤5 samples,
four on exactly 2.** Both datasets are fixed, so no plumbing fix widens the intersection.

**That is the wall: the conditioning signal is not measured in the conditions the phenotype data uses.**
Not "learned methods fail" — the information needed to condition on isn't there at the required resolution.

A second, separate ceiling: on the nitrogen axis **six of thirteen conditions give identical wildtype
growth (0.92593)**, so the model treats those sources as interchangeable and there is almost no
condition-specificity available to find. Some axes are structurally poorer targets than others.

## Corrected map

| regime | status |
|---|---|
| natural-population + zero-shot embedding | **closed negative**, 0-for-5, de-confounded |
| constructed variation → molecular phenotype | **works** — TEM-1 genome-edit path, Spearman 0.761 vs measured ampicillin fitness |
| constructed variation → organism phenotype, **per-condition** | **works** — FBA MCC 0.70–0.74 |
| constructed variation → organism phenotype, **condition-SWITCH** | **open, ~null, bottleneck measured** |

## What would actually move the open cell

**Condition-resolved conditioning data covering the phenotype's conditions** — expression (or any
per-condition signal) measured in the same 28 carbon / 13 nitrogen conditions the Keio fitness data uses,
rather than 11 with a median of a few samples. That is a specific, nameable acquisition target with a
quantified gap, not "we need more data."

Two cheaper things worth trying first, both untested as far as these artifacts show:
1. **Score the continuous knockout growth ratio as a ranking** rather than thresholding at 1% — the
   artifact already shows the deployed cutoff discards signal, and quantifies the ratio.
2. **Pick axes with real dynamic range.** The nitrogen axis has almost none by construction; carbon has
   condition-specificity for a third of PARTIAL_OVERLAP genes. Axis choice is a free lever.

## The lesson worth carrying

Both errors were *strategic summaries asserted from memory* while the grounding artifacts sat on disk. The
session's whole methodological finding — verify the claim against the artifact before publishing it —
applies to my own framing of the project, not just to its measurements.
