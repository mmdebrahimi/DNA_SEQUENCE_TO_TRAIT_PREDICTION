# Track B — sequence→expression, both elements (2026-08-11)

The pre-registered test from `wiki/design_epoch_plan_2026-08-07.md`, run against Kosuri et al. 2013
(PNAS **110**:14024) — 12,563 constructed promoter × RBS combinations with measured DNA, RNA and protein.
Q2 of the design epoch: *will the host express it?*

**Bar set before the data was in hand:** beat protein **R² ≈ 0.82**, baseline re-fit on the training
split only, split **by element**.

## Reproduction gate (run before trusting any new number)

| | reproduced | published |
|---|---|---|
| protein, simple model | **0.7525** | 0.76 |
| RNA, simple model | **0.9238** | 0.92 |
| RNA, full model | **0.9623** | 0.96 |

> **Units trap:** `model.prot.simple` is stored in **log2**, `prot` is raw RFU. Compared in the wrong
> space it returns **R² = −15**, which reads like a broken loader rather than a units mismatch.

## Headline result — a novel part CAN be scored from its sequence

25 repeated `GroupShuffleSplit`s, held out by element, target `log2(protein)`. The held-out element is
never seen in training; the *other* element's identity is supplied (a designer knows their partner part).

| arm | held-out **RBS** | held-out **PROMOTER** |
|---|---|---|
| additive baseline | 0.4991 ± 0.051 | 0.2471 ± 0.061 |
| identity-only model | 0.2678 | 0.0304 |
| other-element-only *(control)* | 0.4991 | 0.2471 |
| sequence-only *(control)* | 0.1429 | 0.1119 |
| **other element + SEQUENCE** | **0.7762 ± 0.034** | **0.4967 ± 0.097** |
| ridge, same features *(comparator)* | 0.0681 ± 0.743 | **−3.3728 ± 2.73** |
| + ΔG *(ORACLE — see below)* | 0.8068 | 0.5521 |
| **per-element mean from sequence alone** | **0.6123** (111 pts) | **0.4165** (112 pts) |

**Both elements generalise from sequence**, and by a wide margin over the baseline (+0.28 RBS,
+0.25 promoter). Feature sets are mechanistic, not learned: k-mers, length, GC, plus the
Shine-Dalgarno core for RBSs and the σ70 **−35 / −10** boxes with their spacer for promoters.

### The promoter is harder than the RBS — the interesting finding

Confound-free, sequence alone: **RBS 0.612 vs promoter 0.417**. The promoter explains *more* of the
protein variance (~54% vs ~30% per the paper's ANOVA) yet is **less predictable from its letters**.

A plausible reading, offered as interpretation rather than result: translation initiation is dominated
by one short, well-understood motif (SD) and its spacing, which simple features capture. Promoter
strength depends on −35/−10 boxes *plus* UP elements, discriminator, spacer geometry, TSS selection and
supercoiling — much of which these features do not represent.

### Two numbers, not one — they answer different questions

| question | number |
|---|---|
| *How strong is this novel part?* (part-level ranking, no partner replication) | **RBS 0.612 · promoter 0.417** |
| *What will this specific construct express?* (novel part × characterised partner panel) | **RBS 0.776 ± 0.034 · promoter 0.497 ± 0.097** |

Reporting only one would mislead in one direction or the other. The per-construct figure is inflated by
partner replication and is conditional on a characterised panel; the per-element figure is confound-free
but measures a narrower thing.

### The comparator is not a strawman — measured, not assumed

The obvious objection is that the GBM only beat a weak additive baseline. It doesn't: **ridge** with
one-hot partner identity + standardised sequence features **collapses** on held-out groups (RBS mean
0.068, p5 **−1.77**; promoter mean **−3.37**, p5 −7.72). Regularised linear models extrapolate badly to
unseen element groups, so the additive baseline is a *strong* comparator.

### ΔG is an ORACLE bound, never a headline

ΔG is **dataset-provided** and spans promoter TSS → +30 of GFP, so it contains promoter-derived
sequence and is **not recomputable at design time**. It is reported only as an explicitly-named upper
bound (`other_plus_sequence_plus_deltaG_ORACLE`) and `sequence_verdict` headlines the **no-ΔG** arm.
`TSS.best` from S1 is excluded from promoter features for exactly the same reason — it was measured by
RNA-seq, not predicted.

> **Correction history.** An earlier version of this memo headlined **0.781** (the ΔG arm) and, before
> that, concluded sequence generalisation "is not demonstrated and this dataset alone cannot answer it."
> Both were wrong: the first smuggled a non-reproducible feature into a "from sequence" claim; the
> second was drawn while only 2 of the 4 supplementary files had been used.

## The composability result (a different question)

Held-out **combination** — both elements seen, the *pairing* is new:
additive **0.795** → GBM **0.893** → +ΔG **0.919**. Clears the 0.82 bar, +0.124 over the fair baseline.

And the falsification that keeps it honest: an **identity**-encoded model given an unseen promoter scores
**−0.014** — below chance and worse than the baseline it beat on combinations. That is a statement about
*encoding*, not about expression being unpredictable; supplied with real sequence the same split reaches
0.497.

## Honest verdict on the pre-registration

**By the stated falsifier — "beat 0.82, split BY ELEMENT" — this FAILS.** Nothing reaches 0.82 on an
element split, the best being 0.776 (RBS, with a characterised promoter panel).

The bar was also **mis-specified for that split**, and the incompatibility went unnoticed until the data
was in hand: 0.82 is a *combination-level, in-sample* number, and an element-strength model has no
strength for an unseen element — the baseline itself reaches only 0.25–0.50 there.

## Out-of-distribution: leave-one-library-out (2026-08-11)

Every part here is **designed** — BIOFAB, BioBrick/Anderson, Salis, cloning vectors — so a
leave-one-*element*-out score can still be interpolation *within a design style*. Provenance survives in
the part names, so a whole library can be held out. That is the closest this data gets to "a part nobody
in this dataset designed".

**The raw LOLO number is not interpretable on its own**, because holding out a library shrinks the
training set at the same time (holding out BIOFAB leaves just **22** promoters to train on). So each
holdout is paired with a **size-matched random control**: same train and test sizes, parts drawn across
all libraries, 20 repeats. The **gap** isolates library shift from data starvation.

| element | held out | n_test | n_train | LOLO | random, same size | **gap** | verdict |
|---|---|---|---|---|---|---|---|
| **RBS** | BIOFAB | 55 | 56 | 0.2524 | 0.6375 ± 0.083 | **−0.3851** | **real shift** |
| RBS | BioBrick/Anderson | 31 | 80 | 0.7143 | 0.5791 ± 0.126 | +0.1352 | no shift |
| RBS | Salis | 13 | 98 | 0.6248 | 0.6596 ± 0.178 | −0.0348 | no shift |
| RBS | vector/other | 12 | 99 | 0.6137 | 0.5947 ± 0.218 | +0.0190 | no shift |
| **PROMOTER** | BIOFAB | 90 | **22** | −0.6547 | −0.0540 ± 0.078 | **−0.6007** | **real shift** |
| PROMOTER | BioBrick/Anderson | 15 | 97 | 0.1478 | 0.2619 ± 0.410 | −0.1141 | inside noise |
| PROMOTER | vector/other | 7 | 105 | 0.5441 | 0.1179 ± 0.544 | +0.4262 | inside noise |

**Neither "OOD is fine" nor "OOD fails" — the answer is specific.** Transfer holds across three of four
RBS library boundaries (0.61–0.71, one of them *above* the within-distribution 0.612) and across the
promoter's small libraries. It degrades materially only when **BIOFAB** is held out — and BIOFAB is
**50% of the RBSs and 80% of the promoters**. The dominant design style is carrying the model; remove it
and performance drops by a gap of 0.39 (RBS) to 0.59 (promoter), both well clear of the control's spread.

**The control earned its place.** Promoter-BIOFAB's raw −0.655 looks catastrophic, but a *random* split
at that same training size scores −0.065 — so most of that collapse is **small-data, not unfamiliarity**.
Reporting the raw number as an OOD failure would have been wrong. Conversely the RBS-BIOFAB gap (−0.385
against a control σ of 0.083) is ~4.6σ and is real.

**Practical reading for a designer:** expect roughly the headline figure for a novel part resembling
these design styles; expect materially less for a part from a genuinely different paradigm. Small holdouts (n=7–15) have very wide control bands (σ up to **0.54**) and should not be over-read in
either direction — a scratch run with a different control seed moved the promoter vector/other gap from
+0.13 to +0.43 while the two BIOFAB verdicts were unchanged. Only the two flagged `real shift` rows,
whose gaps clear 2σ of their own control, are load-bearing.

## Limits

- These remain **designed** parts throughout; nothing here tests truly random sequence.
- Per-element means are **simple** means, not adjusted for partner main effects.
- Per-element means are **simple** means, not adjusted for partner main effects. Adjusted means would
  give a sharper estimate of intrinsic part strength.
- Nothing here is wet-lab validated. It is a prediction about a measured dataset.

## Reproduce

```bash
uv run python scripts/kosuri_expression_validate.py \
  --sd03 <path>/sd03.xls --sd02 <path>/sd02.xls --sd01 <path>/sd01.xls
```

Data is **not committed** (third-party supplementary; PNAS is Cloudflare-gated to scripts).
Sidecar: `wiki/kosuri_expression_2026-08-11.json`. Tests: `tests/test_kosuri_expression.py`.
