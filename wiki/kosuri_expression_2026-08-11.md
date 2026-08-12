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

### Is the headline just identifying the library? — decomposed, and it survives

Every part belongs to a design library (BIOFAB / Anderson / Salis / vectors), and libraries differ in mean
strength. A model that only learned *"this looks like a BIOFAB part"* would score well against the global
mean while ranking nothing. The per-element numbers are therefore decomposed:

| | headline R² | **library identity alone** (no sequence) | **within-library R²** (identity removed) | Spearman | RMSE (log2) |
|---|---|---|---|---|---|
| RBS | 0.6123 | 0.1526 | **0.4951** | 0.808 | **0.872** |
| promoter | 0.4165 | 0.0223 | **0.3522** | 0.563 | **1.465** |

Library identity alone explains little (0.15 / 0.02). After centring **both** truth and prediction within
each library, genuine part-level ranking remains: **0.495 / 0.352**. The claim survives with a haircut,
not a collapse. RMSE is the denominator-free reading — a novel RBS lands within ~0.87 log2 (≈1.8×) and a
novel promoter within ~1.47 log2 (≈2.8×) of its measured mean.

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

> **Correction history.** Four corrections, all of them narrowing the claim:
> 1. An early version concluded sequence generalisation "is not demonstrated and this dataset alone cannot
>    answer it" — drawn while only 2 of the 4 supplementary files had been used.
> 2. A later version headlined **0.781**, the ΔG arm, smuggling a non-reproducible feature into a
>    "from sequence" claim.
> 3. The OOD section stated *"most of that collapse is small-data, not unfamiliarity."* **That was
>    arithmetically false.** Of the promoter's 1.0712 total drop (0.4165 → −0.6547), the scarcity component
>    is 0.4818 (**45.0%**) and the shift component 0.5894 (**55.0%**) — shift is the *larger* share. The
>    section now avoids the decomposition entirely, because at n_train = 22 the fit is degenerate and no
>    clean decomposition is available (see ⚠ above).
> 4. Significance used a `mean − 2σ` rule, which compared one structured point against a control spread
>    while assuming near-normality, using a population σ, with no adjustment across seven tests. Replaced
>    by an empirical randomization percentile over 200 controls. **Both BIOFAB verdicts survive** the
>    stricter test; nothing else was ever significant.

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
training set at the same time (holding out BIOFAB leaves just **22** promoters to train on). Each holdout
is therefore paired with a **size-matched random control**: same train and test sizes, parts drawn across
all libraries, **200** repeats. Significance is the **empirical percentile** of the structured holdout in
that control distribution — the fraction of random splits that scored at or below it.

> **What the control does and does not do.** It compares *structured removal* against *iid same-size
> removal*. It does **not** cleanly isolate shift from scarcity: a random 22-promoter training draw
> usually still **contains** BIOFAB parts (80% of promoters), so the comparison bounds the shift rather
> than decomposing it. An earlier version of this memo described it as a decomposition; that was too strong.

| element | held out | n_test | n_train | LOLO *(vs global mean)* | **within-library R²** | **RMSE** | iid same-size | **percentile** |
|---|---|---|---|---|---|---|---|---|
| **RBS** | BIOFAB | 55 | 56 | 0.2524 | 0.4497 | 1.145 | 0.6266 | **0.005** |
| RBS | BioBrick/Anderson | 31 | 80 | 0.7143 | 0.4953 | 0.737 | 0.6041 | 0.810 |
| RBS | Salis | 13 | 98 | 0.6248 | 0.3799 | 0.876 | 0.6052 | 0.470 |
| RBS | vector/other | 12 | 99 | 0.6137 | 0.6948 | 1.076 | 0.5916 | 0.475 |
| **PROMOTER** | BIOFAB | 90 | **22** | −0.6547 | **0.0000** ⚠ | 2.288 | −0.0653 | **0.005** |
| PROMOTER | BioBrick/Anderson | 15 | 97 | 0.1478 | **−0.3495** | 2.360 | 0.3450 | 0.175 |
| PROMOTER | vector/other | 7 | 105 | 0.5441 | 0.5676 | 1.384 | 0.2182 | 0.655 |

**Only the two BIOFAB rows are load-bearing**, and both survive the randomization test: each was worse
than **199 of 200** same-size random splits (percentile 0.005). Every other row sits between the 17th and
81st percentile — unremarkable. BIOFAB is **50% of the RBSs and 80% of the promoters**; the dominant
design style is carrying the model.

> **⚠ Promoter-BIOFAB is a degenerate fit, and that is the honest report.** At n_train = 22 the regressor
> cannot split (`min_samples_leaf = 20`) and emits **one constant** — 13.1199 — for all 90 held-out
> promoters. Its within-library R² is exactly 0.0000 *by construction*, and its −0.655 is entirely about
> where that constant sits relative to BIOFAB's mean. So: **no clean shift-vs-scarcity decomposition is
> available at this training size**; what can be said is that the structured holdout was worse than 199 of
> 200 iid same-size splits. The artifact flags this as `prediction_is_constant`.

**Two of three promoter libraries have no within-library ranking at all** (0.0000 and −0.3495), which the
global-mean column hides. RBS ranking survives every boundary (0.38–0.69). Reporting only the
offset-inclusive R² would have overstated promoter transfer.

**Practical reading for a designer:** expect roughly the headline figure for a novel RBS resembling these
design styles, and treat novel *promoters* from an unfamiliar library as essentially unranked. Small
holdouts (n = 7–15) sit on very wide control bands and should not be over-read in either direction.

## Limits

- These remain **designed** parts throughout; nothing here tests truly random sequence.
- Per-element means are **simple** means, not adjusted for partner main effects. Adjusted means would
  give a sharper estimate of intrinsic part strength.
- The seven library tests carry **no multiple-comparison adjustment**. At 200 controls the two
  load-bearing rows sit at the resolution floor (0.005 = 1/200); a Bonferroni-style correction across
  seven tests would need more controls to separate them from it.
- The BIOFAB verdicts have not been checked for **model-class sensitivity** — a different regressor
  (especially at n_train = 22, where this one degenerates to a constant) could move them.
- Nothing here is wet-lab validated. It is a prediction about a measured dataset.

## Reproduce

```bash
uv run python scripts/kosuri_expression_validate.py \
  --sd03 <path>/sd03.xls --sd02 <path>/sd02.xls --sd01 <path>/sd01.xls
```

Data is **not committed** (third-party supplementary; PNAS is Cloudflare-gated to scripts).
Sidecar: `wiki/kosuri_expression_2026-08-11.json`. Tests: `tests/test_kosuri_expression.py`.
