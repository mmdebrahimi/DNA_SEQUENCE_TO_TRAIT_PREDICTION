# The stacking question, measured: a third channel DOES add, and the orthogonality screen got the direction backwards on the pair it was applied to

**Date:** 2026-10-26
**Artifact:** `wiki/glm_tile_hybrid_arm_2026-10-26.json` (+ `.md` sidecar)
**Substrate:** GEO GSE144621 peak tiles, leave_peak_out, ceiling 0.941, 10 seeds
**Top-level verdict (the ANSWERED question, unchanged):** `HYBRID_ADDS` — hybrid +0.0339 over max_tx, 10/10 seeds, CI95 [0.0238, 0.0409]
**STACKING verdict (the NEW question):** `WEAK_STACK`

## THE PREDICTION WAS FALSIFIED, AND THAT IS THE RESULT

Recorded in `PREREGISTERED["stacking"]["predicted_outcome"]` **before the run**, quoted verbatim:

> decomposition-as-third is predicted to LOSE, and to lose by MORE than gc did. Measured pairwise orthogonality (10-seed median) puts decomposition <-> max_tx at +0.7812, against gc <-> max_tx's +0.5508, so decomposition is the MORE redundant of the two available third channels -- and gc-as-third already lost (0.4217 against the 2-channel's 0.4254). The expected verdict is therefore NO_STACK for both new combos.

**It is wrong.** Measured, every combo at `n = 10`:

| third channel added | median gain over the 2-channel | seeds winning | CI95 | zero |
|---|---|---|---|---|
| `decomposition` | **+0.0145** | **10/10** | [0.0108, 0.0172] | **excludes zero** |
| `gc` + `decomposition` | +0.0103 | 9/10 | [0.0039, 0.0156] | **excludes zero** |
| `gc` alone | -0.0039 | 4/10 | [-0.0095, 0.0026] | includes zero |

`decomposition`-as-third gains **+0.0145 on 10/10 seeds with a CI that excludes zero**. I
predicted it would lose, and lose by more than `gc` did. It is the only one of the three that wins on
every seed.

**The direction is the opposite of what the orthogonality argument predicts.** Re-measured here:

| pair | median | range over 10 seeds |
|---|---|---|
| `max_tx` ↔ `decomposition` | +0.7812 | +0.7451 … +0.8146 |
| `gc` ↔ `decomposition` | +0.6041 | +0.5828 … +0.6278 |
| `max_tx` ↔ `gc` | +0.5508 | +0.5478 … +0.5754 |
| `conv` ↔ `gc` | +0.4593 | +0.4100 … +0.6557 |
| `conv` ↔ `decomposition` | +0.4213 | +0.3890 … +0.5440 |
| `conv` ↔ `max_tx` | +0.3862 | +0.3590 … +0.4790 |

`decomposition` is the **most** redundant channel with the comparator (`max_tx` ↔ `decomposition`
+0.7812) and it is the one that **adds**. `gc` is markedly less
redundant (+0.5508 with `max_tx`, +0.4593 with
`conv`) and it **subtracts** (-0.0039, 4/10 seeds). Adding `gc` on top of the winning 3-channel
makes it worse, not better (+0.0103 against +0.0145).

So **pairwise Spearman between prediction vectors is not a valid screen for whether a channel
stacks** on this rung. That is a stronger and more useful finding than the negative I expected, and
it could only have been obtained by running the combo the argument was used to decline.

### SCOPE on that finding — and it bounds what this measured

The two channels in the counterexample are **not two modalities.** `scripts/glm_tile_predictions.py`
hands the SAME decomposition dict to both scorers — `tile_features(dec.get(i, []))` builds the full
feature vector, `max_tx_scalar(dec.get(i, []))` returns `[log1p(max Tx_rate)]` — so `max_tx` is a
lossy **scalar reduction** of exactly what the `decomposition` channel sees in full.

Three consequences, and the third is the one that matters:

1. Their **+0.7812** redundancy is mechanically expected, not a
   surprise. Two readings of one computation at different resolutions *should* rank alike.
2. The pair that is most redundant and still stacks is therefore a **reduction/superset pair**, which
   is a special case rather than a general counterexample.
3. What is refuted is prediction-vector Spearman as a screen **across that pair**. It is **NOT**
   shown invalid between genuinely different modality families — none were tested here.

Read the result as *"finer-grained access to a channel already present pays"*, **not** as *"a new
kind of information pays"*. The honest version of the headline claim is narrower than the version
this memo led with, and the narrower version is still worth having: the decline was wrong, and it was
wrong for a reason that generalises less far than first stated.

This fact was already in the document — inside the quoted `negative_scope`, which scopes the **gc
negative**. It was attached to the wrong finding, so the positive headline carried no caveat at all.
That is the correction.

## WHY THE VERDICT IS STILL `WEAK_STACK`, AND WHY IT STAYS THERE

The frozen bar asks for `median >= 0.02` **and** `>= 8` of
10 seeds. `decomposition` clears the seed count outright (10/10) and misses the
margin (+0.0145 against 0.02). Mechanically: `WEAK_STACK`.

**The margin is NOT being revisited, and that is deliberate.** `0.02` was inherited
verbatim from the parent bar and recorded AS an inheritance, not a derivation:

> the parent bar above (margin 0.02, 8 of 10 seeds), itself mirroring scripts/promcalc_decomposition_arm.py PREREGISTERED. Reused VERBATIM rather than derived, so the two 'does more information help?' results stay directly comparable; the inheritance is recorded as an inheritance, not as a derivation.

A consistent gain on every seed with a CI excluding zero is a real signal, and the temptation to
lower a margin that a real signal just missed is exactly the retro-fit this plan was written to
avoid. The bar was frozen before the number was read; it stays frozen after. What changes is the
record: the question is **no longer unmeasured**, and the answer is *"a third channel adds roughly
half the margin, reproducibly"* rather than *"no third channel adds"*.

Anyone who wants to act on +0.0145 should set a NEW bar for a NEW question (how large a gain is worth
a channel's cost?) on its own evidence, not relax this one after the fact.

## WHAT THE ORIGINAL DECLINE COST

`scripts/glm_tile_hybrid_arm.py::PREREGISTERED["precondition"]` carries `not_a_gate`:

> the orthogonality number is REPORTED, not used to skip the run -- a prediction that a combine will
> fail is not a substitute for measuring it

and `tests/test_glm_tile_hybrid_arm.py::test_the_precondition_is_reported_not_used_as_a_gate` pins
it. The run was declined on precisely the reasoning that rule forbids, in the file that carries the
rule. Had it stayed declined, the record would have said *no third channel exists in these tables* —
and the measurement says one does. The rule earned its keep.

## THE 2-CHANNEL QUESTION IS UNMOVED — the regression evidence

Verified against `git show HEAD:wiki/glm_tile_hybrid_arm_2026-10-26.json`, field by field:

- `_paired_hybrid_vs_comparator` +0.0339, 10/10 seeds, CI95 [0.0238, 0.0409] — **byte-identical**
- `hybrid_conv_maxtx` +0.4254, `hybrid_conv_maxtx_gc` +0.4217, `max_tx` +0.3968, `conv` +0.3234, `control_shuffled_conv` +0.2703 — all **identical**
- `_paired_hybrid3_vs_comparator`, `_paired_control_vs_comparator`, `_paired_conv_vs_comparator`, and the whole `orthogonality` block — **identical**
- top-level verdict still `HYBRID_ADDS`

The generalisation added questions and did not perturb the answered one. Without this, the stacking
numbers would not be trustworthy either.

Guards both silent, as expected and not as a pass: `_short_arm` absent and `_skipped_combos` absent,
because all 10 tables carry all four channels. `_short_arm` is therefore
`[inferred]` in production and is exercised only by the mutation test
`test_main_REFUSES_a_SHORT_ARM_end_to_end_under_MUTATION`.

## The statistic that did not previously exist

Before this run the artifact carried `_paired_hybrid3_vs_comparator` = +0.0302
— `conv+max_tx+gc` against **the scalar**. That says the 3-channel beats `max_tx`; it never said
whether it beats the **2-channel hybrid**. It does not: -0.0039 on 4/10 seeds. The two numbers are
both correct and answer different questions, and only the second one is about stacking.

## Honest limits

- **The CI brackets the MEAN while every point estimate here is the MEDIAN, and `z=1.96` is the
  wrong critical value at n=10.** `promcalc_decomposition_arm.paired` builds a normal approximation
  on `fmean(d)` with z, while the headline figures are `median(d)` — two different estimands — and at
  n=10 the correct value is `t(9)=2.262`, so the published intervals are ~15% narrower than their
  nominal level. **Re-derived under t and the conclusion SURVIVES on all three combos**
  (`decomposition` [0.0103, 0.0177] and `gc`+`decomposition` [0.0030, 0.0165] still exclude zero;
  `gc` still includes it), with sign tests corroborating independently (p=0.00098 / 0.0107 / 0.828).
  So nothing above is wrong. It is **deliberately NOT fixed**: `paired()` is shared with the
  published parent arm, so correcting it moves parent numbers and needs its own pre-registration.
- **The verdict vocabulary was TIGHTENED after this run, and that needs disclosing.** As published,
  `WEAK_STACK` fired on `median > 0` alone — no seed condition — so this result (+0.0145
  on 10/10 seeds, CI95 [0.0108, 0.0172] excluding zero) and positive noise on 4 of 10
  with a CI spanning zero received the SAME machine-readable verdict. A fourth key,
  `INDETERMINATE_GAIN_NOT_REPRODUCIBLE`, now separates them, and `WEAK_STACK` additionally requires
  `>= min_seeds_winning` seeds. Editing a pre-registered bar after the result is the retro-fit this
  arm exists to avoid, so two things bound it: **no new threshold** was introduced (the condition
  reuses the already-frozen `min_seeds_winning`), and the change was **outcome-neutral** — every
  measured verdict here is the one it was published with, checked by diffing the regenerated artifact
  against the pre-change copy rather than asserted. Had one moved, the edit was to be reverted.
- **`WEAK_STACK` is an aggregate over more than one combo against ONE bar.** Read the per-combo
  block: `decomposition` +0.0145/10/10, `gc+decomposition` +0.0103/9/10, `gc` -0.0039/4/10.
  Three superset combos against one pre-registered margin is a multiple-comparison exposure; the
  aggregate takes the best of them.
- **One rung, one combine rule, one substrate.** Tile rung, leave-peak-out, ceiling 0.941,
  rank-average with no parameters. A different combine rule (a fitted stack, a weighted average) is
  untested and is a different question.
- **The negative scope recorded before the run still applies to the part that WAS negative** —
  `gc`-as-third. Quoted verbatim: a negative closes THESE FOUR SCORERS on THIS rung, not stacking in general. conv and gc are both composition-family; max_tx and decomposition are both readings of the same Promoter Calculator output, max_tx being a scalar reduction of it. So the four span only TWO modality families, and a channel from a NEW family is untested here. A negative is evidence about the available channels, never about the mechanism.
- **The four channels span two modality families and a channel from a NEW family remains untested.**
  That was true before the run and is still true; what is no longer true is the claim that no
  *available* channel adds.
- **+0.0145 is a RANK improvement.** The combine ranks and does not dose, so this is not a calibrated
  effect size and nothing here licenses a magnitude claim.
- **`max_tx` and `conv` are fitted by ridge on the same split**, so their errors may be correlated
  through the split rather than through the biology — unchanged from the parent arm and inherited by
  every number above.
