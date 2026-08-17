# Nitrogen replicates all three carbon findings — the first non-carbon substrate axis

**Date:** 2026-08-17 · **Verdict:** `ALL_THREE_REPLICATE` · **Determinism gate:** PASS
**Pre-registration:** `wiki/fba_nitrogen_prereg_2026-08-17.md` (frozen before any solve, + a disclosed amendment)
**Artifact:** `wiki/fba_conditional_nitrogen_2026-08-17.json` · **Runner:** `scripts/fba_conditional_nitrogen.py`

## Result

Panel: **13 nitrogen sources** (of 16 in the assay; 2 dipeptides + casamino acids are unmappable to
iML1515 exchanges). 1,339 genes with complete rows → **155 conditionally essential** (two-sided).

| pre-registered prediction | carbon | nitrogen | verdict |
|---|---|---|---|
| **P1** ratio distribution bimodal; threshold band near-empty | 0 of 16,676 in [0.001, 0.05] | **0.0 %** in band | **REPLICATES** |
| **P2** missed-essential cells are mostly FLAT (deletion did nothing) | 76.9 % | **91.1 %** (408/448) | **REPLICATES** |
| **P3** per-cell beats the best-constant null | lift positive | **0.6799 vs 0.5355 (+0.1444)** | **REPLICATES** |

All three were frozen *before* the first solve, each with a number carbon had to supply and nitrogen
could have falsified. None did.

## What this establishes

**P1 is the load-bearing one.** The E-Flux null was explained by the readout being a threshold on a
*bimodal* quantity — and the honest worry was that this was a property of the carbon panel rather than of
the method. It is not: on an independent substrate axis, **exactly zero** cells land in the threshold
band. The readout explanation generalises, and the closed E-Flux negative does **not** need re-opening.

**P2 is stronger on nitrogen than on carbon.** 91.1 % of missed essential cells are deletions that
changed nothing at all (vs 76.9 % on carbon). The conditional deficit is a **model** problem — the
metabolic network has no route by which those knockouts matter — not a threshold-tuning problem. Retuning
the cutoff cannot recover cells whose deletion produced no change.

**P3 is a genuine, if modest, signal.** +0.1444 over the best constant baseline, so the model is doing
real conditional work on nitrogen, not just exploiting class imbalance.

## The determinism gate, and the failure it caught

The pre-registration made determinism mandatory *before* measuring — applying today's retraction lesson.
It earned its keep immediately.

**The first run FAILED it**: 147 of 2,015 cells differed between two identical `processes=1` passes, and
no verdict was reported. Probing showed the largest disagreement anywhere was ~3e-11 — ordinary float64
behaviour for one LP objective divided by another — that **zero** cells crossed the call line, and that
the differing-cell *count* was itself unstable across runs (147, 58, 66), the signature of a threshold
sitting inside float noise. Bit-equality was measuring the wrong thing.

The gate was therefore rebuilt around what the claim actually needs — and the amendment is disclosed in
the pre-registration rather than applied silently, because **any criterion loosened after seeing a
failing run deserves suspicion**:

| | original | amended |
|---|---|---|
| numeric | bitwise equality | **derived** safety margin, not an asserted tolerance |
| calls | *not checked* | **zero flips across the `FRAC=0.01` line** |
| headline metric | *not checked* | **identical across passes** |

The safety factor is `min_margin_to_threshold / max_abs_delta` — how much headroom exists between the
numerical noise and the decision line. It **adapts to the data**: if a future panel's ratios crowd the
threshold, the margin shrinks and the gate fails on its own. A fixed tolerance cannot do that, and
`test_determinism_gate_fails_when_a_cell_sits_near_the_line` pins exactly that behaviour (it fails on
*margin*, with zero call flips — which is the point).

**This run:** 0 call flips · drift 2.9e-11 · nearest cell to threshold 1.0e-2 · **safety factor 3.46e8**
against a bar of 1000 · headline metric **identical** across both passes (0.6799, 0.6799).

**A second control was added because a flag proved insufficient.** The first run wrote
`"deterministic": false` and then printed `per_cell_agreement` and all three verdicts beside it — numbers
a reader could quote in good faith despite the pre-registration saying no verdict may be reported.
`redact_unverified` now **removes** every solve-derived field on a determinism failure. A flag sitting
next to the numbers it invalidates is not a control; removing the numbers is.

## Named caveats

1. **These are not nitrogen-only perturbations.** Alanine, serine, aspartate and glutamine supply carbon
   as well. That matches the real assay (glucose minimal medium + test compound as sole N source), so it
   is biology rather than a modelling error — but the conditions must not be described as pure N swaps.
2. **13 of 16 sources**; the panel is not the full assay.
3. **2 replicates per source**, averaged — thinner than carbon's ~2.5.
4. **`EX_nh4_e` is closed explicitly** in every condition. Without that, residual ammonium would make
   every condition silently score as ammonium — the exact failure `_ALL_CARBON` guards on the carbon axis.

## Next

The three carbon findings now hold on two independent substrate axes. The remaining untouched substrate
in `feba.db` is **55 stress experiments** (a different perturbation *class*, not just a different
nutrient) and the **47-organism Ortholog table** for cross-organism conditional essentiality — the
natural next axes, in that order.
