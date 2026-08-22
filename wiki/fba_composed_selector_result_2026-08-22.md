# The first intervention in this arc that doesn't fail — and it recovers one gene

**Date:** 2026-08-22 · **Verdict:** `H1_WEAKLY_SUPPORTED` (pre-committed, §6)
**Pre-registration:** `wiki/fba_composed_selector_prereg_2026-08-22.md` (frozen at `d08a210`, before scoring)
**Artifact:** `wiki/fba_composed_selector_2026-08-22.json` · **Script:** `scripts/fba_composed_selector.py`
**Tests:** `tests/test_fba_composed_selector.py` (8) · Preflight: `wiki/fba_selector_preflight_2026-08-22.md`

## Result

| endpoint | bar | measured |
|---|---|---|
| **Gate-0** — wildtype unchanged in every condition | hard | **held in all 11, all 3 arms** |
| **Primary** — targets recovered | ≥ 4 of 8 | **1 of 8** (`ilvB`) |
| **Mechanism** — partner isozyme in the gated set | 1.0 | **1/1 at p10, 3/3 at p20** |
| **Guardrail** — false positives | ≤ +20 % | **+2.6 % — held** |
| Recall over the 131 | reported | 0.3376 → 0.3413 |
| Determinism | identical, 2 runs | **passed** |

**This is the first intervention across the whole arc that neither collapses the wildtype nor breaches the
guardrail.** The two before it: expression-gating +334 % FP with recoveries that were pure
wildtype-collapse artifacts; biomass completion +53.6 % FP. Here every safety property held — and the
effect is small.

## Why it's safe: the composition step

The preflight established that both prior failures shared an under-diagnosed cause, and that **two of my
own claims were wrong**:

1. The collapse was **not** mainly about volume — a per-gene threshold marking only **7 genes** still
   collapsed galactose.
2. "Safe by construction" (restricting to genes whose knockout disables nothing) **failed in 2 of 11
   conditions**, because `gpr_disabled_reactions` is a *single-gene* property: gate one member of an
   isozyme pair and the reaction survives; gate **both** and it dies.

The fix is the third step both attempts lacked — verify the selected **set** jointly and drop the
highest-expressed member of each collision until nothing is disabled. It bites exactly where predicted:

| condition | selected (p20) | dropped by the joint check |
|---|---|---|
| D-Ribose | 90 | **13** |
| N-Acetyl-D-Glucosamine | 77 | 9 |
| D-Gluconate | 75 | 8 |
| D-Sorbitol | 71 | 8 |
| D-Glucose | 0 | 0 |

Those are the same conditions that collapsed before. With the joint check, wildtype growth is unchanged
**to five decimal places** in every one.

## Honest size of the effect

| arm | TP gained | FP gained | increment precision |
|---|---|---|---|
| p5 | +1 | +1 | 0.500 |
| **p10 (primary)** | **+2** | **+2** | **0.500** |
| p20 | +4 | +5 | 0.444 |

Baseline model precision is **0.701**. So the increment is *worse than what it was added to* — better than
biomass completion's 0.38, far better than expression-gating, but not a clear win. **Two true positives
for two false positives is a marginal, non-harmful effect, not a fix.**

The pre-registration forbids tuning the percentile to reach the bar, and p20 does not reach it either
(still 1 of 8). Reported in full; the primary stays p10.

## What this changes

The prior conclusion — *"the operator works, the selector doesn't"* — was too pessimistic. Corrected:

> A **composed** selector (per-gene threshold → single-gene eligibility → joint verification) is
> **wildtype-safe and guardrail-clean**, and recovers **1 of 8** masked genes with a perfect mechanism
> check. Expression-based boolean gating is viable; its measured ceiling on this target set is small.

That is a real, if modest, positive after four consecutive negatives — and the mechanism check (every
recovery had its partner gated, 3/3 at p20) says the one recovery happened for exactly the predicted
reason, not by accident.

## Honest limits

1. **11 of 25 conditions scored** — 14 have no matched expression and are excluded, reported, not dropped
   silently. A target cannot be recovered where there is no data.
2. **1 of 8 is below the pre-registered bar.** The verdict is `H1_WEAKLY_SUPPORTED`, not supported.
3. **Increment precision 0.500 < baseline 0.701** — the additions are lower-quality than the existing
   calls, even though the guardrail passes.
4. **The TRN was fetched and verified (8/8 targets covered) but deliberately NOT used here** — folding it
   in would confound two changes. It remains a separate, untried pre-registration.
5. **Strain mismatch** (PRECISE-1K K-12 MG1655 vs `Keio` BW25113) and **low mRNA ≠ absent protein** both
   unchanged.
6. **The drop-the-highest-expressed resolution rule is a choice, not a derivation** — fixed in advance,
   not tuned.

## Next

The honest options, in order of how much they are supported by what was measured:

1. **The TRN variant** — pre-registered separately, already verified to cover 8/8 targets. It is the one
   remaining named in-family lever and it now has a *safe* gating harness to plug into.
2. **Accept the ceiling.** Four of the eight targets are the `ilv` genes, and only `ilvB` moved; the other
   AHAS subunits did not, which suggests the per-condition expression simply does not separate the two
   isozyme complexes in the 11 conditions available.

Not recommended: tuning the percentile, or widening the target set post-hoc. Both were foreclosed in
writing before the run.
