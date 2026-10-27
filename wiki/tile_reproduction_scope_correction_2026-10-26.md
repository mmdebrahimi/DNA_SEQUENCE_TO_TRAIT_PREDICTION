# Scope correction: the tile-table rebuild is a substrate-integrity check, not an independent reproduction

**Date:** 2026-10-26
**Three corrections.** Two to the published record, one to an earlier revision of this memo itself (Correction 3) — written before the stacking run and refuted by it.
**Corrects:** commit `e96b55912e22931dd1da7442db5fee80d4370fef` (short `e96b5591`, pushed to
`origin/main`) and `project_state/glm-umbrella-2026-10-07.md` Action Log **row 17**.
**Correction mechanism is ADDITIVE.** The commit is pushed and cannot be amended, so the over-claim
stays in the git record verbatim. This memo plus the umbrella row that points at it are the
correction; a reader who finds the claim finds this file by name.

## CORRECTION 1 — "a genuine independent reproduction" is an over-claim

**What was published.** Row 17 and the commit body say: *"The substrate was regenerated from scratch
(770 min promcalc) and the conv retrained from scratch, and it lands on the same numbers -- that is a
genuine independent reproduction of the whole tile-rung hybrid result."*

**Why it is wrong, by its own mechanism.** The conv is **deterministic given the substrate**, so a
4-decimal match was guaranteed the moment the substrates matched. Re-grepped at write time, not
copied from the plan:

| location | what it does |
|---|---|
| `dna_decode/glm/encoder.py:354` | `torch.manual_seed(seed)` |
| `dna_decode/glm/encoder.py:361` | `gen = torch.Generator().manual_seed(seed)` |
| `dna_decode/glm/encoder.py:334` | `peak_blocked_val_split(groups, val_fraction=..., seed=seed)` |
| `scripts/glm_tile_predictions.py:390` | `fit_encoder(tr, te, seed=s, device=a.device)` — the per-seed loop value threaded straight through |

Same code path, same seeds, same split. Nothing could have come out differently. An agreement that
is guaranteed by construction cannot corroborate the method it is offered as evidence for — that is
the same error class as a filter that removes nothing being presented as a control.

**What the rebuild DOES establish, and it is not nothing.** Two regenerated datasets are effectively
identical to the originals: the promcalc substrate (rebuilt under a toolchain version that is
**unrecoverable** — the original pin is gone) and the re-fetched MPRA table. That is a genuine
**substrate-integrity check**, and it is the non-trivial half. It licenses every downstream number
computed off `data/promcalc/tile_predictions/`; it does not license the word *reproduction*.

**Correct wording:** *the rebuilt substrate reproduces the published arm exactly, which verifies the
substrate; the method was not independently re-derived.*

## CORRECTION 2 — "redundant by construction" is false as stated

**What was published.** Row 17 says the four scorers span only two modality families, *"biophysics
(max_tx, decomposition, where max_tx is a scalar reduced from decomposition so they are **redundant
by construction**)"*.

**Measured.** Pairwise Spearman between prediction vectors, median over the 10 real seed tables,
re-measured at write time:

| pair | median | range |
|---|---|---|
| `max_tx` ↔ `decomposition` | **+0.7812** | +0.7451 … +0.8146 |
| `gc` ↔ `decomposition` | +0.6041 | +0.5828 … +0.6278 |
| `gc` ↔ `max_tx` | +0.5508 | +0.5478 … +0.5754 |
| `gc` ↔ `conv` | +0.4593 | +0.4100 … +0.6557 |
| `conv` ↔ `decomposition` | +0.4213 | +0.3890 … +0.5440 |
| **`conv` ↔ `max_tx`** — the pair that gained +0.0339 | **+0.3862** | +0.3590 … +0.4790 |

+0.7812 is high. It is not redundant, and "by construction" asserts an identity the data refutes:
a scalar reduction of a 55-feature decomposition discards most of that decomposition, and what it
discards is measurable.

**Self-validating.** `conv ↔ max_tx` recomputes to **+0.3862**, matching the published
`orthogonality.median` in `wiki/glm_tile_hybrid_arm_2026-10-26.json` exactly. The method that
refutes the claim reproduces the one number the artifact already carried.

## CORRECTION 3 — the CONCLUSION does not survive either, and this memo said it did

**An earlier revision of this file, written before the measurement, claimed:** *"the conclusion —
'no viable third channel exists in these four scorers' — had support; it was not the support that
was given,"* citing `gc ↔ conv` at +0.4593 being comparably orthogonal to the working pair at +0.3862
while `gc`-as-third still lost.

**Measured the same day, and it is wrong too** (`wiki/glm_tile_stack_arm_2026-10-26.md`):

| third channel added | gain over the 2-channel | seeds winning | CI95 |
|---|---|---|---|
| `decomposition` | **+0.0145** | **10/10** | [0.0108, 0.0172] — excludes zero |
| `gc` + `decomposition` | +0.0103 | 9/10 | [0.0039, 0.0156] — excludes zero |
| `gc` alone | −0.0039 | 4/10 | [−0.0095, 0.0026] — includes zero |

A third channel **does** add, reproducibly, on every seed. It lands at `WEAK_STACK` only because it
falls under the inherited 0.02 margin — not because it is absent.

**And the orthogonality reasoning is backwards, not merely insufficient.** `decomposition` is the
**most** redundant channel with the comparator (+0.7812) and it is the one that adds; `gc` is
markedly less redundant (+0.5508 / +0.4593) and it subtracts. Pairwise Spearman between prediction
vectors does not order these channels by whether they stack.

**The honest summary of all three corrections:** the reasoning was wrong (Correction 2), the
reasoning's *direction* was wrong (here), and the conclusion it was offered for was wrong (here) —
so the decline cost a real result, not just a weak argument. The only thing that was right was the
rule I overrode: `precondition.not_a_gate`, which says a prediction that a combine will fail is not
a substitute for measuring it.

## The deeper error, named rather than filed around

The decline itself was the failure. `scripts/glm_tile_hybrid_arm.py::PREREGISTERED["precondition"]`
carries `not_a_gate`: *"the orthogonality number is REPORTED, not used to skip the run -- a
prediction that a combine will fail is not a substitute for measuring it"*, and
`tests/test_glm_tile_hybrid_arm.py::test_the_precondition_is_reported_not_used_as_a_gate` pins it.
I declined a run on exactly the reasoning that rule forbids, in the script that carries the rule.

The remedy is the measurement, not an apology — and the measurement overturned the prediction.
See `wiki/glm_tile_stack_arm_2026-10-26.md`: run against a bar pre-registered before the run,
`decomposition`-as-third gains +0.0145 on 10 of 10 seeds with a CI excluding zero, where I had
predicted it would lose by more than `gc` did.

## Honest limits of this memo

- The matrix is 10 seeds on ONE rung (tile, leave-peak-out, ceiling 0.941). It is not a general
  statement about these scorers anywhere else.
- Spearman between prediction vectors is one orthogonality measure. Two scorers can be
  rank-correlated and still carry complementary error structure; a rank-average is only sensitive to
  the rank relationship, which is why this is the measure the arm uses.
- Correction 1 does not claim the rebuild was worthless — the opposite. It reclassifies what it
  licenses.
- The original claim remains in the pushed commit message. Nothing here can change that.
