# Measure The Stacking Question Instead Of Arguing It, And Correct Two Published Over-Claims

## Lens status

Execution mode: sequential (recommended — Steps 2–4 all edit one file).
Repo-index: absent — no `.repo-index/`; source-only reading, no degradation.
Reuse signals: no `graphify-out/`, no conventional util dirs.
Tools: `gh` present; `sentrux` absent (architecture gate `n/a — sentrux not installed`); independent reviewer CLI absent.
Optional context: `wiki/` present; no `project-rules.md`, no `DESIGN.md`, no `AUDIT_REPORT.md`.
Baseline tests: `tests/test_glm_tile_hybrid_arm.py` + `tests/test_glm_tile_predictions.py` — **73 passed, 0 failed**.

## Problem Statement

I declined a run on an orthogonality argument, in a script whose own frozen bar forbids exactly
that, and I published two claims stronger than the evidence.

**The declined run.** `scripts/glm_tile_hybrid_arm.py::PREREGISTERED["precondition"]` carries
`not_a_gate` with the text *"the orthogonality number is REPORTED, not used to skip the run -- a
prediction that a combine will fail is not a substitute for measuring it"*, and
`tests/test_glm_tile_hybrid_arm.py::test_the_precondition_is_reported_not_used_as_a_gate` pins it. I
then told the user no third channel exists in these tables, on the strength of a reasoned
redundancy claim, without running the combo.

**Over-claim 1 — "independent reproduction".** Commit `e96b5591` and umbrella ledger row 17 call the
rebuild *"a genuine independent reproduction of the whole tile-rung hybrid result"*. It is not.
`dna_decode/glm/encoder.py:354` calls `torch.manual_seed(seed)`, `:361` builds a seeded
`torch.Generator`, the validation split is seeded via `peak_blocked_val_split(..., seed=seed)`, and
`scripts/glm_tile_predictions.py:390` passes the per-seed loop value through as
`fit_encoder(tr, te, seed=s, ...)`. The conv is therefore **deterministic given the substrate**, so a
4-decimal match was guaranteed once the substrates matched. What the rebuild actually establishes is
that two regenerated datasets — promcalc (under an unrecoverable toolchain version) and the re-fetched
MPRA table — are effectively identical to the originals. That is a real and useful **substrate-integrity
check**; it is not corroboration of the method, because the code path and seeds were identical and
nothing could have come out differently.

**Over-claim 2 — "redundant by construction".** I said `max_tx` and `decomposition` are redundant
because the former is a scalar reduced from the latter. Measured from the now-existing per-tile
vectors, median over 10 seeds:

| pair | median | range |
|---|---|---|
| max_tx ↔ decomposition | **+0.7812** | 0.745–0.815 |
| gc ↔ decomposition | +0.6041 | 0.583–0.628 |
| gc ↔ max_tx | +0.5508 | 0.548–0.575 |
| gc ↔ conv | +0.4593 | 0.410–0.656 |
| conv ↔ decomposition | +0.4213 | 0.389–0.544 |
| **conv ↔ max_tx** (the pair that gained +0.0339) | **+0.3862** | 0.359–0.479 |

0.7812 is high but not redundant. The method is self-validating: `conv ↔ max_tx` recomputes to
`+0.3862`, matching the published `orthogonality.median` exactly.

**And the matrix undercuts my reasoning while supporting my conclusion.** `gc ↔ conv` at +0.4593 is
comparably orthogonal to the working pair at +0.3862, and `gc`-as-third **still lost** (0.4217 against
the 2-channel's 0.4254). So orthogonality in this range does not buy stacking — which is direct
evidence for "no viable third channel here", and is not the evidence I offered.

Scope: pre-register a stacking bar, make the script able to express the stacking question, add the
decomposition combos, run it, and correct both over-claims. Non-goals: changing the frozen 2-channel
bar; killing any arm; re-running any substrate build; touching `my_skills`.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

- `uv run pytest tests/test_glm_tile_hybrid_arm.py tests/test_glm_tile_predictions.py -q` →
  **73 passed, 0 failed** in 4.05 s.
- `wiki/glm_tile_hybrid_arm_2026-10-26.json`: verdict `HYBRID_ADDS`; `arms` keys are exactly
  `hybrid_conv_maxtx`, `hybrid_conv_maxtx_gc`, `max_tx`, `conv`, `control_shuffled_conv`,
  `_paired_hybrid_vs_comparator`, `_paired_hybrid3_vs_comparator`, `_paired_control_vs_comparator`,
  `_paired_conv_vs_comparator`; every arm reports `n = 10`; `orthogonality` keys are
  `per_seed / median / min / max`.
- Published values this plan must not move: `hybrid_conv_maxtx` 0.4254, `hybrid_conv_maxtx_gc` 0.4217,
  `max_tx` 0.3968, `conv` 0.3234, `control_shuffled_conv` 0.2703,
  `_paired_hybrid_vs_comparator` median 0.0339 CI95 [0.0238, 0.0409],
  `_paired_hybrid3_vs_comparator` median 0.0302 CI95 [0.0226, 0.0351],
  `orthogonality.median` 0.3862.
- Over-claims as currently published: commit `e96b5591` on `origin/main` and umbrella ledger row 17.
  The commit is pushed and **cannot be amended**, so correction is additive.

## Verification Signal

1. **The pre-registered 2-channel verdict is unchanged.** After every step,
   `_paired_hybrid_vs_comparator` is still median 0.0339 CI95 [0.0238, 0.0409] and the verdict is
   still `HYBRID_ADDS`. This plan adds questions; it must not perturb the answered one.
2. **The new stacking statistic exists and is a measurement, not a figure** — every combo in `COMBOS`
   has a paired statistic with an n, a median and a CI, and every combo that is a strict superset of
   `("conv","max_tx")` additionally has one against the 2-channel.
3. **The short-arm guard refuses rather than reporting a short n** — demonstrated by mutation, since
   all 10 real tables carry all four channels and the guard cannot fire on real data.
4. **73 passed stays the floor** — the count rises by the number of tests added and no pre-existing
   test flips.
5. **The stacking verdict reaches the markdown**, not only the JSON (`render_md`'s own stated contract).
6. **Both corrections are reachable from the artifact that carries the over-claim** — the memo names
   commit `e96b5591` so a reader who finds the claim finds the correction.

## Implementation Steps

### Step 1: Pre-register the stacking bar as a new key, with the prediction recorded before any run
Files: scripts/glm_tile_hybrid_arm.py, tests/test_glm_tile_hybrid_arm.py
Depends on: none

**What changes:**
- Add `PREREGISTERED["stacking"]` as a NEW key: `question` ("does a THIRD channel add over the
  measured 2-channel hybrid?"), `comparator: "hybrid_conv_maxtx"`, `margin: 0.02`,
  `min_seeds_winning: 8`, `mirrors` citing the same frozen shape as the parent bar, and a `verdicts`
  map (`STACKS` / `WEAK_STACK` / `NO_STACK`).
- **Do NOT touch `margin`, `min_seeds_winning` or `comparator` at the top level** —
  `test_frozen_bar_mirrors_the_decomposition_arm` pins all three, and changing them would retro-fit
  the answered question's bar.
- **Derivation, not assertion:** the 0.02 / 8-of-10 shape is reused verbatim from the parent bar so
  the two "does more information help?" results stay comparable — the same justification the file
  already records in `mirrors`.
- Record `PREREGISTERED["stacking"]["predicted_outcome"]` BEFORE running: `decomposition ↔ max_tx` is
  +0.7812 against `gc ↔ max_tx`'s +0.5508, and gc-as-third already lost (0.4217 vs 0.4254), so
  decomposition-as-third is predicted to lose by more. Writing the prediction down first is what stops
  a negative being read as a surprise and stops a positive being retro-rationalised.
- Record `PREREGISTERED["stacking"]["negative_scope"]`: a negative closes *these four scorers*, not
  stacking in general — the channels here span two modality families and a new family is untested.

**Test strategy:**
- Unit: new test pins the stacking bar's shape and, in the SAME test, re-asserts
  `margin == 0.02` / `min_seeds_winning == 8` / `comparator == "max_tx"` unchanged — so a future edit
  cannot satisfy the new pin by loosening the old one.
- Unit: new test asserts `predicted_outcome` and `negative_scope` are non-empty, because an unrecorded
  prediction is the thing this step exists to prevent.
- Unit: the 36 pre-existing tests still pass, `test_the_precondition_is_reported_not_used_as_a_gate`
  included.

### Step 2: Generalise the paired statistics and the guards over COMBOS
Files: scripts/glm_tile_hybrid_arm.py, tests/test_glm_tile_hybrid_arm.py
Depends on: Step 1

**What changes:**
- Replace the two hardcoded paired lines with a loop over `COMBOS`, emitting
  `_paired_<combo>_vs_comparator` for each. **Keep `_paired_hybrid_vs_comparator` and
  `_paired_hybrid3_vs_comparator` as the exact names for the two existing combos** — published
  artifacts and tests read those keys, and renaming them would break backward compatibility for a
  cosmetic gain.
- Add `_paired_<combo>_vs_hybrid2` for every combo that is a strict superset of `("conv","max_tx")`.
  This is the statistic that answers the stacking question and which does not currently exist in any
  form.
- Extend the `missing` guard from the literal `("conv","max_tx")` to the **union of every channel named
  in any COMBOS entry**, so a table lacking `decomposition` returns `INDETERMINATE_MISSING_SCORER`
  instead of silently scoring a subset.
- Add a short-arm guard: if any arm's `vals` length differs from `len(tables)`, record
  `arms["_short_arm"]` and return a refusal. Today `if vals` accepts 3 of 10 and records `"n": 3`,
  which a reader comparing medians would not notice.
- Order the new refusal with the existing ones in `verdict()` so a refusal still outranks a headline
  (`test_refusals_outrank_a_headline` pins that property).
- Generalise the per-seed progress print so every combo appears rather than only `hybrid_conv_maxtx`.

**Test strategy:**
- Unit: every COMBOS entry has a `_paired_*_vs_comparator`; every superset combo also has a
  `_paired_*_vs_hybrid2`. Asserted over `COMBOS` itself, so adding a combo later cannot skip its statistic.
- Unit: the missing-scorer guard fires for a synthetic table lacking `decomposition` — **and the test
  derives the expected channel set from `COMBOS`**, so it cannot pass by hardcoding today's four.
- Unit: the short-arm guard fires on a synthetic 3-of-10 arm and the verdict is a refusal, not a headline.
- Unit: re-running the real arm reproduces the frozen published values — the 2-channel at 0.4254 and
  `_paired_hybrid_vs_comparator` at 0.0339 CI95 [0.0238, 0.0409]. Generalising the loop must not move
  the answered question.

### Step 3: Add the decomposition combos to COMBOS
Files: scripts/glm_tile_hybrid_arm.py
Depends on: Step 2

**What changes:**
- Add `hybrid_conv_maxtx_decomp: ("conv","max_tx","decomposition")` — the 3-channel the orthogonality
  argument was used to decline.
- Add `hybrid_conv_maxtx_gc_decomp: ("conv","max_tx","gc","decomposition")` — the 4-channel, included
  because with all four channels on disk it is free, and because it distinguishes "a third channel
  never helps here" from "this particular third channel does not".
- Leave the two existing entries byte-identical, so the published 0.4254 and 0.4217 remain reproducible.
- Depends on Step 2 because the paired statistics and guards must be generic BEFORE a combo relies on
  them — otherwise a missing statistic and a real negative are indistinguishable in the output.

**Test strategy:**
- Unit: `COMBOS` has four entries and the two original keys and tuples are unchanged, asserted literally.
- Unit: the new combos are supersets of `("conv","max_tx")`, so Step 2's superset rule gives each a
  `_paired_*_vs_hybrid2` — checked rather than assumed.
- Unit: the full 73-test floor still passes.

### Step 4: Report the measured pairwise orthogonality matrix in the artifact
Files: scripts/glm_tile_hybrid_arm.py
Depends on: Step 3

**What changes:**
- Add `orthogonality_matrix` to the artifact: for every unordered pair of channels used by any combo,
  the per-seed values plus median/min/max, computed with the existing `orthogonality()` helper.
- **Keep the existing `orthogonality` key byte-identical** (`per_seed / median / min / max` for
  conv↔max_tx) — published artifacts and tests read it, and the matrix is additive.
- Surface the matrix in `render_md` as well, per that function's own stated contract that a caveat
  carried only in the JSON is not a disclosure. Order rows most-redundant-first so the reader sees
  the stacking obstacle before the working pair.
- This is the measurement that refuted the "redundant by construction" claim; putting it in the
  artifact means the next reader does not have to re-derive it from the vectors.

**Test strategy:**
- Unit: `orthogonality["median"]` is still 0.3862 on the real tables and the matrix's conv↔max_tx
  entry equals it — the matrix must agree with the number it generalises.
- Unit: the matrix contains an entry for every unordered pair of channels named in `COMBOS`, derived
  from `COMBOS` rather than a fixed list.
- Unit: `render_md` output contains the matrix.

### Step 5: Run the arm and write the result memo
Files: wiki/glm_tile_stack_arm_2026-10-26.md
Depends on: Step 4
Test-strategy-ignore: tests/test_glm_tile_hybrid_arm.py, wiki/glm_tile_hybrid_arm_2026-10-26.json

**What changes:**
- Run `uv run python scripts/glm_tile_hybrid_arm.py --in-dir data/promcalc/tile_predictions --seeds 10`.
- Write the memo with the stacking verdicts, the orthogonality matrix, and the pre-registered
  prediction quoted from `PREREGISTERED["stacking"]` so prediction and outcome sit side by side.
- State the negative scope explicitly if the result is negative: it closes these four scorers, not
  stacking, and specifically not a channel from a new modality family.
- Record that the 2-channel verdict is unchanged, as the regression evidence.

**Test strategy:**
- Manual: `_paired_hybrid_vs_comparator` median is 0.0339 with CI95 [0.0238, 0.0409] and the top-level
  verdict is still `HYBRID_ADDS` — if either moved, the generalisation changed the answered question
  and the stacking numbers are not trustworthy either.
- Manual: every one of the four combos reports `n = 10`; any short n means Step 2's guard did not fire
  when it should have.
- Manual: the memo quotes `predicted_outcome` verbatim from the artifact rather than restating it.

### Step 6: Correct both published over-claims where a reader of the claim will find them
Files: project_state/glm-umbrella-2026-10-07.md, wiki/tile_reproduction_scope_correction_2026-10-26.md
Depends on: none

**What changes:**
- Write the memo naming **commit `e96b5591`** and **umbrella ledger row 17** as the artifacts carrying
  the over-claims, so the correction is reachable from the claim rather than only filed near it.
- Correction 1: the rebuild is a **substrate-integrity check**, not an independent reproduction, and
  give the mechanism — `encoder.py:354` `torch.manual_seed(seed)`, `:361` the seeded `Generator`, the
  seeded `peak_blocked_val_split`, and `glm_tile_predictions.py:390` threading the loop seed — so the
  conv is deterministic given the substrate and a 4-decimal match was guaranteed. State what it DOES
  establish: two regenerated substrates are effectively identical to the originals, which is
  non-trivial because the promcalc toolchain version is unrecoverable.
- Correction 2: "redundant by construction" is false as stated; `max_tx ↔ decomposition` measures
  **+0.7812** (range 0.745–0.815). Carry the full matrix.
- Record the part that survives, so the correction does not over-swing: `gc ↔ conv` at +0.4593 is
  comparably orthogonal to the working pair at +0.3862 and gc-as-third still lost, so the
  *conclusion* had support — just not the support that was given.
- Append one umbrella Action Log row, class `propose`, pointing at the memo. Commit `e96b5591` is
  pushed and cannot be amended, so the row plus the memo ARE the correction mechanism.

**Test strategy:**
- Manual: `git log -1 --format=%H e96b5591` resolves, so the memo cites a real commit.
- Manual: the four code locations cited in Correction 1 exist at the quoted lines, re-grepped at
  write time rather than copied from this plan.
- Manual: Step 6.5 ledger verification — row count up by exactly one, end marker still appears once.

### Step 7: Record the arm's own executed evidence without killing it
Files: soraya_runs/arms-ledger.json, tests/test_tile_chain_substrate_wall_amendment5.py
Depends on: Step 5

**What changes:**
- Add the Step-5 result to `orthogonal-modality-stack-promcalc-plus-N-channels` as `evidence`. The arm
  now has its OWN executed evidence for the first time, which is the only thing R6 accepts.
- **Status handling depends on the outcome and neither branch is a kill.** If stacking fails: narrow
  the arm's scope to "a channel from a NEW modality family" and keep it `planned`, because the four
  scorers here span only two families and that form is untested. If stacking succeeds: `measure_arm`
  with the metric.
- Edit as raw JSON — `Arm.to_dict()` models five fields while the stored arm carries `name`,
  `cost_tier`, `rationale`, `kill_option_withdrawn`, `verdict_consequence` and the 2026-10-23
  substrate note, so `ArmsLedger.save()` would silently drop them.
- Update `tests/test_tile_chain_substrate_wall_amendment5.py`, which pins `evidence == ""` and will
  fail once evidence exists. Replace the empty-string assertion with one that requires the evidence to
  name the Step-5 artifact, and keep the `blocking()` assertion so the ceiling criterion is still named.

**Test strategy:**
- Unit: the amended amendment5 test passes and `ArmsLedger(cwd).blocking()` still names the ceiling
  criterion — asserted in the same test so the two cannot drift apart.
- Unit: a preserved-key assertion confirms all six unmodelled keys survive the raw-JSON edit.
- Manual: `git diff --numstat` on `arms-ledger.json` is a small additive diff, not a reformat.

## Execution Preview

Waves (topological, from the `Depends on` edges):

| wave | steps | parallelism |
|---|---|---|
| 0 | Step 1, Step 6 | 2 |
| 1 | Step 2 | 1 |
| 2 | Step 3 | 1 |
| 3 | Step 4 | 1 |
| 4 | Step 5 | 1 |
| 5 | Step 7 | 1 |

Total waves: 6. Max parallelism: 2. Critical path: 1 → 2 → 3 → 4 → 5 → 7 (6 of 7 steps).

Intra-wave file overlap: none. Wave 0's two steps touch
`{scripts/glm_tile_hybrid_arm.py, tests/test_glm_tile_hybrid_arm.py}` and
`{project_state/…, wiki/tile_reproduction_scope_correction_2026-10-26.md}`.

**Recommended execution mode: SEQUENTIAL.** The chain is long because Steps 2, 3 and 4 all edit
`scripts/glm_tile_hybrid_arm.py` and are serialised by explicit dependencies. That ordering is
deliberate rather than incidental: the guards and paired statistics must be generic and verified
BEFORE a combo depends on them, or a missing statistic and a real negative produce the same empty
output and cannot be told apart. The 2-way parallelism in wave 0 saves nothing worth splitting a
worktree for.

## Risk Flags

1. **The stacking bar is new, and a new bar on an old question is how results get retro-fitted.** The
   mitigation is structural: it is an additive key, the three frozen parameters are re-asserted
   unchanged in the same test that pins the new bar, and `predicted_outcome` is written before the run.
   The residual is that the 0.02 / 8-of-10 shape is inherited rather than derived from this question's
   own effect sizes — defensible because it keeps the two results comparable, but it is an inheritance.
2. **Steps 2–4 touch the script that produced a published artifact.** Verification Signal 1 exists for
   exactly this: if `_paired_hybrid_vs_comparator` moves off 0.0339 / [0.0238, 0.0409], the
   generalisation changed the answered question and everything downstream is void.
3. **A negative is the predicted outcome, and a predicted negative is easy to under-test.** The
   temptation is to accept the first negative as confirmation of the prediction. Step 5's checks are
   therefore about mechanism (all four combos at `n = 10`, 2-channel unmoved) rather than about the
   sign of the result.
4. **The short-arm guard cannot fire on real data** — all 10 tables carry all four channels, verified
   per-seed. It is tested by mutation only, so its real-world behaviour is `[inferred]`, not
   `[grounded]`.
5. **Pre-existing, now relevant: `verdict()` is single-question by design.** Even after this plan the
   top-level `verdict` field reflects the 2-channel question only, with stacking as a secondary
   verdict. A reader who reads `verdict` alone still learns nothing about stacking. Changing the
   top-level field's meaning would break `test_verdict_*` and the published artifacts, so it is
   deliberately not done — but it means the markdown render is load-bearing for the stacking result.
6. **Commit `e96b5591` cannot be corrected in place.** It is pushed, so the over-claim stays in the
   git record verbatim and the correction is additive only. Anyone reading the commit message without
   the ledger or the memo still sees "independent reproduction".

## Open Questions

1. **Should the 4-channel combo be included, or does it muddy a clean 3-channel test?** Included here
   because all four channels are on disk so it is free, and because a 3-channel negative alone cannot
   distinguish "a third channel never helps here" from "this specific third channel does not". Against:
   two new combos means two new numbers against one pre-registered bar, and the 4-channel has no
   independent rationale beyond being available.
2. **Is inheriting the 0.02 margin right for the stacking question?** The 2-channel gained 0.0339 over
   max_tx, so 0.02 asks a third channel for roughly 60% of that again. A margin derived from this
   question's own effect sizes would be more principled but there are none yet — deriving it after
   seeing the result is exactly the retro-fit this plan is trying not to do.
3. **If stacking fails, does the arm's narrowed form have a candidate at all?** The only declared
   new-family candidate is `promoteratlas`, which carries an unresolved circularity question. So a
   negative here may leave the arm `planned` with no runnable successor — honest, but worth naming
   before it happens rather than after.

## Verification

- `uv run pytest tests/test_glm_tile_hybrid_arm.py tests/test_glm_tile_predictions.py -q` — 73 passed
  rises by the tests added, no pre-existing test flips.
- `uv run python scripts/glm_tile_hybrid_arm.py --in-dir data/promcalc/tile_predictions --seeds 10` —
  top-level verdict `HYBRID_ADDS`, `_paired_hybrid_vs_comparator` median 0.0339 CI95 [0.0238, 0.0409],
  all four combos at `n = 10`.
- `orthogonality.median` still 0.3862 and the matrix's conv↔max_tx entry equals it.
- Every `COMBOS` key has a `_paired_*_vs_comparator`; every superset combo also has a
  `_paired_*_vs_hybrid2`.
- `python -c "...ArmsLedger(cwd).blocking()"` still names the ceiling criterion after Step 7.
- `git diff --quiet` on the frozen five plus `uv.lock` — unchanged; this plan touches none of them.
- `git log -1 e96b5591` resolves, and the correction memo cites it.

## Save-time amendments

Captured at: 2026-10-26
Source: `/save-plan` arguments

**Audit notes only.** `/execute-plan` reads ONLY `## Implementation Steps` for executable
work; nothing in this section is an instruction. Where an amendment describes a STRUCTURAL
change, the executable body still carries the original structure and `/technical-plan` must
be re-run before `/execute-plan` acts on it.

- scope is measure-the-stacking-question plus correct two published over-claims
- the stacking bar is an ADDITIVE PREREGISTERED key rather than an edit, because test_frozen_bar_mirrors_the_decomposition_arm pins margin 0.02 / min_seeds_winning 8 / comparator max_tx and retro-fitting the answered question's bar is the failure being avoided
- predicted_outcome and negative_scope are recorded BEFORE the run since a negative is the prediction
- Steps 2-4 are serialised on one file deliberately so a missing statistic and a real negative cannot produce the same empty output
- the arm is narrowed-not-killed on a negative per R6
- commit e96b5591 is pushed and unamendable so the correction is a memo plus a ledger row that NAME the commit
- whether to include the 4-channel combo, whether inheriting the 0.02 margin is right, and whether the narrowed arm has any runnable successor all left as open questions rather than folded in

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,reuse-missing -->
