# Reconcile The Three Disagreeing Registries And Re-Scope The Stack Arm

## Lens status

Repo-index: absent — no `.repo-index/`, so source-only reading, no degradation.
Reuse signals: no `graphify-out/`, no conventional util dirs; package tree read directly under `dna_decode/`.
Tools: `gh` present; `sentrux` absent (no architecture gate); the independent reviewer CLI absent.
Optional context: `wiki/` present; no `project-rules.md`, no `DESIGN.md`, no `AUDIT_REPORT.md`.
Baseline tests: 14 failed / 7029 passed / 183 skipped (the 14 are a documented absent-data set).

## Problem Statement

Three registries that are supposed to say what has been measured disagree with each other, and in
one session that disagreement produced three successive recommendations to do work that was already
done.

1. **The substrate census says the fragment size lever was never pulled.** `fragments_LB.size_lever`
   reads *"UNDER-USED on this rung -- 30,000 of 321,123 rows (9.3%) is the largest N any measurement
   has used. The size lever here has never been pulled."* and `largest_n_measured` is `30000`. The
   arms ledger records `fragment-rung-321k-rows-size-lever-LB-M9` as **`status: measured`,
   `metric: 0.185`, resolved 2026-10-18**, with a verdict of `DOES_NOT_GENERALISE` at both 20k and
   150k and a margin over GC that **shrinks** with more data (+0.0288 → +0.0191).
2. **The census asserts a figure for a rung whose measurement provenance is ambiguous.**
   `fragments_M9.best_measured_arm` carries `spearman 0.1313` / `fraction_of_ceiling 0.1656`,
   byte-identical to `fragments_LB` including the caveat string. The cited artifact
   (`wiki/glm_genomewide_oracle_2026-10-07.json`) contains exactly ONE fragment arm set under a
   condition-unqualified key `model_arms_random_fragments` (n_train 24,006 + n_test 5,994 = one
   30,000-row subsample) and normalises in its own `q1_detail` against *"the 0.7928 ceiling"*, which
   is LB's. M9's own ceiling is 0.8014, so `0.1313 / 0.8014 = 0.1638`, not the recorded `0.1656`.
3. **An arms-ledger arm reads `planned` although its named combination was published five days
   before the arm was declared.** `orthogonal-modality-stack-promcalc-plus-N-channels` (declared
   2026-10-19, `cost_tier: cheap`) asks whether the 2-channel hybrid stacks to 3+. The combo
   `hybrid_conv_maxtx_gc` is already in `scripts/glm_tile_hybrid_arm.py::COMBOS` and already measured
   in `wiki/glm_tile_hybrid_arm_2026-10-14.json`: median **0.4217** over 10 seeds against the
   2-channel **0.4254** — it did not stack, it cost −0.0037.

A fourth fact makes this a process problem rather than three typos: **a guard for exactly this class
already exists and cannot see it.** `tests/test_mpra_census_claims.py::test_no_rung_INHERITS_another_rungs_ceiling`
compares the `ceiling_for_a_sequence_model` FIELD across rungs. LB 0.7928 ≠ M9 0.8014, so it passes.
The defect lives in the DERIVED `best_measured_arm.fraction_of_ceiling`, which no test relates back
to its own rung's ceiling.

Scope: reconcile the registries, close the guard gap, and re-scope the stack arm on its own executed
evidence. Non-goals: running any new measurement arm; re-running the fragment rung; killing any arm;
touching the frozen five or `uv.lock`.

## Codebase Context

- `wiki/mpra_substrate_census_2026-10-12.json` — four rungs (`fragments_LB`, `fragments_M9`,
  `peak_tiles`, `sigma70_designed`), each with `n`, `single_replicate_agreement`,
  `ceiling_for_a_sequence_model`, `largest_n_measured`, `size_lever`, `best_measured_arm`. Consumed by
  `scripts/project_status.py`, whose output is the project's mandated orientation surface.
- `soraya_runs/arms-ledger.json` — one ceiling criterion with 11 arms: 7 `measured`, 4 `planned`.
  Tracked in git despite living under a gitignored directory. The stored arm objects carry FIVE keys
  the `arms_ledger.Arm` dataclass does not model (`name`, `cost_tier`, `rationale`,
  `kill_option_withdrawn`, `verdict_consequence`, plus the note added 2026-10-23), so
  `ArmsLedger.save()` would silently drop them — this file must be edited as raw JSON.
- `scripts/glm_tile_hybrid_arm.py` — already generic over channel count: `combine(preds, names)`
  takes a name tuple, `orthogonality(a, b)` is implemented, `shuffled(vec, seed)` is the control, and
  `COMBOS` holds both `("conv","max_tx")` and `("conv","max_tx","gc")`. Its `PREREGISTERED`
  explicitly records that the orthogonality number is **`not_a_gate`** — *"a prediction that a
  combine will fail is not a substitute for measuring it"*.
- `wiki/glm_tile_hybrid_arm_2026-10-14.json` — `verdict: HYBRID_ADDS`, 10 seeds, conv↔max_tx
  orthogonality median 0.3862, `_paired_hybrid3_vs_comparator` median +0.0302 CI95 [0.0226, 0.0351]
  against `_paired_hybrid_vs_comparator` +0.0339 CI95 [0.0238, 0.0409].
- `tests/test_tile_chain_substrate_wall_amendment5.py:167-169` pins the stack arm at
  `status == "planned"`, `metric is None`, `evidence == ""`. Step 4 changes those values, so that
  test is a hard dependency of Step 4 and must change in the same step.
- Per-tile prediction tables are NOT present locally: `glm_tile_predictions.DEFAULT_OUT` is
  `D:/dna_decode_cache/promcalc/tile_predictions`, that volume is unreachable, and no local
  `tile_predictions` directory exists. This is why no step here runs an arm.

### Reusable-Code Survey

- **`tests/test_mpra_census_claims.py`** — reuse its `_census()` loader and its non-vacuity idiom
  (`test_the_ceiling_guard_FIRES_on_the_real_census` asserts the guard fires on the shipped file, not
  only on a fixture). Step 1 adds a sibling test in the same file rather than a new file, so the
  census's guards stay in one place.
- **`scripts/glm_tile_hybrid_arm.py::orthogonality` + `combine`** — already N-channel generic. Nothing
  needs writing for a future stack run; only a third channel is missing.
- **`arms_ledger.ArmsLedger.blocking()`** — reused unchanged as the assertion surface in Step 4; the
  criterion must still be named after the edit.
- Searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/utils`, `src/lib`, `src/common`,
  `src/helpers`, `utils`, `lib`, `common`, `helpers` (all absent), `dna_decode/**` (2 levels),
  `scripts/*.py` for `rank_average|rank-average|hybrid`, `tests/*.py` for
  `mpra_substrate_census|glm_tile_hybrid|glm_tile_predictions`.

## Pre-Change Baseline

- Full suite: **14 failed / 7029 passed / 183 skipped** in 489 s. All 14 failures are absent data
  files (cyp4f2/GeT-RM 7, essentiality-ladder 5, pneumo cps 1, README-quickstart 1) and are the
  documented stable baseline. The invariant to preserve is the failure SET, not the count.
- `tests/test_mpra_census_claims.py`: 6 tests, all passing — including the inheritance guard that
  does not catch the defect in scope here.
- `scripts/project_status.py` currently prints, for `fragments_M9`, `under-used:9.4%  30,000 of
  318,457 rows  ceiling 0.8014` and reports its best arm as 0.1313 / 16.6%.
- `ArmsLedger(cwd).blocking()` returns exactly one `CeilingVerdict` with
  `verdict='blocked:untried-arm'`.
- Frozen units unchanged and verified by `git diff --quiet`: `dna_decode/eval/amr_rules.py`,
  `dna_decode/data/calibrated_amr_rules.json`, `dna_decode/data/mic_tiers.py`,
  `dna_decode/data/shipped_decoder_surface.py`, `dna_decode/eval/cohort_manifest.py`, `uv.lock`.

## Verification Signal

1. **The new ratio guard FAILS before Step 2 and PASSES after it.** This is the plan's primary
   signal and it is non-vacuous by construction rather than by simulation: `0.1313 / 0.8014 = 0.1638`
   against a recorded `0.1656` on `fragments_M9`, while `0.1313 / 0.7928 = 0.1656` on `fragments_LB`
   is consistent. A guard written after the fix could not demonstrate this.
2. `ArmsLedger(cwd).blocking()` still names the ceiling criterion after Step 4. Recording an arm's
   own executed evidence must not settle a ceiling — if `blocking()` empties, Step 4 has killed an
   arm by bookkeeping, which is the failure this plan exists to avoid.
3. `scripts/project_status.py` no longer prints a "size lever never pulled" claim for a rung the
   arms ledger records as measured, and its M9 line no longer quotes a figure normalised against
   LB's ceiling.
4. Failure SET identical to baseline; frozen five and `uv.lock` byte-unchanged.

## Implementation Steps

### Step 1: Add the derived-ratio guard the existing ceiling guard cannot see
Files: tests/test_mpra_census_claims.py
Depends on: none

**What changes:**
- Add `test_every_fraction_of_ceiling_is_consistent_with_ITS_OWN_rungs_ceiling`: for every rung
  carrying a numeric `ceiling_for_a_sequence_model` AND a `best_measured_arm` with both `spearman`
  and `fraction_of_ceiling`, assert `abs(spearman / ceiling - fraction_of_ceiling) <= 0.0005`.
- Docstring states why the sibling guard misses it: that one compares the ceiling FIELD across rungs
  and LB 0.7928 ≠ M9 0.8014, so it passes while the derived ratio is wrong.
- Add a companion non-vacuity assertion in the same test: at least two rungs must supply both
  numbers, otherwise the guard has nothing to discriminate and should say so rather than pass.
- **Write and run this step BEFORE Step 2.** It is expected to FAIL on `fragments_M9` at
  authoring time; that failure is the evidence for Step 2.

**Test strategy:**
- Unit: the new test, run alone, must FAIL naming `fragments_M9` and must NOT name `fragments_LB`.
  A failure that names both rungs means the tolerance is too tight and the guard would fire on
  ordinary 4-decimal rounding.
- Unit: the five pre-existing tests in the file still pass.

### Step 2: Correct the census's M9 attribution and its stale size-lever claim
Files: wiki/mpra_substrate_census_2026-10-12.json
Depends on: Step 1

**What changes:**
- `fragments_M9.best_measured_arm`: replace the copied `spearman` / `fraction_of_ceiling` with an
  explicit `attribution` object recording that the cited artifact carries ONE condition-unqualified
  arm set normalised against LB's ceiling, so whether M9 rows entered the 30,000-row subsample is
  **undetermined from the artifact**. Do NOT substitute a re-normalised `0.1638` — that would assert
  an M9 measurement from an LB one, which is the defect in the other direction.
- `fragments_LB.size_lever` and `fragments_M9.size_lever`: replace "the size lever here has never
  been pulled" with the measured outcome, citing the arms-ledger arm by name and date —
  `DOES_NOT_GENERALISE` at 20k and 150k, margin over GC +0.0288 → +0.0191 (it SHRINKS), against a
  pre-registered 0.05 bar never approached.
- `largest_n_measured` on both fragment rungs: **`30000` -> `150000`**, and nothing else. DECIDED
  2026-10-23: the field's name IS its definition ("the largest N any measurement has used"), the arm's
  own evidence names 20k and 150k endpoints with FULL 321k explicitly RAM-blocked, so 150,000 is the
  literal answer and 321,123 would be false. The earlier draft's idea of recording the endpoint PAIR in
  this field is DROPPED: `scripts/project_status.py:188-189` renders it verbatim as "X of Y rows" with
  `verdict` as a SEPARATE field, so a scalar field gets a scalar and the nuance goes where the consumer
  already looks. The arm's honest limits (2 seeds not 10; 2 of 4 replication runs lost when the external
  drive disconnected; `kmer6` at full 321k RAM-blocked so the asymptote is UNMEASURED) go into
  `size_lever` beside the verdict, never into this number.
- Carry across the arm's real sub-finding rather than dropping it: `kmer6` DOES gain from size
  (+0.0749 at 20k → +0.1172 at 150k), so 4096 features were data-starved — it simply never reaches
  `kmer4`'s 20k value and GC rises alongside.
- Write with `indent=2` and `newline=""`; the file is LF and `Path.write_text` would convert it.

**Test strategy:**
- Unit: the Step-1 guard now PASSES, and the transition from fail to pass is the plan's primary
  verification signal.
- Unit: `test_no_rung_INHERITS_another_rungs_ceiling` and the other four census tests still pass —
  the ceilings themselves are not touched.
- Manual: `git diff --numstat` on the census shows an additive-shaped diff, not a whole-file
  reformat; a removal count in the hundreds means the indent or EOL changed.

### Step 3: Mark the shipped candidate row done so it cannot be replayed
Files: project_state/glm-generator-finetune-2026-10-07.md
Depends on: none

**What changes:**
- Candidate row 1, "Build the disjoint train/eval window split (CPU, NOW)", shipped in commit
  `58146e40`; both `scripts/glm_finetune_window_split.py` and the 101,549-byte
  `wiki/glm_finetune_window_split.json` exist. Append the `candidate_replay.done_marker(...)` token
  for that row so `advance_ranker` stops returning it as `next_action`.
- Do not delete or reorder the row; the marker is the documented mechanism and the row is history.

**Test strategy:**
- Manual: `advance_ranker.rank(...)` on the umbrella no longer reports that text as the family's
  `next_action`, and reports no new `degraded-input` for the family.
- Manual: the two shipped artifacts are confirmed present before the marker is written — a marker on
  work that did not actually land is worse than a replay.

### Step 4: Record the stack arm's OWN executed evidence and re-scope it without killing it
Files: soraya_runs/arms-ledger.json, tests/test_tile_chain_substrate_wall_amendment5.py
Depends on: none

**What changes:**
- `orthogonal-modality-stack-promcalc-plus-N-channels` gains `evidence` citing
  `wiki/glm_tile_hybrid_arm_2026-10-14.json`: `hybrid_conv_maxtx_gc` median 0.4217 over 10 seeds vs
  the 2-channel 0.4254, paired-vs-`max_tx` +0.0302 CI95 [0.0226, 0.0351] against the 2-channel's
  +0.0339 CI95 [0.0238, 0.0409] — a third channel made it WORSE by 0.0037.
- **Status stays `planned` and `metric` stays `null`.** This is R6-compliant narrowing, not a kill:
  the arm's stated intent is 3+ **ORTHOGONAL** channels, and `gc` is not orthogonal to `conv` — the
  conv is substantially a GC detector on this rung (it scores 0.3234 against GC's 0.3112 with a
  paired CI including zero). So the cheap form of the arm is answered and the orthogonal form is not.
- Correct `cost_tier` from `cheap` to `medium`: the surviving path needs a genuinely different
  modality, and the only declared candidate is
  `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` (`planned`, `cost_tier: medium`).
- Record that the measured combo PREDATES the arm's own declaration (artifact 2026-10-14, arm
  declared 2026-10-19T12:32:35Z) — that ordering is the registry-drift evidence, not a side note.
- Update `tests/test_tile_chain_substrate_wall_amendment5.py:167-169`, which pins `evidence == ""`.
  Replace the empty-string assertion with: `status == "planned"`, `metric is None`, `evidence` names
  the 2026-10-14 artifact, and `blocking()` still names the criterion. The test's purpose —
  "a sibling question's negative verdict must not settle this ceiling" — is preserved and
  strengthened, because the arm now has evidence AND is still blocking.
- Edit as raw JSON with `indent=2` and `newline=""`; `ArmsLedger.save()` would drop the six
  unmodelled keys.

**Test strategy:**
- Unit: the amended amendment-5 test passes, and `blocking()` still returns the ceiling criterion —
  verified in the same test so the two cannot drift apart.
- Unit: a preserved-key assertion confirms all six unmodelled keys survive the raw-JSON edit.
- Manual: `git diff --numstat` on `arms-ledger.json` shows a small additive diff.

### Step 5: Guard the contradiction REFERENTIALLY, because both string guards were measured broken
Files: wiki/mpra_substrate_census_2026-10-12.json, tests/test_registry_consistency.py
Depends on: Step 2, Step 4

**What changes:**
- **The two string-level guards this step used to specify are WITHDRAWN, and they were withdrawn by
  measurement rather than by argument.** Both were run against the `git HEAD` state they were designed
  to catch. Guard A (census says "never pulled" while a `measured` arm names the rung) fired on
  **NOTHING** — the rung key is `fragments_LB` while the arm is
  `fragment-rung-321k-rows-size-lever-LB-M9`, so the singular/plural mismatch means the token match
  never lands. Guard B (a `planned` arm whose combo is already published) fired on **ALL FOUR** planned
  arms, including `multiresolution-pyramid` and `substrate-change-eukaryotic`, which have nothing to do
  with the published hybrid combos. One vacuous where one fire was needed, one indiscriminate where one
  fire was needed. The open question's case-AGAINST was right, and more right than it knew: not "these
  get loosened later" but "these do not work now".
- **Replace both with ONE referential-integrity check.** Add a structured `measured_by` field to each
  census rung naming the arms-ledger arm that measured it (e.g. `fragments_LB.measured_by =
  "fragment-rung-321k-rows-size-lever-LB-M9"`). The guard asserts that the named arm RESOLVES in
  `soraya_runs/arms-ledger.json` and carries `status: "measured"`. A reference either resolves or it
  does not — there is no stemming to get wrong, no superset matching, and no tolerance to loosen.
- **A null must be a statement, not a gap**, reusing this test file's own established convention
  (`test_a_rung_declaring_no_ceiling_also_declares_WHY`): any rung whose `size_lever` makes a claim
  about measurement MUST carry either a resolving `measured_by` or an explicit `measured_by: null` plus
  a `measured_by_reason`. That closes the one remaining way a referential guard goes vacuous — nobody
  writing the field at all.

**Test strategy:**
- Unit: `test_every_measured_by_resolves_to_a_measured_arm` passes on the real repo after Step 2.
- Unit: non-vacuity by MUTATION rather than by fixture — point a rung's `measured_by` at a nonexistent
  arm name, and then at a `planned` arm, and assert the guard fails in BOTH cases. A referential guard
  that has not been shown to fail on a dangling reference is decoration.
- Unit: `test_a_rung_claiming_measurement_declares_measured_by_or_says_why` fails when the field is
  simply absent, so the convention cannot decay by omission.
- Manual: confirm the guard could not have passed on the pre-Step-2 census by pointing it at the
  `git show HEAD:` copy — the M9 rung there has no `measured_by` to resolve.

### Step 6: Resolve the PromoterAtlas circularity BEFORE any spend decision, or record it as the blocker
Files: soraya_runs/arms-ledger.json, wiki/promoteratlas_circularity_2026-10-23.md
Depends on: Step 4

**What changes:**
- The repo already records an UNANSWERED prior question about this arm, in two places:
  `wiki/mpra_substrate_census_2026-10-12.json:150` and
  `wiki/sigma70_loader_reconciliation_2026-10-13.json:103` both flag that PromoterAtlas's
  `urtecho2018` head may have been fit on GSE108535, which would make it a CIRCULAR comparator.
  Verified still open: the only two occurrences of `urtecho` anywhere in `wiki/` or `project_state/`
  are those two open questions. **Spending the arm before this check risks a contaminated number that
  has to be retracted, and the check costs a model-card read.**
- Attempt the resolution from the published source (the model card's or paper's training-data section)
  and write `wiki/promoteratlas_circularity_2026-10-23.md` with the verdict and its locator.
- **If it cannot be resolved from a published source, that is a RESULT and not a failure:** stamp
  `circularity: "unresolved"` plus the list of what was checked onto the arm itself, so the blocker
  travels with the arm rather than living in a wiki open-question list nobody greps.
- Record the MERIT judgement on the arm as well, explicitly labelled INFERRED: the 2-channel gain
  paired a DIFFERENT PRIOR (thermodynamics, `max_tx` 0.3968) with a learned sequence model (`conv`
  0.3234) at orthogonality 0.3862, so PromoterAtlas would be a THIRD sequence-learned channel beside
  the `conv` already in the stack — the failure mode `gc` has just demonstrated at 0.4217 against
  0.4254. Every sequence-statistics arm on this rung clusters 0.2691–0.3234 while the one
  thermodynamic arm sits at 0.3968, so the next channel should carry a different PRIOR rather than a
  different neural net. The orthogonality measurement that would settle it needs per-tile vectors on
  the unreachable volume, so this stays an inference.
- **Do NOT kill the arm** (R6): it keeps `status: "planned"` so it cannot evaporate. The output of this
  step is a recorded blocker plus a recorded merit judgement, never a retirement.

**Test strategy:**
- Unit: the arm still reads `status: "planned"` and `ArmsLedger(cwd).blocking()` still names the
  ceiling criterion — this step must not become a kill by paperwork.
- Unit: `circularity` is one of `{resolved_non_circular, resolved_circular, unresolved}`, and an
  `unresolved` value must carry a non-empty `checked` list so that "unresolved" cannot silently mean
  "unexamined".
- Manual: the memo cites a locator for whatever it concludes. A circularity verdict with no source is
  the fabrication hazard this step exists to prevent.

### Step 7: Record the episode in the umbrella ledger
Files: project_state/glm-umbrella-2026-10-07.md
Depends on: Step 2, Step 3, Step 4, Step 5, Step 6

**What changes:**
- One Action Log row, class `edit-local-code`, recording: the three registry disagreements; that
  three successive VOI recommendations in one session targeted already-executed work; that the
  existing inheritance guard compares the ceiling FIELD and therefore cannot see a wrong derived
  RATIO; the fragment size lever's real measured verdict including its margin SHRINKING with N; and
  that the stack arm was narrowed on its own evidence rather than killed.
- Record the reusable lesson in the row rather than only the actions: **a registry that records
  "never measured" is a claim with an expiry date, and nothing was checking it against the registry
  that records measurements.**

**Test strategy:**
- Manual: Step 6.5 verification — re-read the ledger, confirm the row count incremented by exactly
  one and the end marker still appears exactly once.
- Manual: the row's class is one of the seven enum values.

## Execution Preview

Waves (topological, from the `Depends on` edges):

| wave | steps | parallelism |
|---|---|---|
| 0 | Step 1, Step 3, Step 4 | 3 |
| 1 | Step 2, Step 6 | 2 |
| 2 | Step 5 | 1 |
| 3 | Step 7 | 1 |

Total waves: 4. Max parallelism: 3. Critical path: 1 → 2 → 5 → 7 (4 of 7 steps).

Intra-wave file overlap: none. Wave 0's three steps touch `tests/test_mpra_census_claims.py`,
`project_state/glm-generator-finetune-2026-10-07.md`, and
`{soraya_runs/arms-ledger.json, tests/test_tile_chain_substrate_wall_amendment5.py}` respectively.
Wave 1's two steps touch `wiki/mpra_substrate_census_2026-10-12.json` and
`{soraya_runs/arms-ledger.json, wiki/promoteratlas_circularity_2026-10-23.md}` — also disjoint, and
Step 6's dependency on Step 4 exists precisely because both write the arms ledger.

**Recommended execution mode: SEQUENTIAL, overriding the 3-way parallelism.** Step 1 is deliberately
authored to FAIL, and a worktree agent that sees a red test it did not expect has no way to
distinguish "this is the designed signal" from "I broke something". The whole plan is six small edits
and two short test runs; the parallelism saves nothing worth that ambiguity.

## Risk Flags

1. **Step 1's tolerance is the only number in this plan that is asserted rather than derived.**
   `0.0005` is chosen to sit above 4-decimal rounding noise (`±0.00005`) and below the real
   discrepancy (`0.0018`), i.e. roughly an order of magnitude clear on each side. If a future rung
   records `fraction_of_ceiling` at 3 decimals the guard will false-fire. Mitigation is in the step:
   the test must name M9 only, and a failure naming both rungs is the signal that the tolerance is
   wrong.
2. **Step 2 deliberately leaves a hole rather than filling it.** Replacing M9's copied figures with
   an explicit "undetermined" is less satisfying than a number, and `scripts/project_status.py` will
   print less for that rung. Filling it would require re-running the fragment arm per condition,
   which is out of scope and needs the substrate on a reachable path.
3. **Step 4 is the step most likely to be mis-read later as a kill.** The arm gains evidence and a
   worse cost tier while staying `planned`. Verification signal 2 exists because of exactly this: if
   `blocking()` empties, the step did the thing it claims not to do.
4. **RESOLVED, and it was a measured finding rather than a risk.** This flag used to warn that
   Step 5's superset-match was `[inferred]` and might prove brittle. It was run: one guard fired on
   nothing and the other on all four planned arms, so both were withdrawn and replaced by the
   referential check now in Step 5. The residual risk MOVES rather than disappearing — a referential
   guard goes vacuous if nobody writes `measured_by`, which is exactly what the null-plus-reason
   requirement and the absent-field test in Step 5 exist to prevent.
5. **The claim that `gc` is not orthogonal to `conv` is `[inferred]`.** The published artifact
   reports orthogonality only for conv↔max_tx (median 0.3862); no gc↔conv number exists, because the
   per-tile prediction tables are on an unreachable volume. The supporting evidence is indirect: the
   conv scores 0.3234 against GC's 0.3112 with a paired CI including zero, and the repo separately
   records that the oracle converts ρ 0.59 into GC at 13%. Step 4's wording must carry the label.
6. **Pre-existing: `UNDECLARED DOCUMENT VALIDITY 694/694`.** No `wiki/` document carries a
   `doc-status` header, so nothing distinguishes a live claim from a historical record. All three
   disagreements in this plan are instances of what that warning predicts; fixing it repo-wide is out
   of scope.

## Open Questions

**All three open questions in the emitted draft were delegated to the executor and RESOLVED on
2026-10-23 (ledger row 14, commit `b8f372b0`). They are recorded here with their resolutions rather
than deleted, because two of them resolved AGAINST the plan as drafted.**

1. **Should Step 5 exist at all? — RESOLVED: yes, but not as written.** Both proposed guards were run
   against the historical state and both failed (0 fires and 4 fires where 1 was needed). Step 5 is
   rewritten to a referential-integrity check. The mechanism is still warranted because the drift cost
   three recommendations in one session; the string implementation was the part that was wrong.
2. **150k or 321k? — RESOLVED: `150000`.** The field's name is its definition and 321,123 would be
   false. The consumer was read rather than guessed: the value renders verbatim beside a separate
   `verdict` field, so the nuance belongs in `size_lever`, not in a scalar.
3. **Is PromoterAtlas the right successor, and worth its medium cost? — RESOLVED: the question is
   premature, and the blocker is cheap rather than expensive.** A recorded, still-unanswered
   circularity question about that arm must be settled first (now Step 6). On merit it also looks like
   the wrong successor — a third sequence-learned channel beside the `conv` already in the stack — but
   that part is `[inferred]`, not measured, and the arm stays declared either way.

**The one genuinely NEW residual this amendment introduces:** Step 6 may be unable to resolve the
circularity from a published source, in which case it records `unresolved` plus what was checked. That
is a legitimate terminal for the step, but it leaves the spend decision parked rather than answered —
and the parking is now visible on the arm instead of in a wiki list.


## Verification

- `uv run pytest tests/test_mpra_census_claims.py -q` — fails on `fragments_M9` after Step 1, passes
  after Step 2.
- `uv run pytest tests/test_tile_chain_substrate_wall_amendment5.py tests/test_registry_consistency.py -q`
  — all pass after Steps 4 and 5.
- `uv run pytest tests/ -q` — failure SET identical to the 14-item baseline; passed count rises by
  the number of tests added.
- `uv run python scripts/project_status.py` — the fragment rung lines no longer claim an unpulled
  size lever, and M9 no longer reports an LB-normalised fraction.
- `python -c "...ArmsLedger(cwd).blocking()"` — still returns the ceiling criterion.
- `git diff --quiet` on each of the frozen five plus `uv.lock` — all unchanged.
- `git diff --numstat` on `wiki/mpra_substrate_census_2026-10-12.json` and
  `soraya_runs/arms-ledger.json` — additive-shaped, low removal counts, LF preserved.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
