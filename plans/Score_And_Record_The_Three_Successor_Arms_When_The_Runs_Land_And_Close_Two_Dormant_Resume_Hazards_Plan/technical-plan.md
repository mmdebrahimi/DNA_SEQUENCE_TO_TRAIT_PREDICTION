# Score And Record The Three Successor Arms When The Runs Land, And Close Two Dormant Resume Hazards

## Lens status
Mode: conversation-driven (no argument)
Repo-index: absent — source-only reading, no degradation
project-rules.md: absent — Project Rules Check omitted
DESIGN.md: absent — no UI in scope
sentrux: present · gh: present · codex: ABSENT
Baseline: measured live this session (see Pre-Change Baseline)

## Problem Statement

Two detached runs are in flight and the plan that launched them has no runnable work left: its Steps 7 and 8 are blocked on artifacts that do not exist yet. This plan covers what happens when they land — and, more usefully, moves everything that does **not** depend on the runs forward now, while the user is still present.

That split is the plan's main contribution, and it corrects a dependency error in the plan it supersedes. The two-windows plan made `collect_unattended_runs.py` depend on Steps 4 and 6 — the runs themselves. But a *reader* of an artifact does not depend on the artifact existing; only *running* it does. Writing the collector and the fine-tune scorer now, against synthetic fixtures covering every state, converts the on-return work from "write two scripts while attention is scarce" into one invocation.

**Three findings from re-reading the code change the shape of the work.**

First, **the fine-tune run artifact does not contain a verdict, and nothing currently can produce one.** `finetune_arm_run.py` emits per-seed, per-epoch validation curves and a frozen-budget value per seed. The frozen bar classifies on a bootstrap interval (`classify(ci_lo, ci_hi)`), and that interval is nowhere in the run artifact. Someone must compute the paired bootstrap of per-seed frozen-budget values against the committed atlas per-seed and run the bar's five gates in order. The two-windows plan's Step 7 said "collect ... + artifacts + tests" and did not name this; it is a scorer, not a collection.

Second, **the fine-tune bar has no committed per-seed test sizes and the heads bar does.** Both arms split the same 44,106 tiles with `leave_peak_out_indices(frac=0.25, seed=s)` over the same ten seeds — verified: both bars record frac 0.25 and n_seeds 10. So the fine-tune's `split_equality` gate should read the **heads** bar's committed `n_test_per_seed`, because that is the committed reference and byte-identical splits are what make the two arms' numbers comparable at all. Inventing a second reference would let the two arms drift apart silently.

Third, **three arms need recording and they carry three different schemas with three different metric paths.** The as-given arm is already finished and *unrecorded*: `wiki/asgiven_arm_result_2026-11-02.json` reads `ASGIVEN_TIES` with its metric at `rung.median_spearman`, while the ledger still lists it `planned`. The heads result will carry `arm_metric` plus `best_primary.gain_bootstrap`; the fine-tune result will carry whatever Step 2 emits. A collector that assumed one shape would silently read `None` for two of three. The as-given artifact existing now also means the collector can be exercised on **real** data for one arm immediately, not only on fixtures.

**And two dormant resume hazards are live, one of which has now recurred.** `.claude/execute-plan-state/Finish_The_Three_Successor_Arms_….json` has Steps 5, 6, 7 and 8 `pending` and is **not** neutralised, while the two-windows plan explicitly took over exactly those four steps; its `wiki/plans-index.md` entry carries no `**Status:**` field at all. This is the same hazard already closed for that plan's own predecessor by renaming the state file to `.json.superseded` and annotating line 906. Having caught it manually twice, this plan converts it into a guard test rather than fixing it by hand a third time.

Non-goals: re-opening either frozen bar (both are frozen and their runs have started), re-running either detached job, and touching the frozen five or `uv.lock`.

## Codebase Context

**The three result schemas and where each arm's metric lives.**

| Arm | Artifact | Schema | Metric path | Interval |
|---|---|---|---|---|
| as-given | `wiki/asgiven_arm_result_2026-11-02.json` (exists, `ASGIVEN_TIES`) | `asgiven-arm-result-v1` | `rung.median_spearman` | `rung.gain_bootstrap` |
| heads | `wiki/expression_heads_result_2026-11-02.json` (pending, run in flight) | `expression-heads-result-v1` | `arm_metric` | `best_primary.gain_bootstrap` |
| fine-tune | `wiki/finetune_arm_run_2026-11-13.json` (pending) | `finetune-arm-run-v1` | `per_seed[s].frozen_budget.value` | **absent — must be computed** |

**The frozen fine-tune bar** (`wiki/finetune_arm_acceptance_bar.json`, frozen 2026-11-03) supplies everything the scorer needs and nothing it does not: `check_order` is `(run_completed → INDETERMINATE_RUN_DIED, budget_reached → INDETERMINATE_BUDGET_EXHAUSTED, split_equality → INDETERMINATE_SPLIT_MISMATCH, shuffled_label_null → INDETERMINATE_NULL_FAIL, primary_decision → FINETUNE_ADDS|TIES|NO_GAIN)`; `inherited_not_chosen.margin.value` 0.02; `inherited_not_chosen.comparator.median_spearman` 0.5076; `multiplicity.correction` `none` with one declared comparison; `frozen_hyperparameters.n_epochs_budget` 7 and `patience` 2.

**Reusable, verified by reading:**
- `scripts/positional_arm_score.py` — `paired_bootstrap(gains, percentiles=, n_boot=, seed=)` (line 165), `median_of` (157), `leave_peak_out_indices(tiles, frac=, seed=)` (108), `ScoreError`.
- `scripts/freeze_finetune_arm_bar.py` — `classify(ci_lo, ci_hi, margin=)`, `OUTCOMES`, `CHECK_ORDER`, `reachable_outcomes`, `BarRefused`. Import them; do not restate.
- `scripts/expression_heads_score.py` — the gate-dispatch and `verdict_from_artifact` pattern to mirror: `GATE_SEQUENCE` derived from the bar's `CHECK_ORDER`, with structural tests that every ordered gate has an implementation and every implementation has a place in the order.
- `scripts/detached_launch.py` — `classify(job)`, `status_all`, the seven states, `ATTENTION_STATES`.
- `scripts/finetune_arm_run.py` — `frozen_budget_value(curve, frozen_epochs=, patience=)` and `diagnostic_tail`. The scorer reads the artifact's already-computed `frozen_budget` block; it must not recompute the rule from curves, or two implementations of one rule can disagree.
- `scripts/positional_arm_memo.py` — `render(res, bar, ceiling=None)` (line 56), the memo-renderer shape.

**The ledger write is a DIRECT dict write, and this is documented rather than folklore.** `scripts/declare_successor_arms.py`'s own module docstring records why: `ArmsLedger.declare_arm` ends in `Arm(**arms[name])`, `Arm` is a seven-field dataclass with no field filter (unlike `kill_arm`, which filters to `Arm.__dataclass_fields__`), so an arm carrying an annotation key — and all three of these carry `rationale` and `declared_from` — raises `TypeError` on the **return** construction *after* having mutated the dict. The declaration path at line 145 therefore writes `arms[name] = {...}` directly, and the recording path must do the same.

`verdict_now()` in that same module recomputes `ceiling_verdict` through the real skill code from disk (`C:/my_skills/skills/soraya/scripts/arms_ledger.py`). Reuse it; do not re-derive the verdict.

`evidence` must stay EMPTY for a `planned` arm — guarded by `test_every_PLANNED_arm_still_has_an_empty_evidence_string`. Recording a `measured` arm is what makes an executed-evidence string legitimate.

**Launcher state is the gate on reading any number at all.** `detached_launch.classify` returns `COMPLETE` / `FAILED` / `REFUSED_LOCK_HELD` / `RUNNING` / `STALLED` / `DEAD_NO_SENTINEL` / `NOT_STARTED`. Only `COMPLETE` licenses reading an artifact as the arm's result. `REFUSED_LOCK_HELD` means the payload never ran; `DEAD_NO_SENTINEL` is this host's twice-recorded unexplained death.

**The heads block refuses a capped run by design.** `expression_heads_score.load_scalars` raises unless the sidecar reads `COMPLETE`, so a `--limit` run cannot be read as the arm's result. The collector must not defeat that by reading the block directly.

### Reusable-Code Survey

Candidates found and read (no graphify report present; no conventional util dirs):
- `scripts/positional_arm_score.py` — bootstrap, median, split indices. **Reuse directly.**
- `scripts/freeze_finetune_arm_bar.py` — outcomes, check order, `classify`. **Import, never restate.**
- `scripts/expression_heads_score.py` — gate-dispatch + verdict-re-derivation pattern. **Mirror.**
- `scripts/declare_successor_arms.py` — direct-dict-write ledger pattern + `verdict_now`. **Reuse.**
- `scripts/positional_arm_memo.py` — memo renderer. **Mirror.**
- `scripts/detached_launch.py` — job-state reader. **Reuse.**

None — searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/utils`, `src/lib`, `src/common`, `src/helpers`, `utils`, `lib`, `common`, `helpers` (all absent); no existing `scripts/collect*` or `scripts/*unattended*` (the collector is genuinely new).

## Pre-Change Baseline

Measured live this session, not recalled:

- **Arms ledger** (`verdict_now()` from disk): 16 arms — **12 measured, 3 planned, 1 killed**; `ceiling_verdict` = `blocked:untried-arm`; `best` = `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` at 0.5076, which is 54% of the 0.941 replicate-derived ceiling. The three `planned` arms are `asgiven-position-resolved-promoteratlas-150bp`, `shipped-expression-heads-zeroshot-150bp`, `finetuned-promoteratlas-encoder-150bp`.
- **The as-given arm is finished and unrecorded**: artifact reads `ASGIVEN_TIES`, ledger reads `planned`, `metric: None`.
- **Test suite**: 12 failed / 7,793 passed / 201 skipped at `098e8ae1`, measured this session; all 12 failures are absent DATA, none a code defect. 282 tests added since across five new files (`test_durable_clone` 24, `test_detached_launch` 76, `test_freeze_finetune_arm_bar` 67, `test_expression_heads_score` 62, `test_finetune_arm_run` 53), each file green in isolation. **The full suite has NOT been re-run since**, because a torch job holds a model and the documented failure mode is an access violation that reads like a test failure.
- **Runs in flight**: `expression_heads` RUNNING, ~3,300/44,106 tiles at last check, ETA ~18:30 on 2026-11-03; `finetune_arm` waiting on it, then 381 epochs / ~129.5 h against a 129.6 h derived ceiling.
- **Uncommitted**: `LESSONS_LEARNED.md`, `TODOS.md`, `FUTURE_FEATURES.md`, `wiki/plans-index.md` — 96 insertions, 10 deletions.
- **Untracked and unrendered**: `wiki/eukaryotic_arm_result_K562_2026-11-01.json`, verdict `SUBSTRATE_LIMITED`, carrying both a `refused` block and a `provisional_verdict_if_complete`.
- **Resume hazards**: `Finish_The_Three_Successor_Arms_….json` has Steps 5–8 `pending`, not neutralised, index entry has no `**Status:**`. `Two_Unattended_Windows_….json` has Steps 7–8 `blocked_on_runs`, which this plan takes over.

## Verification Signal

The plan has worked when all of these hold:

1. `uv run python scripts/collect_unattended_runs.py` reports, for each of the three arms, its launcher state and either a metric or an explicit refusal — and refuses every arm whose job state is not `COMPLETE`.
2. The fine-tune scorer emits a result whose verdict is re-derivable from the artifact plus the frozen bar alone, and whose `reachable_outcomes()` equals the bar's declared `OUTCOMES`.
3. The ledger shows **15 measured / 0 planned / 1 killed**, every recorded arm carrying a non-empty `evidence` string and a `metric`, with `rationale` and `declared_from` preserved byte-for-byte.
4. `verdict_now()` is recomputed and **reported**, not acted on.
5. A guard test fails if any `.claude/execute-plan-state/*.json` corresponds to a `plans-index` entry marked `superseded`.
6. The K562 memo exists and does not present `provisional_verdict_if_complete` as the verdict.
7. Full suite re-run with no new failures beyond the 12 data-absent ones — run only once no torch job holds a model.

## Implementation Steps

### Step 1: Close both resume hazards and convert the recurrence into a guard
Files: .claude/execute-plan-state/Finish_The_Three_Successor_Arms_The_Expression_Heads_And_Fine_Tune_Halves_Re_Grounded_In_What_The_First_Seven_Steps_Measured_Plan.json, wiki/plans-index.md, tests/test_plan_state_hygiene.py
Depends on: none

**What changes:**
- Rename the `Finish_The_Three…` state file to `.json.superseded`, the same remedy already applied to its predecessor. Renaming PRESERVES the record of Steps 1–4 being done; deleting would not.
- Annotate its `plans-index.md` entry (line 918) with `**Status:** superseded` and `**Superseded by:** Two_Unattended_Windows_…/ (for Steps 5-8 only; Steps 1-4 stand as executed)`, matching the wording already used at line 906.
- Mark the two-windows entry (line 929) `**Status:** as-built`, noting that its Steps 7–8 are carried forward by this plan and that execution diverged on the D: volume.
- Add `tests/test_plan_state_hygiene.py`: a resumable `*.json` state file may not coexist with a `plans-index` entry whose `**Status:**` is `superseded`. The hazard has now been caught by hand twice; a third hand-fix is the wrong response.

**Test strategy:**
- Non-vacuity: the guard must FAIL against the current tree before the rename and PASS after. Assert that explicitly by constructing the pre-fix condition in a fixture, rather than trusting that it passes post-fix.
- Assert the superseded file is still READABLE and still records Steps 1–4 as done, so the guard cannot be satisfied by deletion.
- Assert every `plans-index` entry for a plan with a live state file carries some `**Status:**` field, so a missing status is a failure rather than a pass.

### Step 2: The fine-tune arm's scorer
Files: scripts/finetune_arm_score.py, tests/test_finetune_arm_score.py
Depends on: none

**What changes:**
- Read `wiki/finetune_arm_run_2026-11-13.json` plus the frozen bar, run the bar's four indeterminate gates **in its declared order**, and only then classify.
- `run_completed` refuses a `BOUNDED_BY_WALL_CLOCK` status and, separately, a job whose launcher state is `DEAD_NO_SENTINEL` — distinct causes, distinct outcomes.
- `budget_reached` requires every seed to have recorded at least `n_epochs_budget` epochs. Early stopping INSIDE the budget is not a shortfall: the bar records the budget as a MAXIMUM, so a seed that stopped at epoch 5 of 7 with `stopped_early` true passes.
- `split_equality` reads the **heads** bar's committed `n_test_per_seed` against the run's recorded `n_val` per seed, because the fine-tune bar carries no such reference and byte-identical splits are what make the two arms comparable. Refuse if the two bars' `frac` or `n_seeds` disagree rather than comparing across different designs.
- `shuffled_label_null` applies the bar's ratio gate to the run's recorded null values, and reports the untrained baseline beside it so a high null is diagnosable rather than ambiguous.
- Compute the paired bootstrap of per-seed frozen-budget values minus the committed atlas per-seed via `positional_arm_score.paired_bootstrap`, then `freeze_finetune_arm_bar.classify`. Take `frozen_budget.value` from the run artifact; do NOT recompute the early-stopping rule from the curves.
- Emit `wiki/finetune_arm_result_<date>.json` carrying the ineligible diagnostic tail, marked ineligible, and `verdict_is_rederivable_from_this_artifact`.

**Test strategy:**
- Drive all seven outcomes from synthetic run artifacts, including both indeterminates that only this arm has.
- Assert `reachable_outcomes()` equals the bar's `OUTCOMES`, with the padding and trimming non-vacuity proofs.
- Structural: every gate in the bar's order has an implementation and every implementation has a place in the order — a gate named but unimplemented would raise only after the artifact is read.
- Assert a BOUNDED run cannot produce a decision outcome, and that a run missing one seed's budget returns `INDETERMINATE_BUDGET_EXHAUSTED` rather than scoring the seeds it has.
- Assert a seed that stopped early *inside* the budget PASSES `budget_reached`, since the opposite reading would fail every well-converged run.
- Assert the verdict is re-derivable from the emitted artifact plus the bar alone.
- Assert a disagreement between the two bars' `frac`/`n_seeds` is refused rather than silently compared.

### Step 3: The three-schema collector and job-state roll-up
Files: scripts/collect_unattended_runs.py, tests/test_collect_unattended_runs.py
Depends on: none

**What changes:**
- One read-only roll-up over the three arms: launcher state from `detached_launch.classify`, then the arm's result artifact, then its metric via a per-schema extractor (`asgiven-arm-result-v1` → `rung.median_spearman`; `expression-heads-result-v1` → `arm_metric`; `finetune-arm-result-v1` → Step 2's field).
- A per-arm REFUSAL rather than a null metric when the launcher state is not `COMPLETE`, when the artifact is absent, or when its schema is unrecognised. An unrecognised schema must refuse, not fall through to a default extractor — a default would read `None` and report it as a measured absence of signal.
- Exercise it on REAL data for the as-given arm immediately: its artifact exists and its expected metric is known from `rung.median_spearman`.
- Emit `wiki/unattended_runs_collection_<date>.json` plus a short human-readable table, exit non-zero while any arm is still incomplete so it is usable as an on-return readiness check.

**Test strategy:**
- All seven launcher states, including a fresh heartbeat with a dead pid.
- All three schemas, plus an unrecognised schema asserting REFUSAL rather than a null metric.
- A real-data case: the committed as-given artifact, asserting the extracted metric equals `rung.median_spearman` read independently.
- Assert the collector does not read the heads scalar block directly, so the block's own `COMPLETE` refusal cannot be bypassed.
- Assert exit code is non-zero while any arm is incomplete and zero only when all three are collected.

### Step 4: Render the K562 eukaryotic memo
Files: scripts/eukaryotic_arm_memo.py, tests/test_eukaryotic_arm_memo.py
Depends on: none

**What changes:**
- Render `wiki/eukaryotic_arm_result_K562_2026-11-01.json` (verdict `SUBSTRATE_LIMITED`) to a markdown memo, ordered by required reader action.
- Honesty rail, enforced not documented: the memo must headline `verdict`, never `provisional_verdict_if_complete`. A provisional verdict presented as a verdict is the over-claim this artifact is unusually exposed to, because it carries both fields side by side.
- Surface the `refused` block and the `honest_limits` list verbatim rather than summarising them.

**Test strategy:**
- Assert the rendered headline contains `SUBSTRATE_LIMITED` and that `provisional_verdict_if_complete`'s value does not appear in any heading.
- Assert every entry of `honest_limits` appears in the output, so a limit cannot be dropped by a future edit.
- Render the real committed artifact end to end.

### Step 5: Commit the outstanding documentation, extended with this execution's findings
Files: LESSONS_LEARNED.md, TODOS.md, FUTURE_FEATURES.md, wiki/plans-index.md
Depends on: Step 1

**What changes:**
- Commit the 96 insertions / 10 deletions already written by the earlier `/documentation` pass and left uncommitted.
- Add the three findings that pass has not seen, each one a reusable rule rather than a narrative: a volume recorded as attached in a plan is not attached at execution time and killed a launch in one second; a shuffled-label null need not collapse when the representation already encodes the target, so an untrained baseline is what makes a high null diagnosable; and an ordering constraint between two unattended jobs belongs in the launcher, because an instruction to start one after the other requires someone present.
- Depends on Step 1 because both edit `wiki/plans-index.md`.

**Test strategy:**
- No unit tests — documentation. Verification is a diff review plus `git log` showing the four files committed together with the findings attributed to the execution that produced them.
- Assert CLAUDE.md is NOT touched: it sits at 35,459 words of a 35,460 ceiling and any addition needs a comparable removal, which is an Open Question below, not a step.

### Step 6: Score and collect the two finished runs
Files: wiki/finetune_arm_result_2026-11-09.json, wiki/unattended_runs_collection_2026-11-09.json
Depends on: Step 2, Step 3

**What changes:**
- Confirm both launcher states read `COMPLETE` before reading any number; a `FAILED`, `BOUNDED` or `DEAD_NO_SENTINEL` state routes to diagnosis, not to scoring.
- Run the fine-tune scorer and the collector. The heads result is written by its own detached run and needs no further compute.
- If the fine-tune returned `BOUNDED_BY_WALL_CLOCK`, record which seeds finished and relaunch to resume rather than scoring a partial run — finished seeds are read, not recomputed.

**Test strategy:**
- Covered by Steps 2 and 3's tests; this step is execution against real artifacts.
- Manual verification: the collection's three verdicts each trace to an artifact field a reader can open, and the fine-tune verdict re-derives from its own artifact.

### Step 7: Record the three arms and REPORT the recomputed ceiling verdict
Files: scripts/record_successor_arm_results.py, tests/test_record_successor_arm_results.py, soraya_runs/arms-ledger.json
Depends on: Step 6

**What changes:**
- Idempotent, re-runnable recording of all three arms via DIRECT dict writes, never `declare_arm` — whose `Arm(**arms[name])` return raises `TypeError` on an arm carrying `rationale`/`declared_from`, after having already mutated the dict.
- Set `status` to `measured`, `metric` to the collected value, `resolved_at` to a stamp, and `evidence` to an executed-evidence string naming the artifact. `rationale` and `declared_from` are preserved byte-for-byte.
- REFUSE to record an arm the collector refused, and refuse to overwrite an arm already `measured` with a different metric.
- Recompute `ceiling_verdict` with `verdict_now()` and REPORT it. With no `planned` arm left it reads `met` again, which is a statement about the arm set being exhausted and not about the 0.941 ceiling having been reached — the same authority call that surfaced once already. The step reports and stops.

**Test strategy:**
- Assert a round-trip preserves `rationale` and `declared_from` exactly, and that the recording path never calls `declare_arm`.
- Assert re-running is a no-op, and that a conflicting metric on an already-measured arm is refused rather than overwritten.
- Assert the PLANNED-arm empty-evidence guard still passes on the mutated ledger, and that a `measured` arm's evidence is non-empty.
- Assert the recomputed verdict is reported in the output rather than written as a decision.

### Step 8: Render the three-arm memo and journal it
Files: scripts/successor_arms_memo.py, tests/test_successor_arms_memo.py
Depends on: Step 7

**What changes:**
- One memo over the three arms, ordered by required reader action: decisions needed first, then corrections, then findings.
- Carries each arm's pre-declared honest limits rather than restating them: the observed-before-declared hazard on as-given, the windowing approximation and the translation-endpoint mismatch on the heads arm, and the budget-versus-representation distinction on the fine-tune — now settled empirically by the ineligible diagnostics, or explicitly not settled if the run was bounded.
- States the recomputed `ceiling_verdict` and that accepting it is the user's call.
- Reports each arm's metric against the 0.5076 incumbent and the 0.941 replicate ceiling, never as a bare number.

**Test strategy:**
- Assert every arm's `honest_limits` entries appear in the output.
- Assert no ineligible secondary or diagnostic appears as a headline figure.
- Assert the memo states the margin and the comparator alongside every gain.
- Render from the real collected artifacts end to end.

### Step 9: Final suite, architecture gate, and archive
Files: wiki/plans-index.md
Depends on: Step 8

**What changes:**
- Run the full suite, **only** once `detached_launch.py status` shows no job in `RUNNING` or `STALLED`: the documented failure mode of running it while a torch job holds a model is an access violation that reads like a test failure.
- Compare against the 12 data-absent failures by NAME, not by count — a count match can hide one failure replacing another.
- Run the sentrux architecture compare; its output is advisory and never replaces tests.
- Archive both executed plans and update their index statuses.

**Test strategy:**
- The suite itself is the test. Record the failure SET and assert it is a subset of the recorded 12, naming any addition.
- Manual: confirm no job is live before starting, and record the elapsed time against the 8m37s baseline.

## Execution Preview

Waves, from the dependency graph:

- **Wave 0** (max parallelism 4, all `Depends on: none`): Steps 1, 2, 3, 4 — every one doable NOW, while the runs are in flight. No file overlap: Step 1 touches the state file and `plans-index.md`; Steps 2, 3, 4 each touch their own new script and test.
- **Wave 1**: Steps 5, 6 — Step 5 after Step 1 (both edit `plans-index.md`), Step 6 after Steps 2 and 3.
- **Wave 2**: Step 7.
- **Wave 3**: Step 8.
- **Wave 4**: Step 9.

Critical path: Step 2 or 3 → 6 → 7 → 8 → 9, five steps deep.

**The schedule split is the point.** Wave 0 is hours of cheap work available today. Waves 1 onward need the fine-tune artifact, which lands about 2026-11-09. Nothing in Wave 0 competes for the CPU the runs are using, because none of it trains anything.

## Risk Flags

- **The fine-tune run may return `BOUNDED_BY_WALL_CLOCK` and that is a real possibility, not a formality.** The plan is 129.5 h against a 129.6 h ceiling — 0.1 h of slack. Any slowdown bounds the run. Step 2 refuses to score a bounded run, and Step 6 routes it to a resume instead. Resuming costs wall-clock the window may not have, in which case the arm reports an indeterminate with the seeds it finished, which is the honest outcome.
- **A Windows Update reboot inside ten days is plausible** and nothing restarts the run. Recovery is one relaunch, which resumes from the last completed epoch, but it needs someone present.
- **`ceiling_verdict` will read `met` again after Step 7** and that will be the second time this exact statement has surfaced. It means the arm set is exhausted, not that the 0.941 ceiling was reached — `best` would still sit near 0.5076 unless an arm adds. Step 7 reports it; accepting or re-opening is an authority call.
- **The as-given arm's result was observed before its bar's forecast was recorded**, which the artifact itself records under `not_a_blind_pre_registration`. That limit must ride with its number in Steps 7 and 8; a metric quoted without it over-claims.
- **The fine-tune `split_equality` gate reads the HEADS bar.** This is a deliberate coupling and it is the right one, but it means a future edit to the heads bar would move the fine-tune's gate. Step 2 refuses on a `frac`/`n_seeds` disagreement so the coupling cannot drift silently, and the Open Questions below name it for ratification.
- **The full suite cannot run until both jobs are done**, so Step 9 is genuinely last and the 282 new tests are only verified file-by-file until then.
- **D: remains unattached** and every path in this plan is on C:. If D: returns, nothing here should start using it without the content check that caught the last substitution.

## Open Questions

1. **Should the fine-tune's `split_equality` gate inherit the heads bar's committed per-seed sizes?** The plan does this because both arms split the same tiles with the same frac and seeds and byte-identical splits are what make the numbers comparable — but it couples two frozen bars. The alternative is no split gate at all for this arm, which is weaker. Answerable before Step 2 runs.
2. **CLAUDE.md is at 35,459 words of a 35,460 ceiling.** Two findings from this execution are worth adding — the attached-volume correction and the null-vacuity caveat — and each needs a comparable removal. Which to remove, or whether to skip both, is a user decision; Step 5 deliberately does not touch the file.
3. **If the fine-tune returns bounded, is a partial-seed indeterminate acceptable as the arm's recorded result**, or should the arm stay `planned` until a full run completes? The first keeps the ledger honest about what was measured; the second keeps the arm alive for a later window.

## Verification

- Step 1 is verified by the guard test failing against the pre-fix tree and passing after, with the superseded state file still readable and still recording Steps 1–4.
- Step 2 is verified by all seven outcomes being reachable from synthetic artifacts, by `reachable_outcomes()` equalling the bar's declared set under padding and trimming, and by the verdict re-deriving from the emitted artifact alone.
- Step 3 is verified by all seven launcher states, by an unrecognised schema refusing rather than returning a null metric, and by the real as-given artifact extracting the metric an independent read confirms.
- Step 4 is verified by the rendered memo headlining `SUBSTRATE_LIMITED` and by every honest limit appearing.
- Step 5 is verified by diff review; no code changes.
- Step 6 is verified by both launcher states reading `COMPLETE` before any number is read.
- Step 7 is verified by the ledger reaching 15 measured / 0 planned / 1 killed with annotations preserved, by re-running being a no-op, and by the recomputed verdict being reported rather than acted on.
- Step 8 is verified by every arm's limits and comparators appearing beside its number.
- Step 9 is verified by the failure SET being a subset of the recorded 12, named rather than counted.

## Save-time amendments

Captured at: 2026-11-03
Source: `/save-plan` arguments

**Audit notes only.** `/execute-plan` reads executable work from `## Implementation Steps` alone; nothing in this section is an instruction. If an amendment changed a step's files, dependencies or structure, re-validation would fire and the plan would need `/technical-plan` again before execution.

- **No user decision was captured, and none is attributed here.** The `/save-plan` argument arrived as the unfilled template `"Q1 answered: <your decision on the split-gate coupling>"`. A placeholder is not a value, so it was not written into the plan as though Open Question 1 had been answered.
- **Open Question 1 is CLOSED by analysis rather than by decision, because the alternative it offered turns out to be foreclosed.** The question asked whether the fine-tune scorer's `split_equality` gate should inherit the heads bar's committed per-seed test sizes, and offered "no split gate for this arm" as the weaker alternative. That alternative was never available: the frozen bar declares `split_equality` in its `check_order`, `INDETERMINATE_SPLIT_MISMATCH` is in its `outcomes`, and `freeze_finetune_arm_bar.build_bar` RAISES on `reachable_outcomes() != set(OUTCOMES)`. Dropping the gate would make that outcome unreachable and the frozen bar internally inconsistent.
- **Given the gate must exist, the heads bar is the only committed reference and is provably the right one.** The fine-tune bar carries no `n_test_per_seed`; the heads bar does. Both record `frac` 0.25 and `n_seeds` 10, over the same 44,106 tiles through the same `leave_peak_out_indices`, so the test sets are identical by construction. The plan therefore stands exactly as written, and the coupling risk stays mitigated by Step 2's refusal on a `frac`/`n_seeds` disagreement.
- **This is my analysis and remains overridable, at a stated cost.** Making the fine-tune arm independent of the heads bar would require re-freezing the fine-tune bar, which invalidates it and the run already launched against it.
- **Open Questions 2 and 3 remain genuinely open and are user decisions** — the CLAUDE.md word ceiling at 35,459 of 35,460, and whether a bounded fine-tune run may be recorded as a partial-seed indeterminate or should leave the arm `planned` for a later window.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
