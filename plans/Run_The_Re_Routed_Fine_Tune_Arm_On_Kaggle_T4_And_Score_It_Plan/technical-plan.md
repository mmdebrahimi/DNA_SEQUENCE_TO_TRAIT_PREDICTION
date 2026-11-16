# Run The Re-Routed Fine-Tune Arm On Kaggle T4 And Score It — Technical Plan

> Take `finetuned-promoteratlas-encoder-150bp` from `planned` to a graded verdict on the frozen bar, closing the ceiling criterion's last untried arm.

## Lens status
Inputs: conversation
Degradations: codex unavailable (not installed); `repo-index-status.sh` absent on this host, so no repo-index layer was consulted and Step 1 read source directly; `graphify-out/GRAPH_REPORT.md` absent; `project-rules.md` absent, so the Project Rules Check section is omitted; `DESIGN.md` absent and this plan touches no UI; sentrux installed but the change is not structural, so Architecture Impact is omitted.

## Problem Statement

The arm `finetuned-promoteratlas-encoder-150bp` is the last `planned` entry on the ceiling criterion `best achievable prediction of bacterial promoter activity from sequence`, so `ArmsLedger.ceiling_verdict` returns `blocked:untried-arm` and the criterion cannot read `met` no matter what the measured arms say. Everything the arm needs already exists: a frozen bar (`wiki/finetune_arm_acceptance_bar.json`), a trainer (`scripts/finetune_arm_run.py`), a complete scorer (`scripts/finetune_arm_score.py`), and a convergence probe establishing that the curve turns (best epoch 6, `turned: true`).

Its only launch died on 2026-11-04 in an external power loss, not a defect — the third of three unexpected shutdowns in ~42 h (`wiki/finetune_arm_launch_death_2026-11-14.json`). On 2026-11-14 the user re-routed the arm from `local_cpu` to `kaggle_t4`, because the original route's stated ground was retrievability inside a 216 h unattended window that has since ended, not compute. That edit was applied with a diff proving all 19 non-route keys byte-identical, and `scripts/kaggle/build_promoteratlas_finetune_kernel.py` was built and verified against the real cohort.

What remains is to run it and grade it. Step 1 research found one blocking defect in the kernel as built: **a Kaggle script kernel receives no environment variables**, so the kernel's `FT_SMOKE_ONLY`, `FT_WALL_HOURS` and `FT_CKPT_IN` knobs are unreachable on Kaggle and only their defaults can ever apply. Smoke-first sizing and cross-session chaining are therefore both impossible until those values are baked at build time.

## Codebase Context

**The push surface, verified by reading `scripts/kaggle_push_poll.py` rather than inherited from a reference.** Subcommands are `push <code_py> <slug> [--gpu] [--no-internet] [--dataset <slug>]`, `status <slug>`, `pull <slug> <dest>`. `--gpu` sets `machine_shape = "NvidiaTeslaT4"`; the source comments that `enable_gpu` alone provisions a Tesla P100 (CC 6.0) that Kaggle's current torch cannot run. Metadata is `kernel_type: script`, `is_private: True`, `enable_internet` defaulting on. `push` **fails closed**: `kernels_push` does not raise on refusal, so it inspects `.error` and `.ref` and returns 1 — a measured 2026-10-30 failure was `Maximum batch GPU session count of 2 reached` reported as success, after which `status` 404s. `--dataset` attaches at `/kaggle/input/<slug>/`, which is the chaining hook.

**The kernel and its builder.** `scripts/kaggle/build_promoteratlas_finetune_kernel.py` resolves its embed set by AST (21 modules), emits a 189,927-byte kernel, `ast.parse`s the emitted body before writing, and `--self-test` proves that check non-vacuous against a deliberately broken body. The kernel calls `finetune_arm_run.run()`, never `main()` — the route gate lives in `main()` and correctly refuses a non-local route. `PA_REPO` and `MPRA_CACHE` are read by `atlas_embed` and `positional_arm_probe` respectively, so the clone and cohort need no code change; the kernel sets both with `os.environ.setdefault`, which works because the kernel itself sets them before importing.

**The scorer is complete and must not be re-implemented.** `scripts/finetune_arm_score.py` computes the paired bootstrap of per-seed frozen-budget values against the committed atlas per-seed, runs the bar's gates in `CHECK_ORDER` order, derives `ci_percentiles` from the bar's frozen `multiplicity.correction` and refuses if it is anything but `none` with one comparison, and reads the comparator from `data/promoteratlas/local_reconcile_2026-11-01.json` — the PRIMARY — rather than from the heads bar, which is a derived byte-identical copy. It takes `--run-artifact` and `--out`, defaulting to `wiki/finetune_arm_run_2026-11-13.json` and `wiki/finetune_arm_result_2026-11-13.json`. It reads each seed's `frozen_budget.value` from the run artifact rather than re-deriving the early-stopping rule, because two implementations of one rule disagree silently on exactly the seeds that stopped early.

**The phase order is what makes a session cap survivable.** `run()` trains every seed to the frozen 7-epoch budget, then runs the nulls, then extends toward the 38-epoch diagnostic horizon, because the bar makes the shuffled-label null a GATE. A run cut at the wall clock after the first two phases is gradeable; only the ineligible 8-to-38 diagnostic is lost, and the record stamps `BOUNDED_BY_WALL_CLOCK`, never `COMPLETE`.

**Registration is a separate act from measurement, and this is a recorded trap.** `ArmsLedger.measure_arm(criterion, name, metric, evidence)` mutates memory only; `save()` is separate. `wiki/decisions-log.md` records that a first registration attempt was *"a silent no-op whose verdict read unchanged, indistinguishable from a successful write."* `measure_arm` also requires `metric: float`, so an INDETERMINATE outcome has no number to register and must not be given a fabricated one.

### Reusable-Code Survey
- `scripts/kaggle_push_poll.py` — the authenticated push/status/pull helper with the T4 pin and fail-closed push; use as-is, do not write a second pusher.
- `scripts/finetune_arm_score.py` — the complete scorer with the bar's gate chain and paired bootstrap; point it at the pulled artifact via `--run-artifact`.
- `scripts/kaggle/build_promoteratlas_finetune_kernel.py` — the kernel builder with the AST closure and emitted-body validation; extend with build-time flags rather than forking.
- `scripts/kaggle/build_promoteratlas_tile_arm_kernel.py` — the sibling builder whose smoke-then-full ordering and bounded-GEO-retry lessons this kernel already mirrors.
- `scripts/positional_arm_score.py` — owns `gate_split_equality`, `gate_parity`, `gate_shuffled_label_null`; the scorer imports them so a second copy cannot drift.

Searched: `scripts/`, `scripts/kaggle/`, `dna_decode/glm/`, `soraya_runs/arms-ledger.json`, `graphify-out/GRAPH_REPORT.md` (absent), `src/utils`, `src/lib`, `utils`, `lib`, `common`, `helpers` (all absent).

## Pre-Change Baseline

Captured before any push, all re-read from artifacts rather than retyped:

- `ArmsLedger('.').ceiling_verdict(<ceiling criterion>)` returns `blocked:untried-arm`, with 14 arms `measured`, 1 `planned` (`finetuned-promoteratlas-encoder-150bp`), 1 `killed`, and `best = asgiven-position-resolved-promoteratlas-150bp`.
- The arm's ledger entry reads `status: planned`, `metric: null`, `evidence: ""`, `resolved_at: ""`.
- Incumbent comparator: frozen median Spearman **0.5076** over **10** seeds at `split_frac` **0.25**, from `data/promoteratlas/local_reconcile_2026-11-01.json`.
- Bar: `route.chosen = kaggle_t4`, `frozen_hyperparameters.n_epochs_budget = 7`, `horizon.diagnostic_horizon_epochs = 38`, `horizon.n_seeds = 10`, `ceiling.ceiling_hours = 129.6`.
- No run artifact exists: `wiki/finetune_arm_run*.json` is empty and `data/promoteratlas/finetune_checkpoints/` is absent.
- Frozen five + `uv.lock` sha256; prospective lock resolves to `wiki/prospective_lock_manifest_2026-08-31.json`.
- Suite baseline: record `uv run pytest tests/ -q` pass/fail/skip counts, expecting the known 12-14 absent-data failures, so Step 8 can assert a strict subset rather than an absolute.

## Verification Signal

- `wiki/finetune_arm_result_<date>.json` exists and its `verdict` is one member of the bar's `OUTCOMES`, reached through `CHECK_ORDER` with every gate's own record present.
- `ceiling_verdict` no longer returns `blocked:untried-arm`; it returns `met`, `blocked:single-arm` or another verdict computed from the ledger, and the arm's entry carries `status: measured` with a float `metric` and a non-empty `evidence` naming the result artifact.
- The run record's `status` is either `COMPLETE` or `BOUNDED_BY_WALL_CLOCK` and the distinction is carried into the memo verbatim, never collapsed to "done".
- Frozen five, `uv.lock` and the bar are byte-unchanged; `contract_number_audit` still reports `NOTHING_TO_ADJUDICATE`; test failures are a subset of the Step 1 baseline set.
- Negative control: if the push is refused or the kernel dies, no result artifact is written and the arm stays `planned` — an absent verdict, never a defaulted one.

## Implementation Steps

### Step 1: Capture the pre-run baseline and verify the push surface is usable
Files: wiki/finetune_arm_prerun_baseline_2026-11-14.json
Depends on: none

**What changes:**
- New read-only artifact recording every Pre-Change Baseline item above, each re-read from its source file at emit time, with a `rule` field stating that a baseline anchored to retyped numbers is anchored to nothing.
- Resolve and record the active prospective lock and the frozen-five sha256 set so Step 8 can diff rather than assert.
- Verify Kaggle auth non-interactively with `uv run kaggle config view` and record `username` plus `auth_method` (expected `OAUTH`). **Do not** attempt a re-auth here: `kaggle auth login --force` is interactive and belongs in the user's own terminal, never routed through a session.
- Choose a FRESH kernel slug (proposal: `pa-finetune-arm-v1`) and record `status <slug>` output. A 404 confirms it is unused; never reuse a slug that holds retrieved output.
- Record the current count of in-flight batch GPU sessions by polling `status` on recently used slugs, because Kaggle permits only 2 concurrent and a third push is refused.

**Test strategy:**
- Unit: assert the emitted artifact's baseline numbers equal values re-read independently from the bar, the reconcile artifact and the arms ledger (not from the artifact itself, which would be circular).
- Manual: confirm the chosen slug's `status` is a 404/absent before Step 3 pushes to it.

### Step 2: Bake the knobs Kaggle cannot supply, because a script kernel gets no environment
Files: scripts/kaggle/build_promoteratlas_finetune_kernel.py, tests/test_finetune_kernel_builder.py
Depends on: none

**What changes:**
- Add `--smoke`, `--wall-hours H` and `--ckpt-dataset SLUG` build flags that BAKE their values into the emitted kernel's module constants, replacing the current `os.environ.get` defaults for `SMOKE_ONLY`, `WALL_HOURS` and `FT_CKPT_IN`. The env reads stay as a fallback so local dry-runs keep working, but the baked literal wins.
- Keep `WALL_HOURS` default at 10.5 h against Kaggle's 12 h GPU session cap, leaving margin for the clone, the GEO download and the artifact write; record the derivation in the emitted kernel's header rather than leaving the number bare.
- `--ckpt-dataset` sets the baked checkpoint-input path to `/kaggle/input/<slug>` so a chained session seeds finished seeds and reads rather than recomputes them.
- New test file asserting: a `--smoke` build's emitted body contains a baked truthy `SMOKE_ONLY` and a non-smoke build does not; `--wall-hours 3.5` appears as a literal; the emitted body still `ast.parse`s under every flag combination; and the AST closure is unchanged by flag choice.

**Test strategy:**
- Unit: build with each flag combination in a tmp dir and assert on the emitted text and `ast.parse` success; assert `--self-test` still passes.
- Non-vacuity: assert that a build WITHOUT `--smoke` fails the smoke assertion, so the test cannot pass for both.

### Step 3: Push and run the smoke kernel to MEASURE the T4 per-epoch cost
Files: wiki/finetune_arm_kaggle_smoke_2026-11-14.json
Depends on: Step 1, Step 2

**What changes:**
- Emit a `--smoke` kernel, push with `--gpu` to the fresh slug, poll `status`, and `pull` its output.
- The smoke measures one real epoch on the provisioned device in a THROWAWAY checkpoint directory and asserts the real directory is untouched; it scores nothing and writes no run artifact.
- Record in the artifact: the device name and compute capability, the measured seconds/epoch, the projected hours for 38 epochs x 10 seeds at that rate, and whether that projection fits one 10.5 h session.
- If the device report shows CC < 7.0, record that a P100 was provisioned despite the T4 pin and treat the measured rate as a P100 rate, not a T4 rate.

**Test strategy:**
- Manual: `status` reaches `complete`; the pulled log shows a non-zero measured seconds/epoch and the smoke guard's "real ckpt dir untouched" line.
- Unit: assert the emitted artifact's `projected_hours_at_horizon` equals `spe * 38 * 10 / 3600` recomputed independently.

### Step 4: Size and run the full arm, then pull the run artifact
Files: wiki/finetune_arm_run_kaggle_2026-11-14.json
Depends on: Step 3

**What changes:**
- Decide the horizon from the Step 3 measurement rather than from the local figure: if the projection fits one session, build without `--smoke` at the bar's full 38-epoch horizon; if it does not, run the frozen-budget and null phases first and chain the extension with `--ckpt-dataset`, which the phase order makes safe.
- Push, poll to completion, pull, and commit the pulled run artifact under a dated name.
- Record `status` (`COMPLETE` vs `BOUNDED_BY_WALL_CLOCK`) and, if bounded, exactly which phase was reached.
- If chaining is needed, upload the pulled checkpoint directory as a Kaggle dataset and attach it on the next push; the resume property holds by construction because each epoch's shuffle is seeded from `(seed, epoch)` rather than the accumulated global RNG.

**Test strategy:**
- Unit: assert the pulled artifact's per-seed count equals the bar's `horizon.n_seeds` and that every seed carries a `frozen_budget.value`.
- Manual: confirm `status` is one of the two accepted values and that `bounded_by_wall_clock` agrees with it.

### Step 5: Score the pulled artifact against the frozen bar
Files: wiki/finetune_arm_result_2026-11-14.json
Depends on: Step 4

**What changes:**
- Run `uv run python scripts/finetune_arm_score.py --run-artifact wiki/finetune_arm_run_kaggle_2026-11-14.json --out wiki/finetune_arm_result_2026-11-14.json`.
- Change no scoring logic. The bootstrap, the gate chain, the `ci_percentiles` derivation and the comparator source are all already frozen in the scorer; touching them after a result exists is the move the bar exists to prevent.
- Record the verdict and every gate's record, including any `INDETERMINATE_*` and which gate produced it.

**Test strategy:**
- Unit: the existing scorer tests must pass unchanged; assert the emitted result's `comparator` equals the reconcile artifact's values re-read independently.
- Negative control: assert a deliberately truncated run artifact yields `INDETERMINATE_RUN_DIED` rather than a substantive verdict.

### Step 6: Register the arm in the arms ledger, with save() proven to have landed
Files: soraya_runs/arms-ledger.json
Depends on: Step 5

**What changes:**
- Call `measure_arm(criterion, "finetuned-promoteratlas-encoder-150bp", metric, evidence)` then `save()`, with `evidence` naming the result artifact, the run artifact, the verdict and the bounded/complete status.
- **Re-read the ledger from disk afterwards and assert the entry's `status` is `measured` and its `metric` equals what was passed.** `measure_arm` mutates memory only, and the decisions log records a prior registration that was a silent no-op indistinguishable from success.
- Record the new `ceiling_verdict` before and after.
- If the verdict is `INDETERMINATE_*`, do NOT register a metric: `measure_arm` requires a float and fabricating one to clear `blocked:untried-arm` would convert an honest non-result into a fake measurement. Leave the arm `planned` and record why.

**Test strategy:**
- Unit: assert the on-disk ledger changed, by comparing sha256 before and after; assert `ceiling_verdict` recomputed from the re-read ledger.
- Non-vacuity: assert that skipping `save()` leaves the on-disk sha256 unchanged, pinning the trap rather than describing it.

### Step 7: Record the result on both decision surfaces
Files: wiki/finetune_arm_result_2026-11-14.md, project_state/glm-learned-representation-2026-10-07.md, NEXT.md
Depends on: Step 6

**What changes:**
- New memo stating the verdict, the measured interval, the comparator, the device actually provisioned, the `COMPLETE` vs `BOUNDED_BY_WALL_CLOCK` status, and what the arm does NOT establish — in particular that failing to beat 0.5076 leaves the budget-limited and substrate-limited readings open, which is the distinction the bar names.
- Append an Action-Log row to the family ledger, anchored on the `<!-- project-state:end:action-log -->` end marker rather than on the last row, and verify the numeric order afterwards.
- Revise the `## ANSWERED 2026-11-14` block in `NEXT.md` in place to carry the outcome, since that file is explicitly transient.
- Use ASCII minus signs in every number, because `contract_number_audit` extracts decimals with a boundary-guarded regex and U+2212 is invisible to it.

**Test strategy:**
- Unit: `contract_number_audit` reports `NOTHING_TO_ADJUDICATE`; a test asserts every number in the memo appears in the cited artifact.
- Manual: confirm the Action-Log row order reads 26, 27 and that the end marker is still the last line of the block.

### Step 8: Verify nothing frozen moved and close out
Files: wiki/finetune_arm_closeout_2026-11-14.json
Depends on: Step 7

**What changes:**
- Diff the frozen five, `uv.lock` and the bar against the Step 1 sha256 set and assert byte-equality; re-resolve the prospective lock and assert it still returns the 2026-08-31 manifest.
- Run the full suite and assert the failure set is a strict SUBSET of the Step 1 baseline set, not an absolute count.
- Re-run `advance_ranker` and record the family's `next_action` and any warnings, so the portfolio surface reflects the closed arm.

**Test strategy:**
- Unit: the closeout artifact's assertions are executable and the script exits non-zero on any drift.
- Manual: read the closeout artifact and confirm each claim names the file it was checked against.

## Execution Preview

Wave 0 (2 parallel):  Step 1 — pre-run baseline + push surface, Step 2 — bake the Kaggle-unreachable knobs
Wave 1 (1):           Step 3 — smoke push, measure T4 per-epoch cost
Wave 2 (1):           Step 4 — size and run the full arm, pull the artifact
Wave 3 (1):           Step 5 — score against the frozen bar
Wave 4 (1):           Step 6 — register in the arms ledger, prove save() landed
Wave 5 (1):           Step 7 — record on both decision surfaces
Wave 6 (1):           Step 8 — verify nothing frozen moved, close out

Critical path: Step 1 → Step 3 → Step 4 → Step 5 → Step 6 → Step 7 → Step 8 (7 waves)
Max parallelism: 2 agents

(Critical path as computed by the toolkit. Step 3 depends on BOTH Wave-0 steps, so Step 2 is equally
load-bearing for schedule even though the toolkit's path runs through Step 1.)

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — **A Kaggle script kernel receives no environment variables**, so the kernel as built can only ever use its `os.environ.get` defaults on Kaggle: smoke-only sizing and cross-session chaining are both unreachable. Step 2 exists solely to bake them, and until it lands Step 3 cannot measure anything. Verified by reading `scripts/kaggle_push_poll.py`: the metadata it writes carries `id`, `title`, `code_file`, `language`, `kernel_type`, `is_private`, `enable_gpu`, `enable_internet`, `machine_shape` and optional `dataset_sources` — there is no environment field.
- Severity: high — **The T4 per-epoch cost is unmeasured and the local figure does not transfer.** The local 1568.1 s/epoch was measured on 4 CPU threads; carrying it onto a T4 would mis-size the horizon in either direction. Step 3 measures it on the provisioned device before Step 4 commits a session, which is the same discipline that rejected the probe's guessed local rate.
- Severity: medium — **The T4 pin can still yield a P100.** `--gpu` sets `machine_shape: NvidiaTeslaT4`, but the kernel's own device report warns below CC 7.0 because a P100 has been provisioned before. This arm is fp32 either way, so a P100 is slow-but-correct; the risk is a mis-sized horizon, which Step 3's measurement absorbs by measuring whatever device actually arrived.
- Severity: medium — **Only 2 concurrent batch GPU sessions are permitted**, and a third push is REFUSED while reporting success in the raw API. The helper already fails closed on this; Step 1 records in-flight sessions so a refusal is diagnosed as a quota rather than as the recorded wrong-slug trap, which once sent a diagnosis to the wrong place for several minutes.
- Severity: medium — **`measure_arm` mutates memory only and `save()` is separate**, and a missed `save()` is a silent no-op indistinguishable from a successful write; the decisions log records exactly that. Step 6 re-reads the ledger from disk and compares sha256 rather than trusting the call.
- Severity: medium — **A `BOUNDED_BY_WALL_CLOCK` run must never be read as `COMPLETE`.** The phase order makes a bounded run gradeable, which is a strength, but it is also the failure mode: a partial run that grades cleanly invites reporting the arm as finished when only the ineligible 8-to-38 diagnostic was lost. Step 4 and Step 7 both carry the distinction verbatim.
- Severity: medium — **An `INDETERMINATE_*` verdict has no number to register**, and `measure_arm` requires a float. The temptation is to register something to clear `blocked:untried-arm`; Step 6 forbids it explicitly and leaves the arm `planned`.
- Severity: low — File overlaps: none within any wave. Wave 0's two steps touch disjoint files (a new wiki artifact vs the builder plus a new test file). Every later wave holds one step.
- Severity: low — Transitive imports not in `Files:` lists: Step 5 runs the scorer, which imports gates from `scripts/positional_arm_score.py` and reads `data/promoteratlas/local_reconcile_2026-11-01.json`; neither is modified, and the scorer's own `verify()` refuses on drift.
- Severity: low — **Do not smoke into the real checkpoint directory.** The resume contract reads any existing per-seed checkpoint and never recomputes it, so a smoke epoch written there becomes the real run's trajectory prefix. The kernel already guards this and refuses if the real directory changed; the flag work in Step 2 must not weaken that guard.
- Severity: low — External services: Kaggle (auth verified as OAUTH for `emanueleebrahimi`; never route a credential through the session — `kaggle auth login --force` is the user's own terminal), GitHub for the public MIT PromoterAtlas clone, and NCBI GEO for `GSE144621`, whose fetch already carries a bounded 4-attempt retry after a measured `ContentTooShortError`.
- Severity: low — Restructuring applied: Step 2 was given `Depends on: none` after checking its file set is disjoint from Step 1's, which lifts max parallelism from 1 to 2 and shortens the critical path by one wave.
- Severity: low — No code generator or build tool writes any file in this plan; the only generated file is the emitted kernel, which is written to a scratch path and is not committed.

## Open Questions

- If the full run comes back `BOUNDED_BY_WALL_CLOCK` with the frozen-budget and null phases complete, should the arm be registered as `measured` on that gradeable prefix, or held `planned` until the 38-epoch diagnostic finishes in a chained session? The bar makes the prefix gradeable, so registering is defensible, but it closes the ceiling criterion on a run that did not reach its declared horizon.
- `C:\my_skills\.venv` (232 MB, created accidentally by running `uv` with the cwd inside the skills repo) is still present. It cannot be moved, only deleted and recreated, and deletion is a recursive delete that hard-pauses. Confirm deletion, or leave it.
- Should the pulled checkpoint directory be uploaded as a Kaggle dataset unconditionally as cheap insurance against a session loss, or only when Step 3's projection shows chaining is actually required?

## Verification

1. `wiki/finetune_arm_result_2026-11-14.json` exists with a verdict drawn from the bar's `OUTCOMES`, every gate recorded in `CHECK_ORDER` order.
2. `ArmsLedger('.').ceiling_verdict(<criterion>)`, recomputed from the ledger re-read off disk, no longer returns `blocked:untried-arm`; the arm reads `status: measured` with a float metric, or is explicitly still `planned` with a recorded INDETERMINATE reason.
3. The run artifact's `status` and the memo agree on `COMPLETE` versus `BOUNDED_BY_WALL_CLOCK`.
4. Frozen five, `uv.lock` and `wiki/finetune_arm_acceptance_bar.json` byte-identical to the Step 1 baseline; prospective lock still resolves to `wiki/prospective_lock_manifest_2026-08-31.json`.
5. `uv run python scripts/contract_number_audit.py` reports `NOTHING_TO_ADJUDICATE`.
6. Full-suite failures are a strict subset of the Step 1 baseline failure set.
7. `advance_ranker` reports the family's `next_action` with no new warnings beyond the two standing ones (ambition concentration, undeclared document validity).

## Save-time amendments

Captured at: 2026-11-14
Source: `/save-plan` arguments

Audit notes only. `/execute-plan` reads executable work ONLY from `## Implementation Steps`; nothing in this section is an instruction.

- Route re-decided to kaggle_t4 on user authority (neutrality proven, 19 non-route keys byte-identical)
- kernel built and verified against the real cohort
- Step 2 added after finding that Kaggle script kernels receive no environment variables

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,severity-high -->
