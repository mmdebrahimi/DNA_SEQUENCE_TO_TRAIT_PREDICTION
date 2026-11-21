# Can This Instrument Resolve The Generator Question At All — Technical Plan

> Measure the noise floor the committed grid never recorded, decompose the across-seed variance into the part a seed-role split can remove, and freeze a bar before any re-run number exists.

## Lens status
Inputs: conversation (the `/brainstorm` recommendation of 2026-11-21 and its two corrected numbers)
Degradations: repo-index absent (`scripts/repo-index-status.sh` not installed) — Step 1 used direct file reading; independent second-model review unavailable for the preceding `/brainstorm` (`codex` not installed), so the plan it rests on is single-model; sentrux present but the plan is non-structural, so no Architecture Impact section

## Problem Statement

The GLM generator loses to a 3-mer Markov null at 10 of 10 seeds and the project has no measurement that says whether that question is even resolvable with the current instrument. The preceding review recommended one cheap run to find out. Grounding that recommendation against the committed artifacts found two errors in it and one missing number, which this plan is scoped to fix before any compute is spent.

The missing number: `wiki/glm_generator_paired_sweep_2026-11-20.json` grid cells carry exactly `{model, epoch, seed, both, kmer, null{kmer,both}}`. There is **no self-control ceiling field**, so the per-seed margin `null - noise_floor` — the width of the region a generator could possibly win in — cannot be computed from the committed grid at all. Only the single seed-7 ceiling exists, in a different artifact (`wiki/glm_generator_falsifier_REAL_2026-10-07.json`).

The two errors, both mine, both now measured:

- The review stated the kmer-mode winnable region is 0.0543 wide (null 0.5990 minus floor 0.5447). 0.5990 is the seed-7 null, which is the **maximum of the ten**. Across all 10 grid seeds the kmer null averages **0.5302** (min 0.4623, max 0.5990, SD 0.0493) against a floor of 0.5447 — so on average the null sits **below its own noise floor** and there is no winnable region in kmer mode rather than a narrow one.
- The review applied a +/-0.0470 error bar to the model-versus-null gap. That bar is the **finetuned-versus-zeroshot delta's** (SD 0.0658, SEM 0.0208). The model-versus-null gap is a different quantity and is far better resolved: both-mode gap mean **0.1186**, SD **0.0275**, SEM **0.0087**. The "generator loses to the null" claim was never underpowered; the "fine-tuning changed something" claim is.

Scope is the measurement and the instrument change it implies. Explicitly out of scope: giving the generator a ceiling criterion in `soraya_runs/arms-ledger.json` (blocked on a user decision about what that criterion should optimize), and the two alternative directions the review raised (emit an edit against a template; iterative-refinement sampling).

## Codebase Context

Every seed consumer in the generator falsifier lives in one file. `scripts/glm_generator_falsifier.py` threads its single `--seed` (default 7) into five places: `build_natural(a.length, a.limit_natural, a.seed, kind)` at the corpus build, which does `random.Random(seed).shuffle(wins)` and therefore fixes the window order, the `limit_natural` truncation, the fit/eval halves and — because conditional prompts are taken from the eval windows — the prompts themselves; `UniformGenerator(seed=a.seed)`; `MarkovGenerator(order=..., seed=a.seed)`; `HFGenerator(..., seed=a.seed)`; and `run_falsifier(..., seed=a.seed)`, which forwards to `heldout_auroc(..., seed=seed)` for the discriminator's 50/50 held-out split.

`self_control_ceiling` in `dna_decode/glm/falsifier.py` is already seed-independent in its own right — it uses fixed internal seeds `(0,1,2,3,4)` and `random.Random(1000 + s)` — but it is called on the **eval set**, which the corpus seed selects, so the realised ceiling still moves with the corpus seed. That is exactly the per-seed quantity Step 2 measures.

The discriminator is `heldout_auroc`: `build_features` at `mode` in `{kmer, positional, both}`, a 50/50 seeded split, `StandardScaler` plus `LogisticRegression(max_iter=2000, C=1.0)`, scored by `roc_auc_score` on the held-out half. The null is `MarkovGenerator(order=3)` in `dna_decode/glm/baselines.py`, whose `fit` accumulates counts of (3-base context -> next base) with add-one smoothing. That is the 4-mer conditional distribution, and the default discriminator feature set is `kmer_features(k=4)` — normalised 4-mer frequencies. The null's sufficient statistic and the kmer test statistic are the same object, which is the mechanism behind the sub-floor kmer nulls above.

Geometry is already single-sourced. `geometry_from_args` in the falsifier script is the one expression both the `--compare-to` guard and the written artifact read, and `tests/test_glm_geometry_contract.py` pins that correspondence (`test_geometry_from_args_covers_the_declared_field_set_EXACTLY`, `test_the_ARTIFACT_restates_NONE_of_the_frozen_fields_as_its_own_literal`, `test_BOTH_consumers_read_the_same_function`). `COMPARISON_GEOMETRY_FIELDS` in `dna_decode/glm/knob_contract.py` holds nine fields and **none of them is a seed**, so adding seed arguments does not touch the comparison contract.

The artifact record in the same script spreads `**geometry_from_args(a)` and then adds `natural_sets`, `fit_eval_split`, `device`, `knob_contract`, `sampling`, `results`, `honest_limits`. It records **no seed at all**, so an existing falsifier artifact cannot say which seed produced it. That is a provenance gap independent of the variance question.

Local feasibility was verified during planning, not assumed. `data/cache/refseq/GCF_000005845.2/genome.fna` (4,699,745 bytes) and `annotations.gff3` (2,314,238 bytes) are present on C:. A standalone run through `.venv/Scripts/python.exe` loaded the genome and GFF3 in 1.5 s, extracted 4,340 upstream windows in 0.7 s, and computed `self_control_ceiling` on the 600-window eval half in 2.9 s (kmer) and 0.6 s (both), returning **0.5447** and **0.5865** — reproducing the committed ceilings exactly. Ten seeds across two modes is therefore well under two minutes of CPU on this host, with no GPU, no torch, no network and no D: volume, which keeps it inside the standing rule against long CPU-saturating jobs here.

Test infrastructure: pytest under `uv run pytest tests/`. The GLM falsifier surface is already covered by `tests/test_glm_geometry_contract.py`, `test_glm_mode_contract.py`, `test_glm_knob_probe_path.py`, `test_glm_sweep_kernel_attests_the_knob.py` and `test_glm_finetune_kernel_builder.py` — the last two are the pattern for offline-testing a Kaggle kernel's argv without running it.

### Reusable-Code Survey
- `dna_decode/glm/falsifier.py` — `self_control_ceiling`, `heldout_auroc`, `build_features` are pure and already parameterised by mode/k/bins/seed; the probes call them unchanged rather than reimplementing a discriminator.
- `dna_decode/glm/baselines.py` — `MarkovGenerator` is pure Python with a seeded RNG, which is what makes the Step 3 variance decomposition free of GPU and torch.
- `dna_decode/glm/corpus.py` — `load_fasta`, `extract_upstream_windows`, `window_context` already exist; Step 1 moves the shuffle/limit/split wrapper here rather than writing a second copy.
- `dna_decode/glm/knob_contract.py` — `COMPARISON_GEOMETRY_FIELDS`, `geometry_of`, `assert_geometry_comparable` are the existing drift guard the seed split must leave untouched.
- `scripts/kaggle/glm_paired_sweep_kernel.py` — the existing paired-sweep kernel is the template for Step 6's argv construction and its pinned T4 machine shape.
- None additionally — searched: `src/`, `lib/`, `utils/`, `common/`, `helpers/` (absent in this repo), `graphify-out/GRAPH_REPORT.md` (absent).

## Pre-Change Baseline

Captured from committed artifacts before any change, and reproduced live during planning:

- Per-seed kmer Markov null across the 10 grid seeds: 0.5009, 0.5022, 0.4834, 0.5778, 0.5618, 0.5490, 0.4623, 0.5990, 0.5830, 0.4823 — mean **0.5302**, SD **0.0493**.
- Per-seed both-mode Markov null: mean **0.6383**, SD **0.0410**, min 0.5674, max 0.7024.
- Zeroshot-minus-null gap, kmer: mean **0.1785**, SD **0.0662**, model beats null **0 of 10**.
- Zeroshot-minus-null gap, both: mean **0.1186**, SD **0.0275**, SEM **0.0087**, model beats null **0 of 10**.
- Finetuned-minus-zeroshot paired delta, both: mean **0.04576**, SD **0.0658**, SEM **0.0208**, range (the artifact's `spread_delta`) **0.2072**, sign flipped at the anchor seed.
- Seed-7 noise floors, reproduced live: kmer **0.5447** (values 0.5430, 0.5447, 0.4698, 0.4760, 0.5213), both **0.5865** (values 0.5348, 0.5455, 0.4548, 0.5053, 0.5865).
- Per-seed noise floors for the other nine seeds: **do not exist** — the grid cell schema has no ceiling field. This is the gap Step 2 closes.
- Falsifier artifact seed provenance: **absent** — the record contains no seed key.
- Suite state before the change: run `uv run pytest tests/ -q` and record the pass/fail/skip counts. The last recorded figure is 13 failed / 7,238 passed / 180 skipped, of which 12 failures are absent data and 1 is a deliberately-unpinned registry count; treat all 13 as pre-existing.

## Verification Signal

- `wiki/glm_seed_role_ceilings_2026-11-21.json` exists and carries 10 seeds x 2 modes of ceiling, with a per-seed `margin = null - ceiling` joined from the committed grid, and its seed-7 row equals 0.5447 / 0.5865 exactly. The decisive observable is the **sign** of the per-seed kmer margins: if most are at or below zero, kmer mode cannot host a generator win and is retired as a scoreboard by measurement.
- `wiki/glm_seed_role_variance_2026-11-21.json` reports three SDs per mode (all-coupled, generation-only, corpus-only) over 10 replicates each, and the generation-only SD is strictly comparable to the all-coupled SD because both are computed on the same candidate by the same discriminator.
- `scripts/glm_generator_falsifier.py --seed 7` with no other seed flag produces an artifact whose `seeds` block reads `{master: 7, corpus: 7, gen: 7, split: 7, decoupled: false}`, and `tests/test_glm_geometry_contract.py` passes **unmodified** — the proof that the split changed no committed number.
- `uv run pytest tests/ -q` shows no new failures against the Pre-Change Baseline count.
- `wiki/glm_seed_decoupling_bar_2026-11-21.json` exists and is dated before the first re-run artifact, with every threshold traceable to a Step 2 or Step 3 number rather than a chosen literal.

## Implementation Steps

### Step 1: Single-source the corpus builder so the probes cannot drift from the runner
Files: dna_decode/glm/corpus.py, scripts/glm_generator_falsifier.py
Depends on: none

**What changes:**
- Move the body of `build_natural` out of `scripts/glm_generator_falsifier.py` into `dna_decode/glm/corpus.py` as `build_natural_set(reference_dir, *, length, limit, seed, kind) -> (windows, stats)`, keeping the exact current semantics: `load_fasta`, then `random_genomic_windows` for `kind="genomic"` or `parse_gff3` plus `extract_upstream_windows` for `kind="promoter"`, then `random.Random(seed).shuffle(wins)`, then `wins[:limit]`.
- Reduce the script's `build_natural` to a thin delegation so the runner and the new probes read ONE definition. Two parallel descriptions of one corpus build is the precise hazard `geometry_from_args` was created to kill after the two descriptions diverged on a single field; a second copy written for the probes would reintroduce it.
- Export the fit/eval halving too (`half = len(wins) // 2`), since the probes need the same eval half the runner scores on, and that split currently lives inline in `main`.
- No behaviour change. The function moves; the call does the same thing with the same arguments.

**Test strategy:**
- unit: `build_natural_set` on the local reference at `length=150, limit=1200, seed=7, kind="promoter"` returns 1,200 windows whose eval half reproduces the committed seed-7 ceilings 0.5447 / 0.5865 — the same assertion the planning feasibility check already satisfied.
- unit: the script's `build_natural` and `corpus.build_natural_set` return identical window sequences for the same arguments, so the delegation is proved rather than assumed.
- regression: the existing geometry-contract and mode-contract suites pass with their own files UNMODIFIED. Those files are deliberately absent from this step's `Files:` list — not editing them is the contract, and the exact paths plus a `git diff` emptiness check are named in Verification.

### Step 2: Measure the per-seed noise floor the committed grid never recorded
Files: scripts/glm_seed_role_ceilings.py, tests/test_glm_seed_role_ceilings.py
Depends on: Step 1

**What changes:**
- New read-only script. For each of the 10 seeds in `wiki/glm_generator_paired_sweep_2026-11-20.json`, call `corpus.build_natural_set(..., seed=s)` at the frozen geometry (length 150, limit 1200, promoter), take the eval half, and compute `self_control_ceiling` in both `kmer` (k=4) and `both` (k=4, bins=6) modes.
- Join against the committed grid to emit, per seed per mode, `{seed, null, ceiling, margin}` where `margin = null - ceiling`, plus the mean and SD of each column.
- A margin may be **negative** and must be reported as such. A null below its own noise floor is the finding, not an error, so nothing clamps at zero.
- Refuse before writing: assert the seed-7 ceilings equal 0.5447 and 0.5865. A probe that cannot reproduce the one committed ceiling has no standing to report nine new ones.
- Record in the artifact the measured reason this number was missing: the grid cell schema is `{model, epoch, seed, both, kmer, null{kmer,both}}` with no ceiling field.
- Write `wiki/glm_seed_role_ceilings_2026-11-21.json`.

**Test strategy:**
- unit: the margin function is pure; assert `margin = null - ceiling` and that `ceiling > null` yields a negative margin with no clamping.
- unit: the seed-7 reproduction gate is non-vacuous — perturb the ceiling input by 0.01 and assert the script refuses rather than writing.
- unit: the seed list is READ from the committed grid, not restated as a literal, so a grid with different seeds cannot be silently scored against the wrong ten.

### Step 3: Decompose the across-seed variance into corpus-and-split versus generation
Files: scripts/glm_seed_role_variance.py, tests/test_glm_seed_role_variance.py
Depends on: Step 1

**What changes:**
- New read-only script using the **Markov null as the probe candidate**, which is what makes this free: `MarkovGenerator` is pure Python with a seeded RNG, so there is no torch, no GPU, no model download and no D: dependency, and the run stays inside the standing no-long-CPU-job rule for this host.
- Three arms, 10 replicates each, in both feature modes:
  - **A all-coupled** — one seed drives corpus shuffle, null RNG and discriminator split together. Reproduces the status quo.
  - **B generation-only** — corpus seed and split seed pinned at 7; only the null's generation RNG varies.
  - **C corpus-only** — null generation RNG pinned; corpus seed and split seed vary.
- Report per-arm per-mode SD of distinguishability. B against A is the variance that survives decoupling; C against A is what the corpus and split contribute.
- Record the limit in the artifact, not just the memo: this is measured on the Markov candidate, so it bounds the **instrument's** nuisance variance and says nothing about the HF generator's own sampling variance, which may be larger or smaller.
- Write `wiki/glm_seed_role_variance_2026-11-21.json`.

**Test strategy:**
- unit: the arm wiring is a pure seed-plan function returning `(corpus_seed, gen_seed, split_seed)` per replicate; assert arm B holds corpus and split constant while gen varies, arm C the converse, and arm A moves all three together.
- unit: assert the three arms are not accidentally identical — a decomposition whose arms share one seed plan measures nothing, and that degenerate case must fail loudly.
- unit: assert at least two distinct distinguishability values within each arm, so an arm that silently produced one repeated number cannot report an SD of 0 as a finding.

### Step 4: Split the single `--seed` into three named seed roles
Files: scripts/glm_generator_falsifier.py
Depends on: Step 1

**What changes:**
- Keep `--seed` as the master. Add three optional arguments, each defaulting to the master when omitted: `--corpus-seed` (window shuffle, limit truncation, fit/eval halves and therefore the conditional prompts), `--gen-seed` (`HFGenerator`, `MarkovGenerator`, `UniformGenerator`), `--split-seed` (`run_falsifier`, which forwards to the discriminator's held-out split).
- Rewire all five confirmed consumers accordingly. These are all of them: `self_control_ceiling` already uses its own fixed internal seeds, so no sixth site is hiding.
- Add an additive `"seeds": {"master", "corpus", "gen", "split", "decoupled"}` block to the artifact record, where `decoupled` is true iff any derived seed differs from the master. The record currently stores no seed, so this closes a provenance gap that exists whatever Steps 2 and 3 find — which is why this step does not depend on them.
- Do NOT touch `COMPARISON_GEOMETRY_FIELDS`, `geometry_from_args` or `geometry_of`. Seeds are not geometry; a run at a different seed is still comparable to its baseline, and widening the frozen field set would retire every existing `--compare-to` relationship.
- This step overlaps `scripts/glm_generator_falsifier.py` with Step 1, which is why it is sequenced after it rather than run alongside.

**Test strategy:**
- unit: with only `--seed 7`, all three derived seeds equal 7 and `seeds.decoupled` is false — backward compatibility by construction rather than by comparison.
- unit: `geometry_from_args` still returns exactly `COMPARISON_GEOMETRY_FIELDS`; the new arguments must not leak into the geometry contract.
- unit: each of the three flags, set alone, changes only its own derived seed and flips `decoupled` to true.
- regression: the existing geometry-contract, mode-contract and knob-probe suites pass **unmodified** (exact paths in Verification, which also asserts their `git diff` is empty). Modifying any of them to accommodate this step would mean the step changed a committed contract, so they are deliberately not declared in `Files:`.

### Step 5: Freeze the acceptance bar before any re-run number exists
Files: wiki/glm_seed_decoupling_bar_2026-11-21.json, tests/test_glm_seed_decoupling_bar.py
Depends on: Step 2, Step 3

**What changes:**
- A pre-registration whose thresholds DERIVE from Steps 2 and 3 rather than being chosen: the resolvable-effect floor is the measured arm-B SD, and the width it must resolve is the measured per-seed margin.
- Name every outcome before any re-run number exists, ordered, first match wins:
  - `NO_WINNABLE_REGION` — the per-seed margin is at or below zero in a mode, so no generator can beat the null there and that mode is retired as a scoreboard.
  - `DECOUPLING_HELPS` — arm-B SD is below a stated fraction of arm-A SD AND the implied CI half-width at the stated replicate count falls below the measured margin.
  - `DECOUPLING_INERT` — arm-B SD is not materially below arm-A SD, so the seed split buys no power and no re-run is warranted at any replicate count.
  - `INDETERMINATE` — the arms did not separate, or a gate failed.
- State the replicate count the bar requires as a FUNCTION of the measured SD and margin, not a literal.
- Record explicitly that the bar licenses no claim about the HF generator's own sampling variance, carrying Step 3's limit forward so a later reader cannot widen it.
- Record the per-mode asymmetry from the baseline (kmer null mean 0.5302 versus floor 0.5447; both-mode null mean 0.6383 versus floor 0.5865) so the bar cannot be read mode-blind.

**Test strategy:**
- unit: the bar parses as JSON (guarded by the standing `test_every_fba_wiki_artifact_is_parseable_json`-style discipline after a hand-edited artifact was once invalid for a day).
- unit: every named outcome is reachable by some synthetic input, and no two can fire on one input — an unreachable verdict is dead prose and an ambiguous one is not a pre-registration.
- unit: assert the replicate count is stored as a derived expression plus its inputs, and that recomputing it from the cited artifacts reproduces the stored value.

### Step 6: Build the decoupled re-run kernel, offline-tested, execution gated
Files: scripts/kaggle/glm_seed_decoupled_kernel.py, tests/test_glm_seed_decoupled_kernel.py
Depends on: Step 4, Step 5

**What changes:**
- A Kaggle kernel re-running the zeroshot-versus-finetuned paired comparison with `--corpus-seed 7 --split-seed 7` pinned and `--gen-seed` varying over the replicate count the bar derives, at the frozen geometry (length 150, k 4, bins 6, markov-order 3, n 120, limit-natural 1200, modes `kmer,both`, conditional, context 300) and `machine_shape: NvidiaTeslaT4` — the P100 that bare `enable_gpu` provisions cannot run current Kaggle torch in fp16.
- Modelled on `scripts/kaggle/glm_paired_sweep_kernel.py`, and offline-tested the way `tests/test_glm_sweep_kernel_attests_the_knob.py` and `test_glm_finetune_kernel_builder.py` already test kernels: assert the constructed argv, run nothing.
- Requires a FRESH slug. The standing never-re-push list covers `glm-paired-sweep-v1`, `glm-generator-finetune`, `pa-finetune-arm-v1/v2/v3`, `pa-finetune-ckpt-probe/probe2/probe3` and `pa-finetune-chain-s1/s2/s3/s4`.
- **The kernel is the deliverable; running it is a separate, user-gated step.** It costs an outward session and needs `uv run kaggle auth login --force` in the user's own terminal, so execution is not part of this plan — the same shape as the external-cohort arm, where the code plus offline tests were the shipped artifact and the live run was a manual post-build step.

**Test strategy:**
- unit (offline, no network, no torch): assert the argv carries fixed `--corpus-seed 7` and `--split-seed 7` and a DISTINCT `--gen-seed` per replicate; a kernel that varied the corpus seed would silently reproduce the coupled regime it exists to replace.
- unit: assert every frozen-geometry flag is present with the baseline value, field by field against the nine `COMPARISON_GEOMETRY_FIELDS`-bearing arguments.
- unit: assert the slug is not any member of the never-re-push list, and that `machine_shape` is pinned to the T4.

### Step 7: Write the result memo and reconcile the two corrected numbers
Files: wiki/glm_seed_role_decomposition_2026-11-21.md, tests/test_glm_seed_role_memo.py
Depends on: Step 2, Step 3, Step 5

**What changes:**
- The memo, ordered DECISIONS NEEDED, then CORRECTIONS, then result, then limits.
- CORRECTIONS carries both planning-found errors in the review this plan came from: the kmer target width of 0.0543 was computed from the seed-7 null 0.5990, the maximum of ten, where the mean 0.5302 sits BELOW the 0.5447 floor; and the +/-0.0470 error bar belonged to the model-versus-model delta, not the model-versus-null gap, whose both-mode SEM is 0.0087 and which was never underpowered.
- DECISIONS NEEDED carries the two open tradeoffs unchanged: what the generator criterion should optimize, and whether to spend the gated Kaggle session.
- Limits section carries Step 3's Markov-candidate scope, the single-organism scope, the fact that the re-run has not happened, and the standing prohibition on quoting a raw cross-substrate distinguishability level.

**Test strategy:**
- prose-versus-artifact: assert every decimal in the memo appears in one of the artifacts it cites, matched on whole numbers rather than substrings — the matcher-strictness defect that inflated every rate four- to tenfold is on record.
- unit: assert the memo contains no raw cross-substrate level — the substring-absence shape the existing human-falsifier verdict suite already uses for the bacterial `0.6679` — and that it names the re-run as pending rather than done.

## Execution Preview

Wave 0 (1 parallel):  Step 1 — Single-source the corpus builder
Wave 1 (3 parallel):  Step 2 — Per-seed noise floor, Step 3 — Variance decomposition, Step 4 — Seed-role split
Wave 2 (1 parallel):  Step 5 — Freeze the bar
Wave 3 (2 parallel):  Step 6 — Decoupled re-run kernel, Step 7 — Result memo

Critical path: Step 1 -> Step 2 -> Step 5 -> Step 6 (4 waves)
Max parallelism: 3 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — **The plan rests on a single-model review.** `codex` is not installed on this host and the local probe helper is absent, so the `/brainstorm` that produced the recommendation had no independent critic, and two of its numbers were already wrong. Steps 2 and 3 are designed to be cheap and self-refuting for exactly that reason; re-run `/brainstorm` on this plan once Codex is available.
- Severity: high — **A negative margin would retire a scoreboard, not just a mode.** If Step 2 finds most kmer margins at or below zero, then every kmer-mode generator number ever reported here — including the temperature result at T=0.8 that motivated the original proposal — was measured in a mode where the null is already at the noise floor. That is a wide re-reading of past artifacts and must be stated in Step 7 rather than left implicit.
- Severity: medium — **File overlap between Step 1 and Step 4**, both on `scripts/glm_generator_falsifier.py`. Resolved by sequencing Step 4 after Step 1; they are deliberately not placed in the same wave.
- Severity: medium — **Transitive import not visible in the Files lists.** Steps 2 and 3 import `dna_decode/glm/falsifier.py` and `dna_decode/glm/baselines.py` without modifying them, and Step 4 modifies the script that also imports `falsifier.py`. No step edits `falsifier.py` or `baselines.py`, so the shared import is read-only, but an executor that decides to "fix" a signature in `falsifier.py` would break the wave isolation.
- Severity: medium — **Restructuring applied during planning.** The original recommendation was one step ("re-measure the paired-delta SD with seed roles decoupled"). Grounding split it into the missing-ceiling measurement, the variance decomposition, and the code change, because the ceiling number gates whether the other two are worth doing at all and did not exist. The scope is not widened: no new direction is added and the arms-ledger criterion is explicitly excluded.
- Severity: medium — **External tool surface: Kaggle.** Verified surfaces: `machine_shape: NvidiaTeslaT4` in kernel metadata (verified against the in-repo kernels `scripts/kaggle/glm_paired_sweep_kernel.py` and `glm_finetune_kernel.py`, and against the standing recorded finding that `enable_gpu: true` alone provisions a P100 that current Kaggle torch cannot run in fp16); `uv run kaggle auth login --force` for re-auth in the user's own terminal. `[unverified]` — the Kaggle CLI version on this host was not probed during planning and no live `--help` was run, so the push/status/pull flags are inherited from the in-repo helper rather than confirmed. Step 6 ships only the kernel and its offline tests, so an unverified CLI flag cannot corrupt a result; confirm before the gated run.
- Severity: low — **No code generator or build tool writes any file in this plan.** All seven steps touch hand-maintained Python, tests, and `wiki/` artifacts.
- Severity: low — **Environment dependency.** Steps 2 and 3 need `sklearn` plus the local reference at `data/cache/refseq/GCF_000005845.2/`; both were verified present during planning and need no torch, so `uv sync --extra forward` is not required for Waves 0 through 2. Step 6's tests must stay importable without torch.
- Severity: low — **Suite hygiene.** Do not run the full suite while any torch job holds a model; an access violation there reads like a test failure. No torch job is required by Waves 0 through 2.

## Open Questions

- **What should the generator's ceiling criterion optimize?** Distinguishability from natural sequence is built, cheap and has a derived noise floor, but the repo's own evidence says realism and fitness can diverge (the most-fit TrpB variants carry a substitution nearly absent in natural sequence). Oracle-scored design value is what the north star names, but it makes the predictor half load-bearing out of distribution, and that half is validated on natural sequence at 0.6405 of a 0.941 ceiling. This plan deliberately does not add a generator criterion to `soraya_runs/arms-ledger.json` until this is settled.
- **Spend the gated Kaggle session?** Step 6 ships the kernel; running it is outward, costs a session and needs OAuth in the user's own terminal. Worth deciding only after Step 5 returns `DECOUPLING_HELPS` rather than `DECOUPLING_INERT` or `NO_WINNABLE_REGION`.
- **If kmer mode has no winnable region, is it retired or kept as a reported-but-not-scored diagnostic?** Retiring it changes what every future generator artifact reports; keeping it risks another near-miss being read as near-success.
- **Does the seed-role split warrant re-stamping the existing committed falsifier artifacts** with their (now recoverable) seed provenance, or do they stay as-is with the gap recorded? Re-stamping edits committed evidence, which is an authority call rather than an executor task.

## Verification

1. `uv run python scripts/glm_seed_role_ceilings.py` exits 0 and writes `wiki/glm_seed_role_ceilings_2026-11-21.json` with 10 seeds x 2 modes; confirm by hand that the seed-7 row reads 0.5447 (kmer) and 0.5865 (both), and read the sign of the ten kmer margins.
2. `uv run python scripts/glm_seed_role_variance.py` exits 0 and writes `wiki/glm_seed_role_variance_2026-11-21.json` with three arm SDs per mode; confirm arm A's SD is in the neighbourhood of the committed per-seed null SD (0.0493 kmer, 0.0410 both), which is the sanity check that arm A really did reproduce the status quo.
3. `uv run python scripts/glm_generator_falsifier.py --seed 7 --help` resolves, and a no-model run at `--seed 7` writes an artifact whose `seeds` block reads all-7 with `decoupled: false`.
4. `uv run pytest tests/test_glm_geometry_contract.py tests/test_glm_mode_contract.py tests/test_glm_knob_probe_path.py -q` passes with those files unmodified in `git diff`.
5. `uv run pytest tests/ -q` shows no new failures against the Pre-Change Baseline count.
6. `python -c "import json; json.load(open('wiki/glm_seed_decoupling_bar_2026-11-21.json'))"` succeeds, and the bar file's git timestamp precedes any re-run artifact.
7. `git log --oneline -n 8` shows the step commits, and `git status --short --branch` is clean with the branch ahead of origin by the expected count before pushing.

## Save-time amendments

Captured at: 2026-11-21
Source: `/save-plan` arguments

These are audit notes for human readers only. `/execute-plan` reads ONLY the plan body's `## Implementation Steps` for executable work; nothing in this section is an executable instruction.

- kmer-mode null measured below its own noise floor across 10 seeds
- the underpowered quantity is the model-vs-model delta, not the model-vs-null gap
- scope excludes the arms-ledger generator criterion pending the criterion decision
<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high -->
