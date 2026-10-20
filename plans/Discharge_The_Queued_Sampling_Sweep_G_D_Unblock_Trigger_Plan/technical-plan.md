# Technical plan — Discharge the queued sampling sweep (G-D's unblock trigger)

## Lens status
Repo-index: absent (`repo-index-status.sh` not installed) — source-only reading, no degradation claimed.
Dispatch: conversation-driven, revised after review.
Tools: `gh` present; `sentrux` absent (architecture gate `n/a`); `codex` absent (reviews are single-voice).
project-rules.md: absent → Project Rules Check omitted.
wiki-replay: ran. Relevant HIGH-salience priors — `GLM_Ceiling_Arm_Accounting_Correction_Plan` (a plan's
baseline claims are MEASUREMENTS, taken at authoring time, not recalled) and `GLM_G_A_And_G_C_V2_Build_Plan`
(compute a metric's ceiling before setting a bar on it). `plans-index.md:708` explicitly records that
proposing generator work while `features/glm-generator-finetune/` holds a deferred plan "would be the
repo's recorded orientation error #3" — this plan does NOT propose fine-tuning; it discharges that
family's named predecessor.

## Problem Statement

The sampling sweep is not a new idea. It is a **queued, pre-registered action and the named unblock
trigger** for GLM family G-D, and it has never produced a valid artifact.

- `project_state/glm-umbrella-2026-10-07.md:18` — G-D is `ACCEPTED`, `blocked_by: sampling sweep + 8 GB GPU`.
- `project_state/glm-generator-finetune-2026-10-07.md:46` — short-term action 1: "Sampling sweep reports |
  Distinguishability per temperature/top-k setting, so the levers stay separable | **blocked (server)**".
- `features/glm-generator-finetune/technical-plan.md:21` — "**Unblock trigger:** the sampling sweep
  artifact exists AND the 8 GB host is reachable."
- `:5` — "**If the sweep already clears 0.6679, this family's premise weakens materially** before any
  training starts."

Two things blocked it, and only one still does. The hardware blocker (8 GB GPU mid-relocation) applies to
**fine-tuning**, not to the sweep: the sweep is CPU inference and has been runnable here all along. The
real blocker was never named, because nobody knew it existed — **the sampling knobs were inert**
(fixed 2026-10-20, `wiki/glm_temperature_knob_2026-10-20.md`). A sweep run before that fix could not have
varied sampling, so it could not have discharged the trigger.

My own 2026-10-19 attempt was **doubly invalid**: the knob was dead AND the mode was wrong
(`conditional=False` against the aligned promoter set, which `corpus.py:168-175` and the `--conditional`
help both document as not meaningful).

**Success criteria.** A valid sweep artifact — conditional mode, the recorded geometry, knobs verified
live — reporting distinguishability per sampling setting against the comparator already on record, plus
the pre-registered read on G-D's premise.

**Non-goals.** Fine-tuning (that is G-D itself, deferred by design with its own trigger); building a
conditional Markov null; raising `n` above the recorded 120; touching the oracle, the discriminator or
the model.

## Codebase Context

- `scripts/glm_generator_falsifier.py` — the runner. Verified flags: `--conditional`, `--context` (300),
  `--natural-set` (`promoter|genomic|both`), `--n` (400), `--length` (150), `--limit-natural` (2000),
  `--modes` (`kmer|positional|both`, comma-separated), `--markov-order` (3), `--seed` (7), `--temperature`,
  `--top-k`, `--verify-knob`.
- Conditional path (`:166-178`): prompts from `naturals[kinds[0]]["eval_windows"]` via
  `corpus.window_context`; **only the continuation is scored**; breaks at `len(prompts) >= a.n`, so `n` is
  capped by available eval windows.
- `dna_decode/glm/falsifier.py` — `MIN_PER_CLASS = 40`; `MIN_POWER_AUROC = 0.90` on the uniform control.
- `dna_decode/glm/baselines.py:92` — `MarkovGenerator.generate(n, length)` is **unconditional only**; no
  conditional variant exists. Consequence recorded in Risk Flags.
- `dna_decode/glm/knob_contract.py` — the sampling preconditions (`validate_sampling`,
  `assert_knob_wired`, `verify_knob_live`, `TOP_K_KMER_FLOOR = 22`).
- `dna_decode/glm/corpus.py` — `window_context`, `random_genomic_windows`, `load_fasta`.
- Input: `data/cache/refseq/GCF_000005845.2/{genome.fna,annotations.gff3}`, gitignored; the script refuses
  rather than running on nothing.

### Reusable-Code Survey

This plan adds ONE guard plus artifacts. Everything that computes a number already exists.

- `glm_generator_falsifier.py` conditional path + `--verify-knob` — reuse unchanged.
- `glm/falsifier.py` powering + instrument controls — reuse; do not restate thresholds.
- `glm/knob_contract.py` — extend the SAME module with the mode guard rather than opening a second
  contract surface.
- `glm/corpus.py::window_context`, `glm/baselines.py::MarkovGenerator` — reuse unchanged.
- Consulted (not re-derived): `project_state/glm-generator-finetune-2026-10-07.md` (the confirmed
  geometry + comparator), `features/glm-generator-finetune/technical-plan.md` (the trigger + decision
  rule), `project_state/glm-umbrella-2026-10-07.md` (the blocker edge).
- None — searched: `dna_decode/glm/` (12 modules), `tests/test_glm_*.py` (25 files), `scripts/glm_*.py`,
  `wiki/glm_*`, for an existing mode guard or a committed sweep artifact. Neither exists.

## Pre-Change Baseline

Measured at authoring time from the ledger and the committed artifact, not recalled:

| quantity | value | provenance |
|---|---|---|
| GLM distinguishability, mode `both` | **0.7706** | `glm-generator-finetune-2026-10-07.md:20,80` |
| GLM distinguishability, mode `kmer` | **0.7237** | same + `wiki/glm_generator_falsifier_REAL_2026-10-07.json` |
| Markov-k3 null, `both` / `kmer` | **0.6679 / 0.5990** | same |
| verdict, both modes | `WORSE_THAN_MARKOV_NULL` | same |
| uniform control | 0.9408 | committed artifact |
| geometry | window 150 bp, context 300 bp, markov_order 3, n_generated 120, fit/eval 600/600 | ledger, "ALL CONFIRMED" |
| sampling | default; **the knobs were inert until 2026-10-20** | `wiki/glm_temperature_knob_2026-10-20.md` |
| committed sweep artifact | **none exists** | `git ls-files` + `wiki/` listing |
| my 2026-10-19 sweep | `conditional=False`, dead knob — invalid | scratchpad artifacts |

Two mechanism hypotheses are already REFUTED on record and must not be restated: mode collapse (pairwise
identity 0.249 vs ~0.25 random) and coding-like output (stop density 4.83 vs natural 4.67).

## Verification Signal

The comparator and the decision rule are **already pre-registered** — nothing new is invented here.

- Comparator: this run's own Markov null, with the recorded **0.5990 (`kmer`) / 0.6679 (`both`)** as the
  cross-check that the null reproduced.
- G-D premise read (verbatim from `features/glm-generator-finetune/technical-plan.md:5`): if the sweep
  clears the null, "this family's premise weakens materially".

Validity gates — all must hold or the artifact reports INDETERMINATE and no verdict is quoted:

1. `conditional == true` and `natural_set == promoter`.
2. `knob_verified_end_to_end == true` at every non-default grid point.
3. `control_uniform >= 0.90`; `n_candidate >= 40`; verdict not `INSUFFICIENT_SEQUENCES` /
   `INSTRUMENT_UNDERPOWERED`.
4. **The reproduction anchor passed** (Step 3).

Outcome branches, applied per mode:

- **SAMPLING_CLEARS_NULL** — some setting puts GLM distinguishability below the run's null. G-D's premise
  weakens; record that and stop, because deciding what follows is a scope call.
- **SAMPLING_NARROWS_ONLY** — best setting improves on default but stays above the null.
- **SAMPLING_INERT** — no setting beats default beyond the null's own seed spread.
- **INDETERMINATE** — any gate fails.

## Implementation Steps

### Step 1: Refuse the mode mismatch at the contract layer
Files: dna_decode/glm/knob_contract.py, scripts/glm_generator_falsifier.py, tests/test_glm_mode_contract.py
Depends on: none

**What changes:**
- `validate_comparison_mode(natural_set, conditional)` in the existing contract module:
  - `promoter` + not conditional → REFUSE, citing `corpus.py:168-175` and the `--conditional` help.
  - `genomic` + not conditional → REFUSE with the MEASURED reason: that reference has no discriminating
    power (null 0.5086 against a 0.5255 ceiling — below its own noise floor).
  - `genomic` + conditional → REFUSE (conditional sample against unaligned windows).
  - `both` → REFUSE: the runner builds prompts from `kinds[0]` only, so the second set's candidates would
    be prompted from the first set's windows. (Review resolved this against allowing it.)
  - `promoter` + conditional → the one valid pairing.
- Call it in the runner before any corpus build or model load.

**Test strategy:**
- Unit: all five combinations hit their expected accept/refuse.
- Unit (non-vacuity): the `promoter` + unconditional message names the recorded defect, so a reader learns
  why rather than just being blocked.
- Unit: the valid pairing passes, and the refusal fires before any model load (no network, no weights).

### Step 2: Record that the 2026-10-19 sweep was doubly invalid
Files: wiki/glm_sweep_invalidation_2026-10-20.md
Depends on: none

**What changes:**
- Memo stating both defects with evidence: the knobs were inert (so all three grid points used identical
  sampling), AND the mode was `conditional=False` against the aligned promoter set.
- State what the valid baseline IS: the committed `glm_generator_falsifier_REAL_2026-10-07.json` plus the
  ledger's confirmed geometry — and that its direction (`WORSE_THAN_MARKOV_NULL`) is the project's
  recorded premise, not a new finding.
- State explicitly that no committed sweep artifact exists, so G-D's trigger is genuinely undischarged.

**Test strategy:**
- Unit: the memo names both defects and cites the committed artifact (guards against a memo that asserts
  the conclusion without its evidence).

### Step 3: Reproduction anchor — conditional, default sampling, recorded geometry
Files: wiki/glm_sweep_anchor_2026-10-20.json
Depends on: Step 1

**What changes:**
- Run at the recorded geometry with DEFAULT sampling:
  `--model auto --conditional --context 300 --natural-set promoter --n 120 --length 150
  --limit-natural 1200 --modes kmer,both --markov-order 3`
- `HF_HUB_OFFLINE=1`; serialised (no concurrent torch process).
- **Pass condition is the VERDICT and the SIGN of the gap, not a numeric tolerance**, because the recorded
  run's seed is not recorded so exact reproduction is not expected: require `WORSE_THAN_MARKOV_NULL` and
  `GLM > null` in both modes. Report the magnitudes beside the recorded 0.7237/0.5990 and 0.7706/0.6679.
- If the verdict or sign does NOT reproduce, STOP and report: the sweep would be uninterpretable against
  the recorded baseline.

**Test strategy:**
- The run's own validity gates; the artifact must record `conditional: true` before any metric is read.
- Unit: the anchor-check helper returns pass only when verdict AND sign both match.

### Step 4: Sweep sampling downward from the default
Files: wiki/glm_sweep_measurement_2026-10-20.json
Depends on: Step 3

**What changes:**
- Grid `temperature ∈ {0.6, 0.8}` plus the Step 3 default (1.0) as the third point, each with
  `--verify-knob`, all other flags identical to Step 3.
- **Downward only, and that is derived rather than arbitrary:** the uniform control scores 0.9408 (far
  worse than the GLM's 0.7237), so flattening the distribution moves output toward the worst arm on
  record. Sharpening is the only direction with a mechanism for improvement.
- `top_k` is NOT swept in this plan: `TOP_K_KMER_FLOOR = 22` bounds it to k-mer-level truncation whose
  effect on per-base marginals is measured but whose interaction with the discriminator is not, so it
  would add an unattributable second lever. Named as deferred, not forgotten.
- Each grid point stamps `knob_contract.knob_verified_end_to_end`; a point whose knob is unverified is
  excluded from the verdict rather than averaged in.

**Test strategy:**
- Per-point validity gates as in Step 3.
- Unit: the grid assembly refuses a point whose `knob_verified_end_to_end` is false.

### Step 5: Report against the pre-registered comparator
Files: wiki/glm_sweep_result_2026-10-20.md, wiki/glm_sweep_result_2026-10-20.json
Depends on: Step 4

**What changes:**
- Per mode and per grid point: GLM distinguishability, this run's null, the gap, and the branch from
  Verification Signal.
- Record the G-D premise read verbatim against `:5`'s rule, and state that acting on a weakened premise is
  a scope decision rather than a consequence of this plan.
- Honest limits: n=120 as recorded; the Markov null is UNCONDITIONAL while the GLM is conditioned, so the
  comparison handicaps the null and the GLM's loss is if anything understated; one organism; one promoter
  corpus; `top_k` unswept; the recorded run's seed unknown.

**Test strategy:**
- Unit: the result's branch equals the branch function applied to its own recorded inputs.
- Unit: the result cites the pre-registered 0.5990 / 0.6679 and does not invent a band.

### Step 6: Discharge the trigger in G-D's ledger and feature record
Files: project_state/glm-generator-finetune-2026-10-07.md, features/glm-generator-finetune/feature.md
Depends on: Step 5

**What changes:**
- Append the sweep artifact to the ledger's Findings/Evidence with its date and verdict, via
  `/project-state --append-action --class run-tests` (AC9: Action-Log rows take the sanctioned op; the
  Findings/Evidence tables take their own ops, never a hand edit).
- Update the short-term action-1 row's status from `blocked` to reflect that the sweep has reported.
- Do NOT change G-D's `blocked_by` edge for the 8 GB GPU, and do NOT start fine-tuning: the hardware half
  of the trigger is unchanged and writing G-D's plan is that family's own step.

**Test strategy:**
- Unit: the ledger's evidence row cites the new artifact path and the artifact exists.
- Manual: `git diff` shows no edit to any `## Decisions Made` or `## Hypotheses` row (those have their own
  ops and are out of scope here).

## Execution Preview

- **Wave 0** (parallel, 2): Step 1, Step 2
- **Wave 1**: Step 3 (needs the guard)
- **Wave 2**: Step 4 (needs the anchor to pass)
- **Wave 3**: Step 5
- **Wave 4**: Step 6

5 waves, max parallelism 2. Critical path: 1 → 3 → 4 → 5 → 6. No intra-wave file overlap.
Cost: Steps 1-2 minutes; Step 3 ~25-30 min (120 sequences at 150 bp on CPU); Step 4 two more points,
~50-60 min; Steps 5-6 minutes. Total ~1.5-2 h, dominated by generation.

## Risk Flags

- **The Markov null is unconditional while the GLM is conditioned.** `baselines.py:92` has no conditional
  variant. This handicaps the null, so a GLM loss is understated and a GLM win would be the weaker claim.
  Named, reported in Step 5, NOT corrected here (building a conditional null is a separate scope call).
- **The recorded run's seed is not in the ledger**, so Step 3 can only reproduce the verdict and gap sign,
  not the exact value. That is why the anchor is specified on sign, not tolerance.
- **n=120 at length 150 is the dominant cost.** If Step 3's wall-clock exceeds ~45 min, Step 4's two
  points should be run one at a time rather than batched, so a host fault loses one point not three.
- **Intermittent host faults:** Hub-fetch segfault (mitigated by `HF_HUB_OFFLINE=1`) and the
  access-violation-under-concurrent-torch gotcha (mitigated by serialising).
- **`sentrux` absent** → architecture gate `n/a`. **`codex` absent** → any review of this plan is
  single-voice, which is weaker than the usual two and should be stated rather than glossed.
- **Scope discipline:** `plans-index.md:708` records that proposing generator fine-tuning here would be
  orientation error #3. Step 6 explicitly does not start it.

## Open Questions

- **Should `top_k` be swept in a follow-up?** Deferred above to keep levers attributable. The floor (22)
  is measured, but its interaction with the k-mer discriminator is not — a sweep would need its own
  anchor.
- **If the sweep clears the null, does G-D get cancelled or re-scoped?** The pre-registered rule says the
  premise "weakens materially"; it does not say what follows. That is a scope decision for the user, and
  Step 5 records it as one rather than resolving it.
- **Is prompting from EVAL-window context a leakage risk?** The prompt is the 300 bp immediately upstream
  of the scored window. The runner argues no, since only the continuation is scored, and the committed
  2026-10-07 run already lives with this — so it is pre-existing rather than introduced. Worth one
  adversarial pass before any number is published externally.

## Verification

- Steps 1, 2, 5, 6: new unit tests run with the suite; baseline **15 failed / 6688 passed / 185 skipped**,
  and any change to the failure count must be explained before the plan is called done.
- Steps 3-4: the runs' own validity gates; mode and knob-verification fields read BEFORE any metric.
- Frozen five (`amr_rules.py`, `calibrated_amr_rules.json`, `mic_tiers.py`, `shipped_decoder_surface.py`,
  `cohort_manifest.py`) byte-unchanged and the prospective lock re-verified after every step.
- CLAUDE.md word ceiling re-checked if any orientation text is touched (35,435 words, 25 under).

## Save-time amendments

Captured at: 2026-10-20
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads
ONLY `## Implementation Steps` for executable work; these bullets are NOT executable instructions.
If an amendment changes a Step contract (files, dependencies, structure), re-run `/technical-plan`
before `/execute-plan` rather than relying on this section.

- relocated from a new measurement to discharging G-D's queued sampling-sweep trigger
- comparator and decision rule taken from the existing ledger rather than invented
- both/genomic natural-set combinations resolved to REFUSE
- grid restricted downward on the uniform-control argument
- top_k deferred to keep levers attributable

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
