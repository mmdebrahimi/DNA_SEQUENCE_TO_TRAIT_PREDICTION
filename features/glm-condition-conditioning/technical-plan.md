# Technical plan — G-C: condition-conditioned expression oracle

## Lens status
Repo-index: absent — grounding is direct file reading.
Baseline tests: `uv run pytest tests/test_glm_*.py -q` → 59 pass.
sentinel/sentrux: n/a. Mode: sequential (steps share one script).

## Problem Statement

Growth condition carries **real, measured** signal: on 297,599 fragments shared between LB and M9, the
noise-matched cross-condition agreement is **0.4062** against within-condition **0.4676** (gap **−0.0614**),
and the disattenuated true cross-condition correlation is **0.7887** — so roughly **38%** of the real
variance is condition-specific (derived: 1 − 0.7887²).

This is the bacterial analogue of AlphaGenome's cell-state conditioning, and it is the one item from that
critique that was both new and already evidenced. The question is whether **one** model taking
(sequence, medium) beats **two** independent per-medium models.

## Codebase Context

### Reusable-Code Survey
- **`dna_decode/glm/genomewide.py`** — `load_fragments(which="LB"|"M9")`, `position_blocked_split`,
  `noise_ceiling`, `condition_effect`, `kmer_features`, `gc_features`, `spearman`. Everything needed to
  load, split and ceiling-normalise both media already exists and is tested (21 tests).
- **`scripts/glm_genomewide_oracle.py`** — `_fit` is the arm shape to mirror.
- **None — searched:** no existing conditioned/multi-task model anywhere under `dna_decode/` (searched for
  `condition`, `medium`, `multi_task`). New code.

### The substrate's binding constraint
Fragments are **48–499 bp with 359 distinct lengths**, so only **length-invariant** features apply —
`positional_onehot` is deliberately unavailable for this arm and refuses ragged input rather than padding.
The fragment ceiling is low (**LB 0.7928 / M9 0.8014**) because 97% of fragments carry a single barcode.

## Pre-Change Baseline

From `wiki/glm_genomewide_oracle_2026-10-07.json`, position-blocked, 30,000-fragment subsample, LB only:
GC **0.1023** · 4-mer **0.1313** · 6-mer **0.0919**. A size ladder (3k / 10k / 30k → 0.099 / 0.154 / 0.131
for 4-mers) shows this **plateaus rather than climbing**, so the weak absolute numbers are not a
power artifact.

**Per-medium baselines do not exist yet** — Step 1 creates them, and they are the real comparator.

## Verification Signal

**Primary:** Spearman of a single (sequence, medium) model on held-out position-blocked data, against the
**pooled** performance of two independently-fitted per-medium models, margin tested at ≥0.02, each medium
scored against **its own** ceiling (0.7928 / 0.8014 — never a shared normaliser).
**Mandatory null:** a **medium-shuffled** arm. The shared model trains on both media and therefore sees ~2×
the rows, so an unqualified gain is attributable to data volume, not to conditioning. Without this arm the
headline is uninterpretable — this is the experiment, not a nicety.
**Secondary:** a **data-matched** arm — the shared model trained on a half-size sample so its row count
matches a single-medium model. This separates "sharing helps" from "more data helps" directly.

## Implementation Steps

### Step 1: Per-medium baselines and per-medium ceilings
Files: scripts/glm_condition_conditioning.py, tests/test_glm_condition.py
Depends on: none

**What changes:**
- Load LB and M9, intersect to the shared-fragment set (297,599), and build **ONE** position-blocked split
  used identically by every arm — so arms differ only in the model, never in the split.
- Fit two independent per-medium ridges (GC / 4-mer); record each against its own `noise_ceiling`.
- Emit the two-model pooled baseline that the conditioned model must beat.

**Test strategy:** assert the split is identical across arms (same held block ids); assert each medium's
ceiling is computed from that medium's own replicates; refuse if the shared-fragment set is below a floor.

### Step 2: The conditioned model + the three control arms
Files: scripts/glm_condition_conditioning.py, tests/test_glm_condition.py
Depends on: Step 1

**What changes:**
- `shared` arm: features = [sequence features ‖ medium indicator], one model, both media.
- `shuffled` arm: identical, with the medium label permuted. **Must show no gain.**
- `data_matched` arm: `shared` trained on a row-count-matched subsample.
- Pre-registered verdict frozen in the file before any run:
  `CONDITIONING_HELPS` (shared − two_model ≥ 0.02 **AND** shared − shuffled ≥ 0.02 **AND** data_matched also
  clears) / `GAIN_IS_DATA_VOLUME` (shared beats two_model but data_matched does not) /
  `NO_GAIN` / `INDETERMINATE_NULL_NOT_CLEAN` (the shuffled arm itself gains ≥0.02 → the pipeline fabricates
  signal and nothing is graded).

**Test strategy:** the verdict function is pure — unit-test all four branches on synthetic arm dicts,
including `GAIN_IS_DATA_VOLUME`, which is the branch that exists specifically because it is the likely one.

### Step 3: Run, commit the artifact, record the verdict
Files: wiki/glm_condition_conditioning_2026-10-07.json, project_state/glm-condition-conditioning-2026-10-07.md
Depends on: Step 2

**What changes:**
- Execute; commit the artifact whatever the verdict; append an Action Log row naming it.
- Carry the honest limits into the artifact: two conditions only; fragment substrate with a low ceiling;
  conditioning tested as an indicator, not as representation modulation.

**Test strategy:** the MVP predicates are the test.

## Execution Preview
Waves: 3 (Step 1 → Step 2 → Step 3), max parallelism 1 — all three touch one script, so there is nothing to
parallelise and claiming otherwise would be false. Self-budget estimate: ~12 steps. All steps classify
`auto`. No money, no destructive-local, no outward-irreversible.

## Risk Flags
- **The doubled-data confound is the main risk and it is near-certain to appear.** A naive two-vs-one
  comparison will show a gain. `shuffled` + `data_matched` exist to attribute it. If `data_matched` does not
  clear, the honest verdict is `GAIN_IS_DATA_VOLUME`, not a win.
- **n = 2 conditions.** A null result is weak evidence against conditioning in general, and the artifact must
  say so. Flagged in the ledger as a user-facing acceptance question.
- **Low ceiling substrate** (0.79) — absolute numbers will look poor. They must be read as a fraction of
  ceiling, and the two media's ceilings differ, so each arm carries its own.
- **Length-invariant features only** — this arm cannot use the grid's winning positional feature, so it is
  not comparable to the 0.5914 grid number. Stated, not papered over.

## Open Questions
- Whether a medium indicator is a fair test of "conditioning", or whether a negative only rules out the
  weakest form of it (currently: start with the indicator, and say which form was tested).
- Whether the ~38% condition-specific variance is concentrated in regulated/catabolic promoters — a
  follow-on, deliberately out of scope here.

## Verification
`uv run pytest tests/test_glm_condition.py -q` plus the full `tests/test_glm_*.py` set against the 59
passing today, then the live run and the four MVP predicates re-evaluated.
