# Technical plan — G-A: learned representation for the expression oracle

## Lens status
Repo-index: absent (not configured here) — grounding is direct file reading.
Baseline tests: `uv run pytest tests/test_glm_*.py -q` → 59 pass across the three GLM data modules.
sentinel/sentrux: n/a. gh: available. Mode: sequential (steps share files).

## Problem Statement

The expression oracle is ridge regression on hand-specified features. On the designed σ70 grid positional
one-hot reaches Spearman **0.5914** (leave-element-out), but on real genomic peak tiles a single GC-content
number reaches **0.3000** and **beats every learned feature** (one-hot 0.2691 / 4-mer 0.2769 / 6-mer
0.2022) — the best learned arm reaching only **29.4% of the replicate-derived 0.9410 ceiling**. Fixed
features demonstrably do not transfer off the grid.

AlphaGenome's transferable design principle is capacity in the **representation** with a deliberately
trivial head. We have the trivial head and no learned representation at all. This plan builds one and
measures whether it transfers — committing a verdict either way.

## Codebase Context

### Reusable-Code Survey
- **`dna_decode/glm/genomewide.py`** — `load_peak_tiles` (tab-delimited, returns `Tile` with `rep1`/`rep2`),
  `leave_peak_out`, `noise_ceiling`, `positional_onehot` (refuses ragged input), `kmer_features`,
  `gc_features`, `spearman`. The tile substrate and its honest split already exist.
- **`dna_decode/glm/expression.py`** — `load_pairs`, `leave_element_out`, `fit_predict_ridge(...,
  return_pred=True)`. The grid substrate and the ridge baselines already exist.
- **`dna_decode/glm/selection.py`** — `simulate_selection`, `permutation_null`. If the encoder wins, the
  selection measurement can be re-run on it with no new code.
- **`scripts/glm_genomewide_oracle.py`** — the runner whose `_fit` shape the new arms mirror.
- **None — searched:** no existing conv/CNN module anywhere under `dna_decode/` (searched `dna_decode/**`
  for `Conv1d`, `torch.nn`, `Sequential`). This is genuinely new code, not a refactor.

### External surfaces
torch ≥2.2 via `uv sync --extra forward` (already installed: 2.7.1+cu118). CPU is sufficient — 44,106
sequences × 150 bp. No network, no Docker.

## Pre-Change Baseline

Measured and committed in `wiki/glm_genomewide_oracle_2026-10-07.json` (peak tiles, leave-peak-out,
n_train 32,908 / n_test 11,198, designed controls and non-150bp dropped):

| arm | Spearman | fraction of the 0.9410 ceiling |
|---|---|---|
| **GC content (the bar to beat)** | **0.3000** | 0.3188 |
| 4-mer | 0.2769 | 0.2943 |
| positional one-hot | 0.2691 | 0.2860 |
| 6-mer | 0.2022 | 0.2149 |

Grid comparator (`wiki/glm_expression_oracle_result_2026-10-07.md`): GC 0.1414 / 4-mer 0.4137 /
6-mer 0.5872 / one-hot 0.5914.

## Verification Signal

**Primary:** encoder Spearman on genomic peak tiles under `leave_peak_out`, against GC's 0.3000, with the
margin tested at ≥0.05, reported as a fraction of the 0.9410 ceiling.
**Secondary:** the same encoder on the designed grid under `leave_element_out`, against one-hot's 0.5914 —
so a result that transfers and a result that only works on the grid are distinguishable.
**Null:** a label-shuffled encoder run, so a pipeline bug that fabricates signal is detectable.

## Implementation Steps

### Step 1: Linear-headroom diagnostic
Files: scripts/glm_tile_headroom.py, tests/test_glm_tile_headroom.py
Depends on: none

**What changes:**
- Fit ridge on GC, GC+dinucleotide (16), GC+positional-GC-in-10-bins, and dinucleotide+position, all under
  `leave_peak_out` on the same split as the baseline.
- Emit `wiki/glm_tile_headroom_2026-10-07.json` with each arm against the 0.9410 ceiling.
- **This step can CLOSE the family in an afternoon.** If composition-only features already reach ~0.30 and
  adding position adds nothing, GA-H3 is confirmed and a conv encoder is a predicted-negative.

**Test strategy:** unit-test the feature builders for width and length-invariance; assert the GC arm
reproduces 0.3000 ± 0.01 on the committed split (a reconcile gate — if it does not, the split differs).

### Step 2: The encoder module
Files: dna_decode/glm/encoder.py, tests/test_glm_encoder.py
Depends on: none

**What changes:**
- `SequenceEncoder`: one-hot → parallel `Conv1d` branches at **kernel widths 6 / 20 / 75** (motif / spacer /
  whole-promoter scale — this IS the multi-resolution parameterisation, not a separate project) → ReLU →
  global max+mean pool per branch → concat → linear head. Deliberately small (~32 channels per branch).
- `fit_encoder(train, test, *, widths, epochs, seed) -> dict` mirroring `_fit`'s return shape so the runner
  treats it as one more arm.
- Deterministic given a seed; CPU default; **no GPU requirement**.
- **Refusals:** ragged input raises (reusing the `positional_onehot` rationale — no silent padding); an
  empty train set raises rather than returning a degenerate model.

**Test strategy:** a planted-signal test (a motif inserted at a fixed offset must be learnable to ρ>0.8);
a shuffled-label test (must collapse to ≈0); determinism under a fixed seed; the ragged-input refusal;
output-shape pins for each kernel width.

### Step 3: Single-resolution vs multi-resolution arms
Files: scripts/glm_encoder_gate.py, tests/test_glm_encoder_gate.py
Depends on: Step 2

**What changes:**
- Arms: `gc` (baseline), `encoder-w20` (single resolution), `encoder-multi` (6/20/75), `encoder-shuffled`
  (the null).
- Runs on BOTH substrates: tiles (`leave_peak_out`) and grid (`leave_element_out`).
- **Each substrate normalised by its OWN ceiling** — tiles 0.9410, grid none derivable (stated, not faked).
- Pre-registered verdict, frozen in the file before any run:
  `GENERALISES` (encoder − GC ≥ 0.05) / `WEAK` (0 < margin < 0.05) / `DOES_NOT_GENERALISE` (margin ≤ 0) /
  `INDETERMINATE_NULL_NOT_CLEAN` (the shuffled arm exceeds 0.05, meaning the pipeline fabricates signal).

**Test strategy:** the verdict function is a pure function — unit-test all four branches on synthetic arm
dicts, including the null-not-clean branch, and assert it refuses rather than grading when the null is dirty.

### Step 4: Run, commit the artifact, record the verdict
Files: wiki/glm_learned_representation_2026-10-07.json, project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 1, Step 3

**What changes:**
- Execute the gate; commit the artifact whatever the verdict.
- Append an Action Log row naming the verdict via `/project-state --append-action --class run-tests`.
- If `DOES_NOT_GENERALISE`: write the closed-negative memo and mark GA-H1 falsified / GA-H3 confirmed. **The
  family reaches MVP on a negative** — that is the bar as written.

**Test strategy:** the MVP predicates are the test (`file-exists` on the artifact, `test-exit-0` on the
encoder tests, `project-state-row` on the ledger).

## Execution Preview
Waves: 3 (Step 1 ∥ Step 2 → Step 3 → Step 4). Max parallelism 2. Critical path Step 2 → 3 → 4.
Self-budget estimate: ~18 steps. Gate classification: all `auto` (in-cwd writes, local tests). No money, no
destructive-local, no outward-irreversible actions.

## Risk Flags
- **Step 1 may end the family.** That is a feature, not a risk — but it means Steps 2-4 should not start
  before Step 1 reports if time is tight.
- **Framing risk, larger than the architecture risk:** 15,171 active / 31,542 inactive tiles. Spearman over
  a mostly-inactive set may be largely ranking noise. Mitigation: Step 3 additionally reports AUROC on the
  `active` column as a diagnostic. If regression is the wrong frame, that is a finding, not a failure.
- **Overfitting:** 6-mers (4,096 features) scored WORST on tiles, so capacity has already hurt once here.
  Mitigation: deliberately small encoder, and the shuffled-label null must stay clean.
- **No GPU** — accepted; the data is small. If a run exceeds ~30 min on CPU, reduce channels rather than
  reaching for the offline server.

## Open Questions
- Whether to add the classification framing as a full arm or keep it diagnostic (currently: diagnostic).
- Whether a grid-trained / tile-tested cross-substrate arm is in scope (currently: NO, it is a different
  and harder question — candidate 3 in the ledger).

## Verification
`uv run pytest tests/test_glm_encoder.py tests/test_glm_encoder_gate.py tests/test_glm_tile_headroom.py -q`
plus the full `tests/test_glm_*.py` set to confirm no regression against the 59 passing today. Then the
live gate run, and the four MVP predicates re-evaluated.
