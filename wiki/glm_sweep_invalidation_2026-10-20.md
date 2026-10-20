# The 2026-10-19 sampling sweep was DOUBLY invalid — and G-D's trigger is still undischarged

**Date:** 2026-10-20 · Step 2 of `plans/Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan/`

## Two independent defects, either of which voids the run

**1. The sampling knobs were INERT.** `HFGenerator` accepted `temperature` and passed it to
`model.generate`, but the remote `_BPLogitsProcessor` collapses scores to a one-hot BEFORE
transformers appends `TemperatureLogitsWarper` (merge at `generation/utils.py:1257`, temperature at
`:1276`), and dividing a one-hot by any temperature is identity. Measured: T=0.7 vs T=1.5 at the same
seed produced **6/6 byte-identical sequences**, GC identical to six decimals. So the three grid points
at T=0.8 / 1.0 / 1.2 all used **the same (default) sampling**, and the invariance they showed was the
knob being disconnected rather than a property of the generator. Fixed 2026-10-20 —
`wiki/glm_temperature_knob_2026-10-20.md`.

**2. The comparison mode was wrong.** All three artifacts record `conditional: False` with
`natural_set: ['promoter']`. That scores an UNCONDITIONAL sample — drawn from the whole-genome prior,
which is mostly coding — against ALIGNED promoter windows, while the Markov null is fitted on the
promoter FIT half (`glm_generator_falsifier.py:194`). The null holds the target distribution and the
model does not. The repo documented this in two places *before* the run:

- `corpus.py:168-175` — "An unconditional sample cannot be scored against an aligned set fairly, so the
  meaningful task is *given the genomic context, continue into the promoter*."
- the `--conditional` help — "**required** for a meaningful comparison against the ALIGNED promoter set".

So its headline (GLM **0.8141** vs null **0.3876**, `WORSE_THAN_MARKOV_NULL`) does not support a
generator-quality claim. It is close to tautological: an unconditioned sample scores worse against
aligned promoters than a promoter-fitted null.

## What the valid baseline IS

**`wiki/glm_generator_falsifier_REAL_2026-10-07.json`** — committed, `conditional: True`,
`natural_set: promoter`, `n_candidate: 120`, `control_uniform: 0.9408`:

| mode | GLM | Markov-k3 null | verdict |
|---|---|---|---|
| `kmer` | 0.7237 | 0.5990 | `WORSE_THAN_MARKOV_NULL` |
| `both` | 0.7706 | 0.6679 | `WORSE_THAN_MARKOV_NULL` |

Geometry, recorded as "ALL CONFIRMED" in `project_state/glm-generator-finetune-2026-10-07.md:20`:
window 150 bp, context 300 bp, markov_order 3, n_generated 120, fit/eval 600/600.

**That direction is the project's RECORDED PREMISE, not a new finding.**
`features/glm-generator-finetune/feature.md:13` states family G-D exists to determine whether LoRA
fine-tuning brings the generator "below the distinguishability of a 3-mer Markov chain, **which it
currently loses to at default sampling**." Presenting the conditional loss as a discovery would be the
repo's recorded orientation error #3.

Conditioning does close most of the gap — **+0.4265 unconditional → +0.1247 conditional** (`kmer`) —
but it does not close it. Both readings stand together.

## Why the trigger is still undischarged

**No committed sweep artifact exists.** Verified by `git ls-files` and by listing `wiki/`: the only
falsifier artifacts are the two 2026-10-07 runs, both at default sampling. The 2026-10-19 scratchpad
artifacts are not committed and are void for the two reasons above.

`project_state/glm-umbrella-2026-10-07.md:18` has G-D `blocked_by: sampling sweep + 8 GB GPU`, and
`features/glm-generator-finetune/technical-plan.md:21` names "the sampling sweep artifact exists" as
half the unblock trigger. **That half remains open.**

Worth separating, because it was never named before: the **8 GB GPU blocker applies to fine-tuning, not
to the sweep.** The sweep is CPU inference and has been runnable on this host all along. What actually
blocked it was defect 1 — nobody knew the knobs were inert, so any sweep run before 2026-10-20 could
not have varied sampling and therefore could not have discharged the trigger.

## Scope

This memo records invalidity and relocates the baseline. It makes **no** claim that the GLM is good or
bad beyond what the committed 2026-10-07 run already established, and it does **not** start G-D.
Two mechanism hypotheses are already REFUTED on record and must not be restated: mode collapse
(pairwise identity 0.249 vs ~0.25 random) and coding-like output (stop density 4.83 vs natural 4.67).
