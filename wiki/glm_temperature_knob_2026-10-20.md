# The GLM's `temperature` knob was SHIPPED DEAD — found, mechanism read, fixed, verified

**Date:** 2026-10-20 · `dna_decode/glm/generate.py` (NON-frozen) · 6 tests `tests/test_glm_sampling_knob.py`

## Headline

`HFGenerator(temperature=T)` had **no effect whatsoever**. T=0.7 vs T=1.5, same seed, one loaded
model: **6/6 BYTE-IDENTICAL sequences**, GC identical to six decimals (0.516667 both). After the fix:
**0/6 identical**, and GC separates in the mechanistically correct direction.

| | T=0.7 | T=1.5 | identical |
|---|---|---|---|
| before | GC 0.516667 | GC 0.516667 | 6/6 |
| after | GC **0.547222** | GC **0.497222** | **0/6** |

Lower T → sharper distribution → more of the model's preferred (GC-rich) bases; higher T → flatter →
GC drifts toward uniform 0.5. A knob that merely "changed the numbers" would not produce that ordering.

**Sampling itself was never broken** — sequences within a single call always differed. Only the knobs
were dead: `temperature` and `top_k`.

## Mechanism — an ordering inversion, and the upstream docstring asserts the opposite

Three independent lines agree.

1. **Source.** The remote `_BPLogitsProcessor` (in the model's own `modeling_generator.py`)
   marginalises k-mer logits to per-base probabilities, samples with `torch.multinomial`, and returns a
   **ONE-HOT**: `-inf` everywhere, `0.0` at the selected token. Its docstring claims it *"Runs LAST in
   the logits-processor chain so that temperature / top-k / top-p etc. influence the marginal
   distributions before base selection."* **That is backwards.**
2. **HF ordering, read from the installed source** (`transformers/generation/utils.py`, 5.14.1):
   user/custom processors are merged at **line 1257**, and `TemperatureLogitsWarper` is appended at
   **line 1276** — i.e. AFTER. So temperature sees the one-hot, not the marginals.
3. **Arithmetic.** Dividing a one-hot by any temperature is identity: `0/T == 0`, `-inf/T == -inf`.
   Demonstrated model-free and pinned as the keystone test; the same warper on real logits
   `[1,2,3]` at T=0.5 correctly gives `[2,4,6]`, so the keystone is not vacuous.

## The fix

The remote `generate` builds `existing + [bp_proc]` from the caller's `logits_processor`, so anything
we pass there runs **BEFORE** the one-hot collapse. `HFGenerator._sampling_processors()` returns a
`TemperatureLogitsWarper` (when `T != 1.0`) and a `TopKLogitsWarper` (when `top_k`), and both call
sites now pass them via `logits_processor=` and **no longer pass the inert `temperature=` / `top_k=`
kwargs at all**. Nothing in the remote code or the model cache is modified.

`T == 1.0` installs nothing, because a 1.0 warper is identity and would only add chain noise.

**Non-vacuity proven by reverting:** on the pre-fix file 4 of 6 tests FAIL, including the regression
pin asserting `temperature=` is absent from the generate kwargs; all 6 pass after. The 2 that pass
either way are the arithmetic keystones, which test transformers rather than this repo.

## Consequence — what this does and does NOT invalidate

**Does NOT invalidate the measured distinguishability itself.** The 2026-10-19 sweep's
**0.8141** is a valid measurement of this generator under its DEFAULT sampling, and
`WORSE_THAN_MARKOV_NULL` stands.

**DOES void that sweep's claim to have varied sampling.** All three artifacts
(`sweep_t0.8/1.0/1.2.json`) stamp `sampling.temperature` 0.8 / 1.0 / 1.2 and two carry
`"is_recorded_default": false` — but every run used identical sampling, so the invariant 0.8141 was
the knob being dead, not evidence that temperature does not help. **Any GLM artifact dated on or
before 2026-10-19 whose `sampling` block shows a non-1.0 temperature or a non-zero top_k describes
DEFAULT sampling regardless of what it records.**

**It partially REOPENS a closed-sounding question.** "The generator is the bottleneck and sampling
does not rescue it" was never actually tested — the only generator-side tuning knob available was
inert. Whether temperature/top_k move distinguishability is now an OPEN and cheap question (one model
load, ~90 s for 6 sequences). **Not claimed here: that they will help.** This memo establishes only
that the knob now functions and that the earlier invariance had a plumbing explanation.

## Reusable

A knob can be plumbed correctly at the call site and still be dead, when a downstream component
**destroys the quantity the knob acts on**. `do_sample=True` and `temperature=self.temperature` were
both present and correct; the defect was ORDER. Check a knob by varying it and diffing the OUTPUT —
not by reading the call site, which looked right for months.

## The pre-run gate (added the same day, after the user pressed on the real failure)

The diagnosis above is postmortem, and the user's objection was that the check belonged **before**
the run: *"that should have obviously been checked before the runs, not after… we really need to make
sure the run is structured correctly before actually running it, and not just trial and error and hope
for the best and fixing things postmortem."* That is correct, and it is worse than a missed check —
the invariance was **reported as a finding** before the knob was verified. This repo also already
carried the lesson in another vocabulary (two shipped filters that removed nothing, recorded as *"a
filter that removes nothing is not a control"*), and it was not applied to this experiment's own design.

So the answer is a precondition, not a resolution to be careful. `dna_decode/glm/knob_contract.py`:

| guard | when it fires | cost |
|---|---|---|
| `validate_sampling` | `HFGenerator.__init__`, BEFORE the ~80 s load | microseconds |
| `assert_knob_wired` | whenever warpers are built | a list-length check |
| `verify_knob_live` | `--verify-knob`, before any sweep work | one 2-value diff |

- **`validate_sampling`** refuses `temperature <= 0` / NaN, and `top_k` below a measured floor. This
  also closes a defect the fix itself introduced: once the knob was live, `--temperature 0` became a
  `ValueError` raised ~80 s in, after the model load. Now it fails in **0.000 s**.
- **`assert_knob_wired`** refuses a run whose requested knob installed no warper. That is the exact
  state that let three artifacts stamp `is_recorded_default: false` over runs that all used default
  sampling — a false provenance claim no reader could detect.
- **`verify_knob_live`** is the check that was actually missing: generate at two far-apart
  temperatures with a fixed seed and diff the raw output. `--verify-knob` runs it **before** the sweep
  and **exits 5** if the output does not move, because every grid point would otherwise measure the
  same sampling. Verified live on the real model: `knob verification: TEMPERATURE_IS_LIVE (0/1
  identical)`, stamped into the artifact as `knob_verified_end_to_end: true`.
  **Soundness note:** n=1 suffices for liveness ONLY because the pre-fix run established the pipeline
  is deterministic at a fixed seed (6/6 identical). Without that control a single differing sequence
  would prove nothing.

### `top_k` was NOT disabled, and the floor was derived twice

The review's instinct was to disable `top_k`, since it truncates the **4096 six-mers**, not the 4
bases. Measuring instead changed the answer: counting (position, base) cells driven to probability
exactly 0 out of 24, across 5 seeds × logit sharpness 1.0/2.0/4.0 —

| top_k | 4 | 10 | 14–21 | **22** | 24 | 64 | 1024 |
|---|---|---|---|---|---|---|---|
| worst zero-cells | 7–9 | up to 5 | 1 | **0** | 0 | 0 | 0 |

So it is a legitimate graded knob at and above 22, and below it at least one base becomes
**unreachable** at some position — silently greedy there. `top_k=4`, the value someone types meaning
"4 bases", zeroes 7–9 of 24 cells. Disabling would have recreated a dead knob; the contract refuses
the measured-broken region with that explanation and leaves the rest live.

**The floor was first set to 24 and that was wrong.** The scan went 16, 24, 32 — so 17–23 were never
measured and the constant was **padded by 2** because a rounder number felt safer.
`test_the_floor_is_the_SMALLEST_value_with_no_reachability_loss` re-derives the boundary and fails if
the constant is too low **or too high**; it caught the padding, and that is the only reason it reads
22. Padding a threshold is still asserting it.

## Honest limits of the gate

- `--verify-knob` is **opt-in**, so a sweep can still be run without it. What is not optional is the
  artifact field: an unverified run records `knob_verified_end_to_end: false`, so a reader can tell.
- The floor is measured on **synthetic Gaussian logits**, not the real model's. Indicative, not
  exhaustive; the real-model check is the `--verify-knob` diff.
- `assert_knob_wired` checks OUR half of the plumbing. The half living in third-party code and in
  framework ordering is only covered by the end-to-end diff.
- None of this is enforcement — a run can be launched without the flag, and nothing compels reading
  the stamped field. It converts a silent void into a recorded one.
