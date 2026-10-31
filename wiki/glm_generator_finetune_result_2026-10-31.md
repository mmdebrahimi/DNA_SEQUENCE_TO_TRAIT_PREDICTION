# G-D closes NEGATIVE: LoRA fine-tuning made the generator MORE distinguishable, not less

<!-- doc-status: LIVE as of 2026-10-31. This is the measured result of the pre-registered arm in
     wiki/glm_generator_finetune_acceptance_bar.json. Read that bar before quoting anything here;
     the verdict vocabulary and all four outcomes were fixed before any number existed. -->

## The result

**Pre-registered outcome C: `C_NO_IMPROVEMENT`. GD-H1 is FALSIFIED.**

| arm (Tesla T4, same session) | kmer | **both** (the bar's mode) |
|---|---|---|
| 3-mer Markov null — **the bar** | 0.5990 | **0.6679** |
| GENERator zero-shot | 0.7355 | 0.7787 |
| **GENERator + LoRA, 3 epochs** | **0.7559** | **0.8058** |
| uniform control | 0.9408 | 0.9641 |

Lower is better. The fine-tuned generator is **worse by +0.0271** (`both`) and **+0.0204** (`kmer`).
It did not approach the bar; it moved away from it, and away from its own starting point.

`wiki/glm_generator_finetune_verdict.json` · arms in
`wiki/glm_generator_finetune_arm{A,B}_*_2026-10-31.json` · log in
`wiki/fragment_null_batches/glm_generator_finetune_2026-10-31.log`

## Why this is a real measurement and not a plumbing artifact

This project has a recorded trap in exactly this place: the GLM path once scored **zero sequences**
while emitting a complete-looking artifact. So:

| check | ARM A | ARM B |
|---|---|---|
| sequences generated | **120** | **120** |
| alphabet valid / exact length | 1.0 / 1.0 | 1.0 / 1.0 |
| candidate identity in the artifact | `hf:GENERator-v2-prokaryote-1.2b-base` | **`hf:glm_finetuned_merged`** |
| Markov null | 0.6679 | **0.6679** (identical) |
| uniform control | 0.9641 | **0.9641** (identical) |
| natural-vs-natural self-control ceiling | 0.5865 | 0.5865 |

ARM B's candidate name is the load-bearing one: it proves the **merged adapter** was loaded and
scored, not the base model a second time. The null and uniform controls being bit-identical across
arms proves both scored the *same* eval windows.

**Training really happened.** Mean loss fell monotonically across all three epochs —
**7.4361 → 7.1655 → 6.9148** — over 3,737 windows, ~9.6 min of the 15.7 min run. So this is not
"the adapter did nothing": the model learned, and learning made its output easier to tell from real
promoters.

**Disjointness held.** 4,340 windows → 600 fit / 600 eval / 3,140 unused / **3,737 train**, overlap
**0**, asserted on sequence content and re-checked inside the kernel against the committed manifest.
The generator was never scored on anything it trained on, which matters doubly here: had leakage
existed it would have pushed the number *down*, i.e. toward a false pass.

## A mechanistic lead, stated as a lead

GC content moved **away** from natural:

| | GC |
|---|---|
| natural eval windows | 0.4552 |
| zero-shot | 0.4574 (≈ natural) |
| fine-tuned | **0.4346** |

Zero-shot GC already matched natural almost exactly; fine-tuning pulled it off by ~0.02. Since the
discriminator keys on local composition, a composition shift is a plausible route to the worse
score. **This is a lead, not a diagnosis** — it is one summary statistic, not an attribution, and
nothing here tests whether GC drift *causes* the gap.

## Deviation from the pre-registration, declared

**The frozen method says `dtype: float32` for both arms. Both arms actually ran at float16.**

The kernel's `DTYPE` governed only the model it loaded for *training*; each scoring arm goes through
`glm_generator_falsifier.py`, which chooses its own dtype via `recommended_dtype()` — and on a T4
(CC 7.5, `fp16_fast: true`) that is float16. So the declared parameter was one the kernel did not
control. That is a defect in my kernel, not in the falsifier.

**What this does and does not cost.** The *rationale* behind the fp32 clause was that nothing except
the weights may differ between the arms, and that holds exactly: both arms loaded float16, scored the
same eval windows, with identical null and uniform controls. Corroboration that the baseline is not
an fp16 artifact: ARM A's 0.7355 / 0.7787 reproduces the sampling sweep's T=1.0 CUDA anchor
**exactly**, and that is an independently-run artifact. What is lost is the letter of the clause, and
with it the ability to say the arms were run at the dtype the bar named. Training was fp32, scoring
was fp16.

## What this closes, and what it explicitly does NOT

**CLOSED:** GD-H1 ("LoRA fine-tuning brings distinguishability below 0.6679 at mode `both`") is
FALSIFIED at this corpus and this scale. Together with the sampling sweep — which falsified GD-H4,
since no temperature setting beat the null — **both generator levers this family could reach have now
been pulled and neither works.** The generator remains worse than a 3-mer Markov chain.

**NOT CLOSED, and do not let the negative compress into these:**

- **A miss does NOT separate "fine-tuning does not work" from "not enough data."** 3,737 windows is
  ~560 kb for a 1.2B-parameter model, and the ledger itself rates "is that enough" at **LOW**
  confidence. Whole-genome continued pre-training is this family's recorded candidate 2 and remains
  an open **Pending Decision** — a scope call for the user, not an executor task.
- **The frozen method is one point in a space.** r=16, alpha=32, lr 1e-4, 3 epochs. Loss was still
  falling when training stopped. A different setting is a **NEW ARM to be declared**, not a re-read
  of this one — and per the innovation rails it may not be killed for redundancy, only on its own
  executed evidence.
- **GD-H5 is STRENGTHENED, NOT CONFIRMED.** That the weights moved the wrong way is consistent with
  the tokenizer being the limit rather than the weights, but nothing here tests the tokenizer.
- **This is n=120 per arm at one seed**, the recorded geometry — not a powered comparison of two
  generators. The direction is clear and the margin (+0.0271) is not bounded by a CI here.
- **It says nothing about the forward/inverse cells**, which are a different regime (constructed,
  molecular) and are where this project's measured wins are.

## The forecast, scored

The pre-registration put **P(bar met) = 0.2** and this came in at the 0.8. That is the first forecast
I have logged in this project that was *not* falsified — the two in the fragment-null
pre-registration both were (the ~80%-fail prediction, and the claim that clearing the interval would
need ~1,190 seeds). One hit out of three; treat my priors accordingly.

## Cost, measured

**15.7 minutes wall** on a free Kaggle T4 — against the ledger's "+3 days". That estimate was a
GTX-1070-at-fp32 figure (consumer Pascal runs fp16 at ~1/64 of fp32), not a property of the job. The
corpus is ~560 kb and the T4 replaced a ~12.4 s/sequence CPU rate. **Recording this because the
"3 days" figure is exactly what kept this lever un-pulled**, and it was wrong by two orders of
magnitude.
