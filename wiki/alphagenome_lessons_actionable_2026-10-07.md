# AlphaGenome → our GLM: six lessons, one measured as a FORBIDDEN build

**Date:** 2026-10-07 · **Sources:** the Kohli interview (`h5zkzon0gM4`) + a DeepMind research-scientist
technical talk (`rdBtxtcS4nM`, 52,417 chars), both read in full; architecture analysis at
`wiki/alphagenome_integration_analysis_2026-10-07.md`. Citation verified against primary sources: announced
**25 Jun 2025**, *Nature* **28 Jan 2026**, vol 649 pp 1206–1218 (the circulated summary said "Oct 2024 /
Nature 2024" — wrong on both).

## The one-line split

**The ARCHITECTURE does not transfer. The TRAINING METHOD does.** The 1 Mb context, the U-Net
compress/expand, the transformer tower and the 2D contact branch all exist for one purpose — long-range
enhancer–promoter interaction across hundreds of kb — and a bacterial promoter is 100–300 bp. Its eleven
output heads predict splicing, histone marks and chromatin contacts, which a bacterium does not have. Its
licence separately forbids training on its outputs, and the model is human/mouse. So anyone reading the
circulated architecture description is reading the half that is irrelevant to us.

What transfers is how they *train* and how they *validate*.

## The six lessons

| # | Lesson | Status here |
|---|---|---|
| 1 | Validate on the DOWNSTREAM task shape, not the upstream metric | **DONE today** — measured |
| 2 | Trivial head on a LEARNED representation | **Queued** — server Run 3 |
| 3 | Self-distillation on perturbed inputs | **MEASURED DEGENERATE — blocked until #2** |
| 4 | Fine-tune on your own data; don't consume zero-shot | **The un-pulled lever** — server |
| 5 | DNA-only input is a deliberate philosophy, not an omission | **Adopted now, zero cost** |
| 6 | "Prerequisite, not sufficient" | **Honesty rail** on the north star |

---

### 1. Validate downstream, not upstream — DONE, and it paid

Their distillation metric proves nothing on its own, so the claim is validated downstream on real measured
variant-effect benchmarks. Our analogue: the oracle's ρ ≈ 0.59 is an upstream metric, while the loop consumes
a *choice*. Measured today (`wiki/glm_oracle_selection_result_2026-10-07.md`): the oracle captures **79% of
the available selection value at k=10** across 5 held-out element axes, against GC's 13%, verdict
`ORACLE_READY_AS_SCORING_HEAD`. **The loop's scoring half is adequate; the generator is the bottleneck.**

### 2. The trivial head belongs on a LEARNED representation — and ours has no learned half

AlphaGenome puts the capacity in the representation and keeps the heads simple. Our oracle is the **trivial
head without the learned representation**: ridge on fixed features. So its numbers are a FLOOR. Closing that
is server Run 3 (a small CNN on the same 10,898 pairs, split by `leave_element_out`, never randomly).

### 3. Self-distillation — the most transferable idea in the programme, and we cannot use it yet

**The method:** take the pre-trained model, **randomly perturb the input sequence, and train on the model's
own outputs as labels.** Because the labels are predicted rather than measured, the held-out-region
restriction lifts and they train across the whole genome. The problem it solves is exactly ours: measured
data exists only for sequences someone assayed, so there is no direct signal for what an *edit* does — and
our loop will query the oracle on generated sequences off the measured distribution.

**I expected to queue this build. It is forbidden, and that is measured, not argued**
(`scripts/glm_selfdistill_probe.py`, `wiki/glm_selfdistill_probe_2026-10-07.json`):

| arm | ρ vs measured | agreement with teacher |
|---|---|---|
| teacher (ridge, real labels) | +0.5963 | — |
| student_linear (same basis, teacher labels) | +0.5963 | **1.0000** |
| student_mlp (richer basis, teacher labels) | +0.5933 … +0.5992 | 0.999 |
| student_gc (**non-vacuity control**) | +0.1736 | 0.232 |

Verdict `DEGENERATE_DISTILLATION_BLOCKED_UNTIL_LEARNED_REPRESENTATION`, identical at mutation rates 0.02,
0.05 and 0.10.

**Why, exactly:** a ridge teacher's predictions lie *in the span of its own features*, so a linear student
with the same basis recovers the same function no matter where in sequence space you draw the perturbed
samples. Agreement is **1.0000** — not approximately, exactly. There is nothing for the perturbation to
regularise until the student has capacity the teacher lacks.

**Two further things the probe establishes rather than assumes:**
- **The information bound is real.** A properly-trained non-linear student with a richer basis lands on
  **0.5933–0.5992** against the teacher's 0.5963 — it straddles it and cannot exceed it. "Distillation
  cannot create information the teacher lacks" is the caveat that must travel with this method, and here it
  is measured. (The MLP got 400 iterations with early stopping *on purpose*: at 80 it stopped on max_iter
  un-converged, and an under-trained student failing to beat the teacher would be evidence about the
  optimiser, not about the bound.)
- **The probe can detect a difference.** The GC-only student collapses 0.596 → 0.174. Without that control,
  "the student matched the teacher" would be unfalsifiable. Its 0.1736 is *identical* to GC fitted on real
  labels, which is the pipeline confirming itself: a student recovers the best function available in its
  basis, whatever the labels came from.

**Consequence: the sequencing is forced.** Run 3 (learned representation) is a **prerequisite** for the
distillation step, not a parallel improvement. And when it lands, the first half must not ship without the
second: at DeepMind the distilled model is validated on real measured benchmarks, and **adopting
distil-for-smoothness without validate-on-real-data is circular.** We have the validation half already —
`leave_element_out` on the 10,898, plus GSE144621 genome-wide as an independent substrate.

### 4. Fine-tune on your own data — their own advice, and the lever nobody has pulled

Asked by a plant researcher who wanted AlphaGenome as a feature extractor, the author said: not a silly idea,
but it was not trained on plant data, and **the better approach is to fine-tune on your own data using the
open weights rather than consuming the predictions as-is.**

Our measured generator failure is a **zero-shot** failure: GENERator-1.2B at default sampling scores
distinguishability **0.7706** against a 3-mer Markov chain's **0.6679** (lower is better) — a 1.2B model more
detectable than a trivial Markov chain. Two levers follow, and only the first is in flight:
- **sampling** (temperature / top-k) — server Run 1, in progress;
- **fine-tuning on ~4,000 E. coli promoter windows** — never attempted. LoRA at fp32 should fit the 1070's
  8 GB. This is the author's own prescription applied to our exact failure.

**Scope care:** their encouragement was about a eukaryote feature-extractor, and plants have splicing and
nucleosomes, so their case is *far* more favourable than bacteria. The transferable part is the principle,
not the endorsement.

### 5. DNA-only input is a deliberate philosophy — adopt it now, at zero cost

They rejected multimodal inputs **on generality grounds**, as a choice rather than an omission. Recording it
as a standing constraint so a future session does not "improve" the generator by feeding it expression or
chromatin channels. **Independently corroborated by our own measurement:** the orthogonal-modality hybrid
that beats ESM2 on 84–90% of ProteinGym proteins **lost** on CTX-M-14 (0.204 vs 0.352), because the premise
failed — ProSST alone was at chance there. A modality adds where the incumbent is *starved* of signal, not
where it is merely blind in principle.

### 6. "Prerequisite, not sufficient" — the sharpest honesty line in either source

The author is explicit that AlphaGenome does **not** predict disease or height: it predicts processes at the
level of the cell and specifically the nucleus, which is a **prerequisite** for a genetic change having any
effect at all. *Prerequisite* is the precise word — necessary, not sufficient.

That is exactly our north star's step (2) → (3) seam. The expression oracle predicts **expression**, which is
a prerequisite for a **trait**; ρ 0.59 and a 79% selection headroom are not trait claims. `eval/regime.py`
already encodes this boundary — here it is, stated by the authors about their own model.

## Smaller facts worth keeping

- **DNA methylation is not among the 11 outputs** (not available in the needed form at the time).
- **AlphaGenome does not beat Enformer on everything.** Asked about under-performing it on CAGE at 128 bp,
  the author says probably within statistical noise — a useful corrective to the 25-of-26 headline.
- Resolution is **1 bp AND 128 bp**, with 2048 bp contacts — not uniformly base-pair, as I first wrote.
- The multi-species-improves-human claim has **no ablation in the paper**; it is interview-only.

## What this changes about the plan

1. **Run 3 is promoted from "nice improvement" to "prerequisite"** — it gates lesson 3 entirely.
2. **Lesson 4's fine-tuning arm is named and queued** as the un-pulled lever on the generator failure.
3. **The scoring head is no longer the suspect.** Measured at 79% of achievable selection value, so effort
   belongs on the generator.
4. **No distillation code gets written until a learned representation exists** — and when it does, the
   real-data validation half ships with it or the result is circular.
