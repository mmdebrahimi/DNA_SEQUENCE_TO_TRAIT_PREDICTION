# Dog masked reconstruction — model-SCALE test on Kaggle T4 (2026-07-31)

**User question:** "will a stronger GPU give us a better response?" **Honest answer:** the GPU itself
does not (same model = same logits); but the **bigger model** a GPU enables **does** — clearly and
monotonically. I predicted scale would NOT help (citing the project's ESM/Arabidopsis "scale is dead"
lessons). **The data proved me wrong for this task.**

## Result — scale helps, monotonically (leakage-free, per-base marginal NLL)

Same harness as F1′ (disjoint-fit Markov given its best k; NT scored by per-base marginals). Same
1200 bp canFam4 chr1 eval window, 200 masked 6-mer tokens (1200 bases), same Markov baseline
(best k=2, accuracy 0.2933, NLL 1.3561). Kaggle Tesla T4 (free), transformers 4.44.2.

| Model | per-base marginal accuracy | per-base NLL | NLL delta (Markov − NT) | verdict |
|---|---|---|---|---|
| NT-v2-100M | 0.3267 | 1.4635 | **−0.107** | ≈ parity (loses slightly) |
| NT-v2-500M | 0.5083 | 1.0668 | **+0.289** | **beats Markov** |
| NT-2.5B | 0.5525 | 0.9745 | **+0.382** | **beats Markov by more** |

Monotonic in model size on both metrics. Artifact: `wiki/dog_nt_scale_2026-07-31.json`; reproducible
notebook `scripts/kaggle_dog_nt_scale.py`.

## Verification (why this is trustworthy, not a fluke)

- **Cross-environment cross-check:** NT-100M NLL here (Kaggle, transformers 4.44.2, T4) = **1.4635**;
  the local F1′ isolated-env run (transformers 4.30.2, CPU, different window) = **1.4474**. Two
  independent environments, same model, ~same number → the harness is correct.
- **The 100M ≈ parity here reproduces the F1′ finding** (100M does not beat Markov) — so the smaller-model
  result is stable and the ladder is internally consistent.
- **Not a 2.5B-fp16 artifact:** NT-500M ran in fp32 and already beats Markov by +0.289; the win does not
  depend on the 2.5B's fp16 path. (v2 remote code breaks in fp16 at forward → 100M/500M forced fp32; 2.5B
  runs fp16 — a per-model dtype retry handles both.)

## Honest reconciliation with the "scale is dead" prior

The project's prior scale-negatives (ESM2 peaks at 650M then regresses on ProteinGym; Arabidopsis
"do NOT scale embeddings on a bigger GPU") were on **DOWNSTREAM, de-confounded PHENOTYPE / variant-effect**
tasks — where the bottleneck is signal-vs-structure, not model capacity. This test is the model's
**NATIVE self-supervised objective** (masked reconstruction). It is entirely expected that a bigger
masked-LM reconstructs its own training objective better — that is what more parameters + pretraining buy.
**So scale helping HERE does not contradict scale failing on downstream phenotype tasks — they are
different regimes.** My in-conversation prediction wrongly transferred the downstream-task prior to the
native-objective task.

## What this does and does NOT mean

- **DOES:** NT-500M/2.5B learn real dog-sequence structure beyond local composition (they beat a low-order
  Markov chain that only captures composition), and capacity matters for that. The "world model" is not
  vacuous on its native task.
- **Does NOT:** prove usefulness for DECODING (phenotype). A global reconstruction win over a 2-mer Markov
  is expected + is "calibration, not biology" — reconstruction skill ≠ phenotype-prediction skill (the very
  gap the closed embedding-vs-phenotype negatives measured). Whether this translates to a useful decoder
  needs (a) F2 region-stratification (is the win concentrated in conserved/functional regions?) and (b) a
  downstream task with a de-confounded label.

## Scope + next
Smoke scale (1200 bp, one window, no CI). F2: region-stratified full-chromosome sweep (coding/intergenic/
conserved) at NT-500M/2.5B on Kaggle — test whether the scale win is STRUCTURED (the scientific signal)
or flat (mere calibration). Infra pinned: transformers 4.44.2, machine_shape NvidiaTeslaT4, v2 models
fp32 forward, 2.5B fp16. Free (token-only, no dollars).
