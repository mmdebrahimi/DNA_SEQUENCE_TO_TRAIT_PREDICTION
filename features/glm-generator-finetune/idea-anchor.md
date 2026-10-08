# Idea Anchor (DRAFT — parked for ratification) — G-D: Generator fine-tuning (the un-pulled lever)

> **STATUS: DRAFTED, NOT ANCHORED.** `/idea-anchor` output is **user-confirmed by standing directive**, so
> Soraya drafts it and parks. Nothing downstream treats this as anchored until you ratify or redirect it.
> Drafted 2026-10-07 as part of the GLM portfolio decomposition
> (`plans/GLM_Portfolio_Decompose_2026-10-07.md`). Ledger: `project_state/glm-generator-finetune-2026-10-07.md`.

## 1. Formal Rephrase

Determine whether LoRA fine-tuning GENERator-v2-prokaryote-1.2b on real E. coli promoter windows brings its generated sequence below the distinguishability of a 3-mer Markov chain, which it currently loses to at default sampling.

## 2. Fundamental Clarifications

Each carries Soraya's drafted answer, per draft-then-ratify — you RATIFY or REDIRECT, you do not author
these from scratch.

1. If the queued sampling sweep already clears the bar, does this family still run? (Drafted answer: it becomes much lower priority but not void — sampling and weights are different levers and a sweep win does not tell us fine-tuning adds nothing.)
2. Is whole-genome continued pre-training acceptable if 4,340 promoter windows proves too few? (Drafted answer: it is the obvious fallback but materially more compute, so it is a scope call rather than a technical one.)
3. Is k-mer-plus-positional distinguishability the right target at all? (Drafted answer: it is necessary but NOT sufficient — published audits name long-range organisation as the real failure, which these features cannot see.)

## 3. Current Assumptions

- That the failure is ZERO-SHOT and movable by fine-tuning, rather than architectural.
- That LoRA at fp32 fits 8 GB for a 1.2B model. Untested on that card.
- That 4,340 promoter windows is enough to fine-tune a 1.2B model. Low confidence; this is a very small corpus.
- That the custom remote tokenizer is not itself the limiting factor.

## 4. Blunt Opinion

**The most likely of the five to produce a real improvement, and the most likely to die for an unglamorous reason.** 4,340 windows is a tiny corpus against 1.2B parameters; the plausible outcome is memorisation rather than generalisation, and the discriminator will notice.

Two things could make the family moot before it starts. If the sampling sweep clears 0.6679, the premise weakens a lot. And if the real failure is long-range organisation, then a k-mer plus positional discriminator cannot see the thing that matters, so passing its bar would not mean the generator is good — it would mean we chose a measurable target instead of the right one. Worth saying plainly: a PASS here is necessary, not sufficient.

Do not start training before the sweep reports. Two levers in one run cannot be attributed.

## 5. Recommended Next Step

`/probe glm-generator-finetune` — this is a code-touching idea that names real modules and real artifacts, so repo
grounding comes before any plan. The ledger's short-term actions are already written against the existing
loaders, so a probe should verify those entry points rather than re-explore the substrate.

```
/probe glm-generator-finetune
```
