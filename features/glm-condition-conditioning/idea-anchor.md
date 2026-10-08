# Idea Anchor (DRAFT — parked for ratification) — G-C: Condition-conditioned expression oracle

> **STATUS: DRAFTED, NOT ANCHORED.** `/idea-anchor` output is **user-confirmed by standing directive**, so
> Soraya drafts it and parks. Nothing downstream treats this as anchored until you ratify or redirect it.
> Drafted 2026-10-07 as part of the GLM portfolio decomposition
> (`plans/GLM_Portfolio_Decompose_2026-10-07.md`). Ledger: `project_state/glm-condition-conditioning-2026-10-07.md`.

## 1. Formal Rephrase

Determine whether one model taking (sequence, growth medium) predicts expression better than two independent per-medium models, on held-out position-blocked data from the same fragment library measured in LB and M9.

## 2. Fundamental Clarifications

Each carries Soraya's drafted answer, per draft-then-ratify — you RATIFY or REDIRECT, you do not author
these from scratch.

1. Is two conditions enough for a null result to mean anything? (Drafted answer: no — a null at n=2 media is weak evidence against conditioning in general, and the ledger says so. It is still worth running because a POSITIVE at n=2 is informative.)
2. Does a medium INDICATOR count as conditioning, or must conditioning modulate the representation? (Drafted answer: start with the indicator. If it wins, the richer version is justified; if the indicator cannot win, the richer one probably cannot either.)

## 3. Current Assumptions

- That the shared component dominates, so a shared trunk should help. The disattenuated cross-condition correlation 0.789 supports this.
- That the fragment substrate is usable despite its low 0.793 ceiling. It is the ONLY two-condition substrate available, so this is a constraint, not a choice.
- That a gain from sharing is about conditioning rather than about seeing twice the rows. This is the family's main confound and needs its own null.
- That LB vs M9 is a biologically meaningful contrast rather than two near-identical regimes.

## 4. Blunt Opinion

**The strongest family on evidence and the weakest on ambition.** The motivating measurement is already positive and properly noise-matched, so this will probably produce a small real gain. But a small real gain across two growth media is a long way from 'name a trait' — do not let its tidiness make it look more central than it is. Its honest role is a HEDGE: if G-A closes negative, this family still has a result.

The live risk is mundane and will bite if ignored: the shared model trains on both media, so it sees roughly twice the rows. A naive two-vs-one comparison will show a gain that is about DATA VOLUME and not about conditioning at all. The medium-shuffled null is therefore not a nicety, it is the experiment. Without it the headline is uninterpretable.

## 5. Recommended Next Step

`/probe glm-condition-conditioning` — this is a code-touching idea that names real modules and real artifacts, so repo
grounding comes before any plan. The ledger's short-term actions are already written against the existing
loaders, so a probe should verify those entry points rather than re-explore the substrate.

```
/probe glm-condition-conditioning
```
