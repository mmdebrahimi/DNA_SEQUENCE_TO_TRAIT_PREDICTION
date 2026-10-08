# Idea Anchor (DRAFT — parked for ratification) — G-A: Learned representation for the expression oracle

> **STATUS: DRAFTED, NOT ANCHORED.** `/idea-anchor` output is **user-confirmed by standing directive**, so
> Soraya drafts it and parks. Nothing downstream treats this as anchored until you ratify or redirect it.
> Drafted 2026-10-07 as part of the GLM portfolio decomposition
> (`plans/GLM_Portfolio_Decompose_2026-10-07.md`). Ledger: `project_state/glm-learned-representation-2026-10-07.md`.

## 1. Formal Rephrase

Replace the expression oracle's hand-specified sequence features with a learned multi-resolution encoder feeding a trivial head, and determine whether it predicts promoter strength on REAL genomic sequence better than a single GC-content number does.

## 2. Fundamental Clarifications

Each carries Soraya's drafted answer, per draft-then-ratify — you RATIFY or REDIRECT, you do not author
these from scratch.

1. Is a recorded NEGATIVE an acceptable terminal, or must this family stay open until a learned encoder actually wins? (Drafted answer: a negative IS acceptable and the bar is written that way; a bar only a win could clear makes the measured negative unreachable, and the negative is the likelier outcome.)
2. Is the target promoter STRENGTH (regression) or promoter PRESENCE (classification)? The tile set is 15,171 active / 31,542 inactive. (Drafted answer: keep regression for comparability with the 0.3000 baseline, and add classification as a diagnostic, not a substitute.)
3. Does 'multi-resolution' need to be in v1, or is single-resolution enough to answer the prior question? (Drafted answer: single first — it isolates learned-vs-fixed from multi-vs-single-resolution, and the second question is worthless if the first fails.)

## 3. Current Assumptions

- That a conv encoder has capacity ridge-on-fixed-features lacks. Architectural, safe.
- That 44,106 tiles at 150 bp is enough to train a small CNN without the overfitting that took 6-mers (4,096 features) down to 0.2022 — the WORST arm. That is a warning, not a detail.
- That sequence alone determines promoter strength on real genomic tiles at better than composition level. This is the actual bet and it may simply be false.
- That the 0.9410 tile ceiling is the right normaliser. Derived, and the error of using the fragment ceiling instead has already been caught once.
- That CPU is enough. Probably, since the data is small.

## 4. Blunt Opinion

**The prior here is against us, and the family should be run anyway for a different reason than hope.** A single GC number beating 4,096 6-mer features is not a neutral result — it says most of the linearly-available signal at this resolution IS composition, and a CNN is still a composition-plus-position model. Worse, 6-mers were the WORST arm, so added capacity actively hurt; that is the opposite of the pattern you want before committing to more capacity.

The bigger threat is not architecture, it is FRAMING: two thirds of the tiles are inactive, so Spearman over the whole set may largely be ranking noise among inactive tiles. If that is what is happening, no encoder fixes it and the right move is a classification framing.

What justifies the family is that it is cheap, it has a derived ceiling to measure against, and it GATES G-B by a measured degeneracy. Run it to close a question, not to win a bet. And do the linear-headroom diagnostic FIRST — it can close the family in an afternoon.

## 5. Recommended Next Step

`/probe glm-learned-representation` — this is a code-touching idea that names real modules and real artifacts, so repo
grounding comes before any plan. The ledger's short-term actions are already written against the existing
loaders, so a probe should verify those entry points rather than re-explore the substrate.

```
/probe glm-learned-representation
```
