# Idea Anchor (DRAFT — parked for ratification) — G-B: Self-distillation (blocked by G-A)

> **STATUS: DRAFTED, NOT ANCHORED.** `/idea-anchor` output is **user-confirmed by standing directive**, so
> Soraya drafts it and parks. Nothing downstream treats this as anchored until you ratify or redirect it.
> Drafted 2026-10-07 as part of the GLM portfolio decomposition
> (`plans/GLM_Portfolio_Decompose_2026-10-07.md`). Ledger: `project_state/glm-self-distillation-2026-10-07.md`.

## 1. Formal Rephrase

After a higher-capacity oracle exists, determine whether training it on randomly perturbed sequences labelled by its own predictions improves held-out performance on REAL MEASURED data — the AlphaGenome phase-2 method, adapted.

## 2. Fundamental Clarifications

Each carries Soraya's drafted answer, per draft-then-ratify — you RATIFY or REDIRECT, you do not author
these from scratch.

1. Does a G-A negative close this family or merely defer it? (Drafted answer: it depends on WHY. A measured capacity ceiling closes it; an encoder that is merely weak only defers it.)
2. Is off-distribution SMOOTHNESS an acceptable deliverable if accuracy does not improve? (Drafted answer: yes, but it must be stated as smoothness and never as accuracy — and smoothness is hard to validate without labels off-distribution, which is the whole problem.)

## 3. Current Assumptions

- That a conv encoder's capacity makes distillation non-degenerate. UNTESTED for the CNN; this is the premise the block rests on.
- That uniform random point substitution resembles the perturbations that matter. It almost certainly does not — a generator's candidates are not uniform.
- That the information bound (a student cannot exceed its teacher) is the real limit. Already measured for the non-linear case: 0.5933-0.5992 against a teacher at 0.5963.

## 4. Blunt Opinion

**Correctly blocked, and the block is the most valuable thing this family currently owns.** The method is the most exciting idea in the AlphaGenome programme and that is exactly why it needs a guard: the temptation is to build it now. It was measured degenerate on the current oracle at agreement 1.0000 — not approximately, exactly — so early work is predicted-void.

Be clear-eyed about the ceiling even once unblocked. Distillation cannot create information the teacher lacks, and that is measured here, not borrowed from theory. So the realistic upside is self-consistency under perturbation, which is much less than 'manufacture a variant-effect signal' sounds like. If G-A delivers only a marginal encoder, this family should close rather than be attempted on a weak teacher.

## 5. Recommended Next Step

`/probe glm-self-distillation` — this is a code-touching idea that names real modules and real artifacts, so repo
grounding comes before any plan. The ledger's short-term actions are already written against the existing
loaders, so a probe should verify those entry points rather than re-explore the substrate.

```
/probe glm-self-distillation
```
