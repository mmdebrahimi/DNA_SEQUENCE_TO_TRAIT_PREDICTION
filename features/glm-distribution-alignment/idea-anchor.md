# Idea Anchor (DRAFT — parked for ratification) — G-E: Distribution alignment (generator <-> oracle)

> **STATUS: DRAFTED, NOT ANCHORED.** `/idea-anchor` output is **user-confirmed by standing directive**, so
> Soraya drafts it and parks. Nothing downstream treats this as anchored until you ratify or redirect it.
> Drafted 2026-10-07 as part of the GLM portfolio decomposition
> (`plans/GLM_Portfolio_Decompose_2026-10-07.md`). Ledger: `project_state/glm-distribution-alignment-2026-10-07.md`.

## 1. Formal Rephrase

Measure which reference sequence distribution the generator's output falls in — designed-grid-like, genomic-like, or neither — and report an oracle score on that distribution, so the loop's proposal half and scoring half are validated on the same thing.

## 2. Fundamental Clarifications

Each carries Soraya's drafted answer, per draft-then-ratify — you RATIFY or REDIRECT, you do not author
these from scratch.

1. If generated sequence matches NEITHER reference, do we get a new oracle substrate (data work) or constrain the generator (modelling work)? (Drafted answer: constrain the generator first — it is cheaper and conditional prompting is already wired. But this is a direction call.)
2. Should 'distribution' mean k-mer composition, positional structure, or discriminator separability? (Drafted answer: discriminator separability, because that is what the loop's own falsifier already uses — but report all three when they disagree, since disagreement is itself the finding.)

## 3. Current Assumptions

- That the existing naturalness discriminator can separate grid-like from genome-like sequence. It was built for a different axis (generated vs natural), so this is not automatic.
- That 'which distribution' has a reasonably clean answer. It may not — hence 'neither' is an explicitly allowed outcome.
- That matching distributions is sufficient for loop soundness. It is necessary; the oracle also has to be ACCURATE there, which is G-A's job.

## 4. Blunt Opinion

**The most important of the five and the least concrete, which is a dangerous combination.** It is the family that decides whether the loop means anything: today the oracle is validated on grid-like sequence, measured to FAIL on genome-like sequence, and nobody has characterised what the generator emits. Both halves can look like they work while the loop is unsound.

The risk is that it stays a framing exercise. The antidote is already in the ledger: candidate 2 — measure how distinguishable the two REFERENCE distributions are from each other — needs no predecessor and could shrink the whole family to a footnote. If grid and genomic sequence are barely separable, the alignment worry is much smaller than it looks. **Do that first, before waiting on anything.** A family that waits on two predecessors while holding an unblocked, potentially family-closing measurement is mis-sequenced.

## 5. Recommended Next Step

`/probe glm-distribution-alignment` — this is a code-touching idea that names real modules and real artifacts, so repo
grounding comes before any plan. The ledger's short-term actions are already written against the existing
loaders, so a probe should verify those entry points rather than re-explore the substrate.

```
/probe glm-distribution-alignment
```
