# Feature — G-C: Condition-conditioned expression oracle

<!-- feature-record: intake -->
- **Slug:** `glm-condition-conditioning`
- **Status:** `intake` — idea-anchor DRAFTED and parked for ratification; not yet `spec-ready`
- **Created:** 2026-10-07
- **Family:** G-C of the GLM portfolio · ledger `project_state/glm-condition-conditioning-2026-10-07.md`
- **Decomposition:** `plans/GLM_Portfolio_Decompose_2026-10-07.md`
- **Umbrella:** `project_state/glm-umbrella-2026-10-07.md`

## Problem

Determine whether one model taking (sequence, growth medium) predicts expression better than two independent per-medium models, on held-out position-blocked data from the same fragment library measured in LB and M9.

## Primary actor and outcome

- **Actor:** the project itself (this is an internal capability, not a user-facing surface). The eventual
  downstream actor is a user of the generative loop who names a phenotype and receives ranked edits.
- **Outcome that should improve:** whether the GLM's scoring or proposal half can be trusted on the
  sequence distribution it will actually be used on.

## Acceptance

The checkable bar lives in the ledger's `## MVP Criteria` and is **deliberately not duplicated here** —
one definition, one place. Summary: the family is done when its question is ANSWERED with a committed
artifact, which a measured negative satisfies exactly as a positive does.

## In scope

See the ledger's `## Goal Hierarchy → Short-term`. Those actions are written against loaders that already
exist (`dna_decode/glm/`), not against hypothetical infrastructure.

## Explicitly NOT in scope

Carried from the decomposition's anti-scope section, which exists because each item is a mistake already
made once in this repo:

- Comparing a 2-replicate mean against a single replicate. That unmatched form reads 0.5025 and INVERTS this family's founding answer; matched it is 0.4062 vs 0.4676.
- Reporting a two-vs-one gain without a medium-shuffled null. The shared model sees ~2x the rows, so an unqualified gain is attributable to data volume rather than to conditioning.

## Repo grounding

Files read while drafting: `dna_decode/glm/expression.py`, `dna_decode/glm/genomewide.py`,
`dna_decode/glm/selection.py`, `dna_decode/glm/falsifier.py`, `wiki/glm_genomewide_oracle_2026-10-07.json`,
`wiki/glm_selfdistill_probe_2026-10-07.json`, `wiki/glm_generator_falsifier_REAL_2026-10-07.json`,
`wiki/alphagenome_lessons_actionable_2026-10-07.md`.

Key claim: every number in this family's ledger was verified against a committed artifact at init rather
than restated from memory — see the ledger's `## Empirical Concerns` for what was checked and the one
self-contradiction the check found.

## Next step

Ratify or redirect the drafted idea-anchor, then:

```
/probe glm-condition-conditioning
```
