# Feature — G-A: Learned representation for the expression oracle

<!-- feature-record: intake -->
- **Slug:** `glm-learned-representation`
- **Status:** `intake` — idea-anchor DRAFTED and parked for ratification; not yet `spec-ready`
- **Created:** 2026-10-07
- **Family:** G-A of the GLM portfolio · ledger `project_state/glm-learned-representation-2026-10-07.md`
- **Decomposition:** `plans/GLM_Portfolio_Decompose_2026-10-07.md`
- **Umbrella:** `project_state/glm-umbrella-2026-10-07.md`

## Problem

Replace the expression oracle's hand-specified sequence features with a learned multi-resolution encoder feeding a trivial head, and determine whether it predicts promoter strength on REAL genomic sequence better than a single GC-content number does.

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

- Reporting a RANDOM-split number as the headline. Random splits inflate ~3x in this repo (measured four times), and tiles within a peak overlap heavily, so the only honest split here is leave-peak-out.
- Normalising the tile arm by the FRAGMENT noise ceiling. Measured: 0.9410 vs 0.7928, and using the wrong one already flattered a result by 5 points once.

## Repo grounding

### Captured by: brainstorm @ 2026-10-12

Files read: soraya_runs/arms-ledger.json; soraya_runs/negative-results.json;
soraya_runs/precompact-brief.md; dna_decode/glm/encoder.py; dna_decode/glm/genomewide.py;
scripts/glm_encoder_gate.py; tests/test_glm_encoder.py; tests/test_glm_encoder_gate.py;
tests/test_glm_tile_headroom.py; features/glm-learned-representation/technical-plan.md;
features/glm-learned-representation/feature.md; project_state/glm-learned-representation-2026-10-07.md;
project_state/glm-umbrella-2026-10-07.md; plans/GLM_Portfolio_Decompose_2026-10-07.md;
wiki/glm_pyramid_arm_2026-10-11.md; wiki/glm_pyramid_arm_2026-10-11.json;
wiki/glm_learned_representation_2026-10-08.json; wiki/glm_encoder_gate_result_2026-10-08.md;
wiki/chatgpt_architecture_discussion_lessons_2026-10-09.md; CLAUDE.md

Key claims:
- [grounded] The ceiling criterion is ALREADY substrate-bound by its own wording. The `0.941` it names
  is the replicate agreement of GSE144621 150 bp tiles specifically (`noise_ceiling` in
  wiki/glm_pyramid_arm_2026-10-11.json; substrate string "GEO GSE144621 peak tiles, category=tile,
  len=150"), and the adopted G-A goal is likewise bound to real E. coli tiles + `leave_peak_out`
  (project_state/glm-learned-representation-2026-10-07.md, features/.../technical-plan.md). A
  constructed eukaryotic/human panel would have its OWN different ceiling, so it cannot satisfy this
  criterion as written -- the arm is MISFILED rather than the criterion being due for redefinition.
- [grounded] `multiresolution-pyramid-4-decades-UNRUN` cannot execute on this substrate at all, not
  merely at high cost. `dna_decode/glm/encoder.py:148 validate_widths` raises `EncoderDataError` for
  any kernel wider than the sequence; executed check: `validate_widths((10,100,1000,10000), 150)` ->
  "kernel width(s) [1000, 10000] exceed the sequence length 150". Leaving it `planned` therefore holds
  the criterion in a state the repo already knows is impossible.
- [grounded] `arms_ledger` has exactly three statuses (`planned`/`measured`/`killed`) and `Arm` carries
  no revival field (`name`, `cost_tier`, `status`, `metric`, `evidence`, `declared_at`, `resolved_at`).
  So a non-destructive "deferred/blocked-by-substrate" state does not exist, and a kill's revival
  trigger has to live in `soraya_runs/negative-results.json` (which exists).
- [grounded] `kill_arm` refuses only REDUNDANCY reasons (`_REDUNDANCY_REASON` matches bar-met /
  already-cleared / redundant / sibling / not-needed / goal-is-met), so a substrate-infeasibility kill
  is mechanically permitted -- the refusal is not a blocker here.
- [grounded] The 2026-10-11 memo contradicts itself: line 5 says the run "closes one of them with a
  number" while line 70 says "the arm measured here is not either declared arm" and the criterion stays
  `blocked:untried-arm`. A new arm was measured; NEITHER declared arm was closed.
- [grounded] "More scale is not the missing lever" (wiki/glm_pyramid_arm_2026-10-11.md:31) over-reaches
  the measurement, and contradicts that same artifact's own `honest_limits`. The supported claim is
  narrower: within 150 bp tiles, widening 1.10 -> 1.70 decades under the registered small CNN protocol
  did not help.
- [grounded] Splitting the criterion by substrate scope is an ACCEPTANCE-BAR change requiring user
  ratification, not a bookkeeping correction -- because the current criterion's terminal condition is
  defined against the GSE144621-derived 0.941 ceiling.
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
/probe glm-learned-representation
```

### Captured by: brainstorm @ 2026-10-07

Files read: plans/GLM_Portfolio_Decompose_2026-10-07.md, project_state/glm-umbrella-2026-10-07.md, project_state/glm-learned-representation-2026-10-07.md, project_state/glm-condition-conditioning-2026-10-07.md, features/glm-learned-representation/technical-plan.md, features/glm-condition-conditioning/technical-plan.md, dna_decode/glm/genomewide.py, dna_decode/glm/expression.py, dna_decode/glm/selection.py, dna_decode/glm/generate.py, scripts/glm_genomewide_oracle.py, tests/test_glm_genomewide.py, tests/test_glm_expression.py, tests/test_glm_generate.py, pyproject.toml, wiki/glm_genomewide_oracle_2026-10-07.json, wiki/glm_genomewide_oracle_result_2026-10-07.md, wiki/glm_expression_oracle_result_2026-10-07.md

Key claims:
- [grounded] The plan names the single-seed GC comparator `0.3000` four times and mentions
  optimizer / learning rate / seed list / batch size **zero** times; `fit_encoder(train, test, *, widths,
  epochs, seed)` (technical-plan.md:85) exposes no training protocol. With a 0.05 margin on one fixed
  leave-peak-out split that is a live search surface, and the GC baseline was itself fixed on seed 0's split.
- [grounded] The four MVP predicates (module exists / tests exit 0 / artifact exists / ledger row) do not
  require the experiment to have been VALID — a broken run writing a null-result artifact satisfies all four.
  The shuffled-label null that would prevent this lives in the gate script, not in the bar.
- [grounded] torch is absent from default deps (`biopython, pyyaml, scikit-learn, numpy, pandas, h5py,
  pyarrow, requests, pytest`) and present only in the `forward`/`ml` extras, while the plan's verification
  command is plain `uv run pytest`. Verified NOT latent repo-wide: dna_decode/glm/generate.py imports torch
  only inside functions (lines 48/62/135/174/219), so this is a defect the encoder module would introduce.
- [external] torch.nn.Conv1d expects `(N, C, L)` and is documented as potentially nondeterministic on CUDA;
  the plan claims determinism without pinning tensor layout or enforcing CPU in the gate. Latent, since the
  plan is CPU-first.
