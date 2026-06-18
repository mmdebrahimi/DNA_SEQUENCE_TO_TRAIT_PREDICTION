---
title: Rough genome map (evidence-tiered functional annotation)
slug: genome-map
status: intake
created: 2026-06-17
updated: 2026-06-17
related_plans: []
related_tests: []
related_decisions: [wiki/north_star_distance_brainstorm_2026-06-17.md, wiki/genome_map_idea_anchor_prompt_2026-06-17.md]
related_modules: [dna_decode/eval/amr_rules.py, dna_decode/organism_rules/, dna_decode/data/calibrated_amr_rules.json]
---

# Rough genome map (evidence-tiered functional annotation)

## Summary
A personal tool that turns a single microbial genome into an HONEST, evidence-tiered map of what its
parts likely do — molecular function (homology/annotation) + curated-determinant phenotype overlays —
where every feature carries an explicit evidence tier and the unknown rate is reported. The achievable,
label-free form of the dna_decode north star ("DNA → what its parts do"); explicitly NOT a learned
genotype→phenotype predictor (that arm is a closed negative on free data).

## Problem
Understanding "what the parts of a genome do" today means running several external annotation tools and
reading their raw output, which (a) BLURS confident curated function with weak homology guesses, (b)
OVER-STATES phenotype (association/annotation metadata read as causal), and (c) HIDES how much is
genuinely unknown. The project also has validated determinant/AMR cells but no per-genome view that
situates them among the rest of a genome's features. There is no single honest, confidence-aware view to
explore an unfamiliar genome or to sanity-check the existing determinant calls in genomic context.

## Desired Outcome
Point the tool at ONE microbial genome and get a per-feature map where: each feature carries an explicit
evidence tier; PHENOTYPE claims appear ONLY where a validated determinant cell fires; and the unknown rate
is a reported headline. This enables the developer to (1) understand an unfamiliar genome at a glance and
(2) QC / sanity-check the existing AMR/determinant calls in context. Usefulness is judged against a UX/eval
gate (defined in `design`), NOT a learned-model metric.

## Repo Fit
- REUSE: the project's honesty/evidence-tier discipline (SUSPEND/PLUMBING/abstain rails, lineage
  disclosure) + the determinant PHENOTYPE cells — frozen AMR surface (`dna_decode/eval/amr_rules.py` +
  `dna_decode/data/calibrated_amr_rules.json`) and the non-frozen organism cells
  (`dna_decode/organism_rules/`: TB; the fungal cell). NO edit to the frozen AMR surface.
- NEW: an annotation/homology layer driven by FREE external tools (Bakta/Prokka structural annotation;
  Pfam/hmmer + eggNOG/COG homology-function) following the existing Docker / D:-cache patterns.
- `viz/` already exists with a deferred genome browser (out of v1 scope).

## Scope
**In v1 (the scoping spike):** an evidence-tier schema; a per-genome JSON-first tiered map (+ a flat
feature table); the phenotype overlay from existing determinant cells; molecular-function / homology tiers
from free tools; a reported unknown rate; a UX/eval gate; a prototype on 2-3 genomes; a GO/NO-GO verdict on
whether to invest in catalog integration.

**Out of v1:** any learned/embedding phenotype model (closed-negative arm); a visual genome browser
(deferred `viz/`); catalog integration at scale (gated on the spike verdict); non-microbial / eukaryotic
genomes; phenotype claims outside the validated-determinant wall; multi-genome comparison (the
decision-support intent was NOT chosen — this is an exploration/QC tool).

## Actors
- **Primary:** the developer (Farshad) — exploring / QC-ing a single microbial genome. Solo hobby
  project; NO external user. "Useful" therefore means "lets me understand + sanity-check a genome I don't
  already know," not "ships to a microbiologist."

## Workflow / Behavior
_(deferred to `design`)_

## Rules And Decisions
_(deferred to `design`)_

## Edge Cases
_(deferred to `design`)_

## Acceptance Criteria
_(deferred to `design` / `spec`)_

## Test Scenarios
_(deferred to `spec`)_

## Technical Notes
_(deferred to `design` / `technical-plan`)_

## Repo grounding

### Captured by: intake @ 2026-06-17
- Files read: dna_decode/eval/amr_rules.py (referenced — frozen AMR surface), dna_decode/organism_rules/ (TB cell built this session), dna_decode/data/calibrated_amr_rules.json, wiki/north_star_distance_brainstorm_2026-06-17.md, wiki/genome_map_idea_anchor_prompt_2026-06-17.md
- Key claims:
  - The phenotype tier reuses the EXISTING determinant cells (frozen `amr_rules.py` + non-frozen `organism_rules/`); the feature must NOT edit the frozen surface.
  - The learned genotype→phenotype arm is a CLOSED NEGATIVE on free data (3 de-confounded tests) — out of scope by decision, not omission.
  - The honesty discipline (tier labels, phenotype-behind-the-validated-wall, reported unknown rate) is the differentiated contribution; raw annotation tools (Bakta/eggNOG) already exist, so the spike must test differentiation-over-raw-Bakta.

## Changelog
- 2026-06-17 — record created (cascaded from intake; no template present, built from canonical sections)
- 2026-06-17 — intake completed (Summary/Problem/Desired Outcome/Repo Fit/Scope/Actors populated; primary intent ratified = personal exploration/QC tool)
