# Decision — Does NT-frozen-mean-pool embedding have a niche in E. coli AMR? — 2026-06-05

> Strategic decision forced by the decisive Phase-2 falsifier result. Class (e) research/decision memo
> (no build). Grounded in the full 4-drug AMR evidence set. Verdict + next-epoch direction below.

## Verdict: **NO demonstrable niche for frozen-mean-pool NT embeddings in E. coli AMR.**

The architecture should NOT be pursued further for AMR. For AMR phenotypes, the tool should use
**mechanism features (AMRFinder QRDR-POINT / acquired-gene calls)**, which are both stronger and
biologically interpretable. This is evidence-complete, not a judgment call (reasoning below).

## The evidence (all 4 E. coli AMR drugs in scope)

| Drug | Mechanism class | Embedding result | Status |
|---|---|---|---|
| **ciprofloxacin** | concentrated (QRDR point mutations) | NT-XGBoost **0.914** vs QRDR-POINT knowledge baseline **0.943** (−2.9 pp, CI [−9.0,+2.9]); within-lineage NT concordance **0.605 (p=0.365, chance)** vs POINT **1.000 (p<0.001)** | **FAIL vs knowledge baseline** — and the ONLY de-confounded AMR substrate that exists |
| **ceftriaxone** | concentrated (plasmid β-lactamase) | proper N=49 cohort is geography/lineage **CONFOUNDED** (R≈USA, S≈Africa/India, 1 shared lineage) | no honest test buildable; the earlier 0.833 "pass" was a reused cipro-mini smoke cohort, not a real de-confounded test |
| **tetracycline** | distributed (efflux + ribosomal protection) | NT **0.400** anti-predictive vs k-mer 0.722 (EP2 2026-05-17) | **FAIL** |
| **gentamicin** | distributed (aminoglycoside transferases) | cohort 2R/132S — infeasible (2026-05-18 strict-MIC census) | no substrate |

## Why this is decisive (not just "cipro failed")

1. **The cleanest possible test failed.** cipro is the single AMR drug with a de-confounded cohort (passed
   the lineage/geography de-confound gate: 6 shared R/S MLST lineages, country+year non-aliasing). On that
   cohort the embedding lost to the domain-knowledge baseline AND was at chance within-lineage → its
   apparent skill was lineage/genome-content, not the resistance mechanism.
2. **No second de-confounded AMR substrate is buildable.** cef is confounded; gent is infeasible; tet's
   cohort exists but the mechanism is distributed and NT was anti-predictive there. So there is no
   remaining clean AMR test that could overturn the cipro verdict.
3. **The biology explains it.** AMR resistance is localized (QRDR point mutations; acquired genes).
   Mechanism-feature baselines encode that directly; frozen whole-genome mean-pooling dilutes a localized
   signal across the genome. Concentrated mechanisms → the knowledge baseline wins; distributed
   mechanisms → NT was anti-predictive. Neither class leaves room for the embedding.

## Consequences (roadmap)

- **AMR decoder = mechanism-feature tool, not embeddings.** `dna_decode/eval/point_baseline.py` (QRDR-POINT)
  is the AMR-prediction substrate; AMRFinder is the engine. This is also the interpretable, honest-output
  shape the north star wants.
- **Phase-2 (multi-drug AMR embedding validation) is CLOSED** with this finding. The reusable infra it
  produced — de-confound gate, drug-agnostic CI-aware falsifier, POINT baseline, within-lineage diagnostic
  — is the durable value and transfers to any future phenotype.
- **The embedding architecture's remaining honest frontier (Phase 4+, NOT now):** a phenotype with BOTH
  (a) sampling-INDEPENDENT labels (lab measurement, not clinical-site category — per the pathotype/cef
  confound lessons) AND (b) NO curated mechanism/knowledge baseline for the embedding to lose to. That is
  where embeddings could add value k-mer/POINT can't. It is gated on a **de-confounded labeled substrate**
  — the binding constraint that has now blocked the project 3× (pathotype circular labels, cef geography
  confound, cipro within-lineage). Finding that substrate is a research/`/idea-anchor` problem, not a
  modeling one, and is the correct next epoch IF embeddings are pursued at all.

## What is NOT recommended
- Building a new cef-S / gent cohort to retry concentrated-mechanism embeddings — the cipro result already
  shows the knowledge baseline wins on concentrated mechanisms; a new cohort would re-confirm, not overturn.
- Re-tuning NT pooling (mean→mean+max, per-gene windows) for AMR before a substrate exists where the
  embedding even *could* beat a knowledge baseline. Architecture tuning without a winnable test is wasted.

## Provenance
cipro: `wiki/ciprofloxacin_falsifier_2026-06-05.{md,scores.json}` + `wiki/ciprofloxacin_within_lineage_diagnostic_2026-06-05.md`.
cef: `plans/cef_falsifier_brainstorm.md`. tet: `wiki/EP2_cef_tet_verdict_2026-05-17.md`. gent: `wiki/bvbrc_strict_mic_4drug_census_2026-05-18.md`.
Lesson: `~/.claude/.../memory/feedback_embedding_vs_knowledge_baseline_and_within_lineage.md`. Driven via `/soraya` + 3× `/brainstorm`.
