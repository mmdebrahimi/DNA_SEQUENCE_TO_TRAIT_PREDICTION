# North-star distance — strategy reframe (brainstorm, 2026-06-17)

**Question:** how far is the project from the north star — "input a DNA sequence → rough map of what its
parts do (phenotype + properties)"? Adversarially reviewed; the review CORRECTED the first answer.

## The load-bearing correction (verified)
The learned/general decoder is **NOT a gated/untested frontier — it is a CLOSED NEGATIVE on free data.**
Three de-confounded tests, same result every time (genomic-FM embeddings learn lineage / population
structure, NOT mechanism):
1. E. coli cipro AMR — NT AUROC 0.914 beat k-mer but **lost to QRDR knowledge-baseline 0.943**; within-lineage ≈ chance.
2. AMR-suite clonality — raw sens/spec inflated; cluster-weighted collapses (kleb cipro 0.967→0.5).
3. Arabidopsis flowering time — PlantCaduceus on a real T4, 3 seeds, n=1003: **within-group r² −0.13** (all seeds), structure-only +0.038 wins. H2 **falsified in its best-designed test** (`wiki/phase2_arabidopsis_result_2026-06-12.md` + `embedding_niche_cross_domain_synthesis_2026-06-12.md`).

The earlier framing ("pursue the learned decoder via the parked Path-B GPU run or an acquisition") was
**stale** — Path-B already ran and failed. The CLAUDE.md Path-B entry ("GPU run is the only remaining
step") is itself out of date.

**Honest qualifier (don't overclaim):** closed for the *current free-public foundation-embedding
strategy* — NOT logically closed forever. Precise rule: *no learned decoder is justified without a NEW,
named, label-clean substrate + a PRE-REGISTERED classical-baseline / de-confounding gate.*

## The reframed north star (what's actually achievable)
"Rough map of what each part does" = predominantly an **EVIDENCE-TIERED MOLECULAR-FUNCTION / ANNOTATION
map**, NOT a learned phenotype predictor. Three contributors:
- **Structural annotation** (where the genes/features are) — solved by external tools (Bakta/Prokka).
- **Homology/domain MIDDLE LAYER** (the bulk, previously missed): Pfam / eggNOG / COG / orthology /
  gene-neighborhood / structure-PLM transfer → infers molecular function for *uncatalogued* genes,
  **needs NO phenotype labels** (sidesteps the binding label constraint), GTX-860M/free-feasible.
- **Curated-determinant PHENOTYPE cells** (AMR/TB) — bounded validated overlays, only where a calibrated
  determinant + measured-label validation exists.

"Rough phenotype map across arbitrary traits" is a **category error** — most of the map is molecular
function, not phenotype. Phenotype claims attach only where a validated determinant cell exists.

## Distance verdict
- **Rough functional-annotation map: NEAR** — and bigger than first credited (the homology middle layer is
  the unlock; no labels needed).
- **Learned phenotype decoder for the unknown: CLOSED NEGATIVE on free data** — not the path.

## Recommended single highest-leverage next move
A **scoping spike, not a build**: define the **evidence-tier schema** (phenotype-calibrated-determinant /
molecular-function / pathway / homology-only-hypothesis / unknown) + a **hard UX/eval gate**, then
prototype the map on **2-3 genomes** (one E. coli, one TB/C. auris, one harder case) and judge it against
the gate. **Gate (load-bearing — else it drifts into catalog-stacking busywork):** can a user read the
2-3 maps and understand operons/genes/domains/pathways/determinants **+ uncertainty**, WITHOUT mistaking
association metadata for phenotype prediction, with **unknowns kept visible** (not hidden behind weak
annotations)?

## Secondary positions
- **TB cell = product coverage + validation infrastructure**, NOT a scientific decoder frontier. Scientific
  value comes ONLY from the post-2023 INDEPENDENT gold set, never the CRyPTIC knowledge-baseline. Finish it
  as engineering; frame honestly.
- **TransPred / ecoref = parked feasibility CARD**, not a parallel move. Execute only if label-clean +
  provenance-disjoint + de-confoundable + classical-baseline-not-at-ceiling can be shown FIRST.

## Open design questions for the prototype (before adding any source)
- Genome-map unit: whole-genome browser vs gene table vs operon/neighborhood view vs JSON-first report?
- Which 2-3 genomes best test the UX?
- The EXACT evidence tiers (define before adding sources, else "pathway"/"homology"/"determinant" blur).
- Acceptable unknown rate — a useful rough map must PRESERVE unknowns, not hide them.
