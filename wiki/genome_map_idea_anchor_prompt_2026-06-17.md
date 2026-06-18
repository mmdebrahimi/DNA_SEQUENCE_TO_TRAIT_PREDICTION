# `/idea-anchor` prompt — the "rough genome map" (evidence-tiered functional annotation) (2026-06-17)

**Why this exists:** the 2026-06-17 north-star brainstorm (`wiki/north_star_distance_brainstorm_2026-06-17.md`)
reframed the achievable north star. The learned genotype→phenotype decoder is a **closed negative on free
data** (3 de-confounded tests; embeddings learn lineage/structure, not mechanism). The achievable form of
"input DNA → rough map of what its parts do" is an **evidence-tiered molecular-function/annotation map**,
not a learned phenotype predictor. `/idea-anchor` is user-only → this file holds the paste-ready command +
Soraya's drafted framing for you to RATIFY.

---

## Paste-ready command

```
/idea-anchor Build a "rough genome map" — an evidence-tiered molecular-function annotation map for a microbial genome — as the achievable, honest form of the dna_decode north star.

PROBLEM: The north star is "input a DNA sequence -> a rough map of what its parts do (phenotype + properties)." Two hard facts now bound what that can be: (1) the LEARNED genotype->phenotype decoder (genomic-FM embeddings) is a CLOSED NEGATIVE on free data — 3 independent de-confounded tests (E. coli cipro AMR, the AMR-suite clonality audit, Arabidopsis flowering-time on a real T4) all show embeddings capture lineage/population structure, NOT mechanism (synthesis: wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md). (2) The binding constraint across every expansion attempt has been LABELS, not models. So the achievable "rough map of what each part does" is NOT a learned phenotype predictor — it is an EVIDENCE-TIERED MOLECULAR-FUNCTION / ANNOTATION map with three contributors: (a) STRUCTURAL annotation (where the genes/features are — Bakta/Prokka, a solved external step); (b) a HOMOLOGY/DOMAIN MIDDLE LAYER that infers molecular function for UNCATALOGUED genes WITHOUT any phenotype labels (Pfam/HMM, eggNOG/COG, orthology, gene-neighborhood, optionally structure/PLM transfer) — this is the bulk of "what does each part do" and it sidesteps the label bottleneck; (c) the project's existing CURATED-DETERMINANT PHENOTYPE cells (the frozen AMR decoder + the TB cell) as BOUNDED, VALIDATED overlays. "A rough PHENOTYPE map across arbitrary traits" is a CATEGORY ERROR — most of the map is molecular function, not phenotype.

GOAL: A SCOPING SPIKE, not a build. Produce (1) a precise EVIDENCE-TIER schema, (2) a hard UX/EVAL GATE that defines what a useful map is, then (3) a PROTOTYPE map on 2-3 genomes judged against that gate. Decide whether the output actually satisfies the north-star intent BEFORE integrating any catalog at scale.

THE LOAD-BEARING RISK (design around it, do not ignore): without a hard UX/eval gate this degenerates into CATALOG-STACKING BUSYWORK ("more sources added") that never moves the north star. The gate: a user can read the 2-3 maps and understand operons/genes/domains/pathways/determinants PLUS uncertainty, WITHOUT mistaking association/molecular-function metadata for phenotype prediction, and with UNKNOWNS kept VISIBLE (not hidden behind weak annotations). Second risk: tier BLUR — "pathway", "homology-hypothesis", and "phenotype determinant" must be defined BEFORE any source is wired, or they conflate.

CONSTRAINTS: solo hobby project; FREE data + tools only (no money); GTX 860M (4 GiB) — CPU/homology tools, not GPU embedding; Windows host + the Docker/D:-cache patterns already in use. REUSE the project's existing honesty/evidence-tier discipline (the SUSPEND/PLUMBING/abstain rails, lineage disclosure) + the determinant phenotype cells (frozen AMR surface in dna_decode/eval/amr_rules.py + the non-frozen TB cell in dna_decode/organism_rules/). HONESTY (load-bearing): a PHENOTYPE claim attaches ONLY where a validated determinant cell exists; everything else is molecular-function / pathway / homology-hypothesis / unknown, each explicitly tier-labelled; the UNKNOWN rate is REPORTED, never hidden. Do NOT propose a learned/embedding phenotype model for this — it is the closed-negative arm (honest qualifier: closed for the current free-public-embedding strategy, not logically forever; reopening needs a NEW named label-clean substrate + a pre-registered classical-baseline/de-confounding gate).

DELIVERABLE the downstream chain (NOT this idea-anchor) should produce: (1) the evidence-tier schema; (2) the UX/eval-gate spec; (3) a prototype "rough genome map" on 2-3 genomes; (4) a verdict on whether the output satisfies the north-star UX (and therefore whether to invest in catalog integration).
```

---

## Soraya's drafted framing (ratify or redirect)

**Formal rephrase.** Stand up the project's north star in its honest, achievable form: a per-genome,
evidence-tiered molecular-function annotation map (structural annotation + homology/domain function
transfer + curated-determinant phenotype overlays), scoped as a tier-schema + UX-gate + 2-3-genome
prototype spike — explicitly NOT a learned phenotype predictor (that arm is a closed negative).

**Fundamental clarifications (the ~3 the skill will likely ask — drafted answers):**
1. **Output UNIT — genome browser vs gene table vs operon view vs JSON-first report?** → *Draft:
   JSON-first per-feature report + a flat gene/feature table view.* Cheapest, composable, diff-able,
   matches the project's existing `.json` + `.md`-sidecar artifact style; a visual browser (pygenometracks
   — already deferred in `viz/`) is a later nicety, not the spike. *Technical — my call to draft.*
2. **Which 2-3 GENOMES?** → *Draft:* (i) an E. coli strain from the existing cohort (rich AMR-cell +
   well-annotated → tests the phenotype-overlay + curated tiers); (ii) M. tuberculosis H37Rv or a C. auris
   reference (exercises a second organism + the TB/fungal cells); (iii) a homology-heavy / hypothetical-
   protein-rich genome (a less-studied bacterium) to STRESS the middle-layer + unknown tiers — the whole
   point is whether the map is honest when most genes are Tier-homology/unknown. *Technical — my call.*
3. **The exact EVIDENCE TIERS?** → *Draft (define before wiring any source):*
   `phenotype-calibrated-determinant` (a validated AMR/TB cell fires) ▸ `curated-molecular-function`
   (a curated-DB hit: a named gene/enzyme) ▸ `pathway/module` (KEGG/MetaCyc module inference) ▸
   `homology-only-hypothesis` (Pfam/eggNOG/orthology domain hit, no curated identity) ▸ `unknown`
   (no confident annotation). The **unknown rate is a reported headline metric**, not hidden.
   **This carries a hidden authority decision** — the tier boundary between "molecular-function" and
   "phenotype" IS the project's honesty line; ratifying the tiers ratifies that line. *Yours.*

**Current assumptions (flagged for test by the downstream /probe):**
- The homology/annotation tool stack is FREE + CPU-runnable on this host (Bakta/Prokka, Pfam/hmmer,
  eggNOG-mapper, a KEGG/orthology source). *Likely true; the /probe should confirm install + runtime on
  the GTX-860M/Windows+Docker setup before committing.*
- The existing determinant cells can be invoked per-genome to populate the phenotype tier cheaply.
- A JSON-first map on 2-3 genomes is enough to JUDGE the UX gate (no full pipeline needed for the spike).

**Blunt opinion.** This is the right move: it is the honest, label-free, GTX-860M-feasible form of the
actual north star, it REUSES the project's hard-won honesty discipline, and it stops chasing the dead
learned-embedding arm. The ONE thing that makes it fail is skipping the UX/eval gate — without it this
becomes "wire eggNOG, wire KEGG, ship" busywork that never answers "does a rough map of what each part
does actually help?" Define the gate + the tiers FIRST, prototype on 2-3 genomes, and only then decide on
catalog integration. Keep phenotype claims behind the validated-determinant wall.

**Recommended next step.** `/probe genome-map` — this is code-touching + tool-dependent (a new annotation
+ homology stack, a tier schema, a map artifact) and needs the FREE-tool feasibility + the
phenotype-tier-vs-molecular-function honesty boundary grounded before a plan. First probe actions: confirm
the free annotation/homology tools install + run on this host; inventory how the existing determinant
cells expose a per-genome call; pin the exact tier definitions + the unknown-rate metric; pick the 3 test
genomes.

```
/probe genome-map
```
