# Broaden-Horizon Dataset Hunt — decompose & plan (2026-07-02)

**Goal (verbatim):** "broaden our horizon and search for any available dataset with good enough quality on
any living creature."

**Framing (load-bearing):** this is NOT a blind search. The project has spent months discovering exactly what
"good enough quality" means, encoded in `wiki/negative_results_map_2026-06-13.md` (8 rejection gates) + the
embedding-niche 3-part test + the horse-coat joinability lesson. The hunt = **apply that hard-won quality bar
to the tree of life and find substrates that PASS.** The binding constraint has always been LABELS, not
models — so we search for LABELS first, genotype second.

---

## The quality bar (flow-down root — every candidate is scored against this)

A "good enough" genotype→phenotype dataset must clear ALL of:

| # | Gate | Grounded in |
|---|---|---|
| **G1 Accessible** | free + fetchable, no DUA/paywall wall | the LABELS wall (freeze forward-path) |
| **G2 Non-circular label** | phenotype measured independently, NOT derived from a tool the decoder competes with | circular-label gate; HIV win (Stanford PhenoSense ≠ Sierra) |
| **G3 Sampling-independent label** | phenotype not confounded with sampling context (site/source/study/date) | study==class + sampling-defined + surveillance-domination gates |
| **G4 Unit-level joinable** | per-individual genotype ⋈ per-individual phenotype (NOT aggregate/separate) | horse-coat joinability lesson |
| **G5 Provenance-separable** | a leakage-free split exists (temporal / accession / cohort) | provenance-not-separable gate; prospective-lock |
| **G6 Depth** | ≥~100 same-organism units for a LEARNED decoder; ANY N if a curated determinant catalog exists (deterministic path) | embedding-niche depth + dedup-collapses-balance |
| **G7 Genotype fetchable** | actual sequence/variants per unit, not just phenotype | assembly-attrition gate |
| **G8 Label not censored to uselessness** | quantitative labels tierable (not all operator-censored at one bound) | MIC-censoring gate; operator-aware censoring |

**Decoder-type branch (decides which product paradigm applies):**
- **Curated catalog EXISTS** (determinant→phenotype map) → the DETERMINISTIC-rule paradigm — the project's
  VALIDATED product (AMR/TB/HIV/fungal). Depth requirement relaxes; G6 satisfied by the catalog.
- **No catalog + depth≥100 + G2 + G3** → the LEARNED-decoder niche. HIGH bar: embeddings are **0-for-4**
  de-confounded (cipro/pathotype/Arabidopsis/…); a candidate here must beat a domain-knowledge baseline AND
  hold within-lineage, or it's another structure-learner. Enter here only with eyes open.

---

## Decomposition (4 families; NOT yet ledger-spawned — parked for ratification)

Single exploratory effort, 4 phases. Ledgers deliberately NOT created (premature for a hunt; spawn on GO).

### F1 — Quality-bar scorecard (flow-down root; class a/e)
Turn the 8 gates + branch into a machine-checkable candidate scorecard
(`scripts/dataset_candidate_scorecard.py` or a structured `wiki/dataset_hunt_scorecard.md` rubric). Each
candidate → per-gate PASS/FAIL/UNKNOWN + decoder-type verdict + a one-line rejection reason (feeds the
negative-results map). **Mandatory first — nothing ranks without it.**

### F2 — Multi-modal candidate sweep (parallel search axes; class e — /research by hand)
Three blind-to-each-other axes (the multi-modal-sweep pattern):
- **2a by-repository:** ENA/NCBI, Ensembl{Plants,Metazoa,Fungi,Bacteria}, EBI BioStudies, Dryad, Zenodo,
  Figshare, model-organism DBs (SGD, FlyBase, WormBase, MGI, Gramene), GWAS Catalog, AnimalQTLdb.
- **2b by-phenotype-class:** lab-measured quantitative traits (growth/fitness across conditions, drug
  response, metabolite levels, morphometrics), Mendelian visible traits with functional-genotype joins,
  image-paired phenotype (the multimodal north-star).
- **2c by-curated-catalog:** organisms with an existing determinant catalog (ClinVar, dog/cat coat-color
  IDID, more pathogen AMR/antivirals/antiparasitics) → deterministic-path-extensible.

### F3 — Screen + rank (critical path; class e)
Score every F2 candidate with F1's scorecard → ranked shortlist with explicit per-gate verdicts + rejection
reasons. Output: `wiki/dataset_hunt_shortlist_2026-07-02.md`.

### F4 — Pilot-fetch #1 (real-surface-first go/no-go; class b/e)
Fetch a SMALL slice of the top-ranked candidate, confirm the genotype↔phenotype join actually works on real
data, and that a simple determinant/rule beats the K/N null. This is the GO gate for a full build (mirrors
how the TB/HIV cells were de-risked). Real-surface-first: the first test is the real download, not a mock.

**Flow-down / critical path:** F1 → F3 → F4 (F2 runs parallel to F1; both feed F3).

---

## Seed candidates (R2 widen-once — from knowledge, UNVERIFIED, first-pass gate read only)

Concrete starting points so F2/F3 aren't abstract. **Every cell below is a hypothesis to verify in F2/F3, not
a claim.** Bias toward DEEP + measured + free.

| candidate | creature | phenotype | first-pass read | why promising / risk |
|---|---|---|---|---|
| **1002 Yeast Genomes** (Peter 2018) | *S. cerevisiae* | growth/fitness across many conditions | likely G1/G4/G6/G7 pass | 1011 isolates × quantitative lab phenotypes; risk G3 (some traits sampling-linked) |
| **DGRP** (Drosophila Genetic Reference Panel) | fruit fly | ~many quantitative traits | likely G1/G4/G6/G7 pass | 205 inbred lines × deep measured phenome; risk G6 borderline (205 lines) |
| **CeNDR / CaeNDR** | *C. elegans* | quantitative wild-isolate traits | likely G1/G4/G7 pass | wild-isolate NIL panels + measured traits; risk G3/G6 |
| **Mouse Phenome DB / Collaborative Cross** | mouse | huge measured phenome | G1 pass; G7 check | very deep; risk G4 (genotype↔pheno join granularity) |
| **Rice 3000 Genomes + phenotypes** | rice | agronomic quantitative traits | G1/G6/G7 likely pass | 3000+ accessions; risk G3 (field/environment confound) |
| **ClinVar** | human | variant→disease | curated catalog → deterministic | G2 risk (clinical assertions circular-ish); strong for rule-path |
| **NCBI Pathogen Detection (non-E.coli species)** | many bacteria | AMR AST | deterministic-catalog-extensible | extends the validated AMR product to new species |
| **AraPheno (beyond flowering time)** | Arabidopsis | many quantitative traits | flowering-time is a CLOSED embedding negative | deterministic/Mendelian sub-traits only; embeddings foreclosed here |
| **Dog/cat coat-color (IDID)** | dog/cat | Mendelian coat traits | curated catalog → deterministic | like horse-coat; risk G4 (find a JOINED functional-genotype×phenotype table) |

Standouts to verify first: **1002 Yeast Genomes, DGRP, CeNDR** — deep, measured, free, plausibly
sampling-independent. These are the strongest "new kingdom, learned-decoder-eligible" bets IF they clear G2/G3.

---

## Success criterion (for the hunt as a whole)
≥1 candidate PASSES all 8 gates in F3 **and** clears the F4 pilot (real join works + simple rule beats null).
Honest terminal if none pass: a refreshed negative-results map (the gates held across the tree of life) — a
real, publishable-internally finding, not a failure.

## Scope / gates
Research-only (class e) through F3; F4 is a small real fetch (network, `auto`). No money, no destructive ops.
Docker/native tools reused as available. Frozen AMR surface untouched throughout.

## Recommended next step (planning-STOP — awaiting ratification)
1. Ratify or adjust the quality bar (the 8 gates) + the F1→F4 shape.
2. On GO: I build F1 (the scorecard), run F2 (a real `/research`-by-hand sweep — this is where the actual
   web search happens), then F3 rank, then F4 pilot the #1 — as an `--until-mvp` on this plan.
3. Optionally spawn ONE `/project-init` ledger ("cross-kingdom dataset hunt") to track it, rather than 4.
