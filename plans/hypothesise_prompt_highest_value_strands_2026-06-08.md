# /hypothesise prompt — highest-value test strands for the DNA decoder (2026-06-08)

> Paste the block below after `/hypothesise`. It is self-contained (the skill is stateless).
> Purpose: surface RANKED speculative-but-not-impossible strands worth TESTING next, by naming the
> assumptions currently believed binding and asking which can be productively relaxed.

---

**ROLE / GOAL.** You are doing constraint-relaxation discovery for `dna_decode` — an AI DNA→trait decoder
TOOL (north star: a useful tool, not papers; failure-tolerant iteration). Find the highest-value
*strands* (testable directions) to move forward with, where "value" = how much new decoding reach/utility
a strand unlocks per unit of cheap-to-run test. Output ranked grey-zone hypotheses, each with the
cheapest falsifier runnable on OUR substrate.

**GROUNDED STATE — what is already PROVEN, FALSIFIED, and SHIPPED (do not re-derive these):**

- SHIPPED product = **deterministic mechanism-feature decoders** on one shared curated-DB blastn engine:
  8 decoders (amr [bacterial AMRFinder + fungal BLAST-ERG11], pathotype, plasmid, serotype, resfinder,
  pointfinder, disinfinder, mlst) + 3 analyses (concordance, profile/run-all, coloc). Each is
  offline-safe, DB-on-demand, uniform record schema.
- FALSIFIED (0-for-3): the **learned-embedding thesis** (NT-frozen whole-genome mean-pooling). On cipro it
  beat k-mer (0.914) but LOST to the QRDR-POINT knowledge baseline (0.943) AND its within-lineage
  signal = chance → it learned LINEAGE, not MECHANISM. Captured as the **embedding-niche three-part
  test**: a learned decoder only earns its place with (1) a sampling-INDEPENDENT label, (2) NO curated
  catalog already solving it, (3) organism-specific DEPTH (≥~100 same-organism strains). AMR fails (2),
  pathotype fails (1), carbon-util fails (3). "All roads return to AMR" = AMR/BV-BRC is the only deep
  sampling-independent E. coli lab phenotype found so far.
- PROVEN TRANSFER + its BOUNDARIES (deterministic rules):
  - Rules transfer cleanly across **Enterobacterales** (E. coli / Klebsiella / Enterobacter / Salmonella):
    cipro QRDR-POINT, ceftriaxone ESBL-subclass, tet, gentamicin — acc 0.83–1.0.
  - **CONTENT boundary** (Acinetobacter × meropenem): broad CARBAPENEM-class count over-calls because
    intrinsic blaOXA-51-family is in every isolate (spec→0). A literature carbapenemase-STRENGTH refinement
    recovers acc 0.833/spec 1.000 but has a hard ~33% FN ceiling — IS-element/promoter-driven overexpression
    of intrinsic/weak genes is invisible to gene-PRESENCE.
  - **TUNING boundary** (Campylobacter × cipro, different PHYLUM): the QRDR mechanism caller transfers
    PERFECTLY (gyrA T86I = E. coli S83L analog, 15/15 R carry exactly one) but the THRESHOLD doesn't —
    deployed threshold=2 (E. coli double-mutant tuning) calls all R strains S; threshold=1 → 1.000/1.000/1.000.
  - KINGDOM jump (C. auris azoles, EP7): hand-curated ERG11 target-site caller TRANSFERS across the
    kingdom boundary (sens 1.0 across 2 clades); spec is LABEL-limited (reduced-suscept F126L carriers
    below the CDC tentative breakpoint — "suspect the label", not a caller defect).
- The recurring FAILURE SIGNATURE across every boundary: **what gene-presence cannot see** —
  IS-element/promoter-driven overexpression, efflux up-regulation, porin loss, and binary breakpoints
  splitting a reduced-susceptibility continuum. The failures are concentrated, not scattered.

**ASSUMPTIONS CURRENTLY BELIEVED BINDING — your job is to find which to relax (relax ONE per hypothesis,
name it explicitly):**

1. "The decoder predicts binary R/S from PRESENCE of curated determinants." (→ could it read expression
   context — IS-element/promoter insertions upstream of a gene — from the SAME assembly it already has?)
2. "Each organism×drug needs hand-set thresholds/curation." (→ could a small labeled cohort AUTO-CALIBRATE
   the per-organism threshold / content adjustment — a meta-rule learned cheaply, not a deep model?)
3. "Learned embeddings need the full three-part niche to add value." (→ could embeddings earn a NARROW
   role INSIDE the deterministic frame — e.g. only on the FN residue the curated rule provably can't see?)
4. "All roads return to AMR (only deep sampling-independent E. coli lab phenotype)." (→ what OTHER
   sampling-independent, lab-MEASURED phenotype — any organism — has the depth, that we've not checked?)
5. "The kingdom jump = fungi via hand-curated ERG11/FKS1." (→ which OTHER non-bacterial targets have a
   single-locus, curated, target-site resistance mechanism that the existing engine would transfer to —
   antiviral, antiparasitic, herbicide/crop, antimycobacterial?)
6. "Labels come from NCBI/BV-BRC." (→ what independent label source — NARMS, EuSCAPE, CARD prevalence,
   phenotype microarray, challenge assays — would de-confound or deepen a strand?)

**HARD CONSTRAINTS (must hold — do NOT propose strands that violate these):**

- **No money.** No paid APIs, no paid cloud, no purchases. Free public data + local compute only.
- **Compute:** a GTX 860M laptop (4 GB, no bitsandbytes) for everything; ONE Precision-7780 ~12 GB GPU
  job at a time (already committed to the Arabidopsis G2 embedding test). Docker present for bioinformatics
  tools (AMRFinder/Mash/Bakta/BLAST). Prefer strands testable on the LAPTOP with data we already have.
- **Determinism-first product.** A learned component must beat BOTH the naive tool AND the domain-knowledge
  baseline on INDEPENDENT data to earn its place; otherwise the deterministic rule ships.
- **Honesty:** every strand must carry a falsifier that could actually kill it.

**DISCOVERY TARGET (the question to hypothesise on):**

> Given the proven shape above — a working deterministic decoder whose failures are CONCENTRATED in
> "what gene-presence can't see" — what are the highest-value speculative-but-not-impossible strands to
> test next, that would most expand the decoder's reach or close its biggest blind spot, and that are
> cheaply falsifiable on our substrate (laptop + data on hand)? Relax exactly one binding assumption per
> hypothesis.

**OUTPUT CONTRACT — rank the hypotheses by VOI (value-of-information / cost). For EACH:**

- **Claim** (one falsifiable sentence — the speculative-but-not-impossible strand).
- **Relaxed assumption** (which of the 6 above; one only).
- **Why plausible, not impossible** (the grey-zone reasoning — grounded in the state above).
- **Cheapest falsifier** (the smallest concrete test on OUR substrate — name the data already on hand if
  possible: the cached AMRFinder runs for Campylobacter/Acinetobacter/the 4 Enterobacterales, the N=147
  cipro cohort, the C. auris ERG11 cohort, the committed assemblies). State the kill condition.
- **What it unlocks** if it survives (the value).
- **VOI tier** (HIGH / MED / LOW) + a one-line why.

Aim for 6–10 hypotheses spanning DIFFERENT relaxed assumptions (don't cluster them all on one). Flag any
hypothesis whose falsifier needs the workhorse GPU or any non-free resource — those are lower priority.
