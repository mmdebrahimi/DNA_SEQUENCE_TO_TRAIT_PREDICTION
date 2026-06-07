# EP-7 — Fungal AMR (Candida auris azole) — the eukaryotic kingdom jump

> **Status:** STARTED 2026-06-07 (Path A of the ratified eukaryotic cycle). Determinant catalog shipped;
> BLAST caller + cohort + validation are steps 2-4.
> **Anchors on:** `research_outputs/eukaryotic-multimodal-substrate-feasibility-2026-06-07.md` (the survey
> that ranked this #1) + `plans/Trait_Decoding_Roadmap.md` Phase 6.
> **Why this is the chosen eukaryotic entry:** it is the kingdom jump (fungus) that reuses the project's
> PROVEN deterministic-determinant method — same phenotype class (AMR/MIC), documented ERG11/FKS1/TAC1b
> determinants, **NO foundation model / GPU / money** (a BLAST→known-mutation scan, the eukaryotic
> analogue of `amr_rules.py`). It extends the validated decoder across a kingdom boundary at near-zero
> compute cost.

## Terminal claim
`dna-decode amr --drug fluconazole --organism Candida_auris --genome-fasta X.fna` (or a fungal-specific
entry) emits an R/S azole call + the ERG11/FKS1 resistance substitutions driving it + blind-spots, AND on
a de-confoundable C. auris WGS+MIC cohort the deterministic rule clears acc ≥ 0.80 / sens ≥ 0.80 — OR a
documented organism-specific failure mode (efflux/aneuploidy-mediated R is the expected blind spot).

## Substrate (feasibility VERIFIED 2026-06-07)
- **Data:** deep — S. Africa outbreak 188 WGS with near-perfect MIC↔ERG11 (181/188 MIC>32 had ERG11 mut),
  India 350, global 12,644. Source = paper supplementaries + NCBI BioProject assemblies (NOT NCBI Pathogen
  Detection — Candida is not tracked there; cohort assembly is manual supplementary extraction).
- **Labels:** fluconazole MIC (sampling-independent lab measurement). Confound = CLADE structure (analogous
  to bacterial lineage) — must de-confound (within-clade R/S contrast).
- **Determinants:** documented (no AMRFinder-for-fungi). Catalog shipped: `dna_decode/data/fungal_amr.py`.
- **Compute:** NONE (BLAST+ installed at `C:/Users/Farshad/ncbi-blast/bin`; determinant scan only).

## Build steps
1. ✅ **Determinant catalog** — `dna_decode/data/fungal_amr.py` (ERG11 Y132F/K143R/F126T/VF125AL/clade-IV +
   C. albicans hotspots; FKS1 S639; TAC1b/efflux/aneuploidy as `FUNGAL_UNDETECTABLE_MECHANISMS`). 7 tests.
2. ✅ **BLAST caller (machinery, G0)** — `scripts/fungal_erg11_caller.py`: makeblastdb + blastn(CDS-vs-genome)
   + gap-aware codon-translate + catalog match → `call_from_observed_substitutions`. tblastn absent → blastn
   used (C. auris ERG11 intronless = colinear HSP). Offline-safe (INDETERMINATE when BLAST absent).
   **Validated against a planted Y132F via REAL makeblastdb+blastn** (2 tests). 
   **G0-COMPLETION (step 2b, next):** swap the synthetic reference for the real C. auris ERG11 CDS allele
   (NCBI) + validate on a real resistant genome (Y132F/VF125AL) to confirm catalog numbering ↔ reference.
3. **Cohort** — extract a C. auris WGS+MIC table from the S. Africa (188) + India (350) supplementaries;
   fetch assemblies by accession; build a de-confoundable cohort (within-clade fluconazole R/S contrast,
   reuse `cohort_deconfound.py`).
4. **Validation** — run the caller over the cohort; report acc/sens/spec + the efflux/aneuploidy discordance
   (expect FN from non-ERG11 R, the documented multi-locus mechanism). Capstone-style artifact.

## Falsifier
If ERG11-only sensitivity is poor because azole-R is dominantly efflux/aneuploidy-mediated (the multi-locus
finding: clade II had only 3/~38 R with ERG11), document it as the fungal analogue of the tet/Klebsiella
efflux blind spot — genotype-detects-the-target-mechanism, label-validation-limited.

## Honest scope
First eukaryotic substrate; first fungal determinant catalog. The method (target-site mutation scan) is the
proven one; the new infra is the BLAST caller + the hand-curated catalog. No compute/money.
