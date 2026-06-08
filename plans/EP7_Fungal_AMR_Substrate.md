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
2b. ✅ **G0-COMPLETION (real data)** — fetched the real C. auris ERG11 reference (RefSeq `XM_029033208.2`,
   525 aa; numbering confirmed: **Y@132 / K@143 / V@125** match the catalog) + validated `call_erg11`
   end-to-end on real GenBank isolate alleles (study NCCPF470): `PV630306` WT→S (efflux/aneuploidy blind
   spots surfaced), `PV630305` **Y132F→R**, `PV630302` **K143R→R** = 3/3. Numbering-mismatch risk RETIRED.
   Reference + 3 public allele fixtures committed at `data/fungal_ref/`; 3 real-data regression tests
   (5 caller tests green total). The caller is validated against documented mutations, not just a planted one.
3. ✅ **Cohort-validation INFRASTRUCTURE (steps 3+4 machinery)** — `scripts/build_fungal_cohort.py`:
   reads a per-isolate label TSV (`assembly_accession`→download OR `genome_fasta`→local), runs `call_erg11`
   over the cohort, computes overall + **within-clade de-confounded** acc/sens/spec + the **efflux/aneuploidy
   discordance** set (MIC-R but ERG11-S = documented non-target mechanism), emits the G1 verdict
   (PASS / DOCUMENTED_FAILURE_MODE / FAIL) as `.md`+`.json`. CDC tentative breakpoint (fluconazole MIC≥32=R)
   + `mic_to_phenotype` added to `fungal_amr.py`. 4 cohort tests validate the FULL pipeline on the 3 committed
   real alleles via real blastn (16 fungal tests green). **G1 runs in ONE command once the label table lands.**
   - **G1 label-source finding (iron-law):** fluconazole MIC is NOT in NCBI BioSample metadata (53 hits are
     experimental-*treatment* fields, not AST) and BV-BRC has ZERO C. auris AMR rows. The ONLY label source
     is paper supplementaries — exactly as this plan stated. **The one remaining G1 input = a
     `[isolate_id, assembly_accession, fluconazole_mic, clade]` TSV** from the S.Africa(188)/India(350)
     supplementaries (or any published C. auris WGS+MIC table).
4. **Real G1 verdict (BLOCKED on the label TSV above)** — populate the TSV → run `build_fungal_cohort.py`
   → acc/sens/spec + efflux/aneuploidy discordance. Bar: acc≥0.80/sens≥0.80 OR documented failure mode.

## Falsifier
If ERG11-only sensitivity is poor because azole-R is dominantly efflux/aneuploidy-mediated (the multi-locus
finding: clade II had only 3/~38 R with ERG11), document it as the fungal analogue of the tet/Klebsiella
efflux blind spot — genotype-detects-the-target-mechanism, label-validation-limited.

## Honest scope
First eukaryotic substrate; first fungal determinant catalog. The method (target-site mutation scan) is the
proven one; the new infra is the BLAST caller + the hand-curated catalog. No compute/money.
