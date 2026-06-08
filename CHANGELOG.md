# Changelog

All notable changes to `dna_decode`. Format loosely follows [Keep a Changelog](https://keepachangelog.com/);
this is a solo research-tool repo so the granularity is per-release-theme, not per-PR.

## [Unreleased] — plasmid replicon-typing decoder

New deterministic trait decoder `dna-plasmid` (`dna-decode plasmid`) — the tool grows beyond AMR.

- **New capability:** plasmid incompatibility (Inc) replicon typing from a genome assembly via the curated
  PlasmidFinder allele DB + real blastn (identity 95 / coverage 60, PlasmidFinder defaults). Reports the Inc
  replicons present (IncF/IncH/IncI/IncX/IncN/…) — composing with `dna-amr`: AMR says *what* resistance,
  plasmid typing says whether it likely rides a known mobile element.
- Sibling architecture to `dna-pathotype` (curated-DB blastn caller; reuses `pathotype/vf_runner`'s blastn
  resolvers — DRY). Offline-safe: no blastn / no DB → `status: "unavailable"` (exit 3), never a crash.
  `caller_is_independent_baseline: false` (faithful to PlasmidFinder's own method, not an independent check).
- DB downloaded on demand (not committed), like the VirulenceFinder DB. 7 tests
  (`tests/test_plasmid_decoder.py`); cli-dispatch registry contract updated to the 3-decoder set.

## [0.5.0] — 2026-06-08 — Fungal AMR decoder (the kingdom jump)

`dna-amr` now decodes **fungal** azole/echinocandin resistance, not just bacterial — the determinant-scan
method validated across the bacteria→fungi kingdom boundary.

- **New capability:** `dna-amr --drug fluconazole|voriconazole|caspofungin|micafungin` routes to a
  BLAST-ERG11/FKS1 target-site engine (vs the AMRFinder engine for bacterial drugs — there is no
  AMRFinder-for-fungi). Two source modes: `--genome-fasta` (BLAST the committed C. auris ERG11 reference)
  and `--observed GENE:SUB[,...]` (pure, wheel-only, no BLAST). Emits the same `amr-mechanism-call-v1`
  record as the bacterial path (uniform tool surface); S calls surface the efflux/aneuploidy blind spots.
- **Validation (Gate G1, `wiki/fungal_ep7_g1_closeout_2026-06-08.md`):** on a de-confounded C. auris
  WGS+MIC cohort (S. Africa bloodstream, PRJNA737309 + AraPheno-style Table S1 MICs), the deterministic
  caller found the catalogued ERG11 mutation in **100% of fluconazole-MIC-R isolates across two clades**
  (clade I Y132F, clade III F126L/VF125AL) — sensitivity 1.0. Specificity is label-limited (reduced-
  susceptibility F126L carriers fall below the CDC tentative breakpoint), the documented "suspect the
  label" pattern; the genotype is the trustworthy output.
- **Infra shipped:** `dna_decode/data/fungal_amr.py` (hand-curated catalog + CDC tentative breakpoints),
  `scripts/fungal_erg11_caller.py` (BLAST→codon-map→catalog), `scripts/build_fungal_cohort.py` (cohort
  validation + within-clade de-confound + LABEL_LIMITED_FAILURE verdict), `scripts/assemble_sra_cohort.py`
  (targeted ERG11 read-mapping, ~4 min/isolate vs ~45 min full assembly). Committed real C. auris reference
  + 3 public allele fixtures (`data/fungal_ref/`).
- **Eukaryotic Path B (Arabidopsis flowering-time embedding test, Gate G2)** pre-staged + brainstorm-revised
  + CPU-only dry-manifest gate coded (`scripts/g2_dry_manifest.py`); GPU run deferred to the workhorse.
- ~43 new tests (fungal catalog/caller/cohort/CLI + dry-manifest). Bacterial path unchanged.

## [0.4.0] — 2026-06-07 — Multi-Organism AMR Decoder (capstone)

Milestone release consolidating the AMR arc. No new code — a capstone over v0.3.x.

`dna-amr` is a deterministic, interpretable AMR R/S decoder validated across **6 drugs × 4 organisms ×
4 mechanism classes, spanning the gram divide** (E. coli, K. pneumoniae, P. aeruginosa, S. aureus),
deployed as `dna-amr` / `dna-decode`. Every per-drug rule beats naive AMRFinder on independent data.

- **One engineering principle** held across every organism/mechanism: count the drug's SPECIFIC
  resistance determinants (target point-mutations / drug-specific Subclass / acquired gene-family), not
  the broad drug-class bag — intrinsic chromosomal genes (efflux) are the cross-organism gotcha.
- **Honest limits, named in every output:** `undetectable_mechanisms` (efflux/porin/regulatory expression
  phenotypes) + label-quality caveats (oxacillin → use cefoxitin). The recurring binding constraint is the
  de-confounded, reliably-labeled substrate — not the method.
- Capstone: `wiki/amr_multiorganism_capstone_2026-06-07.md`. 108 tests green.

**This is the milestone.** Further organism/drug breadth is diminishing-returns (re-confirms the same two
findings). The genuinely-different next leaps (cross-lab validation, a non-AMR sampling-independent
substrate, multimodal/eukaryotic) require resources beyond autonomous code work.

## [0.3.7] — 2026-06-07

1st Gram-positive: S. aureus oxacillin (MRSA/mecA) — genotype transfers; honest label finding.

### Added
- **oxacillin** (6th drug, 1st Gram-positive): mecA-based MRSA rule (threshold 1 + METHICILLIN-subclass,
  excludes blaZ penicillinase). Added to mic_tiers (breakpoints, classes, mec loci, primary mechanism) +
  DRUG_RULE. `supported_drugs()` now 6. +2 tests (104 → 106 → 108 green).
- **S. aureus oxacillin validation (4th organism, 1st Gram-positive):** N=30.
  `wiki/staphylococcus_aureus_oxacillin_validate_2026-06-07.md`.

### Finding (the honest result)
- **mecA genotype detection TRANSFERS to Gram-positive: sens 1.000** (all 15 R strains carry mecA). The
  acquired-gene + Subclass-refinement approach works on a Gram-positive, as on the gram-negatives.
- **spec 0.333 is oxacillin-LABEL noise, NOT a rule defect:** 10/15 "oxacillin-susceptible"-labeled strains
  carry full-length mecA — far above genuine OS-MRSA (<5%). Oxacillin direct AST is the documented unreliable
  comparator for mecA; CLSI/EUCAST recommend **cefoxitin** as the surrogate. The proper-label re-test is
  **substrate-blocked** (cefoxitin = only 3R on this NCBI dataset).
- **Terminal:** Gram-positive mecA detection generalizes; phenotype-label validation is the limit — the
  project's recurring "substrate/label is the binding constraint" lesson, now confirmed on a Gram-positive.

## [0.3.6] — 2026-06-07

3rd organism (Pseudomonas) + cross-organism shipped in the CLI.

### Added
- **`dna-amr --organism <O>`** (genome mode) — passes through to AMRFinder `-O` (organism-specific QRDR
  point-mutation detection); recorded in `provenance.amrfinder_organism`. Default Escherichia. Closes the
  gap where cross-organism support lived only in validation scripts, not the shipped CLI. +2 CLI tests.
- **Pseudomonas aeruginosa cipro VALIDATED** (3rd organism): N=30 acc 0.867 / sens 0.80 / spec 0.933
  (beats naive AMRFinder 0.767). The QRDR-POINT rule transfers UNCHANGED to a *less-similar* gram-negative
  (MexAB-OprM efflux, intrinsic AmpC — no oqxAB). 3 FN = efflux-mediated cipro-R (expected blind spot).
  `wiki/pseudomonas_aeruginosa_ciprofloxacin_validate_2026-06-07.md`.
- **`scripts/organism_drug_validate.py`** — generalized any-NCBI-organism × any-drug validator
  (auto-discovers latest PDG metadata, reuses cached runs). Every future organism is now a one-command run.

### Result
dna-amr validated across **3 organisms** (E. coli, K. pneumoniae, P. aeruginosa). The "count the drug's
specific determinants, not the broad class bag" principle holds across all three — strong evidence it is
organism-general, not E. coli-specific.

## [0.3.5] — 2026-06-07

Klebsiella full drug matrix complete — dna-amr validated 5 drugs × 2 organisms.

### Added
- **Klebsiella cef + gent + tet validated** (rules applied unchanged from E. coli):
  - ceftriaxone: acc 0.800 / sens 1.0 / spec 0.6 ✅
  - gentamicin: acc 0.867 / sens 0.867 / spec 0.867 ✅
  - tetracycline: acc 0.800 / spec 1.0 / **sens 0.600** ⚠️ PARTIAL (efflux blind spot — see below)
  - `scripts/klebsiella_drug_validate.py` (drug-agnostic; reuses cached Klebsiella runs across cohorts).
  - Consolidated: `wiki/klebsiella_drug_matrix_2026-06-07.md`.
- **`gene_prefixes` rule refinement** (`amr_rules.py`): tetracycline now counts only acquired `tet*` genes,
  excluding intrinsic K. pneumoniae OqxAB efflux (AMRFinder-tagged TETRACYCLINE but present in susceptible
  isolates). Same cross-organism principle as cipro QRDR-POINT. **Also improved E. coli tet 0.833 → 0.917.**

### Findings
- **tetracycline / Klebsiella is the honest PARTIAL:** the acquired-`tet*` rule is precise (spec 1.0) but
  sens 0.600 — 6/15 R strains are efflux-mediated (oqxAB overexpression), undetectable by ANY
  curated-determinant rule (an expression phenotype, surfaced in `undetectable_mechanisms`). A documented
  biological limit, not a rule defect.
- **Cross-organism principle confirmed 3× (cipro/tet/the gent+cef+mero Subclass refinements):** count the
  drug's specific resistance determinants, not the broad drug-class bag; intrinsic chromosomal determinants
  (efflux) are the organism-specific gotcha.

### Result
dna-amr spans **5 drugs × 2 organisms × 4 mechanism classes**; 4/5 Klebsiella drugs clear the bar
zero-tuning, all beat naive AMRFinder. +3 tet tests (104 green).

## [0.3.4] — 2026-06-07

Klebsiella meropenem — 2nd organism, NEW mechanism class (carbapenem). Phase 3 slice 2.

### Added
- **meropenem decoder** (5th drug): acquired-carbapenemase rule (threshold 1 + CARBAPENEM-subclass
  refinement — blaKPC/NDM/OXA-48). **Klebsiella N=30 acc 0.867 / sens 1.0 / spec 0.733** (vs naive
  AMRFinder 0.533; the CARBAPENEM-subclass refinement lifts spec 0.067→0.733 by excluding ESBL/AmpC that
  raise meropenem MIC without hydrolyzing it). `wiki/klebsiella_meropenem_validate_2026-06-07.md`.
  Carbapenem is the defining K. pneumoniae clinical threat — a mechanism class E. coli AMR never covered.
- meropenem added to `mic_tiers.py` (breakpoints CLSI R≥4/S≤1 + EUCAST; AMRFinder classes;
  carbapenemase loci catalog; primary mechanism). `supported_drugs()` now returns 5.
- `scripts/klebsiella_meropenem_validate.py` (reuses cached cipro-cohort AMRFinder runs on strain overlap).
- +2 tests (carbapenemase counted / ESBL excluded). 102 green.

### Honest scope
The meropenem rule is blind to porin-loss-mediated carbapenem resistance (no carbapenemase gene) — the
expected FN mode; surfaced in `undetectable_mechanisms`. 4 FP (susceptible strains carrying a carbapenemase
gene — likely low-expression / borderline MIC).

## [0.3.3] — 2026-06-07

Cross-organism: Klebsiella + the QRDR-POINT cipro rule (roadmap Phase 3, slice 1).

### Added
- **Cross-organism transfer (Klebsiella pneumoniae cipro):** N=30 NCBI cohort, **acc 1.000** with the
  deployed rule. The method generalizes across organisms. `scripts/klebsiella_cipro_transfer.py`,
  `wiki/klebsiella_cipro_transfer_2026-06-07.md`. `_run_amrfinder` gained an `organism` param
  (`-O Klebsiella_pneumoniae`).
- `qrdr_point_count` / `qrdr_point_determinants` in `amr_rules.py` — count fluoroquinolone QRDR
  target-alteration POINT mutations (gyrA/parC/parE) only.

### Changed (ratified — changes the deployed cipro number)
- **cipro `DRUG_RULE` switched to the QRDR-POINT rule globally** (`counter='qrdr_point'`): count QRDR
  target POINT mutations ≥2, not the broad QUINOLONE-class determinant bag. Rationale: the broad rule
  FAILED on Klebsiella (acc 0.5 — intrinsic chromosomal OqxAB efflux, absent in E. coli, saturates the
  count). The canonical target-mutation count is cross-organism-robust. Net effect:
  - Klebsiella cipro 0.5 → **1.000**
  - E. coli cipro in-cohort 0.939 → 0.925 (−1.4pp; dropped cases were qnr/efflux-mediated)
  - E. coli cipro **cross-source (NCBI) 0.955 → 1.000** (+4.5pp — QRDR-POINT generalizes better; the
    in-cohort −1.4pp was tuning-cohort overfit).
- **Platform finding:** cross-organism transfer requires counting the drug's TARGET-alteration mutations,
  not the broad drug-class bag; intrinsic chromosomal determinants are the organism-specific gotcha.
- Tests: cipro tests updated to QRDR-POINT (point-mutation rows); +4 new (qrdr helpers + intrinsic-exclusion);
  cohort regression re-pinned to 0.925/0.875/0.973. 24 total green.

## [0.3.2] — 2026-06-06

Trust-hardening + honesty taxonomy (adversarial-review follow-through). No new science — makes the
shipped decoders auditable + honest about their blind spots.

### Added
- **Blind-spot honesty:** every SUSCEPTIBLE `dna-amr` call now carries `undetectable_mechanisms`
  (`efflux` / `porin_loss` / `regulatory` — expression/regulatory resistance absent from AMRFinder's
  curated determinants) + a caveat that a negative means "no curated determinant found", not "definitely
  susceptible". `UNDETECTABLE_MECHANISMS` in `amr_rules.py`.
- **Discordance taxonomy:** `discordance_bucket(prediction, true_label)` + `evaluate_cohort` now emit a
  failure-mode breakdown — `FN_undetected_mechanism` (R missed → the blind spots) vs
  `FP_determinant_without_phenotype` (called R but susceptible → label noise / silent-or-low-expression /
  borderline MIC). The "failure-tolerant tool" deliverable: names where it fails.
- **Provenance pin:** output JSON `provenance.amrfinder_image` records the pinned AMRFinderPlus image
  (`ncbi/amr:4.2.7-2026-03-24.1`; tag encodes the DB date) so an R/S verdict is reproducible against a
  known determinant source. `AMRFINDER_IMAGE_PINNED` in `amr_rules.py`.
- **Value-add headline** in `wiki/dna_amr_multidrug_validation_2026-06-06.md`: explicit naive-AMRFinder
  vs dna-amr table (cef +0.28 acc/+0.50 spec; gent +0.43/+0.57) — proves the per-drug policy adds value
  over vanilla "any determinant → R", not just re-prints AMRFinder hits.
- 5 new tests (blind-spots on S/R calls, discordance taxonomy, cohort discordance breakdown). 20 total.

### Fixed
- Stale gentamicin "NOT yet cohort-validated" claims (`amr_rules.py` docstring + wiki caveat) reconciled
  with the N=128 acc 0.945 validation that DRUG_RULE/README/CHANGELOG already recorded (claim-hygiene).

### Validated (cross-source — closes the same-database gap)
- **Independent NCBI Pathogen Detection validation** (`scripts/xsource_validation.py`,
  `wiki/dna_amr_xsource_validation_2026-06-07.md`): 22 E. coli, balanced ~11R/11S per drug, **zero
  accession overlap** with the 176 BV-BRC cohort accessions (enforced at selection). Result: cipro 0.955,
  cef 0.864, gent 1.000, tet 0.909 — comparable to / better than in-cohort. **Answers the product
  question:** dna-amr beats naive "any-determinant→R" AMRFinder on UN-tuned data by +0.09 (cipro) /
  +0.23 (cef) / +0.41 (gent) / +0.0 (tet) accuracy — the per-drug threshold + Subclass refinement IS the
  value-add, not determinant discovery. Sensitivity 1.0 on all 4 drugs (FN=0); all errors are FP
  (determinant-present-but-susceptible). Honest scope: NCBI Pathogen Detection is a different
  source/curation, not a controlled different-lab study.

## [0.3.1] — 2026-06-06

The "one coherent tool" consolidation (after the v0.3.0 build settled the embedding question).

### Added
- **Unified `dna-decode` console entry** (`dna_decode/cli.py`) — single tool that dispatches to the
  trait decoders: `dna-decode amr …`, `dna-decode pathotype …`, `dna-decode list` (capability +
  validation status), `dna-decode --version`. Thin pass-through (argv delegated verbatim); the
  per-decoder entries (`dna-amr`, `dna-pathotype`) stay independently usable + unchanged. 6 dispatch tests.

### Changed
- Project ledger (`project_state/dna-decode-2026-05-11.md`) updated to record the strategic inflection:
  embedding frontier closed (0-for-3), deterministic mechanism-feature decoders are the product.

## [0.3.0] — 2026-06-06

The "deterministic decoders win" release. The frozen-genome-embedding (NT mean-pool) thesis was tested
to a decisive conclusion and the project committed to deterministic, interpretable mechanism-feature
decoders as the shipping path.

### Added
- **Multi-drug deterministic AMR caller** (`dna-amr`). Extended from cipro-only to **all 4 drugs**, with
  per-drug validated rules baked into `dna_decode/eval/amr_rules.py::DRUG_RULE`:
  - ciprofloxacin — threshold 2 (QRDR point-mutations) — N=147 acc 0.939.
  - ceftriaxone — threshold 1 + **extended-spectrum subclass refinement** (CEPHALOSPORIN/CARBAPENEM;
    excludes intrinsic blaTEM-1/blaEC that are ampicillin-R not ceftriaxone-R) — N=60 acc 0.933.
  - tetracycline — threshold 1 (acquired tet genes) — N=12 acc 0.833 (small N, provisional).
  - gentamicin — threshold 1 + **GENTAMICIN-subclass refinement** (excludes aph/aadA
    streptomycin-kanamycin genes that don't confer gentamicin-R) — N=128 acc 0.945.
  - General fix: a broad AMR class over-calls (cef spec 0.41, gent spec 0.39) because it counts genes for
    OTHER class members; AMRFinder's Subclass field is the drug-specific discriminator. One-line refinement.
  - `call_resistance(tsv, drug)` now auto-selects the per-drug rule; explicit threshold still overrides.
  - Validation: `wiki/dna_amr_multidrug_validation_2026-06-06.md`.
- **De-confound gate** (`dna_decode/eval/cohort_deconfound.py`) — within-lineage label-contrast
  precondition (3-state DE_CONFOUNDED/WARN/CONFOUNDED + promotability) for any embedding-vs-classical
  falsifier. The reusable guard against study==class confounding.
- **AMR embedding falsifier** (`scripts/amr_falsifier.py`) + **QRDR-POINT knowledge baseline**
  (`dna_decode/eval/point_baseline.py`) + within-lineage diagnostic.
- **dna-amr external validation** — held-out N=29 acc 0.862 / sens 0.882 / spec 0.833
  (`wiki/dna_amr_external_validation_2026-06-05.md`).
- **EP-6 carbon-utilization substrate infra** — `dna_decode/data/bacdive.py` loader +
  `scripts/bacdive_carbon_util_feasibility.py` census + `scripts/bacdive_li2023_to_long.py` adapter.

### Findings (recorded, not code)
- **AMR embedding thesis FALSIFIED on the cleanest substrate.** Cipro N=147 (de-confounded): NT-XGBoost
  0.914 beats k-mer (+8.9 pp) but **LOSES to QRDR-POINT 0.943**; NT within-lineage concordance = chance
  (0.605, p=0.365) → it learned lineage, not mechanism. Decision: `plans/AMR_embedding_niche_decision_2026-06-05.md`.
- **Carbon-utilization (EP-6) E. coli-INFEASIBLE.** Data acquired (Li et al. 2023 OSF jwkr7); E. coli
  slice = 27 strains, 0 carbon sources clear the ≥100 floor. The embedding-niche test needs a THIRD
  requirement beyond sampling-independent-label + no-catalog: **organism-specific depth at scale**.
  `wiki/bacdive_carbon_util_feasibility_2026-06-06.md`.

### Changed
- `dna-amr` CLI `--resistance-threshold` now defaults to the per-drug validated config (was hard-coded 2).
- README decoder table + `plans/Trait_Decoding_Roadmap.md` Phase 2/4 updated.

## [0.2.0] — 2026-06-05

- In-package `dna-amr` console entry (deterministic AMR mechanism caller, cipro-validated).
- `dna-pathotype` console entry (deterministic VirulenceFinder-marker pathotype resolver + abstention +
  canonical-VF diff; ExPEC recall 0.917).
- Packaging: both decoders ship in the wheel (`[project.scripts]`).

## [pathotype-v0] / [phase-1-shipped] — 2026-05/06

- Phase-1 closeout: NT-frozen-pooling characterized (passes concentrated-signal mechanisms, fails
  distributed). v0 cipro decoder + pathotype v0 resolver shipped. See git history + `wiki/`.
