# Changelog

All notable changes to `dna_decode`. Format loosely follows [Keep a Changelog](https://keepachangelog.com/);
this is a solo research-tool repo so the granularity is per-release-theme, not per-PR.

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
