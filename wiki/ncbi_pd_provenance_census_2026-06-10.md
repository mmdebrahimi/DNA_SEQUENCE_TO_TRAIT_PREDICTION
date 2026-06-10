# NCBI-PD provenance census (Stage-1, no model) — the strand REOPENS — 2026-06-10

> The free move a `/brainstorm` pass surfaced (and the prior census/sweep overclaimed away): the NCBI
> Pathogen Detection metadata TSV the validator already streams carries submitter columns
> (`bioproject_center` / `collected_by` / `sra_center`). Filtering AST to submitters OUTSIDE the integrated
> public-health surveillance networks (NARMS / GenomeTrakr / PulseNet / CDC / FDA / USDA) yields a FREE,
> genome-linked, **provenance-disjoint (different-submitter-lab)** validation subset the decoder has never
> been run on. This Stage-1 census (no modelling, no genome download) asks only: **is that subset powered?**
> Script: `scripts/ncbi_pd_provenance_census.py`. Powering bar: ≥20 isolates per class, both classes.

## VERDICT: POWERED — the strand is NOT closed. A free provenance-disjoint validation exists.

| Organism × drug | with AST + downloadable genome | ecosystem (NARMS/CDC/FDA/…) | **OTHER (provenance-disjoint)** | powered? |
|---|---|---|---|---|
| Campylobacter × ciprofloxacin | 3,984 | 791R / 3030S | **33R / 130S** | **YES** |
| Klebsiella × ciprofloxacin | 816 | 109R / 26S | **404R / 277S** | **YES** |
| Salmonella × ciprofloxacin | 9,349 | 602R / 8656S | **4R / 87S** | **NO** (4R) |

**Organism-dependent:** the move works where AST is clinical-isolate-heavy (Klebsiella, Campylobacter) and
FAILS where it is surveillance-dominated (Salmonella — 602R/8656S are NARMS/CDC; the few non-ecosystem
submitters are veterinary, S-skewed → only 4R). So the provenance-disjoint subset must be checked per
organism; it is not universal. Klebsiella cipro is the clear Stage-2 lead (404R/277S, diverse clinical labs).

**Top non-ecosystem (provenance-disjoint) submitters found:**
- Campylobacter: The University of Melbourne (163).
- Klebsiella: Brigham & Women's Hospital (193), Walter Reed Army Institute of Research (93), CHU Habib
  Bourguiba (92), Day Zero Diagnostics (80), JCVI (57), DS-I Africa (42), Pelé Pequeno Príncipe (35),
  Broad Institute (24).

These are genuinely independent clinical/academic labs (different submitters, different countries) — not the
BV-BRC/NARMS/CDC ecosystem the decoder was tuned + cross-source-validated on. 2/2 organisms tested clear the
powering bar comfortably; Klebsiella in particular has a large, diverse, well-balanced (404R/277S) disjoint set.

## What this establishes (and its honest tier)

- **Establishes:** a free, genome-linked, **provenance-disjoint** validation subset EXISTS and is powered for
  ≥2 decoder substrates. The earlier "strand closed / do not pay" conclusion was an overclaim — corrected in
  the banners on `independent_phenotype_source_sweep_2026-06-10.md` + `independent_phenotype_label_census_2026-06-10.md`.
- **Independence tier (do NOT inflate):** this is **provenance-disjoint** (different submitter / lab /
  country), which is *stronger* than the same-ecosystem re-test but *weaker* than **methodology-disjoint**
  (most NCBI submitters still use CLSI broth microdilution) and weaker still than a fully-external clinical
  validation. The right headline for any Stage-2 result is *"holds on submitter/provenance-disjoint isolate
  phenotype not used in tuning"* — a stress-test against provenance leakage, NOT "external clinical validation."
- **Caveats:** (a) some "other" submitters may still deposit to NCBI via shared downstream normalization;
  (b) submitter strings are heuristic (substring match) — a Stage-2 run should pin exact `bioproject_acc`
  exclusion lists; (c) genome download + AMRFinder for the disjoint set is the Stage-2 cost (~95s/strain).

## Stage-2 — DONE for Klebsiella cipro (2026-06-10): decoder HOLDS, leakage-free

`scripts/provenance_disjoint_validate.py` scored the deployed `call_resistance(organism=Klebsiella,
drug=ciprofloxacin)` on **60 FRESH provenance-disjoint strains** (30R/30S, excluding all 97 prior-cohort
accessions — verified 0 overlap with tuning/prior-validation):

| metric | value |
|---|---|
| n scored | 60 (TP 29 · FP 1 · TN 29 · FN 1; 0 abstain) |
| accuracy / sensitivity / specificity | **0.967 / 0.967 / 0.967** |

**The deployed Klebsiella-cipro decoder generalizes** to a different-submitter-lab cohort it never saw —
0.967 across the board, vs the in-cohort LOO bal-acc 1.0 (the honest OOS number, 1 FP + 1 FN). Tier reminder:
provenance-disjoint, NOT methodology-independent. Artifact:
`wiki/provenance_disjoint_validation_klebsiella_cipro_2026-06-10.{md,json}`.

**LEAKAGE GOTCHA (caught live):** a first pass "preferred cached" accessions → pulled the 30-strain
calibration cohort into the "disjoint" set (leaked 0.967→0.983). Fix: the selector EXCLUDES every accession
in any prior `<slug>_*` cohort. Any future organism×drug Stage-2 must keep that exclusion.

Remaining Stage-2 candidates (free, same method): Campylobacter cipro (33R/130S disjoint pool). Salmonella is
NOT powered (4R). Per-organism, not universal.

## Bottom line
The `/brainstorm` adversarial pass paid off: the "no independent validation possible / spend nothing / closed"
conclusion was wrong. A **free, powered, provenance-disjoint** validation subset exists in data already on
hand. **No payment needed** (the original do-not-pay survives — but for the right reason: the free move exists,
not because validation is impossible). Stage-2 is the next deliberate build.
