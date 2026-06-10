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
| Salmonella × ciprofloxacin | (streaming — pending) | — | — | — |

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

## Stage-2 (the now-unblocked, free next action)
Score the deployed decoder (`call_resistance` / `evaluate_cohort`) on the provenance-disjoint subset per
organism×drug: download those `asm_acc` genomes, run AMRFinder (cached), compute acc/sens/spec, and report
strictly as "provenance-disjoint NCBI-PD validation" with the tier caveat. Start with Klebsiella cipro
(largest, best-balanced disjoint set: 404R/277S → subsample to a balanced ~50–75/class) then Campylobacter.
This is a real, free, higher-independence validation — the genuine continuation of Anchor 3.

## Bottom line
The `/brainstorm` adversarial pass paid off: the "no independent validation possible / spend nothing / closed"
conclusion was wrong. A **free, powered, provenance-disjoint** validation subset exists in data already on
hand. **No payment needed** (the original do-not-pay survives — but for the right reason: the free move exists,
not because validation is impossible). Stage-2 is the next deliberate build.
