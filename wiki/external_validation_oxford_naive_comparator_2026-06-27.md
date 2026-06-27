# Oxford external validation — the curated layer BEATS naive AMRFinder use (2026-06-27)

Closes the wrapper-vs-underlying-tool gap on the 2026-06-15 Oxford independent validation. That
run reported the FROZEN decoder number alone; the committed rail (*a policy layer over a curated-DB
tool must beat NAIVE use of that tool on INDEPENDENT data, else the in-cohort number only proves the
tool works, not that the layer adds value*) requires a paired baseline. This is that baseline, on the
SAME measured-MIC labels, same Oxford cohort (Lipworth et al., `github.com/samlipworth/ecoli_mic_arg`,
PRJNA604975 + PRJNA1007570; 2897 isolates, broth-microdilution MIC + the cohort's own AMRFinder).

## Definitions
- **naive** = the non-expert use of the raw tool: R iff ANY AMRFinder determinant whose `Class` is in
  `mic_tiers.amrfinder_classes_for(drug)` is present — broad class match, NO subclass/point/threshold
  refinement, no abstain.
- **frozen** = the shipped `call_resistance` DRUG_RULE (subclass_any / qrdr_point / threshold curation).
- **metric** = balanced accuracy `(sens+spec)/2` — the honest net: a naive baseline games one axis
  (call everything R → sens≈1, spec≈0), so a per-axis comparison flatters it; balacc nets the over-call out.

## Result — curated beats naive on ALL 3 drugs

| drug | frozen sens/spec | naive sens/spec | frozen balacc | naive balacc | Δ balacc | verdict |
|---|---|---|---|---|---|---|
| ciprofloxacin | 0.935 / 0.963 | 0.962 / 0.571 | 0.949 | 0.767 | **+0.18** | CURATED_LAYER_ADDS_VALUE |
| ceftriaxone | 0.945 / 0.709 | 1.000 / **0.003** | 0.827 | 0.502 | **+0.33** | CURATED_LAYER_ADDS_VALUE |
| gentamicin | 0.922 / 0.995 | 0.969 / 0.459 | 0.959 | 0.714 | **+0.24** | CURATED_LAYER_ADDS_VALUE |

n scored: cipro 2841, cef 2868, gent 2873 (binary measured-MIC labels, frozen and naive on the identical set).

## Reading it
- **The curated layer demonstrably adds value over the raw tool on fully independent data** (+0.18 to
  +0.33 balanced accuracy). The gain is entirely on SPECIFICITY (+0.39 / +0.71 / +0.54), at a trivial
  sens cost (−0.03 to −0.06) that is just the naive baseline's over-call (it "wins" sens only by calling
  nearly everything R).
- **cef is the sharpest demonstration:** naive spec = **0.003** — near-universal intrinsic chromosomal
  AmpC (`blaEC`, `Class=BETA-LACTAM`) makes "any β-lactam determinant → R" call essentially every E. coli
  R. The curated `subclass_any` extended-spectrum refinement lifts spec to 0.709 even on the cohort's OLD
  AMRFinder genotype (the deployed v4.2.7 cef number is higher still — see the 2026-06-15 confound note;
  blaEC reclassifies to `BETA-LACTAM/BETA-LACTAM` and drops out). So the curated rule's value-add is
  understated here, not overstated.
- This is the first time the decoder's value-over-baseline is shown on an INDEPENDENT cohort (the prior
  k-mer / domain-knowledge baselines were in-cohort). It complements, does not replace, the 2026-06-15
  cipro/gent independence headline.

## Honesty + scope
- The naive baseline uses the SAME (older) cohort AMRFinder genotype the frozen rule was scored on, so the
  comparison isolates the RULE's contribution cleanly (same input, different policy layer).
- This validates the curated rule's value-add; it is NOT a claim that the deployed Docker pipeline was run
  on Oxford (Oxford is reads-only / 0 GCA assemblies → the genome-download revalidation arm cannot fetch
  it; the pure-join path is the correct and only executable route for this cohort).
- **FROZEN AMR surface byte-unchanged** (`amr_rules.py` + `calibrated_amr_rules.json`); this is a new
  standalone comparator script, READ-only over the frozen rule.

## Artifacts
`wiki/external_validation_oxford_naive_comparator_2026-06-27.json`; scorer `scripts/oxford_naive_baseline.py`
(reuses `scripts/oxford_score.py` helpers + `mic_tiers.amrfinder_classes_for`); data (gitignored)
`data/raw/oxford/`. Companion: `wiki/oxford_external_validation_result_2026-06-15.md` (the frozen-number headline).
