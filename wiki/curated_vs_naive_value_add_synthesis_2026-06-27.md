# The curated layer adds value over naive AMRFinder use — cross-surface synthesis (2026-06-27)

Closes the **wrapper-vs-underlying-tool** rail across BOTH validation surfaces. The committed lesson:
a curated policy layer over a curated-DB tool must BEAT *naive use of that tool* on independent data —
else the headline number only proves the tool works, not that the layer adds value. The report card and
the Oxford result reported the frozen `call_resistance` number ALONE; this is the paired naive baseline.

- **naive** = the non-expert use of AMRFinder: R iff ANY determinant whose `Class` is in
  `mic_tiers.amrfinder_classes_for(drug)` is present — broad class match, NO subclass/point/threshold
  refinement, no abstain.
- **frozen** = the shipped `call_resistance` DRUG_RULE.
- **metric** = balanced accuracy `(sens+spec)/2` (a naive baseline games one axis: call everything R).

## Surface 1 — Oxford independent measured-MIC cohort (2897 isolates)
`scripts/oxford_naive_baseline.py` → `wiki/external_validation_oxford_naive_comparator_2026-06-27.json`

| drug | frozen balacc | naive balacc | Δ | verdict |
|---|---|---|---|---|
| ciprofloxacin | 0.949 | 0.767 | **+0.18** | CURATED_LAYER_ADDS_VALUE |
| ceftriaxone | 0.827 | 0.502 | **+0.33** | CURATED_LAYER_ADDS_VALUE |
| gentamicin | 0.959 | 0.714 | **+0.24** | CURATED_LAYER_ADDS_VALUE |

## Surface 2 — the 10-cell provenance-disjoint trust surface (NCBI-PD)
`scripts/naive_baseline_provdisjoint.py` → `wiki/provdisjoint_naive_comparator_2026-06-27.json`.
Reuses each cell's ALREADY-CACHED per-isolate AMRFinder `main.tsv` (no Docker, no re-download). **M4
reconciliation guard:** the frozen confusion matrix recomputed from cache MUST equal the committed
provdisjoint artifact before any delta is trusted.

| cell | frozen balacc | naive balacc | Δ | verdict |
|---|---|---|---|---|
| Campylobacter cipro | 1.000 | 1.000 | +0.00 | NAIVE_TIES_CURATED |
| E. coli cipro | 0.817 | 0.650 | **+0.17** | CURATED_LAYER_ADDS_VALUE |
| E. coli ceftriaxone | 0.967 | 0.734 | **+0.23** | CURATED_LAYER_ADDS_VALUE |
| E. coli gentamicin | 0.950 | 0.633 | **+0.32** | CURATED_LAYER_ADDS_VALUE |
| E. coli tetracycline | 0.933 | 0.850 | **+0.08** | CURATED_LAYER_ADDS_VALUE |
| Klebsiella cipro | 0.967 | 0.700 | **+0.27** | CURATED_LAYER_ADDS_VALUE |
| Klebsiella meropenem | 0.683 | 0.500 | **+0.18** | CURATED_LAYER_ADDS_VALUE |
| Klebsiella ceftriaxone | — | — | — | RECONCILE_MISMATCH (cache 54/60) |
| Klebsiella gentamicin | — | — | — | RECONCILE_MISMATCH (cache 3/60) |
| Klebsiella tetracycline | — | — | — | RECONCILE_MISMATCH (cache 33/60) |

**7/10 cells reconcile; 6/7 CURATED_LAYER_ADDS_VALUE.**

## Reading it
- **Across both surfaces (3 + 7 = 10 reconciled cells), the curated layer adds value on 9 and ties on 1**
  — the first systematic demonstration of value-over-baseline on INDEPENDENT / provenance-disjoint data
  (prior baselines were in-cohort k-mer / domain-knowledge controls).
- The gain is almost entirely on **specificity**: naive over-calls because AMRFinder reports benign /
  intrinsic class determinants. The curated `subclass_any` / `qrdr_point` / `gene_prefixes` refinement is
  exactly what removes them. Two sharpest cases:
  - **ceftriaxone** (both surfaces): naive calls nearly every E. coli R on the near-universal intrinsic
    chromosomal AmpC `blaEC` (`Class=BETA-LACTAM`); the extended-spectrum refinement is the whole value.
  - **Klebsiella ciprofloxacin** (+0.27): the intrinsic `OqxAB` efflux pump (`Class` matches the broad
    fluoroquinolone set) makes naive over-call; the curated acquired/strength-specific rule fixes it —
    the exact pattern the `intrinsic-genes-break-broad-AMR-class-rules` lesson predicted.
- **The one tie (Campylobacter cipro) is honest, not a weakness:** Campylobacter cipro resistance is the
  `gyrA` T86I point mutation, and AMRFinder reports no benign quinolone-class determinants in this set —
  so the only class-matching determinant IS the resistance one, and there is nothing for curation to
  refine. The rail's scope is "curation beats naive WHERE the tool reports confounders"; this correctly
  shows where it doesn't.

## Honesty + scope
- **M4 guard is load-bearing:** the 3 Klebsiella RECONCILE_MISMATCH cells have incomplete local AMRFinder
  caches (6/57/27 isolates missing vs the committed scored set) — so NO delta is reported for them rather
  than a misleading partial number. Completing them is a heavy-Docker follow-up (~90 AMRFinder re-runs +
  genome fetches); the Klebsiella value-add pattern is already established by the cipro (+0.27) +
  meropenem (+0.18) reconciled cells.
- **FROZEN AMR surface byte-unchanged** (`amr_rules.py` + `calibrated_amr_rules.json`); both comparators
  are standalone scripts, READ-only over the frozen rule. Leak guard 9/9.

## Artifacts
`scripts/oxford_naive_baseline.py` + `scripts/naive_baseline_provdisjoint.py`;
`wiki/external_validation_oxford_naive_comparator_2026-06-27.{json}` +
`wiki/provdisjoint_naive_comparator_2026-06-27.json`; pure-helper tests `tests/test_naive_comparator.py`
(8 tests). Companion headlines: `wiki/oxford_external_validation_result_2026-06-15.md`,
`wiki/decoder_validation_report_card.md`.
