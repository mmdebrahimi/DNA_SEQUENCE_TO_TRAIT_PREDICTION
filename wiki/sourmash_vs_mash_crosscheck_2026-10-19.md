# sourmash vs Mash crosscheck — INSUFFICIENT_DATA

- bar: `wiki\sourmash_vs_mash_crosscheck_bar.json`
- sourmash: ksize=21 scaled=1000
- refseq cache: `C:\Users\Laptop Portable\pythonprojects\dna_decode\data\refseq_cache`
- pairs observed: 0 / 20

## INSUFFICIENT_DATA — no agreement number

A cohort whose assemblies do not all resolve produces LOWER-BOUND effective-lineage counts. Reporting an agreement rate over a partial cohort would compare a lower bound against an exact-equality clause, so no agreement number is emitted. INSUFFICIENT_DATA is a valid result.

| cohort | reason |
|---|---|
| campylobacter_provdisjoint_ciprofloxacin | no selected.tsv |
| escherichia_coli_shigella_provdisjoint_ceftriaxone | no selected.tsv |
| escherichia_coli_shigella_provdisjoint_ciprofloxacin | no selected.tsv |
| escherichia_coli_shigella_provdisjoint_gentamicin | no selected.tsv |
| escherichia_coli_shigella_provdisjoint_tetracycline | no selected.tsv |
| klebsiella_provdisjoint_ceftriaxone | no selected.tsv |
| klebsiella_provdisjoint_ciprofloxacin | no selected.tsv |
| klebsiella_provdisjoint_gentamicin | no selected.tsv |
| klebsiella_provdisjoint_meropenem | no selected.tsv |
| klebsiella_provdisjoint_tetracycline | no selected.tsv |

## Clauses

| # | clause | met | detail |
|---|---|---|---|
| 1 | verdict_identity | — | not evaluated -- INSUFFICIENT_DATA |
| 2 | exact_minority_count_near_the_bar | — | not evaluated -- INSUFFICIENT_DATA |
| 3 | the_trip_must_reproduce | — | not evaluated -- INSUFFICIENT_DATA |

## Honest limits

- FracMinHash is not Mash's bottom-sketch MinHash, so agreement on these cohorts does not make the two estimators interchangeable in general.
- The thresholds (0.001, 0.005) were calibrated against MASH distances. Passing this bar shows sourmash reproduces Mash's PARTITION VERDICTS at those thresholds on these cohorts; it does not re-derive the thresholds.
- Clause 3 rests on a SINGLE TRIP pair. One negative is enough to make the bar non-vacuous and far too few to characterise sourmash's behaviour near the bar in general.
- These are provenance-disjoint AMR cohorts, not CABBAGE. Accepting the distance source for CABBAGE is an extrapolation across cohorts.
