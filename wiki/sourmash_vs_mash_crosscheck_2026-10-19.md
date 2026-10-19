# sourmash vs Mash crosscheck — OK

- bar: `wiki\sourmash_vs_mash_crosscheck_bar.json`
- sourmash: ksize=21 scaled=1000
- refseq cache: `C:\Users\Laptop Portable\pythonprojects\dna_decode\data\refseq_cache`
- pairs observed: 20 / 20

## Clauses

| # | clause | met | detail |
|---|---|---|---|
| 1 | verdict_identity | yes | 20/20 verdicts identical (NECESSARY BUT NEARLY VACUOUS -- 19/20 PASS by margin alone) |
| 2 | exact_minority_count_near_the_bar | yes | 3/3 near-bar pairs matched exactly (DISCRIMINATING) |
| 3 | the_trip_must_reproduce | yes | 1/1 TRIP pair(s) reproduced (implied by clause 1; retained as a severity pre-commitment) |

## Honest limits

- FracMinHash is not Mash's bottom-sketch MinHash, so agreement on these cohorts does not make the two estimators interchangeable in general.
- The thresholds (0.001, 0.005) were calibrated against MASH distances. Passing this bar shows sourmash reproduces Mash's PARTITION VERDICTS at those thresholds on these cohorts; it does not re-derive the thresholds.
- Clause 3 rests on a SINGLE TRIP pair. One negative is enough to make the bar non-vacuous and far too few to characterise sourmash's behaviour near the bar in general.
- These are provenance-disjoint AMR cohorts, not CABBAGE. Accepting the distance source for CABBAGE is an extrapolation across cohorts.
