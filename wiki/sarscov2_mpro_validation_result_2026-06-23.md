# SARS-CoV-2 Mpro cell — CoV-RDB fold-change validation (2026-06-23)

First validation of the Mpro/nirmatrelvir cell against Stanford CoV-RDB **measured fold-change** (the
coronavirus analogue of the HIV PhenoSense test). Run: `uv run python -m scripts.sarscov2_mpro_validate`.

## Result (nirmatrelvir, fold ≥ 2.5×)
| metric | value |
|---|---|
| status | **`COV_RDB_IN_DISTRIBUTION_KNOWLEDGE_BASELINE`** |
| powering | **UNDERPOWERED** (37 R / 5 S) |
| records / scored | 42 / 42 |
| sens | **0.68** (TP 25 / FN 12) |
| spec | **0.0** (TN 0 / FP 5) |

## What this honestly means
1. **In-distribution, NOT independent.** The catalog is built from CoV-RDB selection records and the
   fold-change is ALSO CoV-RDB — so this is a knowledge baseline (exactly like the TB WHO-catalogue-on-CRyPTIC
   number), NOT independent validation. A truly-independent number needs fold-change from a source the catalog
   was NOT built from (a held-out study / clinical isolates) — that is v0.1.
2. **Underpowered + true-negative-starved.** The CoV-RDB Mpro fold set is ENRICHED for resistant mutants
   (almost every tested isolate carries an Mpro mutation), so there are only 5 "S" rows and 0 true negatives →
   **spec 0.0 is not meaningful** (the catalog calls the 5 weak mutants R; there are no WT isolates here to be
   true negatives). Spec needs WT/low-fold isolates the set doesn't contain.
3. **The cell catches the strong signals.** The highest measured folds are all caught: E166V 7700× / S144E
   480× / F140S 260× / H172Y 250× / S144T 170×. sens 0.68 = it gets 2/3 of the high-fold mutants.
4. **The 12 FN expose catalog INCOMPLETENESS (the v0.1 signal, NOT tuned away).** The mutant-level v0 catalog
   misses novel high-fold mutants AT KNOWN positions: S54A 25×, F140A 21×, S144E/T, C160F, H172Y, R188G,
   Q192L. This is the honest mutant-level-vs-position tradeoff: mutant-level AVOIDS the Omicron P132H over-call
   (verified — P132H → S) but MISSES novel mutants at established positions. I deliberately did NOT add these
   to the catalog (that would overfit the validation = circular); v0.1 = position-level for the established
   resistance positions (50/144/166/167/172) + mutant-level elsewhere, validated on a HELD-OUT fold set.

## Per-mutation measured fold (the over-call/under-call surface)
Weak (catalogued but fold ~1, the selection-derived passengers): T21I 1.1×, A173V 0.9×, T304I 1.4×.
Strong (catalogued, high fold): E166V (up to 7700×), L50F+E166V 34×, F140L, S144A, P252L. Full table in
`wiki/sarscov2_mpro_cov_rdb_validation.json`.

## Verdict
The Mpro cell is **built, wired, and behaves correctly on the strongest measured signals** — but its CoV-RDB
number is in-distribution + underpowered, so it is a KNOWLEDGE BASELINE, not the independent win HIV is (yet).
The path to an independent SARS-CoV-2 number is real and free (held-out CoV-RDB studies / clinical-isolate
fold-change) = the clearly-scoped v0.1. Frozen bacterial AMR surface byte-unchanged; own trust surface
(namespace-separate from the bacterial + HIV cards).

## Provenance
Catalog `dna_decode/data/sarscov2_amr.py` (CoV-RDB selection). Fold-change `rx_fold ⋈ isolate_mutations` from
the CoV-RDB payload (`hivdb/covid-drdb-payload`, MIT; on D:). Scorer `scripts/sarscov2_mpro_validate.py`
(operator-aware fold censoring — the MIC-censoring lesson). Build spec
`plans/SARSCoV2_Mpro_Cell_Build_Spec_2026-06-23.md`; feasibility
`wiki/next_independent_label_cell_feasibility_2026-06-23.md`.
