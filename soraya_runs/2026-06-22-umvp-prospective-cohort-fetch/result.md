# result — prospective-cohort-fetch
**Verdict: mvp-reached.** scripts/fetch_prospective_cohort.py finished, tested (8 offline tests), and
live-run end-to-end: queries the SCORED organisms (taxa 194/562/573) for measured AST + assemblies via
BV-BRC, resolves authoritative NCBI release dates, applies the prospective-lock temporal gate, and writes
the cohort TSV to D:\dna_decode_cache\data files donwload\prospective_cohort.tsv.

Live result today: 0 eligible post-lock rows (header-only TSV) — EXPECTED + honest, not a failure:
BV-BRC's newest E. coli ingestion is ~2026-04-19, ~2 months BEFORE the 2026-06-13 lock, and the prospective
window is only days old. The cohort ACCRUES over time; re-run periodically.

Recovery rounds: 1 (BV-BRC gt(date_inserted) needed the full ISO timestamp with URL-encoded colons + plain
&-clauses, not the and() wrapper — diagnosed via probe, fixed, re-verified).
Frozen AMR surface byte-unchanged; full suite 1588 passed.
