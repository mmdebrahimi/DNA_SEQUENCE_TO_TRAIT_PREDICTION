<!-- intent-contract: 1.0 -->
# Soraya --until-mvp — prospective-cohort-fetch (IMMUTABLE)
- Directive: finish the fetch script for the prospective-lock accrual cohort; download to D:\dna_decode_cache\data files donwload
- MVP bar (draft-then-ratify, execute-mode):
  1. scripts/fetch_prospective_cohort.py exists + pure logic unit-tested (test-exit-0)
  2. funnel: BV-BRC measured AST (ne Computational) x SCORED organisms/drugs, JOIN genome (assembly + date_inserted>lock pre-filter) -> NCBI Datasets release_date as the authoritative first_public -> prospective_lock.is_prospective_eligible (strictly-after 2026-06-13, fail-closed)
  3. emits cohort TSV (biosample, first_public_date, gca, drug, label) to D:\dna_decode_cache\data files donwload\
  4. live run produces the TSV OR honestly reports 0 eligible accrued; frozen surface byte-unchanged; suite green
- Gates: BV-BRC/NCBI reads = auto; in-cwd writes + D: output = auto (free public APIs, no money/dep-install)
