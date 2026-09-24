# PGx allele-frequency panel audit (2026-09-24)

Every asserted population frequency on the PGx surface, **derived from the 1000G panel on disk** rather than compared against a hand-entered value. Tolerance +/-0.1.

**CLAIMS_OUT_OF_BAND** — 6 of 6 claims measurable, 1 out of band, 0 unmeasurable (no local panel).

| gene | allele | pop | asserted | measured | delta | verdict |
|---|---|---|---|---|---|---|
| CYP4F2 | `*3` | EUR | 0.29 | 0.2773 | -0.0127 | ok |
| CYP4F2 | `*3` | EAS | 0.79 | 0.2128 | -0.5772 | **OUT OF BAND** — matches 1−ALT better than ALT |
| UGT1A1 | `*80` | EUR | 0.3 | 0.3009 | 0.0009 | ok |
| UGT1A1 | `*6` | EAS | 0.14 | 0.1368 | -0.0032 | ok |
| CYP2B6 | `*6-proxy` | ALL | 0.32 | 0.3196 | -0.0004 | ok |
| ABCG2 | `141K` | EUR | 0.09 | 0.0908 | 0.0008 | ok |
| ABCG2 | `141K` | EAS | 0.29 | 0.2872 | -0.0028 | ok |
| NUDT15 | `*3` | EAS | 0.1 | 0.0957 | -0.0043 | ok |

## Where each claim is asserted

- **CYP4F2 `*3`** (rs2108622) — dna_decode/pgx/cyp4f2.py docstring + cell_registry CYP4F2 contract
- **UGT1A1 `*80`** (rs887829) — dna_decode/pgx/runner.py UGT1A1 block + pgx_af_corroboration.py
- **UGT1A1 `*6`** (rs4148323) — scripts/pgx_af_corroboration.py
- **CYP2B6 `*6-proxy`** (rs3745274) — dna_decode/pgx/cyp2b6_catalog.py PROVENANCE block
- **ABCG2 `141K`** (rs2231142) — dna_decode/pgx/abcg2.py docstring + cell_registry ABCG2 contract
- **NUDT15 `*3`** (rs116855232) — dna_decode/pgx/nudt15_catalog.py + runner.py NUDT15 block

_bounds whether a shipped FREQUENCY claim matches this panel; NOT a phenotype validation._
