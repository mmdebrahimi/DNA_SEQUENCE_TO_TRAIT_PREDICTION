# Prospective-lock accrual — first sweep (2026-07-10)

**Verdict: `ACCRUING` — a GENUINE zero.** 4 weeks after the lock (`2026-06-13`), no isolate is yet
provably first-public after it. `overall_status=OK`, so this zero is a real signal and not an outage.

## Funnel (NCBI Pathogen Detection, `latest_snps/Metadata/PDG*.metadata.tsv`)

| stage | count |
|---|---:|
| PD metadata rows (Campylobacter + E. coli/Shigella + Klebsiella) | 901,520 |
| with measured AST **and** a downloadable assembly | 13,423 |
| `sra_release_date` valid and **pre**-lock → skipped | 9,754 |
| `sra_release_date` valid and **post**-lock | **0** |
| `sra_release_date` undatable (`NULL`) → resolved by assembly `release_date` | 3,669 |
| resolved via NCBI Datasets | 3,669 |
| **prospective-eligible** | **0** |

Every one of the 3,669 `NULL`-dated rows resolved to an assembly released **2013–2015** (sampled: 8/8),
so they are old records with a missing SRA date, not new arrivals. The binding factor is **PD/NCBI
ingestion lag**, not the decoder.

## Why this run is trustworthy (two defects fixed to get it)

1. **A dead source must never be able to emit the "nothing accrued" signal.** BV-BRC (the original
   source) was fully down and answers an outage with **HTTP 200 + a `{"status":500,...}` error
   envelope**, so the first sweep printed `0 recent genomes`. Now `SourceUnavailable` is raised at the
   fetch boundary, per-scope status is recorded, and the script **refuses to write a cohort TSV** unless
   every scope queried cleanly (exit 2 / 1). A `--row-cap` smoke run likewise cannot claim ACCRUING.
2. **`"NULL" > "2026-06-13"`.** PD writes the literal string `NULL` for a missing `sra_release_date`, and
   a string `<=` pre-filter therefore admitted every NULL row as post-lock (inflating the funnel to
   3,669 "post-lock candidates"). Only the downstream fail-closed assembly-date check caught it.
   `is_iso_date()` now guards every date comparison, and the three date states are counted separately.

Also corrected: PD's `AST_phenotypes` is **comma-separated and quoted** (`"ampicillin=ND,ceftriaxone=R"`),
not semicolon-separated — the `;` form appears only in the derived `candidates.tsv`.

## Re-run

```bash
uv run python -m scripts.fetch_prospective_cohort --source ncbi_pd
```
Exit 0 = OK (a real ACCRUING zero, or a written cohort). Exit 2 = source unavailable. Exit 1 = partial
outage or truncated smoke run. Score a non-empty cohort with
`scripts/prospective_lock_validate.py --cohort-tsv`.
