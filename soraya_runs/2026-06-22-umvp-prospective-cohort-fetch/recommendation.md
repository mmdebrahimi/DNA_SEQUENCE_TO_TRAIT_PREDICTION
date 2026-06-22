# recommendation — prospective-cohort-fetch
The standing accrual pipeline is shipped + correct. Remaining gap (0 rows today) is EXTERNAL (data
ingestion lag + a days-old window), not code-closable by more building on this script.

Two forward options (user's call):
- WAIT: re-run `python -m scripts.fetch_prospective_cohort` periodically; BV-BRC will surface post-lock
  isolates as it ingests them (months out). No further code needed.
- FRESHER SOURCE (code-closable, if sooner accrual matters): add NCBI Pathogen Detection (continuous
  ingestion + antibiogram) as a second source behind the same funnel/eligibility gate. Larger build;
  only worth it if you want prospective rows before BV-BRC catches up.
