# Execution Log — Oxford_MIC_Ingester_Plan
Date: 2026-06-15
Waves: 7 (max parallelism 3; executed SEQUENTIALLY — 4 shipped arm modules modified + 5-frozen-file invariant)
Files changed: dna_decode/eval/biosample_resolver.py, dna_decode/data/external_mic_labels.py, dna_decode/data/external_crosswalk.py, dna_decode/data/external_mic_ingest.py, scripts/oxford_w0_probe.py, scripts/build_oxford_labels.py, scripts/external_cohort_preflight.py, scripts/external_cohort_revalidate.py, scripts/build_external_validation_report.py, scripts/run_oxford_revalidation.py, 10 tests/test_*.py, CLAUDE.md, README.md, wiki/oxford_revalidation_runbook.md
Sentrux verdict: n/a (sentrux not installed)
Commit: 2edad52 (direct-to-main; no PR. Step commits 8273ba8..2edad52)
Tests: baseline 1161 -> final 1238 (+77 offline unit tests; 0 regressions)
Frozen invariant: HELD — amr_rules.py / mic_tiers.py / cohort_manifest.py / build_validation_report_card.py / compute_lineage_metrics.py byte-unchanged
Deferred (manual, network+Docker): fetch the Oxford MIC table -> W0 probe -> run_oxford_revalidation.py (the empirical go/no-go)
