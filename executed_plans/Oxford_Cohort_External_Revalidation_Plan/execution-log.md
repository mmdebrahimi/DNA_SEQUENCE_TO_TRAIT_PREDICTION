# Execution Log — Oxford_Cohort_External_Revalidation_Plan
Date: 2026-06-15
Waves: 5 (max parallelism 2; executed SEQUENTIALLY — forced to protect the frozen-file invariant)
Files changed: dna_decode/eval/biosample_resolver.py, dna_decode/data/external_mic_labels.py, dna_decode/data/external_cohort_genomes.py, scripts/external_cohort_preflight.py, scripts/external_cohort_revalidate.py, scripts/build_external_validation_report.py, tests/test_biosample_resolver.py, tests/test_external_cohort_preflight.py, tests/test_external_mic_labels.py, tests/test_external_cohort_genomes.py, tests/test_external_cohort_revalidate.py, tests/test_build_external_validation_report.py, CLAUDE.md, README.md, LESSONS_LEARNED.md, TODOS.md, wiki/decisions-log.md
Sentrux verdict: n/a (sentrux not installed)
Commit: 262220f (direct-to-main; no PR — project syncs via main. Step commits 0b1fba0..262220f)
Tests: baseline 1056 -> final 1142 (+86 offline unit tests across 6 new modules; 0 regressions)
Frozen invariant: HELD — amr_rules.py / mic_tiers.py / cohort_manifest.py / build_validation_report_card.py / compute_lineage_metrics.py byte-unchanged
Deferred (manual, network+Docker): Gate-0 preflight on PRJNA604975 -> per-drug scorer -> roll-up (the empirical go/no-go run)
