# Execution Log — Clonality_Disclosure_Layer_Plan
Date: 2026-06-11
Waves: 4 (max parallelism 1 — linear data flow 1→2→3→4)
Files changed: dna_decode/eval/clonality.py, dna_decode/data/cell_key.py, scripts/compute_lineage_metrics.py, scripts/build_validation_report_card.py, tests/test_clonality.py, tests/test_compute_lineage_metrics.py, tests/test_build_validation_report_card.py, CLAUDE.md, README.md, docs/ARCHITECTURE.md, LESSONS_LEARNED.md, wiki/decisions-log.md, .claude/testable-modules.md
Sentrux verdict: n/a (sentrux absent on host)
Commit: 6837836 (sequential commits on main; no PR — solo workflow)

Notes:
- All 4 steps fully implemented; none skipped/blocked.
- Tests: baseline 1000 → final 1053 passed (43 net new), 0 regressions. `tests/test_models_foundation.py` ignored throughout (pre-existing torch/paging-file OSError WinError 1455 on this host, unrelated to this plan).
- Manual verification DEFERRED: the Docker Mash run (`scripts/compute_lineage_metrics.py`) was not executed — it requires Docker Desktop. All math validated via pure-helper unit tests + an M4 reconcile spot-check against the real `klebsiella_provdisjoint_ciprofloxacin` cohort (recomputed 29/1/29/1 == committed artifact). Run `uv run python scripts/compute_lineage_metrics.py` on a Docker-equipped host, then rebuild the report card, to populate `wiki/provdisjoint_lineage_metrics.json` with real lineage data.
