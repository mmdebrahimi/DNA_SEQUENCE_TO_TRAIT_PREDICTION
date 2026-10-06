# Execution Log — Cross_Organism_Essentiality_Transfer_Ladder_V2_Coverage_Primary_Plan

Date: 2026-10-06
Waves: 9 (sequential — max parallelism 1; no remote-parallel mode, steps 2→6 are a serialising schema→parser→runner chain)
Files changed: dna_decode/eval/transfer_ladder.py, dna_decode/eval/regime.py, dna_decode/essentiality/annotation_join.py, dna_decode/essentiality/phrasing_floor.py, dna_decode/essentiality/label_sources.py, dna_decode/essentiality/cli.py, scripts/essentiality_ladder_w0_probe.py, scripts/essentiality_transfer_ladder.py, scripts/essentiality_missed_vocabulary.py, scripts/build_essentiality_report_card.py, tests/test_transfer_ladder.py, tests/test_essentiality_ladder_w0_probe.py, tests/test_essentiality_annotation_join.py, tests/test_essentiality_phrasing_floor.py, tests/test_essentiality_label_sources.py, tests/test_essentiality_transfer_ladder_runner.py, tests/test_essentiality_transfer_ladder_verdict.py, tests/test_essentiality_report_card_ladder.py, tests/test_regime_boundary.py, wiki/essentiality_transfer_ladder_2026-10-06.{md,json}, wiki/essentiality_ladder_w0_probe_2026-10-06.json, wiki/essentiality_transfer_ladder_memo_template.md, wiki/essentiality_report_card.{md,json}, wiki/learned_regime_map.json, wiki/decisions-log.md, CHANGELOG.md, CLAUDE.md, README.md, LESSONS_LEARNED.md, TODOS.md, FUTURE_FEATURES.md
Sentrux verdict: n/a — sentrux not installed on this host
Commit: 1777b9d (sequential, committed directly to `main`; run spans 5290311..1777b9d)

## Notes

- **Frozen-surface invariant held.** The five sha256-pinned files were byte-unchanged for the whole run
  and `verify_lock` re-verified OK afterwards; `core_decoder.py` is byte-identical
  (`92907f774ca6df29cd2c7352faa529664badf269c0b47b76279c79b87d7251ee`), so the ladder applies the decoder
  UNCHANGED rather than a variant of it.
- **The reconcile gate earned its place.** Step 6's first run raised `ReconcileMismatch` (0.3122 vs the
  committed 0.3125) and the cause was a duplicate gene symbol (`mrcB`), not a metric bug. Both dedup rules
  were measured against the committed targets rather than the target being tuned to the code.
- **One verdict retracted pre-publication.** `PHRASING_GAP_DOMINATES` was withdrawn on measurement; the
  shipped verdict is the one the numbers support.
- **State-file deviation recorded, not hidden.** The per-step state updates required by sequential-mode.md
  were missed mid-run (only step 1 was recorded); corrected at the retrospective with a
  `state_file_maintenance_note` rather than back-filling silently.
- **Four open user-authority decisions were surfaced, not taken** (step-3 scope; a label file for the three
  walled rungs; `on_duplicate="last"` vs `"first"`; a 4th gene for the supervised-blind-spot arm).
