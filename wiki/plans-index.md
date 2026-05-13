# Plans Index
<!-- Auto-maintained by /save-plan. Do not edit manually. -->

## [plan_file: Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md] 2026-05-13
**Summary:** Scope-reduced delta from `Audit_Calibration_NT_AutoModel_Plan.md` after `/review` synthesis — 5 steps / 4 waves → 2 commits, 2 waves. Drops dual-verdict columns (institutionalization risk), drops wiki update (no text to replace), splits NT refactor into a separate gated commit (equivalence test required).
**Key decisions:**
- Asymmetric warning banner replaces dual-verdict columns (D1)
- Drop wiki/GATE_B_REPORT.md update entirely (D2)
- Split NT refactor into a separate, gated commit (D3)
- Keep default-semantics test as the regression lock (D4)
- thresholds_block(rules) helper, not inline string list (D5)

---

## [plan_file: Audit_Calibration_NT_AutoModel_Plan.md] 2026-05-13
**Summary:** Fix a credibility bug in the just-shipped audit cohort generator (`scripts/audit_cohort.py`) — the "GO" verdict was emitted under silently-relaxed thresholds; defaults produce "WARN" — AND simultaneously replace `NucleotideTransformerModel`'s `AutoModelForMaskedLM` with `AutoModel` to eliminate the `output_hidden_states=True` workaround.
**Key decisions:**
- Audit report header MUST surface threshold values (D1)
- Two verdict columns — Phase 1 production + Gate B infra-only (D2)
- Pin default semantics in tests (D3)
- NT switches to `AutoModel`, not `AutoModelForMaskedLM` (D4)
- Pooling-strategy tag stays "single_seq_mean" (D5)

---

## [plan_file: BVBRC_Genome_Metadata_Adapter_Plan.md] 2026-05-12
**Summary:** Wire `BVBRC_genome.csv` (BV-BRC Genomes-tab export) into the cohort path as a new adapter module, bypassing the wrong-contract `pilot.fetch_ncbi_assembly_quality` scaffold and feeding the existing `--assembly-metadata` wire that `cohort.candidates_from_bvbrc_ast` already accepts.
**Key decisions:**
- Bypass the scaffold instead of implementing it (D1)
- New CLI flag rather than overloading existing `--assembly-metadata` (D2)
- Coverage-log line surfaces ID-namespace mismatches early (D3)
- `fetch_ncbi_assembly_quality` stays scaffolded (D4)

---

## [plan_file: Ecoli_G2P_Phase1_Closeout_Plan.md] 2026-05-12
**Summary:** Wrap up the stalled `/execute-plan` epilogue for `Ecoli_G2P_Phase1_Ship_Path_Plan.md` — toolchain restore, doc reconciliation, first authoritative test pass, archive, state cleanup, push, final report.
**Key decisions:**
- Selective expansion over hold scope (D1)
- Real-data validation = Phase 2 entry criterion, not Phase 1 closeout (D2 — pending)
- Toolchain restore approach: uv vs pip (D3 — pending)
- Archival convention: status-header + git tag recommended (D4 — pending)
- `/documentation` before commit (D5)
- Retrospective re-derivation, not skip (D6)
- Test outcome recorded, not gated (D7)
- Delete both state files at end (D8)
- `.claude/execute-plan-state/` added to `.gitignore` (D9)

---

## [plan_file: Ecoli_G2P_Phase1_Ship_Path_Plan.md] 2026-05-12
**Summary:** Contracted path to ship Phase 1 of `Ecoli_G2P_Platform_Technical_Plan.md`. Captures the `/review` synthesis verdict (HOLD scope + selective contraction within remaining steps) plus the deferred Wave 3.5 hardening fixes from the post-Wave-3 `/brainstorm`. Estimated remaining work: ~700 LOC across 5 implementation steps + 4 hardening edits.
**Key decisions:**
- HOLD scope, do not expand (D1)
- Reorder — Step 15 (smoke + fixtures) BEFORE Step 14 (CLI) (D2)
- Step 14 collapses to one `scripts/pipeline.py` with subcommands (D3)
- Step 13 visualization uses matplotlib + TSV export, NOT pygenometracks (D4)
- Step 17 leaderboard collapses to a shell loop over `pipeline.py train` (D5)
- Step 16 docs trimmed to README + ARCHITECTURE.md only (D6)
- Apply Wave 3.5 hardening BEFORE Step 14 wiring fires (D7)
- Add quantization-fidelity micro-step (selective addition) (D8)

---
