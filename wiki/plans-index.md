# Plans Index
<!-- Auto-maintained by /save-plan. Do not edit manually. -->

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
