# Intent contract — minister run `2026-06-03-1309-ep4-v01-expec-recall`

**Mission:** ep4-v01-expec-recall (maiden real `/soraya minister` run, BOUNDED sub-mission).
**Mission artifact:** `project_state/ep4-v01-expec-recall/big-idea.md` (validated).
**Gap (1, frozen):** `improve-expec-recall` — strategy-budget 2.
**Endpoints (gap retires when BOTH live-MET):**
1. `test-exit-0 python -m pytest tests/test_pathotype_expec_recall.py` (NEW test: ExPEC recall ≥ 0.85 on the
   24-genome H4 cohort AND confident-supported-call precision still 1.0).
2. `project-state-row project_state/ecoli-pathotype-prediction-cli-2026-05-26.md:ExPEC recall hardened`.

**Finiteness:** Φ = Σ(B_g − U_g) = 2 at freeze; strictly decreasing, never refunded. ≤ 2 proposals total.

**Gate model (T1=(c) — attended AUDIT-EVIDENCE, NOT airtight):** promotion needs 2 `/interrogate-me`
receipts. NOT high-stakes (local code + a test + a ledger row; no shipped skill / global-state / migration /
auth) → `/technical-plan` + `/brainstorm` NOT required.

**Money gate:** ONLY hard gate; armed via lease for the duration of each `run_to_stop`. Mission is
pure-Python CPU-only vs the coverage cache → no money/destructive action expected.

**Named residuals (not gated):** forgeable receipts/lease under bypassPermissions; PreToolUse hook blind to
post-admission MCP activity; guard self-protection deferred. OT1: lifecycle SEQUENCING is code
(`run_minister.py` + `minister_driver.py`); model supplies the generate / until-mvp / project-init seams.
