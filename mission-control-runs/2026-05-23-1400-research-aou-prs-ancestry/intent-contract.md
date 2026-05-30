<!-- run-format: 1.0 -->
<!-- intent-contract: 1.0 -->

# Intent Contract — aou-prs-ancestry — 2026-05-23 14:00 UTC

## Verbatim Input

```
All of Us research program polygenic risk score ancestry-stratified accuracy 2025
```

User context: Phase 2b L1 validation continuation. Second `/research` run testing harness-layer discoverability (vs 2026-05-22 manual execution). Topic extends prior PRS cross-ancestry portability research (`polygenic-risk-score-cross-ancestry-portability-2025-2026.md`) — fills the "No All of Us-specific cohort breakdown" honest gap from that run.

## Restated Decomposition

Standard `/research` v0.4 workflow: web search + synthesis → intake validation → followup queue update. Slash-command invocation via harness Skill tool (vs manual execution last run).

Sub-tasks:
1. Web research (Steps 1-3 produce `<slug>.raw.md`)
2. Intake validation (Step 4 produces `<slug>.md` + `<slug>_unsupported.md`)
3. Followup queue update (Step 5 updates `_followup_queue.md`)

## Method Used

Implicit decomposition by `/research` orchestrator pipeline (v0.4). Harness-layer slash-command invocation (Skill tool dispatch) — first end-to-end harness-discoverability validation.

## Verification Criteria

| Sub-task | Verification criterion | Evidence shape |
|---|---|---|
| Web research | ≥5 audit-grade rows OR honest-gap declared | `<slug>.raw.md` exists with V1 13-column table |
| Intake validation | Rows pass audit floor + mapping floor + banned-phrase scan + cite-token scan | `<slug>.md` (supported memo) + `<slug>_unsupported.md` |
| Followup queue update | `_followup_queue.md` modified-date ≥ run-start | mtime of `_followup_queue.md` ≥ run-start timestamp |

## Out-of-Scope Boundary

- Will NOT write outside `dna_decode/research_outputs/` or `dna_decode/mission-control-runs/<run-id>/` (plus global ledger at `~/.claude/mission-control-runs-ledger.md`)
- Will NOT modify rules YAML, ADSP units, wiki/* files
- Will NOT fabricate sources to hit row count targets
- Will NOT auto-promote candidates into rules/wiki
- Will NOT spend more than 15% of daily token budget (SOFT TARGET: aim for 5-8 quality sources to conserve)
- Will NOT run longer than 30 minutes wall-clock
- Will NOT make more than 100 tool calls

## Escalation Conditions

- Zero useful web search results after 3 attempts
- Intake rejects ≥80% of rows
- Slug validation fails
- Followup queue file write fails
- Any budget cap breached
- Unresolved-uncertainty count ≥5

## Run Metadata

- **Run ID:** 2026-05-23-1400-research-aou-prs-ancestry
- **Started:** 2026-05-23T14:00:00Z
- **Autonomy Level:** L1
- **Department(s):** Research
- **Skill invoked:** /research v0.4 (harness slash-command via Skill tool)
- **Caller:** user (Phase 2b harness-layer validation test)
- **Budget caps:**
  - Token: 15% daily
  - Wall-clock: 30 minutes
  - Tool-calls: 100
  - Unresolved-uncertainty: 5

## Contract version

v1.0
