<!-- run-format: 1.0 -->
<!-- intent-contract: 1.0 -->

# Intent Contract — prs-cross-ancestry — 2026-05-22 21:15 UTC

## Verbatim Input

```
Polygenic risk score cross-ancestry portability (2025-2026)
```

User context (Phase 2b L1 validation gate, Mission Control v1 Roadmap): user picked this topic from a 4-option AskUserQuestion menu after approving the L0→L3 autonomy plan. Topic chosen because it's directly relevant to user's DNA decoder project (verdict v2 flagged ancestry-bias as critical technical-buyer concern) and exercises audit-rich peer-reviewed sources.

## Restated Decomposition

Standard `/research` v0.3 workflow: web search + synthesis → intake validation → followup queue update.

Sub-tasks:
1. Web research (Steps 1-3 produce `<slug>.raw.md`)
2. Intake validation (Step 4 produces `<slug>.md` + `<slug>_unsupported.md`)
3. Followup queue update (Step 5 updates `_followup_queue.md`)

## Method Used

Implicit decomposition by `/research` orchestrator pipeline (v0.3). MANUALLY EXECUTED in this session (not via slash-command invocation) because `/research` skill not in session-start-cached harness skill list. Validates the L1 contract DESIGN even though it bypasses slash-command discovery layer.

## Verification Criteria

| Sub-task | Verification criterion | Evidence shape |
|---|---|---|
| Web research | ≥5 audit-grade rows OR honest-gap declared | `<slug>.raw.md` exists with V1 13-column table |
| Intake validation | Rows pass audit floor + mapping floor + banned-phrase scan + cite-token scan | `<slug>.md` (supported memo, ≥1 row OR explicit empty-memo marker) + `<slug>_unsupported.md` |
| Followup queue update | `_followup_queue.md` modified-date ≥ run-start | mtime of `_followup_queue.md` ≥ run-start timestamp |

## Out-of-Scope Boundary

This run explicitly will NOT:
- Write outside `dna_decode/research_outputs/` or `dna_decode/mission-control-runs/<run-id>/` (plus the global ledger at `~/.claude/mission-control-runs-ledger.md`)
- Modify rules YAML, ADSP units, wiki/* files
- Fabricate sources to hit row count targets
- Auto-promote candidates into rules/wiki (Promotion Gate remains manual)
- Spend more than 15% of daily token budget (PHASE 2B SOFT TARGET: aim for 5-8 quality sources rather than full 15-25 to conserve)
- Run longer than 30 minutes wall-clock
- Make more than 100 tool calls (WebSearch + WebFetch + Skill + Write)

## Escalation Conditions

- Zero useful web search results after 3 attempts → halt + ask user to refine topic
- Intake rejects ≥80% of rows → halt + surface likely topic-shape failure
- Followup queue file write fails → halt; supported memo still preserved
- Any budget cap breached → halt + ask user to continue/stop/replan
- Unresolved-uncertainty count ≥5 → halt + surface confused state

## Run Metadata

- **Run ID:** 2026-05-22-2115-research-prs-cross-ancestry
- **Started:** 2026-05-22T21:15:00Z
- **Autonomy Level:** L1
- **Department(s):** Research
- **Skill invoked:** /research (manual execution — slash command not yet in harness skill list this session)
- **Caller:** user (via L0_to_L3_Autonomy_Path_Plan Phase 2b)
- **Budget caps:**
  - Token: 15% daily (soft target: lower to conserve session budget)
  - Wall-clock: 30 minutes
  - Tool-calls: 100
  - Unresolved-uncertainty: 5

## Contract version

v1.0
