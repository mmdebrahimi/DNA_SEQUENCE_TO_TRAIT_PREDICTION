<!-- run-format: 1.0 -->
<!-- intent-contract: 1.0 -->

# Intent Contract — ecoli-pathotype-substrate — 2026-05-27T01:40Z

## Verbatim Input

```
E. coli pathotype labeled-genome substrate survey for an open-source pathotype prediction CLI. Focus on per-record label provenance metadata schemas (clinical / outbreak-investigation / lab-pathotype-assay / curated-literature vs gene-rule-derived) across these candidate sources: (1) EnteroBase E. coli database — what label fields, how labels were assigned, what fraction is gene-rule-derived vs independent annotation; (2) NCBI Pathogen Detection isolate browser — pathotype-label fields and provenance; (3) GenomeTrakr — pathotype labeling discipline; (4) EcoCyc DEC reference strains panel — curated reference isolates with independent labels; (5) published E. coli outbreak isolate sets (BioProject-level deposits with curated pathotype labels). For each source: (a) document the label-provenance metadata schema, (b) estimate the independent-label fraction (target floor: ≥70% labels independent of v0 marker rules), (c) report per-class isolate-count availability against pathotype floors N≥50 for EHEC/STEC, EPEC, ETEC, EAEC, UPEC/ExPEC + N≥75 commensal/low-marker. Audit-quality output: cite primary sources; flag where source documentation is silent on label provenance; rank sources by H1-passing likelihood. Project context: this gates the architecture-fork lock for a deterministic multilabel cluster resolver + abstention v0 CLI (per project_state/ecoli-pathotype-prediction-cli-2026-05-26.md).
```

## Restated Decomposition

Standard /research v0.3 workflow: web search + synthesis → intake validation → followup queue update.

Sub-tasks:
1. Web research (Steps 1-3 produce `<slug>.raw.md`) — survey 5 candidate substrate sources for E. coli pathotype labels, extract per-source label-provenance metadata + per-class isolate counts + independent-label fraction estimates
2. Intake validation (Step 4 produces `<slug>.md` + `<slug>_unsupported.md`)
3. Followup queue update (Step 5 updates `_followup_queue.md`)

## Method Used

Implicit decomposition by `/research` orchestrator pipeline (v0.3) — five candidate sources are the natural axes; per-source ≥3 rows targeted (provenance schema + label-derivation method + per-class counts).

## Verification Criteria

| Sub-task | Verification criterion | Evidence shape |
|---|---|---|
| Web research | ≥5 audit-grade rows OR honest-gap declared | `<slug>.raw.md` exists with V1 13-column table; row count + honest-gap status visible |
| Intake validation | Rows pass audit floor + mapping floor + banned-phrase scan + cite-token scan + source-identity advisory | `<slug>.md` (supported memo, ≥1 row OR explicit empty-memo marker) + `<slug>_unsupported.md` (rejected rows) exist |
| Followup queue update | `_followup_queue.md` exists + has been touched within this run | mtime of `_followup_queue.md` is ≥ run-start timestamp |

## Out-of-Scope Boundary

This run explicitly will NOT:

- Write outside `<cwd>/research_outputs/` or `<cwd>/mission-control-runs/2026-05-27-0140-research-ecoli-pathotype-substrate/`
- Modify rules YAML, ADSP units, wiki/* files
- Modify `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md` (the project ledger is updated separately, not as part of this research run)
- Fabricate sources to hit row count targets
- Auto-promote candidates into rules/wiki (Promotion Gate remains manual)
- Spend more than 15% of daily token budget
- Run longer than 30 minutes wall-clock
- Make more than 100 tool calls

## Escalation Conditions

- Zero useful web search results after 3 attempts → halt + ask user to refine topic
- Intake rejects ≥80% of rows → halt + surface likely topic-shape failure
- Slug validation fails (exotic topic with no alphanumerics) → halt + ask user for manual slug
- Followup queue file write fails → halt; supported memo still preserved
- Any budget cap breached → halt + ask user to continue/stop/replan
- Unresolved-uncertainty count ≥5 → halt + surface confused state

## Run Metadata

- **Run ID:** 2026-05-27-0140-research-ecoli-pathotype-substrate
- **Started:** 2026-05-27T01:40Z
- **Autonomy Level:** L1
- **Department(s):** Research
- **Skill invoked:** /research
- **Caller:** user
- **Parent project:** ecoli-pathotype-prediction-cli-2026-05-26 (`project_state/ecoli-pathotype-prediction-cli-2026-05-26.md`)
- **Budget caps:**
  - Token: 15% daily
  - Wall-clock: 30 minutes
  - Tool-calls: 100
  - Unresolved-uncertainty: 5

## Contract version

v1.0
