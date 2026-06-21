<!-- intent-contract: 1.0 -->
# Intent Contract — viral-antiviral-resistance-gp-datasets

- **Run ID:** 2026-06-21-1314-research-viral-antiviral-resistance-gp-datasets
- **Timestamp:** 2026-06-21T13:14Z
- **Level:** L1
- **Departments:** Research
- **Skill:** /research
- **Caller:** user (via /soraya Wave B de-risk)

## Verbatim Input
free, redistributable, isolate-level laboratory genotype–phenotype susceptibility datasets for antiviral resistance interpretation (HIV reverse-transcriptase/protease/integrase, HCV, SARS-CoV-2 protease), screened against the project's 8 rejection gates + license terms

## Decomposition
Standard /research v0.5 workflow: web search + synthesis → intake validation → followup queue update. Purpose: GO/NO-GO on whether a FREE, de-confounded viral genotype↔phenotype substrate exists to validate the wired-but-unvalidated viral determinant decoder (the "virus" half of the bacteria/virus phenotype→trait tool).

## Verification criteria
| Sub-task | Criterion | Evidence |
|---|---|---|
| Web research | ≥5 audit rows OR honest-gap declared | raw.md with V1 table |
| Intake validation | rows pass floors + scans | <slug>.md + _unsupported.md |
| Followup queue | _followup_queue.md touched | mtime ≥ run-start |

## Out of scope
- No writes outside research_outputs/ + this run dir
- No fabrication to hit row counts
- No auto-promotion into rules/wiki/code

## Escalation conditions
- Zero useful results after 3 attempts → halt
- Budget cap breached → halt
- (OBSERVED) WebSearch usage-policy filter blocks broad multi-pathogen queries → degrade to narrow domain-restricted queries + honest-gap the blocked sub-topics; do NOT hammer

## Budget caps (L1)
token ≤15% daily / wall-clock ≤30m / tool-calls ≤100 / uncertainty ≤5
