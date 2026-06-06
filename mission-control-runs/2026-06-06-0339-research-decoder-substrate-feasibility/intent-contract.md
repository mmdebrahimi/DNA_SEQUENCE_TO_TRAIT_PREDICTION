<!-- intent-contract: 1.0 -->
# Intent Contract — decoder-substrate-feasibility

- **run-id:** 2026-06-06-0339-research-decoder-substrate-feasibility
- **timestamp:** 2026-06-06T03:39Z
- **level:** L1 · **departments:** Research · **skill:** /research · **caller:** user (via /soraya rec #2)

## Verbatim Input
Candidate genotype-to-phenotype prediction targets beyond AMR and pathotype where a DNA-sequence model could add value: which bacterial (and adjacent) phenotypes have (a) sampling-INDEPENDENT labels — a lab measurement or assay, NOT a clinical-site/isolation-source category that causes study==class confounding — AND (b) public paired genome+phenotype cohorts large enough (>=100 strains) to build a de-confounded cohort where the two classes co-occur within lineages? For each candidate, note whether a curated knowledge/mechanism-feature baseline already exists; the open niche for a learned/embedding decoder is phenotypes with sampling-independent labels but NO curated catalog. Focus on E. coli first, then other bacteria. Goal: a ranked shortlist of the most feasible next decoder substrate.

## Decomposition
Standard /research v0.5 workflow: web search + synthesis → intake validation → followup queue update.

## Verification criteria
| Sub-task | Criterion | Evidence |
|---|---|---|
| Web research | ≥5 audit-grade rows OR honest-gap | `<slug>.raw.md` 13-col table |
| Intake | rows pass audit/mapping/banned/cite floors | `<slug>.md` + `_unsupported.md` |
| Followup | queue touched ≥ run-start | `_followup_queue.md` mtime |

## Out of scope
- No writes outside research_outputs/ + this run dir. No rules/wiki/code mods. No fabrication. No auto-promotion.

## Escalation
- 0 useful results ×3 → halt. Intake ≥80% reject → halt. Budget cap (15% tok / 30 min / 100 calls / 5 uncertainty) → halt.

## Budget caps: token 15% · wall-clock 30m · tool-calls 100 · uncertainty 5
