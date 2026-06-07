<!-- intent-contract: 1.0 -->
# Intent Contract — eukaryotic-multimodal-substrate

- **run-id:** 2026-06-07-1500-research-eukaryotic-multimodal-substrate
- **timestamp:** 2026-06-07T15:00Z
- **level:** L1 · **departments:** Research · **skill:** /research · **caller:** user (via /soraya — Phase 5/6 entry)

## Verbatim Input
eukaryotic and multimodal genotype-to-phenotype substrate feasibility for a solo developer beyond bacterial AMR: rank candidate entry points by (a) public paired genome+phenotype dataset existing at depth (>=100 samples), (b) sampling-independent de-confoundable labels, (c) compute requirement (which DNA foundation model + GPU/VRAM, or whether a deterministic mechanism-feature approach like the shipped AMR decoder transfers), (d) whether a curated mechanism/determinant catalog already exists. Specific candidates: Candida auris/albicans azole resistance with WGS+MIC (fungal AMR, ERG11/TAC1); Arabidopsis 1001 Genomes flowering-time GWAS; human T2D polygenic risk (UK Biobank/All of Us); a multimodal bacterial-colony-image + WGS dataset. Goal: ranked shortlist of the single most feasible next substrate + its compute + data prerequisites.

## Decomposition
Standard /research v0.5 workflow: web search + synthesis → intake validation → followup queue update.

## Verification criteria
| Sub-task | Criterion | Evidence |
|---|---|---|
| Web research | ≥5 audit-grade rows OR honest-gap | `<slug>.raw.md` 13-col table |
| Intake | rows pass audit/mapping/banned/cite floors | `<slug>.md` + `_unsupported.md` |
| Followup | queue touched ≥ run-start | `_followup_queue.md` mtime |

## Out of scope
- No writes outside research_outputs/ + this run dir. No code/rules/wiki mods. No fabrication. No auto-promotion. **No compute provisioning (money gate).**

## Escalation
- 0 useful results ×3 → halt. Intake ≥80% reject → halt. Budget cap (15% tok / 30 min / 100 calls / 5 uncertainty) → halt.

## Budget caps: token 15% · wall-clock 30m · tool-calls 100 · uncertainty 5
