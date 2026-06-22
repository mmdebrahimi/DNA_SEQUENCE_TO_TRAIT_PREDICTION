<!-- intent-contract: 1.0 -->
# Intent Contract — noncircular-label-sources (IMMUTABLE)

- **Run ID:** 2026-06-22-1000-research-noncircular-label-sources
- **Level:** L1 · **Departments:** Research · **Skill:** /research · **Caller:** /soraya (label-acquisition epoch)
- **Timestamp:** 2026-06-22

## Verbatim Input
Non-circular, sampling-independent phenotype label sources for bacterial (especially E. coli /
Enterobacterales) genotype-to-phenotype decoding. For each candidate: (a) wet-lab/clinical MEASUREMENT vs
genome-tool-derived? (b) ACCESS PATH (free download / public API / registration / MTA / paid)? (c) rough
SCALE — same-organism isolates with the label AND downloadable assemblies? (d) outside the
NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA surveillance ecosystem? Verify+rank: NARMS raw MIC, EUCAST MIC
distributions, Pfizer ATLAS / Paratek KEYSTONE / Venatorx, SENTRY, CRyPTIC TB MIC, EBI EGA / dbGaP,
Whittam DECA, von Mentzer 2021 ETEC, BacDive, growth/fitness assays. Ranked shortlist by acquirability
(free first), access path + realistic N-after-assembly-filter each.

## Decomposition
Standard /research v0.3: web search + synthesis -> intake validation -> followup queue update.
1. Web research (Steps 1-3 -> <slug>.raw.md)
2. Intake validation (Step 4 -> <slug>.md + <slug>_unsupported.md)
3. Followup queue update (Step 5)

## Verification criteria
| Sub-task | Criterion | Evidence |
|---|---|---|
| Web research | >=5 audit rows OR honest gap | <slug>.raw.md with 13-col table |
| Intake | rows pass floors | <slug>.md + <slug>_unsupported.md |
| Followup | queue touched this run | _followup_queue.md mtime >= run start |

## Out of scope
- No writes outside research_outputs/ + this run dir. No rules/wiki/code mods. No fabrication to hit row
  counts. No auto-promotion. Money spend = hard gate (paid sources flagged, not purchased).

## Escalation
- 0 useful results x3 -> halt; intake rejects >=80% -> halt; budget cap breach -> halt. Caps: token 15%
  daily / wall-clock 30m / tool-calls 100 / uncertainty 5.
