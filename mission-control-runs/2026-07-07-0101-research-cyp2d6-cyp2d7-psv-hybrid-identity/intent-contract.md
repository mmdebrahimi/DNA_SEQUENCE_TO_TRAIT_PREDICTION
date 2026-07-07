<!-- intent-contract: 1.0 -->
# Intent Contract — cyp2d6-cyp2d7-psv-hybrid-identity

- **run-id:** 2026-07-07-0101-research-cyp2d6-cyp2d7-psv-hybrid-identity
- **timestamp:** 2026-07-07T01:01:00Z
- **level:** L1  · **departments:** Research · **skill:** /research · **caller:** user

## Verbatim Input
`fresh research arc (PSV curation + read-level pileup)  hybrid identity (which of *13/*36/*68 via CYP2D6-vs-CYP2D7 PSV analysis)`

## Decomposition
Standard /research v0.5 workflow: web search + synthesis → intake validation → followup queue update.
Goal: curate the CYP2D6-vs-CYP2D7 paralogous-sequence-variant (PSV) positions + the read-level pileup /
copy-number method the field's tools (Cyrius / Aldy / StellarPGx / PharmVar) use to resolve CYP2D6-CYP2D7
HYBRID allele IDENTITY (*13 / *36 / *68 / *4N / *61 / *63) — the input needed to lift dna_decode's shipped
hybrid-PRESENCE detector (sens 0.62 / spec 1.0) to hybrid IDENTITY.

Sub-tasks: (1) web research → `<slug>.raw.md`; (2) intake validation → `<slug>.md` + `_unsupported.md`;
(3) followup queue update.

## Verification criteria
| Sub-task | Criterion | Evidence |
|---|---|---|
| Web research | ≥5 audit-grade rows OR honest-gap declared | `<slug>.raw.md` V1 13-col table |
| Intake validation | audit floor + mapping + banned-phrase + cite-token scans | `<slug>.md` + `<slug>_unsupported.md` |
| Followup queue | `_followup_queue.md` touched ≥ run-start | mtime |

## Out of scope
- No writes outside `research_outputs/` or this run dir; no rules/wiki/code edits; no fabricated sources;
  no auto-promotion; budget ≤15% token / ≤30 min / ≤100 tool-calls / ≤5 uncertainty.

## Escalation conditions
- 0 useful results after 3 attempts · intake rejects ≥80% · slug-validation fail · queue-write fail ·
  any budget cap breached · uncertainty ≥5.

## Budget caps
token 15% · wall-clock 30m · tool-calls 100 · uncertainty 5.
