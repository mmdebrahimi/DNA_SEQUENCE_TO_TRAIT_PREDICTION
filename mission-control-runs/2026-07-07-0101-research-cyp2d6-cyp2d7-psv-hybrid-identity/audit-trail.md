<!-- audit-trail: 1.0 -->
# Audit Trail — cyp2d6-cyp2d7-psv-hybrid-identity

- **Verdict:** COMPLETED
- **run-id:** 2026-07-07-0101-research-cyp2d6-cyp2d7-psv-hybrid-identity
- **level:** L1 · **departments:** Research · **skill:** /research · **caller:** user
- **summary:** Curated the CYP2D6-vs-CYP2D7 PSV method (Cyrius/DRAGEN/PharmVar) for resolving hybrid
  IDENTITY (*13/*36/*68) — the input to lift dna_decode's shipped hybrid-PRESENCE detector to identity.

## Departments invoked
| Sub-task | Duration | Outcome |
|---|---|---|
| Web research | ~6 tool calls | PASS — 16 audit-grade rows, 3 honest gaps flagged |
| Intake validation (by hand — /research-intake not model-invocable) | inline | PASS — 16 supported / 0 unsupported (all pass audit floor) |
| Followup queue update | inline | PASS — +5 candidates appended |

## Skills called
1. WebSearch ×2 (Cyrius/PSV; hybrid *13/*36/*68 detection)
2. WebFetch ×4 (Cyrius PMC7997805; DRAGEN CYP2D6 caller; Gaedigk PMC6556886; Aldy genome.cshlp — redirect, skipped)
3. /research-intake protocol executed by hand (audit-floor + banned-phrase + cite-token scans)
4. /research-followup protocol executed by hand (queue append)

## Verification results
| Sub-task | Result | Evidence |
|---|---|---|
| Web research (≥5 rows OR honest gap) | PASS | 16 rows in `.raw.md` |
| Intake (audit floor + scans) | PASS | `.md` (16 supported) + `_unsupported.md` (0) |
| Followup queue touched ≥ run-start | PASS | `_followup_queue.md` appended 2026-07-07 |

## Budget
- token: within cap · wall-clock: ~10 min (cap 30) · tool-calls: ~10 (cap 100) · uncertainty: 3 (cap 5, honest gaps)

## Escalations
- none

## Adversarial review
- none (this skill does not use /brainstorm internally)

## Result path
- research_outputs/cyp2d6-cyp2d7-psv-hybrid-identity-2026-07-07.md

## Artifacts
- research_outputs/cyp2d6-cyp2d7-psv-hybrid-identity-2026-07-07.raw.md
- research_outputs/cyp2d6-cyp2d7-psv-hybrid-identity-2026-07-07.md
- research_outputs/cyp2d6-cyp2d7-psv-hybrid-identity-2026-07-07_unsupported.md
- research_outputs/_followup_queue.md (+5)
- mission-control-runs/2026-07-07-0101-research-cyp2d6-cyp2d7-psv-hybrid-identity/intent-contract.md
- mission-control-runs/2026-07-07-0101-research-cyp2d6-cyp2d7-psv-hybrid-identity/audit-trail.md

## Lessons
- The 117 Cyrius PSV coordinates are NOT in the literature body — they live in the Cyrius GitHub `CYP2D6.json`
  config. PSV *curation* is a repo-config lift, not a fetch; that is the true first build step.
- Hybrid identity is (breakpoint × direction × CN-profile) = an algorithm, not a single-SNP lookup —
  confirms the v0.3 scoping that this is genuinely Cyrius-class, not a tag-SNP cell.
