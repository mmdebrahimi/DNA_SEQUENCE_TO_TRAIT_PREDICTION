<!-- audit-trail.md for 2026-06-03-0000-ep4-v0-cohort-calibrate -->
# Audit Trail — 2026-06-03-0000-ep4-v0-cohort-calibrate (UNATTENDED)

Persona: Soraya `--advance`, overnight autonomous (user asleep). Lock acquired. Disk checked (44 GB free) before any heavy op.

| # | action | gate | outcome |
|---|---|---|---|
| 0 | check disk + VF tooling feasibility | auto | 44 GB free; no BLAST/KMA/conda/pip-aligner -> VF diff not unattended-installable |
| 1 | write `scripts/pathotype_v0_cohort_eval.py` (coverage cache) | auto | created |
| 2 | run cohort eval (24 genomes) — the long run | auto | ~5 min; baseline H4: EPEC recall 1.0, ExPEC 0.25, abstain 0.58, conf-precision 1.0 |
| 3 | preserve baseline JSON; inspect ExPEC profiles | auto | identified 6 genomes = 1 strong + 2 support -> AMBIGUOUS |
| 4 | calibrate resolve.py (RULE_EXPEC_001 LOW_CONFIDENCE) + markers + 2 tests | auto | 20/20 tests pass |
| 5 | re-run cohort eval from cache (instant) | auto | calibrated H4: ExPEC recall 0.75, abstain 0.08, conf-precision 1.0 (both targets MET) |
| 6 | regenerate JSIS demo with calibration | auto | ExPEC_COMPATIBLE [LOW_CONFIDENCE] |
| 7 | write VF-feasibility doc (rec 3 gated) | auto | research_outputs/pathotype_vf_sidebyside_feasibility_2026-06-03.md |
| 8 | result/recommendation/audit + commit + push + handoff + retrospective | mixed | run; push per user mandate |

## Gate events (unattended discipline)
- Money: none. Deletes: none.
- **Dep-install (VF aligner) DECLINED** — judgment widened the gate past the money-only default; an unattended BLAST/KMA install could break the env and the user can't approve. Documented + handed to user. Correct fail-safe direction.
- No analyst override needed; calibration auto-verdict (H4 targets met) matched the substance.
