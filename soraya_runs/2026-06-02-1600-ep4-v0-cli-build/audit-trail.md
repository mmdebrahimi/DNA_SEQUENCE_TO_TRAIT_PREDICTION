<!-- audit-trail.md for 2026-06-02-1600-ep4-v0-cli-build -->
# Audit Trail — 2026-06-02-1600-ep4-v0-cli-build

Persona: Soraya `--until-mvp`. Lock acquired; bounded by MVP bar (not action cap).

| step | action | gate | outcome |
|---|---|---|---|
| 1 | read ledger v0 Output Contract (23 clusters, 11-class table, provenance schema, v5 reframings) | auto | grounded the build in the locked spec |
| 2 | write `dna_decode/pathotype/{markers,resolve,detect,cli,__init__,__main__}.py` | auto | package created |
| 3 | write `tests/test_pathotype_resolver.py` (18 tests) | auto | created |
| 4 | run decision-table tests | auto (test) | **18/18 PASS** — MVP predicate 2 MET |
| 5 | smoke CLI on EPEC AIEY01 | auto | tEPEC_COMPATIBLE + provenance JSON — predicate 1+3 MET |
| 6 | smoke CLI on ExPEC JSIS01 | auto | AMBIGUOUS (on-spec) + provenance JSON |
| 7 | validate provenance JSON schema completeness | auto | both schema-complete (7 keys + full hit fields + db sha256) |
| 8 | write result/recommendation/audit + dogfood; emit `/project-state`; commit + push | mixed | run; push per user intent this turn |

## MVP loop
3 frozen predicates. Evaluated at stop: all 3 live-MET → `resolve_stop` = `mvp-reached`. No recovery rounds needed (clean build). No money/deletes/dep-installs.

## Honesty note (no override needed, but flagged)
ExPEC JSIS01 → AMBIGUOUS is on-spec (conservative UPEC ≥2-strong rule), not a failure. Recorded as a v0.1 calibration item rather than silently "fixed", to preserve the contract's abstention-first design and let ledger H4 drive the threshold decision with data.
