<!-- audit-trail.md for 2026-06-02-1500-ep4-interpret-confound -->
# Audit Trail — 2026-06-02-1500-ep4-interpret-confound

Persona: Soraya `--advance` (money-only gate). Lock acquired via `run_dir.acquire_lock`. Dynamic cap 8–10; 4 steps planned/executed.

| # | action | gate verdict | gate decision | outcome |
|---|---|---|---|---|
| 0 | dogfood `action_gate.classify_action` on planned batch | — | — | commit=auto, push=irreversible, script run/write=auto, `rm -rf`=irreversible(un-gated under money-only) |
| 1 | `git commit` EP-4 bake-off work (3 scripts + 3 result JSONs + handoff) | auto | run | commit `2d878bb` on `main` (local; not pushed) |
| 2 | write `scripts/pathotype_model_interpret_confound.py` | auto | run | created |
| 3 | run interpretation/confound probe | auto | run | `research_outputs/pathotype_model_interpret_confound_2026-06-02.json` |
| — | `git push origin main` | irreversible | **HELD** (judgment widened gate) | emitted as recommendation, not executed (user weekly-sync preference + unanswered commit-vs-push Q) |
| 4 | record finding + analyst override | mixed | run local; EMIT `/project-state` | result/recommendation/dogfood written; `/project-state` emitted for user |

## Gate events
- No `money` verdicts encountered. No deletes. No dep-installs.
- One HELD irreversible (push) — conservative override, fail-safe direction.
- Fail-safe never triggered (no auto-classified action suspected of spending money).

## Key execution note
Step 3's script emitted an auto-verdict (`SPARSE → hopeful`) from concentration math ALONE. Soraya applied an analyst override after inspecting k-mer IDENTITY (all top-10 are rare CTAG-motif/restriction-site 8-mers with partial presence → batch-artifact-suspected, not biology). Override recorded in `result.md`. This gap (verdict logic ignores feature identity) is logged as a dogfood finding.
