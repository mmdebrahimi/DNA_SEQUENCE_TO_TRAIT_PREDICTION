<!-- audit-trail.md for 2026-06-02-1530-ep4-v0-marker-resolver -->
# Audit Trail — 2026-06-02-1530-ep4-v0-marker-resolver

Persona: Soraya `--advance` (money-only gate). Lock acquired; dynamic cap 8–10; 4 steps.

| # | action | gate | decision | outcome |
|---|---|---|---|---|
| 0 | de-risk dependency: fetch VirulenceFinder E. coli DB | auto | run | OK — 7.75 MB, 4942 alleles, 680 genes; all needed markers present |
| 1 | cache DB → `data/virulencefinder_db/virulence_ecoli.fsa` | auto | run | cached (gitignored data/) |
| 2 | write `scripts/pathotype_v0_marker_screen.py` | auto | run | created (pure-Python k=15 seed presence, both strands) |
| 3 | screen 24 genomes + confound-immune biology test | auto | run | `research_outputs/pathotype_v0_marker_screen_2026-06-02.json` |
| 4 | result/recommendation/audit + dogfood + emit `/project-state` | mixed | run local; EMIT | written; `/project-state` emitted for user |
| — | `git push origin main` | irreversible | HELD | per weekly-sync; not executed |

## Gate events
- No money. No deletes. No dep-install (DB is data; detection is pure-Python — deliberately avoided BLAST/KMA to keep zero-dep + laptop-local).
- Dependency de-risked BEFORE building (step 0) — confirmed the marker DB was fetchable + complete before writing the screen.

## Outcome
Strong positive: known virulence markers separate EPEC/ExPEC perfectly (eae AUROC 1.0) and interpretably → v0 resolver validated, confound-immune. No analyst override needed this run (the auto-verdict BIOLOGY_TRACKS_LABELS matched the substance: eae 12/0 + textbook ExPEC arsenal).
