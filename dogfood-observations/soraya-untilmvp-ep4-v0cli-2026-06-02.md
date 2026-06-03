# Dogfood — Soraya `--until-mvp` on the v0 CLI build (2026-06-02, run 3)

Run: `soraya_runs/2026-06-02-1600-ep4-v0-cli-build/`. First `--until-mvp` flight (runs 1–2 were `--advance`). Goal: ship the ledger-locked v0 CLI to a checkable MVP bar.

## What worked
- **Freezing a checkable MVP bar up front kept the build honest.** Three predicates (cli.py exists / decision tests exit 0 / valid provenance JSON on a real genome) gave an unambiguous "done" — no scope creep, no vibes-based "looks finished." The `test-exit-0` + `file-exists` predicate kinds mapped cleanly onto real build artifacts.
- **Grounding in the locked contract first** (read the ledger v0 Output Contract before writing code) meant the 23 clusters / 11-class table / provenance schema / v5 reframings landed as specified, not reinvented.
- **The standalone-runnable test file** (pytest OR `python tests/...`) made the `test-exit-0` gate trivially satisfiable without depending on pytest being installed — good fit for the gated-runner allowlist.
- **No recovery rounds needed** — clean build, all predicates MET first pass. The dead-end detector / attempt budget never engaged (nothing to recover from).

## Friction / honest gaps
- 🟡 **The MVP bar didn't encode QUALITY, only existence+pass.** All 3 predicates were MET, but the ExPEC genome abstained to AMBIGUOUS (on-spec, but under-sensitive). A checkable bar can certify "it runs + tests pass" without certifying "it calls ExPEC well." That's correct for an MVP, but the loop has no notion of "MET but weak" — the quality judgment stayed with the model (flagged it as a v0.1 item rather than letting `mvp-reached` imply the science is done). Lesson: `--until-mvp` certifies the bar you wrote, nothing more — write bars that include a quality predicate when quality matters, or keep the model-level honesty pass explicit (as here).
- 🟡 Audit-trail appended at the end again (model-discipline), not per-step — consistent across all 3 runs this session. For `--until-mvp`'s potentially-long loops this is the most likely place for the OT1 caveat (skipped lifecycle steps) to bite; a code-owned per-step append would close it.

## Net (session, 3 runs)
`--advance` ×2 (confound discovery → confound-immune validation) then `--until-mvp` ×1 (ship the CLI) took the project from "learned bake-off on confounded data" to "tested, interpretable v0 resolver CLI that emits FASTA→pathotype+provenance." Soraya's safety pieces (gate, lock, de-risk-first, MVP predicates) all behaved. The two durable persona lessons: (1) auto-metric/auto-MVP verdicts certify the literal predicate — keep a model-level substance/quality check explicit; (2) per-step code-owned audit append would harden the OT1 seam.
