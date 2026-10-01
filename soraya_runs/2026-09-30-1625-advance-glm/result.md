# Result — 2026-09-30-1625-advance-glm

**Terminal reason:** authority fork surfaced + the picked family's next action executed to completion.
**Budget:** ~20 planned, ~20 executed, no extension. **Spillover:** not needed — the picked family's own
next action was the highest-VOI reversible move after hygiene.

## Executed
1. Ledger hygiene across 2 families (3 rows retired non-destructively, 2 genuine candidates added,
   both mirrored). Moved `learned-narrow` LOW -> HIGH, which changed what the gate picked.
2. `scripts/glm_headroom_allrt.py` — the picked family's next action. **7 of 8 genuine RT drugs pass
   leave-one-STUDY-out vs 0 of 8 zero-shot; median supervised 0.8141 vs 0.5121; 6/6 anchors bind.**
3. Caught and fixed the nominal-vs-genuine pass-count trap mid-run (DOR, n=37, burden 0.783).

## Verification
65 targeted tests pass; frozen AMR surface byte-unchanged; full suite launched at hand-back.

## Gates
None triggered. One authority fork left for the user — see `recommendation.md`.

## No resume
This run is closed. The next move is in `recommendation.md`; it needs no state from this run.

## Spillover (continued past the first family's action)

**Terminal reason (FINAL):** genuine plateau + an authority fork. Step-8 checklist items (2)+(3) hold.
**Steps:** ~20 planned, ~85 executed across three self-extensions (all announced in the audit trail).
**Pivots:** 1 (into the HIV deployability provenance defect, reached from `recommendation.md`'s next move).

4. Measured leave-one-study-out blind-spot deployability **per gene**: NNRTI 0.8102 (anchor reproduces the
   committed artifact exactly), **INSTI 0.8923 — the second gene, passing the same bar**, PI **withheld**.
   This closes the one-GENE limit on the morning's GLM result and replaces two unsourced shipped numbers.
5. Fixed the generator so a rebuild cannot reintroduce a literal; corrected four restatements.
6. Fixed one of my own guards that the commit itself broke, caught by the full suite.

## Commits
`b67ee2f` (all-RT generalization) · `cd469d0` (deployability correction) · `16966c8` (guard baseline fix).

## Verification
**5136 passed / 11 skipped / 0 failed** on the final committed state.

## Next (no state from this run needed)
`learned-narrow` row 5 — **does CONTEXT beat the linear floor?** Both its prereqs are now met (7 of 8 RT
drugs; 3 genes). 0.8102/0.8923 come from logistic regression over one-hot tokens, the weakest supervised
family member, so they are a FLOOR and say nothing about attention or epistasis. Same data, same bar, same
leave-study-out split — swap the estimator. `scripts/hiv_epistasis.py` already exists for the pairwise term.

## Spillover round 2 + final terminal

7. **Answered candidate row 5** — does context beat the linear floor? **No: it IS the ceiling.** All 4
   context arms lose (2 significantly), against a bar frozen and committed before the run, with the linear
   arm reproducing the floor to +0.0000. Corroborates the independent 2026-07-11 epistasis negative in a
   different framing.
8. Updated `eval/regime.py` + `CLAUDE.md` (both still called the linear model a floor).

## Commits (5)
`b67ee2f` all-RT generalization · `cd469d0` deployability correction · `16966c8` guard baseline fix ·
`2cf723d` frozen bar · `ebbc781` context-vs-floor · `2a46b30` orientation update.

## Verification (final)
**5147 passed / 11 skipped / 0 failed.**

## What needs the user (nothing else is blocked on me)
1. **AUTHORITY FORK** — restrict the gentamicin `rmt` rescue to E. coli? Edits the frozen surface and
   invalidates the v2 lock. Drafted recommendation in `recommendation.md`: **do not restrict yet.**
2. **USER-GATED** — the one-VIRUS limit needs a new LABEL, and the standing record puts label acquisition
   with the user (public-label expansion is a recorded closure). Not an executor task.
