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
