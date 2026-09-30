# Audit trail — 2026-09-30-1625-advance-glm

Mode `--advance`, spillover ON, money-only gate posture. Pre-planned estimate ~20 steps; no extension
needed. Every step classified before execution; nothing hit money / destructive-local /
genuinely-irreversible-outward.

## Step 4 — code-owned eligibility gate
`advance_ranker.rank` -> 5 eligible. Coarse order: label-acquisition(med,d0) > catalog-curation(med,d1) >
doubt-layer / evidence-surface / learned-narrow (low). Ranker emitted replay-skip warnings on FOUR of the
five families.

## Steps executed

1. **Adjudicated the ranker's advisory warnings against the real record** (not auto-skipped):
   - `label-acquisition` "Re-run the prospective accrual sweep" (100% term match) — WARNING CORRECT: the
     v2-clock accrual landed 2026-09-28 (Campylobacter x cipro, N=426). Re-running two days later is
     low-value given PD ingestion lag. Left in place, not re-run.
   - `catalog-curation` row 1 — **USER AUTHORITY** (restrict the gentamicin rescue to E. coli; edits the
     frozen surface + invalidates the v2 lock). SURFACED, never self-executed.
   - `evidence-surface` #1 and `learned-narrow` #2 — JOURNALED done, still pending.

2. **LEDGER HYGIENE (edit-local-code).** Verified each `[done:<id>]` marker really was in the journal
   BEFORE touching anything, then retired the rows NON-DESTRUCTIVELY (re-appended verbatim under a
   `### Retired candidates` heading placed AFTER the end-marker, so the ranker cannot see them but the
   content survives). `learned-narrow` needed TWO retired, not one: the row beneath the flagged one was a
   correction NOTE, not an action, so the ranker would have picked that instead. Added 2 genuine
   candidates. Both mirrored via `/project-state --append-action`; all 6 end-markers still unique on both
   ledgers. EFFECT ON THE GATE: `learned-narrow` moved bucket LOW -> **HIGH** with a real next action,
   which changed what `--advance` picked.

3. **Executed the picked family's next action (run-tests).** `scripts/glm_headroom_allrt.py` —
   generalized the GLM headroom result off ONE drug across all 11 RT drugs x 2 dataset variants, judged by
   the SAME bar (imported, not restated) against the committed zero-shot all-RT baseline.
   **Result: 7 of 8 genuine drugs PASS under leave-one-STUDY-out vs 0 of 8 zero-shot; median supervised
   0.8141 vs 0.5121. 6 of 6 per-drug anchors bind.**

4. **Caught the documented nominal-vs-genuine trap mid-run.** First print read "zero-shot PASSES 1/6";
   that pass is DOR (>=0.65 on a 37-isolate subset, burden 0.783) — the exact nominal-not-genuine case the
   committed artifact already flags. Applied the same split to both arms + named every excluded drug with
   its n; zero-shot -> 0/5, independently reproducing the committed `genuine: 0`.

5. Artifact + memo + 10 tests + 4 LESSONS_LEARNED entries + 2 ledger journal rows.

## Verification
65 targeted tests pass (both GLM suites + regime consumers + complement suites). Frozen AMR surface
byte-unchanged. Full suite launched.

## Gates
None triggered. One authority fork surfaced and left for the user (catalog-curation row 1).
