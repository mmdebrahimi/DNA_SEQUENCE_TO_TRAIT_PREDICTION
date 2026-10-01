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

## SPILLOVER (step 8, default-on) — continued past the picked family's action

6. **Re-aimed to the next high-VOI reversible move** named in `recommendation.md`: the one-GENE limit. The
   per-class non-vacuity artifact already reported a `leave_study_out_blindspot_auroc` for INSTI on gene
   IN, so before building anything I checked whether that number was MEASURED or CITED.
   **It was a hardcoded literal.** `build_hiv_complement_model.py` carried `leave_study_out` as config
   literals and called all three "deployment-validated"; `hiv_supervised_deployability.py` hardcodes
   `DRUG="EFV"` + the NNRTI dataset, so only RT was ever measured.

7. **Built + ran `scripts/hiv_deployability_per_gene.py`** (run-tests). RT anchor reproduces the committed
   artifact EXACTLY (same AUROC/n/R; deployed comparator agrees 4222/4222). **NNRTI 0.8102 · INSTI 0.8923
   (2nd gene, passes the same bar, 69 studies) · PI WITHHELD** (its catalog calls 74.2% positive, leaving
   3 R of 910 — unmeasurable, not merely unsourced).

8. **Closed the root cause, not the values** (edit-local-code): `deployability_block()` now reads the
   measured artifact and RAISES if it is absent or if the anchor missed, so a rebuild cannot restore a
   literal. The correction then propagated on its own into the derived non-vacuity artifact. Corrected the
   builder docstring, a `doubt.py` comment and two memos that restated the figure.

9. **Ledger hygiene round 2**: retired row 4 (its work is complete + journaled), updated row 5 (both
   prereqs now met), added row 6 (the one-VIRUS limit — a LABEL question, not an estimator one).

10. **A guard of my own broke on commit** and the full suite caught it: the weights-unchanged test compared
    against `HEAD`, which became the corrected state once committed. Re-pinned to the pre-correction
    commit + a non-vacuity assertion on that baseline; proven both directions.

## Verification (final, post-commit)
Full suite on the committed fixed state: **5136 passed, 11 skipped, 0 failed (exit 0)**. The preceding run
on the pre-fix state showed exactly ONE failure — the guard above, and nothing else. Frozen AMR surface
byte-unchanged throughout; prospective lock v2 re-verified `ok=True`; every CLI call unchanged.

## SPILLOVER round 2 — candidate row 5 answered

11. **Froze and COMMITTED the acceptance bar first** (`2cf723d`,
    `wiki/hiv_context_vs_linear_floor_acceptance_bar.json`), with the prior stated up front (the
    2026-07-11 epistasis negative) and an explicit argument for why this is not a duplicate of it.
12. **Ran it** (run-tests): all 4 context arms LOSE, 2 with a CI entirely below zero; anchor reproduces
    the floor to +0.0000. Verdict by the frozen rule: `LINEAR_FLOOR_IS_THE_CEILING`.
13. Named the limitation that decides which arm carries the claim (pairs chosen by PREVALENCE are NRTI
    TAMs in an NNRTI cohort, so the general claim rests on the non-linear arm) + a test that fails if it
    is dropped. Updated `eval/regime.py`, whose caveat had said this "says nothing about context".
14. Retired row 5 (answered), added row 7 recording the closure + the one remaining form of the question.
15. Updated `CLAUDE.md` — the auto-loaded orientation still called the linear model a floor and the result
    one-drug/one-gene.

## Terminal (step-8 checklist, assessed not assumed)
Checklist items **(2) + (3)** hold. Every remaining move is one of: an **authority fork**
(catalog-curation, `user_only=True`), a **user-gated label decision** (learned-narrow row 6 — the standing
record puts label ACQUISITION with the user, and public-label expansion is a recorded closure), a move
**structurally blocked by its own named reason** (doubt-layer: TN-starved/no free source; evidence-surface:
unverifiable BY measurement), or **low-value re-confirmation** (label-acquisition's sweep fired 2 days ago;
PD ingestion lag). ~85 of the ~100-step self-budget spent. Banking rather than manufacturing another
artifact — that is the motion-not-signal trap this checklist exists to catch.

## Verification (final)
**5147 passed / 11 skipped / 0 failed (exit 0)** on the final committed state. Frozen AMR surface
byte-unchanged across all five commits; prospective lock v2 `ok=True`; every CLI call unchanged.
