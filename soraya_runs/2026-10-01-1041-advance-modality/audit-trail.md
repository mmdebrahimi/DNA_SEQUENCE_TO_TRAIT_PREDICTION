# Audit trail — 2026-10-01-1041-advance-modality

Mode `--advance`, spillover ON, money-only. Pre-planned ~35 steps; ~22 executed (the probe answered the
question earlier than a build would have). Every step classified before execution; no gate triggered.

## Step 4 — code-owned eligibility gate
5 eligible, all med/low: label-acquisition(med) · learned-narrow(med) · catalog-curation(med, user_only)
· doubt-layer(low) · evidence-surface(low). Picked **learned-narrow** — T2 of the gLM plan belongs to it.

## Steps executed

1. **Answered two mid-run user questions** (one-hot; orthogonality) grounded in this repo's own measured
   numbers rather than abstractly — the ESM2+GEMME case (a WEAKER model improving a stronger one) vs the
   ESM2+AlphaMissense non-lift is the orthogonality demonstration, and both bear directly on step 2.

2. **Did NOT build the channel first — probed for room** (`hiv_onehot_unseen_position_headroom.py`,
   run-tests). The verdict rule was frozen in-file before the first run and `self_check` drives all three
   verdicts, so `HEADROOM_EXISTS` was reachable. Floor reproduces the committed 0.8102 / 0.8923 **exactly**
   (a miss exits 3, because strata measured against a non-reproducing floor are uninterpretable).

   - **Blindness CONFIRMED and measured:** max |coef| over every zero-training-signal column = **exactly
     0.0**, both genes. The mechanism is as hypothesised.
   - **Room ABSENT:** 7.5% of RT blind-spot isolates exposed (bar 20%), and the model scores **better**
     where blind — 0.8229 vs 0.8163, gap **−0.0066, wrong direction**.
   - **Verdict `NO_HEADROOM_ONEHOT_ALREADY_COVERS_IT`.** INSTI exposed stratum = 8 R → AUROC **withheld**,
     cell reads `UNDERPOWERED`.

3. **Measured the CAUSE rather than asserting it** (this family has had three asserted causes refuted).
   Per-isolate redundancy: 12.44 active columns, median blind fraction **1 in 12**. Honest exception kept
   visible: RT max blind fraction **1.0** — one isolate entirely blind (n=1).

4. **Reconciled the two modality results**: the forward cell's +0.056/90.5% lift was on SINGLE variants;
   this blind spot is multi-substitution isolates with redundant signal. Consistent with (not proven by)
   the 2026-07-03 single-mutant head at 0.691.

5. **Recorded the closure on all four standing surfaces** (edit-local-code) — candidate row 8, the hybrid
   plan's Phase 5, `eval/regime.py` + its regenerated map artifact, and `CLAUDE.md` — so the lever is not
   re-proposed. 4 LESSONS_LEARNED entries.

## Verification
10 new tests; 60 pass across regime consumers + all three HIV gLM suites; 93 pass on doc-sync +
contract-number audit after the CLAUDE.md edit. Frozen AMR surface byte-unchanged; prospective lock v2
re-verified `ok=True`. Full suite launched at hand-back.

## Gates
None triggered. The two user-authority items (gentamicin fork; label/arm selection) were NOT touched.
