# Recommendation — 2026-09-30-1625-advance-glm

## The one thing needing the user

**`catalog-curation` row 1 is an AUTHORITY fork, not a technical question.** Restricting the gentamicin
`rmt` rescue to E. coli edits the FROZEN surface and therefore **invalidates the 2026-08-31 v2 prospective
lock** — the same class of decision as deploying the rescue was. Surfaced, not executed.

Drafted position (ratify or redirect): **do not restrict it yet.** The measured Klebsiella over-call
(PPV 0.475) fails this repo's own source-diversity bar — 98.4% of its susceptible carriers come from one
study — and NCBI-PD says 53R/0S across 12 BioProjects, which CLEARS that bar. Under contradictory evidence
the L2 `organism_scope` warning already ships the caveat at call time, which is the cheap half of the
protection; spending a lock retirement on the expensive half is premature. Revisit if a second
source-diverse Klebsiella cohort lands.

## Next high-VOI reversible move (no authority needed)

`learned-narrow` row 5: **does CONTEXT beat the linear floor?** Today's 0.8141 comes from logistic
regression over one-hot tokens — the weakest supervised family member — so it is a floor, not a ceiling,
and says nothing about attention. `scripts/glm_alphabet_headroom.py` + `glm_headroom_allrt.py` are the
harness; it is an estimator swap, not a new pipeline.

Its two prereqs are now unevenly met: **"more than one drug" is DONE** (7 of 8 RT drugs), but **"more than
one gene" is NOT** — every drug here shares RT. A second gene is the cheaper and more informative of the
two remaining prereqs, and the per-class non-vacuity work already says which genes are worth it: INSTI
(integrase) has a real blind spot; **PI (protease) does not** (2 resistant isolates in 614), so protease
would be a wasted run.

## Not worth re-running
`label-acquisition` "re-run the prospective accrual sweep" — fired 2026-09-28, N=426 Campylobacter cipro.
PD ingestion lag means two days buys nothing.
