# CYP2D6 hybrid-identity plan — brainstorm review (2026-07-07)

Adversarial review of the plan to lift the shipped hybrid-PRESENCE detector to hybrid-IDENTITY (*13/*36/*68).
All findings accepted + verified; the plan was sharpened on 3 fronts before any code.

## Critical issues (fixed in the technical plan)
1. **The per-PSV base-count is the real build — and D6-coordinate-only is wrong.** No mpileup/base-count
   contract existed (only coarse `samtools depth`). MUST count BOTH `pos_CYP2D6` AND `pos_CYP2D7` per PSV
   (D6-only measures aligner placement, not biology). Fix: build a read-only PSV evidence table first that
   reproduces sane profiles before any identity call. `[grounded]` → **DONE in Phase A (falsifier GO).**
2. **PSV source misidentified.** No `CYP2D6.json` in the Cyrius repo (verified) — it ships
   `data/CYP2D6_SNP_38.txt`, GRCh38-native (117 PSVs), paired-coord schema. "Curate from CYP2D6.json + lift
   to GRCh38" is wrong/needless. Fix: consume `CYP2D6_SNP_38.txt` directly, no liftover. `[grounded-verified]`
3. **≥90%-on-3-anchors bar is not statistically meaningful.** Committed set = 13 hybrid rows (*68=4, *36=8,
   *13=1 [NA19785 only]); anchors NA24631/HG01161 aren't in the data. Fix: per-allele honest denominators;
   *13 single-sample-UNPOWERED; anchors are additional fetch targets. `[grounded-verified]`

## Medium issues (in the plan)
- Anchor-reproduction is necessary-but-not-sufficient → add a cheap decisive falsifier on a tiny mixed panel
  (normal / pure-dup / *5 / *68 / *36 / NA19785 *13). **DONE — falsifier GO.**
- Identity needs copy-normalized / read-linked evidence, not per-site fractions alone (*36x2/mixed → over-call).
  Default to `hybrid_present_identity_unresolved`.
- Phase B outputs THREE levels: resolved / identity-unresolved / evidence-not-callable. Never force a call.

## Deferred
- *36 embedded-exon-9 gene-conversion blind spot (no CN change; needs read-level fraction) — documented residual.
- *36x2 / multi-copy explicit CN modelling — until read-level linkage is strong.

## What's solid
- The PSV-based direction (Cyrius-validated, 96.5%/99.3% ceiling) — not a tag-SNP problem.
- Reusing the proven remote-CRAM tooling; the abstaining + honest-residual discipline; the benchmark anchors.

## Open tradeoffs
- mpileup flag contract (permissive vs MAPQ/baseQ-filtered) — derive empirically on the panel (Phase A used
  permissive `-B -q0 -Q0`; the falsifier passed).
- Aligner-portability — scope to GRCh38-aligned CRAMs with the benchmark mapper, or invest in robustness.

## Recommended next step (executed)
Phase A only, gated by the falsifier — **DONE + GO.** Next: draft Phase B (this folder's technical-plan.md),
then ratify before building the classifier.
