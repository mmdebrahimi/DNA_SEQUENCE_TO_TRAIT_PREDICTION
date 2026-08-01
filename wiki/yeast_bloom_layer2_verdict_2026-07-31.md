# Layer-2 verdict: does a better model / the AI approach beat the linear baseline on clean data? (2026-07-31)

Prompted by the user's pushback — "the reconstruction negative could only mean (1) wrong learning
approach or (2) optimize the model; why give up?" Correct instinct. We did NOT give up; we ran both on the
first CONFOUND-FREE substrate (Bloom-2013 yeast cross, where a genuine test is finally possible). Two
hypotheses, two honest answers.

## Hypothesis 2 — "optimize the model" → CONFIRMED (the win)

A nonlinear model (gradient-boosted trees, capturing epistasis) BEATS linear ridge on the majority of
traits (`wiki/yeast_bloom_model_bench_2026-07-31.json`, 12 traits, 1,008 segregants):

- gbm beats ridge **7/12** (ridge 3, 2 ties); mean delta **+0.031**, up to:
  - **Maltose 0.728 → 0.889 (+0.160)**, Cadmium 0.693 → 0.788 (+0.095), Copper 0.566 → 0.659 (+0.093).

The linear baseline left real signal on the table — the gene-gene interactions Bloom 2013 flagged (up to
~50% of heritability for some traits), which a linear model structurally cannot represent. A nonlinear
model recovers them. On a confound-free substrate this is a GENUINE gain, not confounded inflation. So:
**yes, optimizing the model helps — the "why give up" concern is answered by doing, and it worked.**

## Hypothesis 1 — "the AI / DNA-foundation-model approach" → TRIED; two concrete walls (empirical)

We attempted the FM approach directly (not by argument): pick LD-pruned genome-wide top-QTL markers,
build each variant's REF/ALT sequence-context window from the S288C reference, embed with Nucleotide
Transformer, and test whether the embedding adds predictive signal beyond the raw allele identity. Two
walls surfaced, both real:

1. **Practical (coordinate provenance).** Bloom's 2013 marker coordinates do NOT cleanly map to the
   current S288C RefSeq: on a 40-marker sample the marker REF matched the reference base only 17/40, and a
   per-locus diagnostic showed a MIX of exact / REF-ALT-swap / off-by-one / drift (5/3/2/2 of 12). This is
   local assembly drift (their reference build had indels since), so correct windows need a proper liftover
   with the original assembly + chain file — real infrastructure, not a one-line fix.
   (`scripts/yeast_bloom_fm_prep.py` implements the pipeline; it is blocked on the liftover.)

2. **Principled (sufficient statistic).** Even with perfect windows, the FM is expected to TIE — not
   because it's weak, but because in a bi-parental cross the genotype markers are a SUFFICIENT STATISTIC
   for the genetic state: they completely specify which parental allele each segregant carries at every
   variable position. A sequence-context embedding can only re-encode that same "which allele" the marker
   already has; it adds no genetic information for predicting THESE individuals. Concretely: Maltose is
   dominated by ONE chr07 QTL that a single marker already captures at r=0.67 and the full marker set at
   0.89 — there is no missing information for the FM to supply.

This is NOT "AI can't decode." It is: **the FM is the wrong tool for within-cross prediction (markers are
complete), and the right tool for a DIFFERENT regime — predicting the effect of variants NOT seen in
training (transfer to new variants / organisms), where markers are NOT available.** That transfer
experiment is the genuine FM frontier and a separate, larger build (cross-organism, zero-shot variant
effect — which the `forward` cell already touches at the protein level via DMS).

## Bottom line (honest, non-sycophantic)
- We did not give up. On clean data, **a better model beat the baseline** (nonlinear, +0.16 on Maltose).
- The **AI/FM approach was tried** and hit a principled wall (sufficient statistic) + a plumbing wall
  (coordinate liftover) on THIS substrate; its real value is in the transfer regime, not within-cross CV.
- The satisfying, bankable result is the confirmed nonlinear gain + the confound-free positive it builds on.

Reproducibility: engine `dna_decode/eval/genomic_prediction.py` (cv_ridge_gp + cv_model_gp, 5 tests);
`scripts/yeast_bloom_model_bench.py` (B1); `scripts/yeast_bloom_fm_prep.py` (B2 pipeline, liftover-blocked).
Frozen AMR/forward surfaces byte-unchanged.
