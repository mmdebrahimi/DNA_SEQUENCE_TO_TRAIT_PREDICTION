# Generality test: does the yeast finding hold on a SECOND organism's confound-free cross? (2026-08-02)

Second confound-free cross = the **mouse BXD recombinant-inbred panel** (C57BL/6J × DBA/2J; R/qtl2
mirror, free/no-DUA) — a MAMMAL, maximally different from the yeast segregants. Same pipeline:
Layer-1 ridge genomic prediction (+ permutation null) + Layer-2 gbm-vs-ridge. 7,320 markers × 198
strains; quantitative traits at coverage n~80-90 (the BXD phenotype matrix is sparse — each trait is a
small study). Marker-subsampled to 2,440.

## Layer 1 — "genomic prediction decodes" → GENERALIZES across kingdoms

6/12 traits beat the label-permutation null; the heritable ones decode clearly:
- **Brain weight r=0.574** (+ variants 0.53 / 0.48 / 0.40), **Body weight r=0.301**, Trait_22 0.273.

So on a mammalian confound-free cross, genomic prediction decodes quantitative traits — the SAME core
result as the yeast segregants (which decoded 12/12). **The robust finding generalizes fungi → mammal.**
The 6 that don't beat null are lower-heritability / noisier traits at this small n (eye/heart weight).

## Layer 2 — "nonlinear beats linear" → does NOT replicate in mouse BXD

gbm beats ridge only 3/12, **mean delta −0.060** (linear ridge is generally better here) — the OPPOSITE
of yeast (gbm won 26/46, +0.023). So the nonlinear/epistasis advantage is **NOT a universal property.**

**Cause — OPEN (honestly not settled).** Two candidates, not disentangled:
1. **Genetic architecture.** Mouse morphological traits (brain/body weight) are classically highly
   ADDITIVE + polygenic with modest epistasis; yeast growth traits (Maltose etc.) carry strong epistasis
   (Bloom flagged up to ~50% for some). A linear model is near-optimal for additive traits, so there is
   simply less nonlinear signal for gbm to recover in mouse.
2. **Sample size.** BXD n~85 vs yeast n=1008 → less power to fit interactions.

**A sample-size kill-test was INCONCLUSIVE, so I do NOT claim (2) as the cause.** Subsampling yeast (where
gbm won) to n=85: gbm STILL beat ridge on 2/3 traits (Maltose delta +0.49, Cadmium +0.37 — larger, not
smaller), only Copper flipped (−0.06). That test was also partly confounded (fixed vs tuned ridge alpha),
so it neither confirms nor cleanly refutes small-n — it just rules out "small n straightforwardly kills
gbm." Architecture (1) is the more plausible primary driver, but it stays a labelled hypothesis.

## Verdict

- **The core finding is GENERAL:** genomic prediction decodes quantitative traits on a confound-free cross,
  across kingdoms (yeast fungus → mouse mammal). This is the robust, bankable generality result.
- **The nonlinear/epistasis advantage is CONTEXT-DEPENDENT, not universal:** it appears where epistasis is
  strong + n is large (yeast growth) and not on more-additive mammalian morphological traits at low n. The
  exact driver (architecture vs power) is an open, honestly-unresolved question.

## Honest scope
BXD per-trait n~85 is low-powered (yeast had 1008) — brain weight is the cleanest signal; the weaker
traits are underpowered. A denser mammalian cross (or more BXD strains) would sharpen Layer 2.
Reproducibility: `scripts/bxd_gp_arm.py` (loader + pipeline); data `D:/dna_decode_cache/bxd/` (free R/qtl2
mirror; not committed). Fixed a real bug: `cv_ridge_gp(n_perm=0)` crashed on the empty null (guarded +
tested). Frozen AMR/forward surfaces byte-unchanged.
