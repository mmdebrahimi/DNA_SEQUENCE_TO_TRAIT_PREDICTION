# The fair test of the alphabet idea — it fails on the gap task, for a documented reason

**2026-09-29.** Follow-up to `wiki/pear_genotype_alphabet_2026-09-29.md`, which I conceded was a strawman:
it tested a *context-free* table on *one* protein. This is the fair version — with pretrained context, on a
real clinical substrate, on the one task the headroom argument cannot reach. It still fails, and the reason
is mechanistic rather than a power problem.

## Why this test and not the last one

The real objection to a GLM here was never the CTX-M-14 run. It is **headroom**: on HIV the additive
bag-of-tokens model reaches **Spearman 0.911 / R² 0.844** (n=2,168, EFV), and `scripts/hiv_epistasis.py`
(2026-07-11) showed pairwise interactions add nothing — **2 of 24 drug-cells CI-positive against a 0.5 bar,
verdict `FAIL_ADDITIVE_SUFFICES`**. The DMS side agrees independently: *"the additive null remains the robust
deployed default."* A transformer cannot add much where a linear model explains 88% of the variance.

That epistasis result **already existed and I nearly rebuilt it** — the second time in one session that
checking the artifact directory first would have saved the work.

**But the headroom argument has one precise blind spot, and it is exactly the alphabet idea.** Every model in
that comparison is keyed on token **identity**. A mutation absent from training is an all-zero row and scores
**exactly zero by construction** — it cannot be wrong, it simply cannot answer. That is why the catalog's
only measured failure mode is *completeness*, and why both confirmed gaps (rmtE1, V179F) were found by
accident when new labels arrived.

A model keyed on token **features** still has a row for an unseen token, so it *can* answer. That is what an
alphabet actually buys, and the 0.911 ceiling does not bound it.

## Design

Stage 1: fit the incumbent additive identity model → one effect coefficient per mutation token.
Stage 2: **leave-one-token-out** — predict a held-out token's coefficient from its features alone, having
never seen its data.

**Stage 1 was validated against known biology before anything downstream was read.** Top-12 coefficients
are dominated by textbook NNRTI DRMs; **K103N ranks 4/421 (+1.231, 529 carriers)**, **Y188L ranks 2
(+1.584)**, and model-free marginals agree (103N Δ+1.231, 188L Δ+1.292, 181C Δ+0.646). Y181C's low fitted
coefficient (+0.185) is the known co-occurrence effect already recorded in this repo, not a defect.

Features: position, carrier count, mutant-AA and wild-type-AA one-hots, BLOSUM62, and **ESM2-650M
masked-marginal conservation** (a legitimate pretrained feature — it sees no fold-change data).

## Result

| | R² | weighted R² | Spearman |
|---|---|---|---|
| **ALL 421 tokens** | | | |
| feature model (alphabet + ESM2) | +0.050 | −0.019 | +0.105 |
| mean-coefficient baseline | −0.005 | −0.006 | — |
| "is this one of the 8 catalogued positions" | **+0.373** | **+0.517** | +0.301 (full fit) |
| **GAP SUBSET — the 398 non-catalog tokens (94.5%)** | | | |
| feature model (alphabet + ESM2) | **−0.079** | −0.056 | **+0.016** |
| mean-coefficient baseline | −0.005 | −0.032 | — |
| catalog baseline | 0.000 — constant here by construction | | |

**B1 (beats the mean) passed on all tokens and FAILS on the gap subset.** The +0.050 was carried entirely by
the 23 catalog tokens. Where gap-prediction actually lives, **the alphabet is worse than a constant.**

**B2 (beats the catalog-position baseline) FAILS decisively on all tokens** — one binary feature explains
R² 0.373 (weighted 0.517) against the alphabet's 0.050 (weighted **−0.019**).

## Why it fails — position 179, which sits outside the catalog

| token | fitted effect | carriers |
|---|---|---|
| 179D | **+0.675** | 44 |
| 179E | **+0.449** | 23 |
| 179F | **+0.261** | 15 |
| 179V | +0.004 | 94 |
| 179I | **−0.098** | **289** |
| 179T | **−0.279** | 13 |

One position carries strong drivers *and* neutral/protective variants, and the **most common substitution
there is the neutral one**. Position and chemistry cannot separate them: from wild-type V, both D and I are
non-conservative, yet their effects are opposite in sign.

This is the same mechanism already recorded for the ESM-vs-catalog negative — *resistance is reached via
chemically conservative substitutions at averagely-conserved sites, so exchangeability scorers call the
resistant residue benign*. Measured here directly on the gap task. **Independent corroboration:** 179F is
precisely the token the doubt layer flagged (15/15 carriers resistant, p = 8.8e-06), reached by a completely
different route.

## Two defects in my own method, both found and fixed before reporting

1. **A cross-class contaminated control.** The catalog baseline initially scanned `hiv_amr` for any integer
   collection and unioned them, giving **17** positions: the 8 NNRTI sites *plus* the 9 NRTI sites. Both are
   in RT, which is why it looked plausible — but M184V and the 41/210/215 TAMs are common in
   treatment-experienced isolates and **NNRTI-neutral**, so against an EFV label the contaminated baseline
   read Spearman **−0.617**. **A non-vacuity guard ("≥ 5 positions") passed the whole time: it checked the
   baseline was not empty, never that it was the right one.** The giveaway was the sign being biologically
   backwards, not the guard. Now class-matched, with an explicit cross-class refusal.
2. **A leave-one-out artifact on a group-mean predictor.** The catalog baseline's LOO Spearman reads
   **−0.743** while the same quantity from a single full fit is **+0.301** (catalog group mean +0.593 vs
   +0.008 elsewhere — biologically correct). Mechanism: dropping a high-coefficient catalog token lowers the
   remaining catalog mean, so that token gets a *lower* prediction — an anti-correlation induced purely by
   excluding a point from the mean it is scored against. With ~23 catalog tokens it dominates. **R² is
   near-immune; Spearman is not.** R² is now the metric of record and the artifact is documented in source.

## Honest limits

- **This is still a linear model over alphabet features, not a transformer.** It is a fair test of *whether
  those features carry generalisable signal* — and they do not — but a transformer was not run. The
  mechanism finding is what makes this more than "we didn't try hard enough": the features themselves are
  uninformative for the gap task, and no architecture recovers signal absent from its inputs.
- One drug, one gene, one class. ~96% subtype B, which removes the ancestry confound but says nothing about
  cross-subtype transfer.
- The target is a **fitted coefficient** with estimation noise and residual collinearity between
  co-occurring tokens; leave-one-token-out does not remove that coupling.
- 421 tokens, only **23** at catalog positions — the all-token comparison rests on a small positive class.
- **NOT ruled out:** features the alphabet does not contain — structural proximity to the drug-binding
  pocket, co-evolution couplings, explicit 3D context. Those are not "alphabet" features, and pursuing them
  is a different proposal from the one tested here.

## What is now measured rather than argued

The incumbent for gap-prediction is **not zero and not weak**: one bit of information ("is this position in
the 8-position catalog") explains R² 0.373 / weighted 0.517 across all tokens. Any learned approach has to
beat that — and on the gap subset it has to beat a constant, which this one does not.
