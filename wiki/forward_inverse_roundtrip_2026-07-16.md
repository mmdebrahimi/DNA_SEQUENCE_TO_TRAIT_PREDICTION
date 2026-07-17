# Oracle-in-the-loop molecular INVERSE — decisive falsifier (2026-07-16)

**Verdict: `PASS_INVERSE_BEATS_BASELINES`** — margin **+52.9%** vs the better real baseline (material bar 25%); splits separated: `True`.

*The handoff (2026-07-16) said: run this before any build. This is the go/no-go.*

## The claim under test

given a desired molecular-effect target T, the forward oracle can PROPOSE the edit achieving it (label-free inverse design) — i.e. inverse design **using the DMS-validated forward cell as
label-free ground truth**, which is the move that dodges this project's binding constraint (labels,
not models).

## Result

blaTEM (`BLAT_ECOLX_Stiffler_2015`), **1546 single-nt-accessible variants** across 263 positions · 39 targets × 6 position-splits.

Mean |measured − target| on the **wet-lab label** (lower is better):

| method | top-1 | best-of-5 | best-of-5 across splits |
|---|---:|---:|---|
| **esm2 (the oracle)** | 0.4750 | **0.1164** | [0.0903 .. 0.1435] |
| blosum62 (real baseline) | 0.8819 | 0.2474 | [0.2234 .. 0.3012] |
| empirical null (no oracle) | 1.1823 | 0.2833 | [0.2798 .. 0.2901] |

The worst ESM split beats the best baseline split, so the margin is not split luck.

## What this does NOT license (the scope is the finding)

**1. It licenses `propose 5, assay 5, keep the best` — NOT `propose 1 and trust it`.**
Top-1 error is 0.475; best-of-5 is 0.116. The 4× gap IS the honest cost of the loop: the metric assumes you can measure the 5 proposals. A single-shot inverse is roughly 4× worse than the headline. (The null is scored the same way — expected min of 5 draws — so the comparison is fair.)

**2. SELECTION works; calibrated MAGNITUDE does not. Informative intervals: `0/6` splits.**
This is the handoff's own warning, now measured: *0.76 is a RANK correlation — the inverse can claim
direction/rank far more than calibrated magnitude.* Confirmed. The conformal interval brackets the
target in every split **and that proves nothing** — coverage holds even for a useless model (J2's
Family-B rail, restated in `forward/dosage.py`). The honest number is
`interval_halfwidth_over_effect_span`, and it says the interval spans >50% of the effect range.
**So: the inverse may say *this edit lands near your target*; it may NOT certify the dose.**

**3. Regime B (molecular fitness) only.** This is blaTEM enzyme fitness, not clinical resistance —
where the same scorer class is *below chance* (ESM2 0.454 vs the catalogue's 0.926).

## Two things the first run got wrong (kept as method evidence)

- **A thin grid overstated the margin.** v1 (9 deciles × 1 split) gave **+71.6%**; the robust
  version (39 targets × 6 splits) gives **+52.9%**.
- **n=1 per cell manufactured a fake failure mode.** v1's single unlucky pick at target −1.774
  (|err| 1.101) read as a *mid-range failure*. Two hypotheses for it — distribution sparsity and
  calibrator-expressibility gaps — were both tested and **both falsified** (that target has *more*
  neighbours than the extreme it nails: 203 vs 168; and the closest expressible prediction is 0.011
  away). It was noise. The fix was the error bar, not a story.

## Design (why it is not circular)

- **non_circular**: calibrator fit on CALIBRATION positions; selection among HELD-OUT positions; grading against the proposed variant's MEASURED wet-lab DMS value
- **split**: by POSITION, interleaved with a per-split offset -- not by variant, not contiguous
- **targets**: quantiles of the real measured-effect distribution
- **why_multi_split**: v1 used 9 deciles x 1 split = n=1 per cell; a single unlucky pick read as a mid-range failure mode. The across-split spread is the error bar.

## Recommended next step

The gate is PASSED, so a build is licensed — but build the **ranking/selection** inverse, not a
dose-certifying one, and carry the `propose-k, assay-k` framing into its interface. The
`cooc-multiedit-inverse` extension stays unproven and adjacent to a closed negative
(`FAIL_ADDITIVE_SUFFICES`); it needs its own falsifier before any build.
