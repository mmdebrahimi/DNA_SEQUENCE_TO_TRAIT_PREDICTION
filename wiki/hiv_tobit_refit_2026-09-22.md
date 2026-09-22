# Tobit refit: the frozen verdict is `IMPLEMENTATION_SUSPECT`, and the control — not the code — is why

**Verdict as frozen: `IMPLEMENTATION_SUSPECT`.** It is recorded unchanged. The substantive result below
is therefore reported as **PROVISIONAL and is NOT claimed**.

Run: `uv run python scripts/hiv_tobit_refit.py` (offline, **CPU, no GPU needed or used** — a
censored-regression MLE on ~2,000 rows; ~3 min for 19 drugs).
Bar: `wiki/hiv_tobit_refit_acceptance_bar.json` (frozen after the mechanism analysis, before any fit).

## The question

`wiki/hiv_fold_censoring_audit_2026-09-22.md` measured that the Stanford fold tables are right-censored
at 100 and named this as the decidable follow-on it did **not** answer. The NRTI/PI/INSTI v0.1 catalogs
were selected by thresholding a multivariate-OLS log10-fold coefficient at ≥ log10(1.5), fit on that
censored response. Does refitting with a censoring-aware likelihood change **which mutants are
selected**? (It asks nothing else — no CV is re-run, so **no new balanced-accuracy gain figures** are
produced.)

## What happened

The bar's negative control fired: Tobit must reduce to OLS as censoring → 0, and I required max drift
< 0.01 on the <2%-censoring drugs. **stavudine (0.2% censored) drifted 0.0101 and cabotegravir (1.6%)
drifted 0.0204.** By the frozen rule that is `IMPLEMENTATION_SUSPECT` — interpret nothing.

**The implementation is not broken, and two independent lines of evidence say so.**

**1. The zero-censoring limit is exact.** The two drugs with *literally* 0% censoring both drift
**exactly 0.0000**:

| drug | censored | max drift |
|---|---|---|
| tenofovir | 0.0% | **0.0000** |
| bictegravir | 0.0% | **0.0000** |
| didanosine | 0.1% | 0.0032 |
| stavudine | 0.2% | 0.0101 |
| abacavir | 0.3% | 0.0097 |
| dolutegravir | 0.8% | 0.0077 |
| cabotegravir | 1.6% | 0.0204 |

That is a smooth gradient, and my fixed cut sliced through the middle of it.

**2. Synthetic recovery, where the truth is known** (now a permanent test). Data generated from known
coefficients, then censored:

| censoring | OLS mean abs error | Tobit mean abs error | OLS estimate of a true **0.900** |
|---|---|---|---|
| 2.4% | 0.0196 | 0.0164 | 0.872 |
| 9.9% | 0.0373 | 0.0160 | 0.774 |
| 26.9% | 0.1007 | **0.0200** | **0.559** |

Tobit is 5× more accurate at 27% censoring, and **the attenuation is worst for the largest
coefficient** (38% for a true 0.9 vs 11% for a true 0.1) — the audit's stated mechanism, confirmed
where the answer is known.

## The control was mis-specified — the third such bar this session

I set an **absolute** drift threshold (0.01) to police a quantity that is **continuous in the censored
fraction and scales inversely with n**. Cabotegravir has n=64 and *one* censored observation; a 0.02
coefficient shift there is noise against a selection threshold of 0.176. Stavudine failed by 0.0001.

This is the same error class as the two bars overridden earlier this session — a summary statistic
policing a failure mode it cannot resolve — and it is the **third**. The pattern is now unmistakable:
**I keep freezing a scalar where the mechanism is graded.** The right control was available and I did
not use it: *recover known coefficients from synthetic censored data*, which tests the implementation
directly instead of proxying it.

**The verdict is NOT rewritten.** Two overrides are already on the record; a third would make
pre-registration decorative. Re-freezing this bar with a scale-aware control is a user-authority call.

## PROVISIONAL result — not claimed under the frozen verdict

With that caveat stated plainly, the numbers:

| drug | censored | flips (knife-edge marked \*) |
|---|---|---|
| **lamivudine** | **45.2%** | +K70R 0.144→0.192, **+K70T −0.202→1.312**, +T215X 0.171→0.207\*, +V75L 0.088→0.328 |
| zidovudine | 13.7% | +T215E 0.107→0.188, −T215V 0.205→0.152 |
| atazanavir | 11.8% | +L33F 0.161→0.194, −V82A 0.181→0.139\*, −V82C 0.197→0.143 |
| lopinavir | 14.4% | +V32I 0.173→0.200\* |
| nelfinavir | 9.2% | −G48M 0.186→0.118\*, −M46V 0.194→0.120, +N88D 0.149→0.207 |
| darunavir | 6.2% | −V82C 0.240→0.162, +V82I 0.107→0.182 |
| indinavir | 5.5% | +L33I 0.166→0.191 |
| fosamprenavir | 5.4% | −G48M 0.199→0.162, −L33M 0.196→0.126 |
| didanosine | **0.1%** | +K70N 0.174→0.177\* |

**19 flips: 14 decisive, 5 knife-edge.** The knife-edge column is load-bearing — didanosine's `K70N`
moves 0.174 → 0.177 at 0.1% censoring, which is a coin toss on a coefficient that already sat on the
threshold. Unseparated it would read identically to lamivudine's `K70T` going **−0.202 → +1.312**, a
sign flip from severe attenuation at 45% censoring.

The pre-registered directional prediction **holds in aggregate**: 11 gains vs 8 losses. It does **not**
hold per-drug — fosamprenavir loses two entries and gains none. Losses are not paradoxical (Tobit
reweights rather than shifting everything up), but the prediction was stated per-drug and is reported
as failing there rather than explained away.

## Honest limits

- **Produces no new gain figures.** The 5-fold CV was not re-run; this asks only whether the
  *selection decision* moves. Re-deriving the published gains is a separate job.
- **Not a decoder change under any verdict.** No v0.1 catalog is deployed — `hiv_amr.py` routes
  NRTI/PI/INSTI through position-based classes. `hiv_amr.py` is unmodified.
- The censoring point is **inferred** from a pile-up with an empty tail above it, not read off assay
  documentation.
- Tobit assumes normal homoscedastic latent errors on the log10 scale — the same assumption OLS already
  makes here, not extra evidence for it.
- Frozen AMR surface untouched.

## Reusable

**An anti-overfitting control can fail in the same way a bar can: by being a scalar where the mechanism
is graded.** Three times this session a frozen threshold was the wrong instrument — a net count for a
directional failure, an unachievable target, and now an absolute drift cut across a continuous
gradient. When the thing you are policing varies smoothly with a measurable quantity, the control has
to vary with it too, or be replaced by a direct test.

**Prefer a control that tests the property directly over one that proxies it.** "Does the estimator
recover known coefficients?" was available, cheap, and decisive. "Does it drift less than 0.01 on
drugs that happen to have little censoring?" was neither.
