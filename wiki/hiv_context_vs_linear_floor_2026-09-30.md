# The linear floor IS the ceiling — context does not beat additive on the HIV blind spot

**2026-09-30.** The shipped supervised blind-spot complement is logistic regression over one-hot
substitution tokens — an **additive** model, the weakest member of the supervised family. Its measured
leave-one-study-out blind-spot AUROC (**RT 0.8102, IN 0.8923**) was therefore a **floor**, not a ceiling,
and said nothing about whether interactions help. Same data, same split, same subset, same bar — swap the
estimator.

**The bar was frozen and committed (`2cf723d`) before any number existed**
(`wiki/hiv_context_vs_linear_floor_acceptance_bar.json`), and the script **reads** it rather than restating
it. The prior was stated up front, not discovered afterwards: `wiki/hiv_epistasis_result_2026-07-11.json`
returned `FAIL_ADDITIVE_SUFFICES` (pairwise interactions beat additive on **2 of 24** powered cells).

## Result

| gene | blind spot n | R | linear floor | pairwise | non-linear | Δ pairwise | Δ non-linear |
|---|---|---|---|---|---|---|---|
| RT (EFV) | 1,883 | 201 | **0.8102** | 0.7796 | 0.8034 | **−0.0305** `[−0.050, −0.011]` | −0.0068 `[−0.027, +0.012]` |
| IN (RAL) | 844 | 59 | **0.8923** | 0.8790 | 0.8362 | −0.0133 `[−0.043, +0.011]` | **−0.0562** `[−0.097, −0.013]` |

**Verdict by the frozen rule: `LINEAR_FLOOR_IS_THE_CEILING`.** All four context arms lose, and **two lose
significantly** (paired-bootstrap CI entirely below zero). This is stronger than "no improvement" — added
capacity actively **hurts**. It hurts on overall AUROC too (RT 0.9545 → 0.9407/0.9398; IN 0.9668 →
0.9611/0.9335), so it is not blind-spot-specific.

The anchor reproduces the measured floor **exactly, +0.0000 on both genes** — the linear arm here is
bit-identical to the deployability measurement, so the deltas are a real comparison rather than two
unrelated numbers.

## Why the negative is strong rather than weak

`build_pairs` ranks pairs by co-occurrence, and on RT it recovered **canonical biology**: `41L:215Y`,
`41L:210W`, `215Y:210W`, `184V:41L` — the **TAM cluster** (M41L / L210W / T215Y) plus M184V. The
construction found genuine, well-known epistatic pairs and they **still did not help**. The failure is not
"we only found noise."

## Named limitation — and it decides which arm carries the claim

**Pair selection is by PREVALENCE, which handicaps the pairwise arm.** The top-40 most prevalent tokens in
an NNRTI cohort of treatment-experienced patients are dominated by **NRTI** TAMs — real epistasis, but
**not efavirenz drivers**. So the pairwise delta is a statement about prevalent-token pairs, not about the
most *predictive* pairs.

**The general claim therefore rests on the non-linear arm**, which sees every feature and can form any
interaction — and which also lost on both genes (−0.0068 RT, −0.0562 IN). A cheap unexplored follow-up, if
anyone wants to push: select pairs by univariate association rather than prevalence. **Not run here.**

## Scope — what this does and does not establish

- **This does NOT show "HIV resistance is additive."** With ~1,063 (RT) / 405 (IN) features against
  1,883 / 844 blind-spot isolates, capacity may simply be **unaffordable out-of-distribution** under
  leave-one-study-out. That is a different statement, and the earlier epistasis memo was careful about the
  same distinction.
- **It does NOT test attention over raw sequence.** "Context" here means interactions among substitution
  tokens; a sequence model is a different and costlier question.
- IN-DISTRIBUTION to Stanford, ~96% subtype B. Leave-one-study-out is the OOD split, **not** an independent
  cohort.
- Two genes of **one virus**. PI is excluded — its blind spot holds 3 R of 910 and is unscoreable.
- Fold labels are right-censored at 100.0, but the binary label is a threshold at fold ≥ 3 far below the
  ceiling and **every arm sees identical labels**, so censoring cannot favour one.

## What this licenses

Exactly what the frozen bar said a negative would: **shipping the linear complement as-is**, with the floor
documented as the measured ceiling for this feature space. It **closes candidate row 5** rather than
leaving it open indefinitely. Two unrelated framings now agree — continuous regression under KFold
(2 of 24) and blind-spot classification under leave-study-out (0 of 4) — which is the useful part.

## Reproduce

```
uv run python scripts/hiv_context_vs_linear_floor.py --self-check   # bar + floor + prior, no data
uv run python scripts/hiv_context_vs_linear_floor.py               # 2 genes x 3 estimators
```

Artifact `wiki/hiv_context_vs_linear_floor_2026-09-30.json`; 11 tests
`tests/test_hiv_context_vs_linear_floor.py`. Exits 3 if the floor does not reproduce. Frozen AMR surface
byte-unchanged (read-only).
