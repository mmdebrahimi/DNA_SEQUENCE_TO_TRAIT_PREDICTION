# Cross-organism conditional essentiality — PRE-REGISTRATION (frozen before any solve)

**Date:** 2026-08-20 · **Status:** frozen before execution
**Axis:** the first NON-E.-coli test. Different organism, different genome-scale model, different gene set.

## Why this is the strongest replication test so far

The three findings have now replicated across a **substrate** axis (carbon → nitrogen), but everything so
far shares one organism and one model. If they are properties of the *method*, they must survive changing
both. If they are properties of *iML1515*, they will not.

The stress axis just closed as a structural NO-GO (a stoichiometric model cannot represent a poison), so
cross-organism is the remaining named axis.

## Feasibility, probed before freezing this (per R2 — derive, don't assert)

`feba.db` has **48 organisms**. Exactly **two** have a genome-scale model in the repo's BiGG registry:
E. coli (`iML1515`, already done) and **Pseudomonas putida KT2440** (`iJN1463`). So the axis is n=1 new
organism — but unlike the stress axis's n=2 conditions, that is fine: the replication unit here is the
*panel*, and P. putida has a full one.

| check | result |
|---|---|
| `iJN1463` loads | 2,927 reactions / 1,462 genes / 348 exchanges |
| wildtype growth (default medium) | 0.586118 |
| carbon sources in the assay | 47 |
| **mapped to iJN1463 exchanges** | **13 (28 %)** |
| fitness data | 4,778 genes |
| **model genes with fitness data** | **1,440 of 1,462 (98.5 %)** |

The last row is the one that makes this viable. P. putida locus IDs (`PP_0026`) match between the model
and the Fitness Browser **directly** — no identifier-mapping layer, and therefore none of the 0 %-overlap
failure mode documented for gene symbols in `CLAUDE.md`.

13 conditions is exactly the size of the nitrogen panel that replicated all three findings.

## Pre-registered predictions (E. coli → P. putida)

| # | prediction | E. coli carbon | E. coli nitrogen | P. putida bar |
|---|---|---|---|---|
| **P1** | ratio distribution bimodal; threshold band near-empty | 0 of 16,676 | 0.0 % | **< 1 %** in [0.001, 0.05] |
| **P2** | missed-essential cells are mostly FLAT | 76.9 % | 91.1 % | **≥ 50 %** |
| **P3** | per-cell beats the best-constant null | positive | +0.1444 | **> 0** |

**P1 is again load-bearing.** It is the claim that underwrites the closed E-Flux negative. Two organisms
and two models agreeing would make it a property of the *readout*; a failure here would confine it to
iML1515 and re-open that question.

**A failure of P3 is NOT a failure of the run.** P. putida's model is less curated than iML1515 and its
gold standard is different. If conditional essentiality is simply not predictable there, that bounds the
method's generality — which is the point of testing it.

## Determinism — mandatory, checked BEFORE any verdict is read

Same claim-level gate as nitrogen (`dna_decode.fba.nitrogen.determinism_verdict`), reused unchanged:

1. **Zero call flips** across `FRAC = 0.01` between two full passes.
2. **Headline metric identical** between passes.
3. **Derived safety margin** `min_margin_to_threshold / max_abs_delta ≥ 1000` — computed from the data,
   so a panel whose ratios crowd the threshold fails on its own.

On failure, `redact_unverified` **removes** every solve-derived number rather than printing it beside a
false flag. `processes=1` is pinned.

## Method (frozen)

- Model `iJN1463`; labels Fitness Browser RB-TnSeq `orgId=Putida`, `expGroup='carbon source'`,
  threshold `fit < -2.0`; `FRAC = 0.01`; NaN growth → essential (the shipped genuine-essentiality coding).
- `apply_carbon_condition` is reused **unchanged** — it is already organism-agnostic (it takes a model and
  an exchange id), and closes every other candidate carbon exchange so a residual uptake cannot make every
  condition silently score as glucose.
- Two-sided subset recomputed, not inherited.

## Named caveats (recorded now)

1. **n = 1 new organism.** This is a replication, not a survey. Two organisms agreeing is not "the method
   generalises across bacteria".
2. **13 of 47 sources.** P. putida's assay is rich in compounds iJN1463 does not carry (diols, alcohols,
   fatty acids); the panel is a biased subset toward central metabolism.
3. **Glucose is over-represented in the assay** (29 of 128 carbon experiments); replicate averaging makes
   its condition better determined than the others.
4. **iJN1463 is less validated than iML1515.** A weaker result here may reflect model quality rather than
   the method.
