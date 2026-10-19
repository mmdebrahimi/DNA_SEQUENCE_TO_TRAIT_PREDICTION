# CABBAGE scored cell — *Salmonella* × ciprofloxacin, N=7,834, and why `threshold=1` does not transfer

**2026-10-11.** The frozen cipro rule (`counter=qrdr_point`, `threshold=1`) is near-perfect on
*Campylobacter* and over-calls badly on *Salmonella*. This cell measures the gap and locates its cause.
Artifacts: `wiki/cabbage_scored_cell_salmonella_ciprofloxacin_2026-10-11.json`,
`wiki/cabbage_fn_anatomy_salmonella_ciprofloxacin_2026-10-11.json`,
`wiki/cabbage_qrdr_multihit_2026-10-11.json`.
Reproduce: `uv run python scripts/cabbage_score_cell.py --species "Salmonella enterica" --drug ciprofloxacin --organism Salmonella --max-isolates 8000`,
then `scripts/cabbage_fn_anatomy.py` and `scripts/cabbage_qrdr_multihit.py` with the same cell.

## Result

| | |
|---|---|
| strict-labelled | 21,280 (196 R / 21,084 S) |
| seeded-random subset | 8,000 (`sample_seed=20261010`) |
| **scored** | **7,834** |
| confusion | tp 55 / fn 12 / tn 7,421 / **fp 346** |
| accuracy | 0.9543 |
| **sensitivity** | **0.8209** (55/67) |
| **specificity** | **0.9555** (7,421/7,767) |
| fetch failures | 166 (2.1%) |

Set beside the two cells scored the same way, the specificity is the outlier: *Campylobacter* × cipro
0.9966 and *E. coli* × gentamicin 0.9951, against **0.9555** here. 346 false positives is not noise.

This is the first POWERED *Salmonella* AMR measurement in this repo — the provenance-disjoint arm has
only 4–5 R for Salmonella cipro/gent/tet and the report card reads UNDERPOWERED. **Scope care: this is a
separate, namespace-separate arm, so it does NOT flip that provdisjoint powering verdict.** It supplies
a powered number where that arm could not.

## The over-call is one substitution, and the cause is the threshold

Per-determinant PPV, worst and best (full table in the anatomy artifact):

| determinant | carriers | R | S | PPV | Wilson CI |
|---|---|---|---|---|---|
| **`gyrA_D87Y`** | **169** | **2** | **167** | **0.012** | [0.003, 0.042] |
| `gyrA_D87N` | 93 | 36 | 57 | 0.387 | [0.294, 0.489] |
| `gyrA_S83F` | 86 | 46 | 40 | 0.535 | [0.430, 0.637] |
| `gyrA_D87G` | 29 | 10 | 19 | 0.345 | [0.199, 0.527] |
| `gyrA_S83Y` | 23 | 4 | 19 | 0.174 | [0.070, 0.371] |
| `gyrB_S464F` | 14 | 0 | 14 | 0.000 | [0.0, 0.215] |
| … | | | | | |
| `parC_S80R` | 28 | 27 | 1 | 0.964 | [0.823, 0.994] |
| `parC_S80I` | 18 | 18 | 0 | 1.000 | [0.824, 1.0] |
| `aac(6')-Ib-cr5` | 24 | 21 | 3 | 0.875 | [0.690, 0.957] |

`gyrA_D87Y` alone — 167 susceptible carriers — accounts for roughly half the false positives. Every
single-hit *gyrA* substitution sits below 0.54; every *parC* substitution sits at or above 0.96.

**Tested, not asserted** (`scripts/cabbage_qrdr_multihit.py`). The textbook account is that one *gyrA*
QRDR hit gives *Salmonella* only reduced susceptibility — historically nalidixic-acid resistance, cipro
MIC still under the clinical breakpoint — while clinical cipro resistance needs several hits. Two
predictions were written down before looking:

**P1′ — the count carries what presence does not.** R-rate at exactly one counted determinant
**0.0233 [0.012, 0.045]** (n=344) against two or more **0.8246 [0.706, 0.902]** (n=57): a **35.4×**
difference with **disjoint** confidence intervals. A `threshold=1` rule over-calls here by construction.

**P2 — *parC* marks the multi-hit genotype rather than predicting alone.** **All 46 *parC* carriers also
carry *gyrA* — fraction 1.000, no exceptions.** So *parC*'s 0.96–1.00 PPV is not *parC* being
informative on its own; it is *parC* appearing only on top of a *gyrA* hit. A rule that counted only
*parC* would get the right answer for the wrong reason and would break on any population where *parC*
arises first.

**The cross-organism contrast is the cleanest statement of all of it.** Same frozen rule, same drug, at
exactly one counted determinant: ***Salmonella* 0.0233 vs *Campylobacter* 0.9891 — 42.5×.** One QRDR
hit predicts resistance in *Campylobacter* (`gyrA_T86I`, PPV 0.995 on 1,000 carriers) and does not in
*Salmonella*. That is where the specificity gap between the two cells comes from, and it is a property
of the organism's resistance genetics, not of the implementation.

*Campylobacter* has **zero** *parC* carriers in its cohort, so the *parC* story is Salmonella-specific.

## My own pre-registered test failed, and it is reported as it fell

P1 was originally specified as a **strictly non-decreasing** R-rate across determinant-count strata.
**Salmonella FAILS it**: the 4-determinant stratum is 0.857 (n=14) against the 3-determinant stratum's
1.000 (n=17), with heavily overlapping Wilson intervals. The premise was "R-rate rises with count";
strict monotonicity across small strata was a brittle way to operationalise it, and one dip on two
isolates is enough to break it. The failed test is kept in the artifact next to the contrast that
replaces it, rather than swapped out for the version that passes — same discipline as the
`salmserovar` threshold episode, where following a badly-specified pre-registered rule would have
regressed the cell.

The underlying progression is steep and unambiguous: R-rate **0.0016 → 0.0233 → 0.3333 → 1.0000** across
0, 1, 2 and 3 counted determinants.

**All 12 false negatives carry ZERO counted determinants** — the declared blind spot again.

## Honest limits

- **DIAGNOSES an over-call; authorises nothing.** Changing the rule's threshold or its per-substitution
  catalog edits the sha256-pinned frozen surface and retires the active 2026-08-31 lock — a **USER
  AUTHORITY** call. Frozen five byte-unchanged; `verify_lock` ok=True, no drift.
- **SOURCE-CONCENTRATED.** 95.4% of the strict-labelled cohort is one country. Unlike the *E. coli*
  gentamicin cell, this one **fails** this repo's own diversity bar, so it cannot reveal a determinant
  family it does not contain.
- **Determinant count is OBSERVATIONAL and MARGINAL, not a causal dose.** Carriers of two determinants
  differ from carriers of one in more than the count, and a carrier may hold determinants the rule's
  class filter does not admit — so the count is the rule's count, not a complete genotype.
- **Reduced-susceptibility isolates sit near the breakpoint**, which is exactly the population in
  question, so some S labels here are genuinely borderline. STRICT tiers reduce that but do not remove
  it, and this is the one caveat that could partly explain `gyrA_D87Y`'s 167 susceptibles without any
  rule defect.
- **The 2+ stratum is small (n=57)** and the 2-determinant stratum is n=12; the 35.4× contrast rests on
  disjoint CIs, not on a large sample.
- **NOT clonality-corrected** (`cabbage-clonality-g8-blocked-no-mash`). Salmonella surveillance is
  clonally structured and a single-determinant stratum may be a few clones.
- **The *Campylobacter* arm is a CONTRAST, not a matched control** — different cohort, size and source
  composition.
- Subset is a **seeded random** sample of 8,000 of 21,280, not a prefix (an accession-ordered prefix is
  measurably easier). 166 fetch failures (2.1%) excluded and counted. Coverage degraded on later
  accessions — zero failures through the first ~5,800, rising to 166 — a datapoint toward
  `cabbage-genomes-tree-coverage-uncensused`, not a census.
- Not temporally locked: these isolates are **not** prospective relative to the active lock.
