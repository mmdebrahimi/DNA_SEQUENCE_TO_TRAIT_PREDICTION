# CABBAGE scored cell — *Salmonella* × ceftriaxone, N=7,832, and the control that exonerates the cohort

**2026-10-11.** Scored the same day as *Salmonella* × ciprofloxacin, on the same organism and the same
cohort source, under a different rule shape. It is the best-performing of the four CABBAGE cells and it
settles where the cipro over-call comes from. Artifacts:
`wiki/cabbage_scored_cell_salmonella_ceftriaxone_2026-10-11.json`,
`wiki/cabbage_fn_anatomy_salmonella_ceftriaxone_2026-10-11.json`,
`wiki/cabbage_qrdr_multihit_2026-10-11.json`.
Reproduce: `uv run python scripts/cabbage_score_cell.py --species "Salmonella enterica" --drug ceftriaxone --organism Salmonella --max-isolates 8000`.

## Result

| | |
|---|---|
| strict-labelled | 22,785 (1,865 R / 20,920 S) |
| seeded-random subset | 8,000 |
| **scored** | **7,832** |
| confusion | tp 591 / fn 36 / tn 7,190 / fp 15 |
| **accuracy** | **0.9935** |
| **sensitivity** | **0.9426** (591/627) |
| **specificity** | **0.9979** (7,190/7,205) |
| fetch failures | 168 (2.1%) |

**627 scored R** — by a wide margin the best-powered resistant class of any cell measured here, against
4–5 R for Salmonella in the provenance-disjoint arm. Determinants are individually excellent:
`blaCMY-2` PPV **0.987** on 608 carriers, `blaCTX-M-65` **0.984** on 184, and every other `blaCTX-M`
variant at 1.000. Only two sit low and both are tiny — `blaOXA-10` 0.400 (n=5), `blaDHA-1` 0.500 (n=2).
All 36 false negatives carry zero counted determinants: the declared blind spot again.

## The within-organism control, and what it does *not* show

Same organism, same cohort source, different rule. R-rate at **exactly one** counted determinant:

| cell | rule | R-rate at 1 determinant | n |
|---|---|---|---|
| *Salmonella* × ceftriaxone | `subclass_any` (acquired β-lactamase) | **0.9756** [0.963, 0.984] | 820 |
| *Salmonella* × ciprofloxacin | `qrdr_point` (target-site) | **0.0183** [0.010, 0.034] | 493 |

**53× apart in the same organism on the same cohort.** That exonerates the data: the Salmonella cipro
over-call cannot be blamed on Salmonella label quality, on the portal's AMRFinder run, or on this
cohort's composition, because all three are shared with a cell that reaches spec 0.9979. Whatever is
wrong is in the rule.

**The tempting generalisation is REFUTED, and the third cell is what kills it.** The obvious reading —
`threshold=1` suits acquired-enzyme mechanisms and mismatches target-site ones — is wrong.
*Campylobacter* × ciprofloxacin is **also** target-site (`gyrA_T86I`) and a single hit predicts at
**0.9891**. So mechanism *class* does not decide it. What decides it is whether **one determinant is
clinically sufficient for that organism–drug pair**: `gyrA_T86I` alone clears the Campylobacter cipro
breakpoint, a single *gyrA* hit alone does not clear the Salmonella cipro breakpoint. `threshold=1` is
a sufficiency claim per organism–drug pair, and it is a per-pair empirical question — not something
inferable from the mechanism's class.

Three cells, three rule-shape/organism combinations:

| cell | mechanism | R-rate at 1 determinant | `threshold=1` appropriate? |
|---|---|---|---|
| *Campylobacter* × cipro | target-site | 0.9891 | yes, measured |
| *Salmonella* × ceftriaxone | acquired enzyme | 0.9756 | yes, measured |
| *Salmonella* × cipro | target-site | **0.0183** | **no, measured** |

Supporting detail from the same probe: R-rate is **strictly monotone** in determinant count for
ceftriaxone (0.0048 → 0.9756 → 1.0) and **not** for Salmonella cipro, where it fails on one dip; and
**all 64 *parC* carriers in the Salmonella cipro cohort also carry *gyrA*** (fraction 1.000, up from
46/46 before the cache grew), so *parC*'s high PPV there is *parC* marking the multi-hit genotype.

## Honest limits

- **Authorises nothing.** Everything here is a diagnosis. Changing `threshold` or the per-substitution
  catalog edits the sha256-pinned frozen surface and retires the active 2026-08-31 lock — a **USER
  AUTHORITY** call. Frozen five byte-unchanged; `verify_lock` ok=True, no drift.
- **SOURCE-CONCENTRATED**, like the cipro cell: 97.6% of the strict-labelled cohort is one country, so
  this cell **fails** this repo's own diversity bar. The *E. coli* gentamicin cell is the only CABBAGE
  cell that passes it. A concentrated cohort cannot reveal a determinant family it does not contain.
- **The control is strong on the data question and silent on prevalence.** It shows the cohort and
  labels support a near-perfect call under one rule; it does not show the cipro label set is free of
  borderline isolates. Reduced-susceptibility isolates sit near the cipro breakpoint by definition,
  which remains the one caveat that could partly explain `gyrA_D87Y`'s 167 susceptibles without any
  rule defect — and the ceftriaxone control cannot address that, because ceftriaxone has no analogous
  reduced-susceptibility stratum.
- **Determinant-count strata are OBSERVATIONAL and MARGINAL**, not a causal dose, and the probe scores
  whatever is in the local cache — so its denominators grow as the cache grows (the Salmonella cipro
  one-determinant stratum read n=344 before the ceftriaxone run and n=493 after). Compare rates, not
  raw counts, across runs of it.
- **NOT clonality-corrected** (`cabbage-clonality-g8-blocked-no-mash`).
- **The three cells are CONTRASTS, not matched controls** — different sizes and source compositions.
  The Salmonella cipro/ceftriaxone pair is the closest to matched, sharing organism and cohort source.
- Subset is a seeded random sample of 8,000 of 22,785; 168 fetch failures (2.1%) excluded and counted,
  again degrading on later accessions. Not temporally locked — not prospective relative to the active
  lock.
