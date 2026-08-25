# Systematic rule-gap screen — the rmt gap is the ONLY actionable one (2026-08-25)

**Verdict: the screen works, and it found NO new actionable rule gap.** `_MEASURED_GAPS` correctly stays
at one row. This bounds the disclosure layer's scope with evidence instead of leaving it as "we found one,
there are probably more."

## Why the screen was built

The 2026-08-24 prospective cohort exposed one gap **by hand**: E. coli × gentamicin missed 24 of 28
resistant isolates carrying a 16S rRNA methyltransferase the frozen `Subclass=GENTAMICIN` rule cannot see.
`wiki/determinant_blindness_atlas.json` suggested the shape might be widespread — across all 12 NCBI-PD
external-validation cells, **every** invisible resistant isolate is `rule_limited` (a determinant IS
present, the rule just doesn't count it) and **zero** are `truly_invisible`. But the atlas reports an
aggregate fraction and never names *which* determinant. `scripts/rule_gap_screen.py` names them.

**The statistic** is the one that made the rmt finding convincing, not the raw count: a token frequent
among missed-R isolates proves nothing (passengers are frequent everywhere). What separated `rmt` was that
it was frequent in missed-R and **absent from the susceptible set**. So each candidate is scored on both
`miss_rate` and `s_rate` and ranked by the gap.

## Result: 2 candidates from 12 cells, and neither is actionable

| cell | R | missed | verdict |
|---|---|---|---|
| N. gonorrhoeae × azithromycin | 110 | 110 | **0 candidates — correct** |
| N. gonorrhoeae × tetracycline | 34 | 23 | `porB1b_G120K` gap 0.57 → already-disclosed ceiling |
| N. gonorrhoeae × cefixime | 19 | 4 | `porB1b_A121D` gap 0.42 → **underpowered** (4 missed) |
| 9 other cells | — | 0–1 | no missed-R, or underpowered |

**The azithromycin zero is the screen's validation, not a miss.** That cell misses 110 of 110 resistant
isolates, so a naive pattern-matcher would emit its most common tokens. The screen emits nothing, because
no token clears the passenger filter — correct: gonococcal azithromycin resistance is mtr-efflux driven and
there is no gene to count. A screen that produces a confident answer there would be worthless.

## Why the one real hit is NOT a new gap

`porB1b` mutations recover a genuine mechanism — which is exactly why this needed checking rather than
celebrating. Three findings, in order of decisiveness:

**1. The rule already documents it.** `neisseria_amr.call_ng_tetracycline`'s own docstring:

> chromosomal low-level tet-R (rpsJ + mtrR + **penB** cumulative, MIC 2-4) is NOT cleanly
> determinant-separable from tet-S (which carries the same near-universal markers) → ~68% of tet-R is not
> catchable from AMRFinder determinants. **This is a genuine multi-locus-cumulative ceiling, not a rule bug.**

*(`penB` is standard nomenclature for the porB1b G120/A121 porin mutations — that identity is my
**[inferred]** reading, not stated anywhere in this repo. Supporting it: the docstring names three
chromosomal markers and the screen's top three discriminators for this exact cell are `rpsJ_V57M`, `mtrR`,
`porB1b_G120K` — two exact name matches plus one positional.)*

**2. The data confirms the ceiling empirically.** Among this cell's 23 missed-R and 26 S:

| determinant | missed-R | S | reading |
|---|---|---|---|
| `mtrR` | 23/23 | **26/26** | completely non-discriminating |
| `rpsJ_V57M` | 23/23 | **10/26** | the v0.1 over-call, measured |
| `porB1b_G120K` | 14/23 | 1/26 | most discriminating, but catches only 61% of the missed |

**3. Counting it would repeat a measured failure, on the tuning set.** The v0.1 rule promoted the
chromosomal markers to primary and scored **spec 0.0** on this cohort — an all-R over-call. v0.2 narrowed
to plasmid `tet(M)` for that reason. And proposing `porB1b` now would be tuning on the **same NCBI-PD
cohort v0.2 was validated against** — circular by construction.

## Clonal-confound check (run before any of the above was believed)

The atlas records that a burden-based *rescue* of this ceiling died as a clonal artifact (pooled +6.2
collapsing to −1.0 within SNP clusters), so the enrichment was checked against the cohort's SNP clusters:

- tetracycline: the 23 missed-R span **18 distinct SNP clusters** (largest 3/23); the 14 token carriers
  span 10 clusters (largest 3/14). The enrichment is **not** one clone.
- cefixime: 4 missed-R across 4 clusters — too few to interpret either way.

So the tetracycline signal is real co-occurrence, not clonality. It is still not an actionable gap, for
the three reasons above. **Caveat on this check:** the cohort's `NULL`-clustered isolates are scored as
unique singletons by the scorer but were lumped into one bucket here, so the within-cluster contrast is
weaker than the across-cluster spread; the spread (18 clusters) is the load-bearing number.

## What this licenses

- `dna_decode/amr/uncounted.py::_MEASURED_GAPS` **stays at one row**. That is now an evidence-backed
  scope statement, not an unexamined default.
- A high missed-R fraction is **expected and honest** for determinant-invisible mechanisms (gono
  azithromycin efflux; gono tetracycline multi-locus-cumulative). The right response there is the atlas's:
  disclose the blindness, do not manufacture a call.
- **DESCRIPTIVE, hypothesis-generating only.** Nothing from this screen may enter `_MEASURED_GAPS` without
  its own validation on a cohort it was not derived from.

Artifacts: `scripts/rule_gap_screen.py` · `wiki/rule_gap_screen_2026-08-25.json`.
Frozen AMR surface byte-unchanged (READ-only screen).
