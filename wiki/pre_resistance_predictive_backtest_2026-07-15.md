# Pre-resistance PREDICTIVE backtest — the real test (2026-07-15)

The `/idea-validation-council` named ONE test as decisive for the escape/pre-resistance forecaster idea:
*"did one-nt-away flags computed on catalogue v(N) predict the DRMs added in v(N+1), above the base rate of
all adjacent codons? This is the only thing separating a codon-table lookup from a forecaster."* This runs it.

**Substrate:** the WHO *M. tuberculosis* mutation catalogue, v1 (2021) → v2 (2023), with the version diff
ENCODED in the committed v2 master file (`CHANGES vs ver1`: `New AwR`/`New AwRI` + `UP from Uncertain to
AwR|AwRI` = mutations newly graded resistance in v2; the v1-R set = the seed the forecaster "knew" at v1).
Deterministic, offline (no network), no GPU. Frozen decoder surface byte-unchanged (`verify_lock OK`).
`scripts/pre_resistance_predictive_backtest.py`. Scope: protein single-substitutions (codon concept defined);
non-coding/rRNA excluded.

## Result

| metric | value | reading |
|---|---:|---|
| protein substitutions parsed | 23,065 | catalogue coverage |
| v1-R positions (seed) | 145 | known resistance codons at v1 |
| v2-added-R protein substitutions | 97 | the prediction target |
| **CEILING** — v2-additions at a known v1-R position | **43/97 = 44%** | the MAX a one-nt-from-known-codon forecaster could catch |
| **base rate** (background one-nt-adjacency at v1-R positions) | **0.558** | one-nt-away is the MAJORITY state already |
| test rate (v2-added one-nt-adjacency) | 0.721 | 31/43 |
| **enrichment ratio** | **1.29×** | weak |
| Fisher p (one-sided, position-matched) | **0.022** | real but small |
| **verdict** (material = ≥2× AND p<0.05) | **FAIL** | |

## The finding — this confirms the council's dominant read

**The one-nt-away flag does NOT materially predict catalogue growth.** Two independent reasons, both decisive:

1. **Structural ceiling (56% miss).** More than half of v2-added DRMs occur at *new* positions/mechanisms
   (rpoB→new gene, a promoter, a new codon). A "one nt from a KNOWN resistance codon" flag cannot fire there
   by construction — so the forecaster is *blind to the majority of real catalogue growth*, regardless of
   codon resolution. This is the exact Contrarian/First-Principles prediction: "new DRMs appear at positions/
   mechanisms, not just adjacent codons."

2. **Near-universal base rate (weak discrimination).** Even restricted to known resistance codons — the
   forecaster's home turf — the base rate of one-nt-adjacency is already **0.558**. Being flagged is the
   *majority* state, so it carries little information. v2-additions are enriched to 0.721 (1.29×, p=0.022):
   a **statistically real but operationally weak** signal, far below the ≥2× material bar. The p=0.022 says
   new alleles *do* cluster slightly at reachable neighbors (real biology); the 1.29× says that clustering is
   nowhere near strong enough to make "one nt away" a forecast rather than a codon-table restatement.

## Honest scope

- **Robust to the aa-vs-codon caveat.** Adjacency is aa-level (generous upper bound; census method). Codon-
  exact (H37Rv reference codons) would only *tighten* the flagged set — it will not lift a 1.29× to ≥2×, and
  the 56% ceiling is codon-resolution-independent. The negative is robust.
- **Pre-registered bar.** "Material = ≥2× enrichment AND Fisher p<0.05" was set in the run framing before
  results; 1.29× is unambiguously below it. Even discarding the threshold, the 56% structural miss alone
  disqualifies a *forecaster* claim.
- **Not zero signal — honestly weak signal.** p=0.022 is reported transparently; the flag is not noise, it is
  a real-but-immaterial 29% relative lift on top of a majority base rate.

## Combined verdict across BOTH council-named tests

| test | question | result |
|---|---|---|
| **1 census** (`pre_resistance_base_rate_census`) | is the substrate non-empty (do genomes carry intermediates)? | **GO** — 231 real isolates one nt from resistance |
| **2 backtest** (this) | does one-nt-away *materially predict* new DRMs above base rate? | **FAIL** — 1.29× weak, 56% structural miss |

**The decisive test (2) fails.** The substrate exists (things *are* one nt away — test 1), but being one nt
away has **no material forecasting power** for what actually becomes resistance (test 2). Combined with the
council's untouched novelty (genetic-barrier = 15-yr published metric; "pre-resistance" = named 2021 paper)
and actionability objections, the recommendation is clear: **do NOT build the pre-resistance forecaster as a
predictive tool** — the predictive claim, the one thing that would separate it from a codon-table lookup,
does not hold empirically. A cheap, decision-grade negative (the north star's "failure-tolerant iteration").

Artifact: `wiki/pre_resistance_predictive_backtest_2026-07-15.json`. Run:
`uv run python scripts/pre_resistance_predictive_backtest.py`.
