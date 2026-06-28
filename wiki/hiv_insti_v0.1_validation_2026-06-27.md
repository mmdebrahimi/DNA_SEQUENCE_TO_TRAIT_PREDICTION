# HIV INSTI v0.1 - mutant-specific catalog (data-derived, CV held-out) (2026-06-27)

Method: data-derived resistant mutants = MULTIVARIATE OLS log10-fold coefficient >= log10(1.5) (independent >=1.5x effect, controlling for co-occurrence -> deconfounds accessory/polymorphic riders), >=5 carriers; 5-fold CROSS-VALIDATED held-out (derive on train, eval on test) -> NOT in-sample.

**Cutoff (delta-honest):** DELTA-HONEST: both v0 (position-based) and v0.1 (mutant-specific) scored at the SAME UNIFORM illustrative fold>=3 (no per-drug INSTI clinical cutoff sourced in-repo); the headline is the v0->v0.1 balanced-accuracy GAIN at that fixed cutoff, NOT an absolute-calibration claim.

Label = Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra).

Position-based v0 over-calls accessory/polymorphic riders at the major integrase positions; v0.1 keeps only fold-associated mutants. Both metrics are 5-fold CV held-out (out-of-sample). **5/5 INSTI drugs improve-or-hold; mean balacc gain 0.087.**

| Drug | n | prev R | v0 pos-based sens/spec/**balacc** | v0.1 mutant sens/spec/**balacc** | balacc gain |
|---|---|---|---|---|---|
| raltegravir | 753 | 0.413 | 0.939/0.814/**0.877** | 0.923/0.95/**0.937** | 0.06 |
| elvitegravir | 754 | 0.5 | 0.87/0.857/**0.863** | 0.828/0.923/**0.875** | 0.012 |
| dolutegravir | 370 | 0.195 | 1.0/0.534/**0.767** | 1.0/0.664/**0.832** | 0.065 |
| bictegravir | 287 | 0.178 | 1.0/0.419/**0.71** | 0.941/0.547/**0.744** | 0.034 |
| cabotegravir | 64 | 0.812 | 1.0/0.333/**0.667** | 0.865/1.0/**0.933** | 0.266 |

## Deliverable catalog (derived on all data)
- **raltegravir** (15): E92Q, G118R, G140A, G140C, G140S, N155H, Q148H, Q148K, Q148R, Y143A, Y143C, Y143G, Y143H, Y143R, Y143S
- **elvitegravir** (18): E92G, E92Q, G118R, G140A, G140C, G140S, N155H, Q148H, Q148K, Q148R, R263K, S147G, T66I, Y143A, Y143C, Y143G, Y143R, Y143S
- **dolutegravir** (12): E138K, E92Q, G118R, G140A, G140S, N155H, Q148H, Q148K, Q148R, R263K, S147G, Y143R
- **bictegravir** (13): E138A, E138K, E92Q, G118R, G140A, G140C, G140S, N155H, Q148K, Q148R, R263K, S147G, Y143R
- **cabotegravir** (5): E138K, G140A, G140S, Q148H, Q148R

## Honest caveats
- DELTA-HONEST at the uniform fold>=3 cutoff (no per-drug INSTI clinical breakpoint sourced); the GAIN is the signal, not the absolute balacc
- deconfounded (multivariate-OLS-coefficient) derivation, mirroring the shipped NRTI/PI v0.1 arc; 5-fold CV held-out -> out-of-sample, not optimistic in-sample
- in-distribution vs HIVDB-PhenoSense (NOT provenance-disjoint external validation)
- INSTI headroom is MIXED, not uniform like PI: RAL/DTG have positive OLS-vs-catalog deltas, EVG/BIC are ~flat (position call near-optimal), CAB's large delta is at n=64 (UNSTABLE). v0.1 is expected to improve-or-hold positive-delta drugs + hold the flat ones; a per-drug result.
- cabotegravir n is small (~64) -> its v0.1 number is unstable; treat as indicative only

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; method Stanford DRMcv.R.