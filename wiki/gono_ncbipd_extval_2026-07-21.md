# N. gonorrhoeae external validation on NCBI-PD — cipro + cefixime v0.1 CONFIRMED (2026-07-21)

**Status:** ✅ **ciprofloxacin + cefixime externally CONFIRMED** on 160+ independent isolates (the fresh
AR-Bank cells + the cefixime v0.1 rule generalize) · 4 other gonococcal cells honestly characterized ·
**Cell:** `neisseria_amr` (NON-FROZEN) · **Frozen surface:** byte-unchanged.

## Why this run

The AR-Bank Neisseria validation (2026-07-21) endorsed cipro + cefixime on **N=20**, but the cefixime v0.1
mosaic-penA-34 rule was **derived + validated on that same 20-isolate cohort** — an honest tension (a
cohort-fit rule is not yet externally confirmed). This scores the SAME cell on **170 independent NCBI
Pathogen Detection isolates**, provenance-disjoint from the AR-Bank set (the 21 AR-Bank BioSamples excluded).

**Feasibility win — no compute wall.** NCBI-PD publishes its own AMRFinderPlus calls in the `AMR_genotypes`
field *including point mutations* (`gyrA_S91F=POINT`, `penA_I312M=POINT`), and this cohort's isolates have
downloadable assemblies (NCBI ran AMRFinder for us). So this is a pure **metadata join + score**, offline,
in seconds — no Kaggle, no Docker, no assembly download. The penA numbering was verified compatible between
the two AMRFinder sources ({I312M,V316T,N512Y,G545S} present in both; `penA_G542S` is a separate real
allele, not a renumbered G545S).

## Result (170-isolate provenance-disjoint cohort; powers 5 drugs)

| drug | n | R/S | sens | spec | acc | verdict |
|---|---|---|---|---|---|---|
| **ciprofloxacin** | 163 | 94R/69S | **0.989** | **0.986** | **0.988** | **SCORED_ENDORSED** ✅ |
| **cefixime** (v0.1) | 166 | 19R/147S | 0.789 | **0.905** | 0.892 | **SCORED_ENDORSED** ✅ |
| ceftriaxone | 169 | 2R/167S | — (2R) | 0.946 | 0.935 | UNDERPOWERED (spec confirmed) |
| azithromycin | 156 | 110R/46S | 0.0 | 1.0 | 0.295 | **DEGENERATE** (data gap — see below) |
| penicillin | 31 | 14R/17S | 1.0 | 0.0 | 0.452 | **DEGENERATE** (over-call) |
| tetracycline | 60 | 34R/26S | 1.0 | 0.0 | 0.567 | **DEGENERATE** (over-call) |

## The two confirmations (the headline)

- **ciprofloxacin — externally CONFIRMED, acc 0.988** (163 isolates). The gyrA/parC QRDR rule is a
  *literature* rule (not fit to any cohort); its near-perfect score on 163 independent isolates is a genuine
  strong external confirmation (N=20 AR-Bank → N=163 NCBI-PD).
- **cefixime v0.1 — externally CONFIRMED, spec 0.905** (147 independent S isolates). The mosaic-penA-34 rule
  I *derived* on the 20 AR-Bank isolates **generalizes**: it holds specificity ≥0.85 on 147 unseen S
  isolates. This directly closes the honest tension — the rule is now cohort-derived **AND** externally
  confirmed. (sens 0.789 = 15/19 R, slightly below the AR-Bank 0.917 on independent data — honest.)

## The four honest negatives (external validation doing its job)

The other 4 gonococcal cells were NEVER AR-Bank-validated (that run scored only cipro + cefixime); this run
is their first test, and it exposes that they are **not reliable** on this substrate:

- **azithromycin — DEGENERATE (sens 0.0): a DATA GAP, not a cell bug.** The cell keys on 23S rRNA mutations
  (A2045G/C2611T…), but NCBI-PD's published AMRFinder calls **do not include 23S rRNA point mutations at all**
  (0/110 azithro-R isolates carry any 23S marker). The determinant the cell needs is absent from the input.
  Validating azithromycin requires a 23S-aware caller on the assemblies (a separate build) — not fixable
  from NCBI-PD's published calls.
- **penicillin + tetracycline — DEGENERATE (spec 0.0, all-R over-call).** Both cells' own docstrings flagged
  an over-call risk (penicillin: chromosomal penA/mtrR promoted → near-universal; tetracycline: rpsJ V57M
  near-universal). This run **confirms the over-call on independent data** — spec 0.0. They need
  acquired-marker narrowing (like the cefixime v0.1 fix) before they can be endorsed.
- **ceftriaxone — spec confirmed (0.946), sens underpowered (2R).** The A501-specific v0.1 rule holds
  specificity on 167 S isolates; sens is untestable (only 2 R, both non-A501).

## Scorer hardening (verify-in-batch)

The first run mislabeled azithromycin `SCORED_ENDORSED` — a cell predicting all-S gets spec 1.0 + passes a
spec-only endorsement check while being useless (sens 0.0). Added a **DEGENERATE guard**: a cell that
predicts all-one-class (sens 0 with powered R, or spec 0 with powered S) is `DEGENERATE_NOT_ENDORSED`
regardless of the other metric. Same class of integrity fix as the empty-assembly-as-S gate — a degenerate
output must never read as a validated one.

## Honest caveats

- **Provenance-disjoint** (different isolates than AR-Bank; the 21 AR-Bank BioSamples excluded) but **NOT
  methodology-independent** (same AMRFinderPlus + same `neisseria_amr` cell) — the standard caveat.
- **RAW sens/spec is CLONALITY-INFLATED** (gonococci are clonally structured). This is a provenance-disjoint
  stress test; the lineage-collapsed number (Mash on the assemblies) is the follow-up and would be lower.
- NON-FROZEN cell; the frozen decoder surface is byte-unchanged (`verify_lock` OK).
