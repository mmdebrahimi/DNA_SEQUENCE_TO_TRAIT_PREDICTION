# Staphylococcus aureus oxacillin (MRSA/mecA) — 1st Gram-positive — 2026-06-07

> Roadmap Phase 3, Gram-positive generality test. Does the deterministic acquired-gene approach transfer
> beyond gram-negatives? Substrate: S. aureus oxacillin (mecA = MRSA, the signature Gram-positive
> mechanism). Cohort: 30 NCBI Pathogen Detection (15R/15S), `-O Staphylococcus_aureus`.
> NOTE: curated finding (the script's auto-template reports the headline only).

## VERDICT: GENOTYPE TRANSFERS (sens 1.0) — but NOT validatable against oxacillin labels (spec 0.333)

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr oxacillin (mecA, METHICILLIN-subclass)** | 30 | 0.667 | **1.000** | 0.333 |
| naive AMRFinder (any beta-lactam/methicillin determinant) | 30 | 0.600 | 1.0 | 0.2 |

## The finding: the rule is sound; the oxacillin LABEL is the weak axis

- **All 15 R strains carry mecA → sens 1.0.** mecA detection transfers to a Gram-positive perfectly; the
  acquired-gene + Subclass-refinement approach (same shape as cef/gent/meropenem) works on S. aureus.
- **10 of 15 "oxacillin-susceptible"-labeled strains ALSO carry full-length mecA** (EXACTX, 100% id) → the
  false positives → spec 0.333. **67% mecA-carriage among susceptible-labeled strains is implausible for
  genuine OS-MRSA** (oxacillin-susceptible mecA+ S. aureus is typically <5%).
- **Root cause = oxacillin AST label unreliability, NOT a rule defect.** Oxacillin direct susceptibility
  testing is method-dependent and well-known to misclassify mecA+ strains as susceptible — which is exactly
  why **CLSI/EUCAST recommend CEFOXITIN as the mecA surrogate**, not oxacillin. The genotype (mecA+) is the
  more reliable indicator of clinical methicillin resistance than the oxacillin phenotype label in this cohort.

## Interpretation (for the long-term vision)

This is a **label-provenance finding** — the same binding constraint the project has hit repeatedly
(pathotype sampling confound; AMR label noise). The deterministic method's Gram-positive genotype transfer
is GOOD (sens 1.0, mecA correctly detected); it cannot be *scored against oxacillin labels* because the
oxacillin phenotype is a noisy comparator for mecA. The honest scientific verdict is: **mecA-based MRSA
detection generalizes to Gram-positive; the oxacillin AST label does not validate it — use cefoxitin.**

## Disposition

- oxacillin support SHIPPED (mic_tiers + DRUG_RULE: mecA / METHICILLIN-subclass, threshold 1). mecA
  detection is sound + clinically the correct MRSA call; the rule is kept.
- Verdict on the *oxacillin label* validation: NOT VALIDATED (label-confounded), distinct from a rule
  failure. Honest, not papered over.
- **Follow-up attempted + BLOCKED by substrate:** censused cefoxitin (the proper CLSI/EUCAST mecA
  surrogate) on the same NCBI dataset → only **3 cefoxitin-R strains** (cannot build a balanced cohort).
  So the proper-label validation is **substrate-blocked** on this data source. Terminal: mecA genotype
  detection transfers to Gram-positive (sens 1.0); phenotype-label validation is the limit (oxacillin
  noisy, cefoxitin sparse) — the project's recurring "substrate/label is the binding constraint" lesson.

## Honest scope

1 organism, 1 drug, N=30, NCBI Pathogen Detection labels. The label-unreliability conclusion is grounded in
the mecA-carriage rate among susceptible-labeled strains (10/15) + the established cefoxitin-surrogate
guidance, not just the accuracy number.
