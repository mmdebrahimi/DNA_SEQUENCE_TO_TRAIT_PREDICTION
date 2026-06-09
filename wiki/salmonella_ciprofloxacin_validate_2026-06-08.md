# Salmonella ciprofloxacin — cross-organism validation — 2026-06-08

> Deployed dna-amr ciprofloxacin rule applied UNCHANGED. AMRFinder `-O Salmonella`.
- NCBI group `Salmonella`; cohort 30 (15R/15S), 30 runs; `ncbi/amr:4.2.7-2026-03-24.1`

## VERDICT: FAILS_BAR (COMBINED TUNING + CONTENT — and the same design choice that wins on Klebsiella)

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (cipro, QRDR-POINT counter, threshold 2)** | 30 | **0.567** | **0.133** | 1.000 |
| naive (≥1 broad quinolone determinant, threshold 1) | 30 | **1.000** | 1.000 | 1.000 |

The deployed rule catches only 2 of 15 R (sens 0.133, FN=13). naive is perfect.

## Root cause — TWO mechanisms the QRDR-point-threshold-2 rule misses

Per-strain dissection of the 15 R:

| QRDR-point count | # R strains | mechanism | caught by rule (≥2)? |
|---:|---:|---|---|
| 3 | 2 | gyrA + parC double/triple | **yes** (the only 2 caught) |
| 1 | 7 | single gyrA (S83F / D87Y) | no — **TUNING** (threshold too high) |
| 0 | 6 | plasmid **qnr** (qnrB19/S1/S2), no QRDR point | no — **CONTENT** (counter excludes qnr) |

Dominant determinants in R: **qnrB19 (10 strains)**, gyrA_D87Y (5), gyrA_S83F (4), qnrS2/S1 (4), parC_S80I (2).

1. **TUNING** (7 strains): single gyrA S83F/D87Y confers clinically relevant cipro/low-level-FQ resistance
   in Salmonella, but the E. coli-tuned threshold=2 (high-level needs gyrA+parC) misses them — exactly the
   Campylobacter pattern.
2. **CONTENT** (6 strains): plasmid-mediated **qnr** genes are a major Salmonella FQ-R mechanism, but the
   cipro rule uses the **QRDR-point-only counter** (which counts gyrA/parC/parE point mutations and nothing
   else) — so qnr is invisible. These 6 R strains have zero QRDR point mutations and are called S.

## The sharp lesson — the same design choice wins on Klebsiella and loses on Salmonella

The cipro rule deliberately uses a **QRDR-point-only counter** instead of broad quinolone-determinant
counting. WHY: on Klebsiella that choice EXCLUDES the intrinsic OqxAB efflux that broad counting sweeps in
→ Klebsiella cipro acc **1.000** (naive 0.500). The exact same choice on Salmonella EXCLUDES the legitimate
qnr signal → Salmonella acc **0.567** (naive 1.000). **The counter choice is itself organism-specific:**
QRDR-point-only is right where the broad-class noise is intrinsic efflux (Klebsiella), wrong where the
broad class carries real acquired determinants (Salmonella qnr). There is no single counter+threshold that
is correct for both — organism-specific calibration is unavoidable.

## Strengthens the calibrate_organism hypothesis (H3)

naive (broad counter, threshold 1) = acc 1.000 on Salmonella. An auto-calibration over (counter, threshold)
on a ≥15R/15S Salmonella cohort would recover broad-counter-threshold-1 here and QRDR-point-threshold-2 on
Klebsiella — exactly what `calibrate_organism(cohort)` (ledger H3) must do. This is now the 3rd organism
(after Campylobacter thr=1, Klebsiella thr=2) showing auto-calibration would pick the right config, AND the
first showing the COUNTER (not just the threshold) needs to be organism-selected.

## Boundary-type map (updated)

| organism | drug | failure | boundary |
|---|---|---|---|
| Klebsiella / E. coli | cipro/cef/tet/gent | — | transfers |
| Acinetobacter / Pseudomonas | meropenem | spec→0 over-call | CONTENT |
| Campylobacter | cipro | sens→0 (single gyrA vs thr 2) | TUNING |
| Enterobacter cloacae | ceftriaxone | sens→0.375 (derepressed AmpC) | EXPRESSION |
| **Salmonella** | **cipro** | **sens→0.133 (single gyrA + qnr)** | **TUNING + CONTENT** |

## Honest scope / caveats
- 1 organism, 1 drug, N=30, NCBI labels (different source/curation, not a different-lab study).
- The mechanism split (7 single-gyrA / 6 qnr / 2 double) is the robust finding; exact acc is N-bounded.
- No refinement wired — documented as a calibrate_organism candidate (broad counter + threshold 1 for
  Salmonella), pending an independent cohort.
