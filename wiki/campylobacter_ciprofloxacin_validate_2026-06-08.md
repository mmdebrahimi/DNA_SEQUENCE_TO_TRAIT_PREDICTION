# Campylobacter ciprofloxacin — cross-organism validation — 2026-06-08

> Deployed dna-amr ciprofloxacin rule applied UNCHANGED. AMRFinder `-O Campylobacter`.
- NCBI group `Campylobacter`; cohort 30 (15R/15S), 30 runs; `ncbi/amr:4.2.7-2026-03-24.1`
- First organism in a DIFFERENT PHYLUM (Campylobacterota / Epsilonproteobacteria) — all prior
  validation was gram-negative Pseudomonadota (Enterobacterales + Acinetobacter).

## VERDICT: FAILS_BAR (the THRESHOLD doesn't transfer — but the MECHANISM does, perfectly)

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (cipro, deployed QRDR-POINT threshold=2)** | 30 | **0.500** | **0.000** | 1.000 |
| naive (≥1 QRDR/quinolone determinant) | 30 | **1.000** | 1.000 | 1.000 |
| **candidate: QRDR-POINT threshold=1** | 30 | **1.000** | 1.000 | 1.000 |

The deployed rule calls **every strain S** (sens 0.0) → no better than always-S. This is the **first
FAILS_BAR caused by a TUNING parameter, not a mechanism gap or label noise.**

## Root cause — single-mutation resistance vs an E. coli-tuned threshold

The cipro rule counts **QRDR target-POINT mutations** (gyrA/parC/parE) and requires **≥2** (threshold 2).
That threshold was validated on E. coli (N=147, acc 0.925), where **high-level FQ-R typically needs a
gyrA + parC double mutation** — single-mutants are often low-level/borderline, so threshold 2 raises
precision there.

**Campylobacter is different biology:** a **single gyrA T86I** confers full clinical fluoroquinolone
resistance (it has no parC-equivalent contribution in the same way; T86I alone is the dominant,
well-documented C. jejuni/coli FQ-R mechanism). The per-strain structure is perfectly clean:

| label | n_QRDR_point distribution | determinant |
|---|---|---|
| R (15) | **all 15 have exactly 1** | `gyrA_T86I` (100%) |
| S (15) | **all 15 have 0** | — |

Threshold sweep on this cohort:
- **threshold=1 → acc 1.000, sens 1.000, spec 1.000** (TP15 TN15 FP0 FN0) — perfect separation.
- threshold=2 (deployed) → acc 0.500, sens 0.000, spec 1.000 (TP0 TN15 FP0 **FN15**).

So the **mechanism caller transfers perfectly across the phylum boundary** — gyrA T86I in Campylobacter
is the exact analog of gyrA S83L in E. coli, and AMRFinder `-O Campylobacter` reports it cleanly. Only
the **count threshold** is organism-specific.

## Candidate refinement — per-organism QRDR-POINT threshold

Make the cipro threshold organism-aware: **2 for Enterobacterales (E. coli/Klebsiella/Salmonella), 1 for
Campylobacter.** Campylobacter threshold=1 is **biologically grounded** (single gyrA T86I = clinical FQ-R
is textbook), NOT merely cohort-fit — distinct from the Acinetobacter strength-tier candidate, which
needs external validation. This is the cleanest transfer-boundary the wider-AMR thread has found:
- The Acinetobacter case = a CONTENT gap (intrinsic OXA-51 over-counted) needing organism curation.
- The Campylobacter case = a TUNING gap (one threshold), with a principled per-organism value.

**NOT wired into the deployed rule yet** — DRUG_RULE currently has one threshold per drug. Parameterizing
it by organism is a scoped follow-on (touches the E. coli/Klebsiella-validated numbers if done wrong; add
regression pins for those before changing the default). Documented here as the validated candidate.

## Honest scope / caveats
- 1 organism, 1 drug, N=30, NCBI labels (different source/curation, not a different-lab study).
- The threshold=1 number is in-cohort, but the single-gyrA-T86I → FQ-R biology is independently
  established, so confidence the refinement generalizes is HIGH (higher than the Acinetobacter candidate).
- naive ≥1-determinant also scores 1.000 here **because** the only quinolone determinant AMRFinder reports
  for these strains IS the QRDR point mutation — there are no confounding acquired/efflux quinolone genes
  in this Campylobacter cohort (unlike Enterobacterales, where naive ≥1 over-calls via intrinsic efflux).

## Takeaway — updated transferability map
| organism (phylum) | drug | deployed rule | transfers? | boundary type |
|---|---|---|---|---|
| E. coli / Klebsiella / Enterobacter / Salmonella (Pseudomonadota) | cipro/cef/tet/gent | class/QRDR rules | **yes** | — |
| Acinetobacter (Pseudomonadota) | meropenem | CARBAPENEM-class | **no** | CONTENT (intrinsic OXA-51); strength-tier candidate |
| **Campylobacter (Campylobacterota)** | **cipro** | **QRDR-POINT thr=2** | **no** | **TUNING (single-mutation R); thr=1 candidate (perfect, biologically grounded)** |

The mechanism abstraction (QRDR target-POINT mutation) is phylum-robust; the resistance **threshold** is
not — it encodes organism-specific genetics of how many hits = clinical resistance.
