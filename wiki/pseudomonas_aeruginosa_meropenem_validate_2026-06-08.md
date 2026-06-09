# Pseudomonas_aeruginosa meropenem — cross-organism validation — 2026-06-08

> Deployed dna-amr meropenem rule applied UNCHANGED. AMRFinder `-O Pseudomonas_aeruginosa`.
- NCBI group `Pseudomonas_aeruginosa`; cohort 30 (15R/15S), 30 runs; `ncbi/amr:4.2.7-2026-03-24.1`

## VERDICT: FAILS_BAR (CONTENT boundary — 2nd carbapenem over-call, like Acinetobacter)

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (meropenem, CARBAPENEM-class rule)** | 30 | **0.500** | 1.000 | **0.000** |
| naive (≥1 carbapenem determinant) | 30 | 0.500 | 1.000 | 0.000 |

All 15 S strains called R (FP=15, FN=0) — pure all-R over-call, the same shape as Acinetobacter × meropenem.

## Root cause — efflux-regulator + porin point-mutations counted as determinants

The over-call is driven by determinants AMRFinder reports for *P. aeruginosa* that are **present in
susceptible strains too** (sampled from the directly-inspectable subset):

- **nalC** point mutations (nalC_G71E, nalC_S209R) — nalC is the **repressor of the MexAB-OprM efflux
  pump**; these variants are common and only *conditionally* up-regulate efflux. Counting them = over-call.
- **oprD** variants (oprD_V359L, oprD_G151PfsTer39) — oprD is the **carbapenem porin**; loss/frameshift
  *can* confer R, but missense variants are common in S strains and a frameshift alone isn't always
  sufficient. Counting all oprD variation = over-call.
- (one true acquired carbapenemase, blaVIM-11, also present — genuinely R.)

So unlike Acinetobacter (intrinsic *gene presence*, OXA-51), the Pseudomonas over-call is from **intrinsic
regulatory + porin POINT MUTATIONS** that AMRFinder flags under a carbapenem-relevant class but that don't
reliably confer the phenotype. Same CONTENT failure mode (counting determinants that don't predict R),
different determinant type.

## Why it's hard — the deeper irony

*P. aeruginosa* meropenem-R is genuinely **multifactorial**: oprD loss + MexAB-OprM/MexXY efflux
overexpression + AmpC (blaPDC) derepression ± acquired carbapenemase, often in combination, all
**expression/regulation-driven**. Gene-presence (even point-mutation-aware) cannot weigh these — it sees
the variants but not their net regulatory effect. This is the **CONTENT over-call AND the EXPRESSION
blind-spot at once**: the determinants it counts are conditional, and the ones that matter are expression
levels it can't read. P. aeruginosa carbapenem is arguably the hardest gene-presence target seen so far.

## Boundary-type map (updated — 4 failing organisms, 3 flavors)

| organism | drug | failure | boundary type |
|---|---|---|---|
| Klebsiella / E. coli / Salmonella(?) | cipro/cef/tet/gent | — | transfers (Enterobacterales) |
| Acinetobacter | meropenem | spec→0 (intrinsic OXA-51 gene presence) | CONTENT |
| **Pseudomonas** | **meropenem** | **spec→0 (nalC/oprD conditional variants)** | **CONTENT (+EXPRESSION)** |
| Campylobacter | cipro | sens→0 (single gyrA T86I vs thr 2) | TUNING |
| Enterobacter cloacae | ceftriaxone | sens→0.375 (derepressed AmpC) | EXPRESSION |

## Honest scope / caveats
- 1 organism, 1 drug, N=30, NCBI labels (different source/curation, not a different-lab study).
- Per-strain determinant identities sampled from the directly-inspectable subset; the all-R verdict
  (FP=15) is the validator's authoritative full-30 result.
- No refinement attempted — a P. aeruginosa carbapenem rule would need expression-level inference (efflux/
  porin regulation), which gene-presence fundamentally cannot do. This is a NEGATIVE result for the
  presence-based approach on this organism×drug, honestly recorded.
