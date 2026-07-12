# non-determinant → determinant cross-axis — lineage de-confound (leave-one-clade-out) (2026-07-12)

**Verdict: CROSS_AXIS_GENERALIZES_BEYOND_LINEAGE** — median clade-grouped AUC = 0.922 (vs naive 0.97); 40/45 AMR determinants still predicted at AUC >= 0.7 when their clade is held out.

Mash: 240 E. coli genomes → 85 clades at threshold 0.005 (largest clade 0.342).

Does the E. coli non-determinant -> determinant cross-axis signal survive leave-one-clade-out CV (real beyond lineage), or collapse (lineage-mediated)?

| AMR determinant | n | AUC naive | **AUC clade-grouped** | drop | generalizes |
|---|---|---|---|---|---|
| mrx(A) | 98 | 0.996 | **0.997** | -0.0 | YES |
| mph(A) | 98 | 0.993 | **0.994** | -0.0 | YES |
| cmlA1 | 13 | 0.996 | **0.992** | 0.005 | YES |
| aadA5 | 90 | 0.997 | **0.991** | 0.005 | YES |
| aac(6')-Ib-cr5 | 69 | 0.994 | **0.989** | 0.005 | YES |
| aph(3'')-Ib | 100 | 0.999 | **0.982** | 0.016 | YES |
| gyrA_D87N | 140 | 0.997 | **0.98** | 0.017 | YES |
| blaOXA-1 | 69 | 0.987 | **0.979** | 0.008 | YES |
| parC_S80I | 141 | 0.998 | **0.979** | 0.019 | YES |
| dfrA17 | 94 | 0.988 | **0.978** | 0.01 | YES |
| dfrA12 | 27 | 0.974 | **0.97** | 0.004 | YES |
| gyrA_S83L | 159 | 0.964 | **0.964** | 0.0 | YES |
| catB3 | 69 | 0.988 | **0.963** | 0.025 | YES |
| ptsI_V25I | 94 | 0.998 | **0.961** | 0.037 | YES |
| aadA2 | 31 | 0.969 | **0.957** | 0.013 | YES |
| sul2 | 102 | 0.98 | **0.955** | 0.026 | YES |
| tet(M) | 11 | 0.99 | **0.952** | 0.039 | YES |
| sat2 | 10 | 0.97 | **0.949** | 0.022 | YES |
| sul1 | 107 | 0.955 | **0.926** | 0.029 | YES |
| aadA1 | 40 | 0.956 | **0.924** | 0.032 | YES |
| estX-3 | 12 | 0.986 | **0.924** | 0.062 | YES |
| sul3 | 25 | 0.996 | **0.924** | 0.072 | YES |
| aph(6)-Id | 86 | 0.982 | **0.92** | 0.062 | YES |
| blaCTX-M-55 | 15 | 0.86 | **0.917** | -0.057 | YES |
| uhpT_E350Q | 111 | 0.969 | **0.913** | 0.055 | YES |
| blaCTX-M-15 | 98 | 0.957 | **0.905** | 0.051 | YES |
| blaCTX-M-27 | 22 | 0.933 | **0.897** | 0.036 | YES |
| aac(3)-IIe | 54 | 0.973 | **0.881** | 0.092 | YES |
| parE_I529L | 88 | 0.997 | **0.848** | 0.148 | YES |
| dfrA14 | 27 | 0.917 | **0.843** | 0.073 | YES |
| qnrS1 | 15 | 0.816 | **0.825** | -0.009 | YES |
| catA1 | 23 | 0.89 | **0.811** | 0.079 | YES |
| blaCMY-2 | 15 | 0.805 | **0.784** | 0.021 | YES |
| dfrA1 | 15 | 0.941 | **0.758** | 0.183 | YES |
| blaTEM-1 | 80 | 0.862 | **0.741** | 0.121 | YES |
| aph(3')-Ia | 23 | 0.935 | **0.737** | 0.198 | YES |
| aac(3)-IId | 17 | 0.899 | **0.734** | 0.164 | YES |
| cyaA_S352T | 31 | 0.893 | **0.727** | 0.166 | YES |
| floR | 19 | 0.829 | **0.715** | 0.113 | YES |
| parE_S458A | 33 | 0.978 | **0.711** | 0.267 | YES |
| tet(B) | 57 | 0.872 | **0.656** | 0.216 | no (lineage) |
| blaTEM | 10 | 0.692 | **0.592** | 0.1 | no (lineage) |
| tet(A) | 109 | 0.851 | **0.59** | 0.261 | no (lineage) |
| parE_L416F | 11 | 0.97 | **0.198** | 0.772 | no (lineage) |
| parC_E84V | 77 | 1.0 | **clade-concentrated** | None | no (clade-concentrated) |

## Literature / mechanism
- The C-core finding (PASS_LINKAGE_STRUCTURE) showed AMR determinants co-occur in cassette blocks (integron sul1->aadA/dfrA, QRDR clusters) under a determinant-profile-DEDUP clonality proxy. This applies the STRONGER phylogenetic (Mash-clade leave-one-out) control: an integron cassette is itself a MOBILE unit, so co-cassette determinants should predict each other even across held-out clades (real cassette linkage); a collapse would mean the co-occurrence is clade-signature, not cassette-intrinsic.

## Honest caveats
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- E. coli/Shigella only (the determinant axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE of the target feature (a common feature present across many clades has cross-clade signal to learn; a clade-restricted accessory feature does not).