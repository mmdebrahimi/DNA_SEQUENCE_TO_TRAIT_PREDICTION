# non-determinant → determinant cross-axis — lineage de-confound (leave-one-clade-out) (2026-07-12)

**Verdict: CROSS_AXIS_GENERALIZES_BEYOND_LINEAGE** — median clade-grouped AUC = 0.908 (vs naive 0.975); 39/45 AMR determinants still predicted at AUC >= 0.7 when their clade is held out.

Mash: 232 escherichia_coli_shigella genomes → 80 clades at threshold 0.005 (largest clade 0.353).

Does the E. coli non-determinant -> determinant cross-axis signal survive leave-one-clade-out CV (real beyond lineage), or collapse (lineage-mediated)?

| AMR determinant | n | AUC naive | **AUC clade-grouped** | drop | generalizes |
|---|---|---|---|---|---|
| mrx(A) | 98 | 0.997 | **0.997** | 0.0 | YES |
| cmlA1 | 13 | 0.996 | **0.993** | 0.003 | YES |
| mph(A) | 98 | 0.993 | **0.992** | 0.001 | YES |
| aac(6')-Ib-cr5 | 69 | 0.994 | **0.99** | 0.005 | YES |
| aadA5 | 90 | 0.996 | **0.988** | 0.008 | YES |
| blaOXA-1 | 69 | 0.987 | **0.98** | 0.007 | YES |
| aph(3'')-Ib | 100 | 0.999 | **0.978** | 0.02 | YES |
| dfrA17 | 94 | 0.989 | **0.977** | 0.012 | YES |
| sat2 | 10 | 0.981 | **0.977** | 0.004 | YES |
| gyrA_D87N | 140 | 0.997 | **0.976** | 0.02 | YES |
| parC_S80I | 141 | 0.997 | **0.976** | 0.021 | YES |
| tet(M) | 11 | 0.984 | **0.972** | 0.012 | YES |
| dfrA12 | 27 | 0.993 | **0.97** | 0.022 | YES |
| gyrA_S83L | 159 | 0.968 | **0.962** | 0.006 | YES |
| sul2 | 102 | 0.978 | **0.957** | 0.021 | YES |
| aadA2 | 31 | 0.964 | **0.955** | 0.009 | YES |
| catB3 | 69 | 0.984 | **0.953** | 0.031 | YES |
| ptsI_V25I | 94 | 0.997 | **0.94** | 0.058 | YES |
| sul1 | 107 | 0.957 | **0.931** | 0.025 | YES |
| uhpT_E350Q | 111 | 0.968 | **0.915** | 0.053 | YES |
| estX-3 | 12 | 0.982 | **0.91** | 0.072 | YES |
| aadA1 | 40 | 0.952 | **0.909** | 0.043 | YES |
| aph(6)-Id | 86 | 0.981 | **0.907** | 0.074 | YES |
| sul3 | 25 | 0.997 | **0.907** | 0.09 | YES |
| blaCTX-M-15 | 98 | 0.963 | **0.899** | 0.064 | YES |
| blaCTX-M-55 | 15 | 0.854 | **0.898** | -0.044 | YES |
| blaCTX-M-27 | 22 | 0.947 | **0.896** | 0.051 | YES |
| aac(3)-IIe | 54 | 0.975 | **0.868** | 0.107 | YES |
| parE_I529L | 88 | 0.998 | **0.852** | 0.146 | YES |
| qnrS1 | 15 | 0.851 | **0.824** | 0.027 | YES |
| blaCMY-2 | 15 | 0.81 | **0.819** | -0.009 | YES |
| catA1 | 23 | 0.892 | **0.817** | 0.076 | YES |
| dfrA14 | 27 | 0.933 | **0.811** | 0.121 | YES |
| aph(3')-Ia | 23 | 0.922 | **0.809** | 0.113 | YES |
| blaTEM-1 | 80 | 0.855 | **0.744** | 0.111 | YES |
| floR | 19 | 0.84 | **0.736** | 0.105 | YES |
| aac(3)-IId | 17 | 0.896 | **0.72** | 0.176 | YES |
| tet(B) | 57 | 0.883 | **0.715** | 0.168 | YES |
| dfrA1 | 15 | 0.901 | **0.708** | 0.193 | YES |
| parE_S458A | 33 | 0.968 | **0.682** | 0.285 | no (lineage) |
| blaTEM | 10 | 0.745 | **0.605** | 0.14 | no (lineage) |
| tet(A) | 109 | 0.844 | **0.579** | 0.265 | no (lineage) |
| cyaA_S352T | 31 | 0.877 | **0.56** | 0.317 | no (lineage) |
| parE_L416F | 11 | 0.987 | **0.245** | 0.742 | no (lineage) |
| parC_E84V | 77 | 1.0 | **clade-concentrated** | None | no (clade-concentrated) |

## Literature / mechanism
- The C-core finding (PASS_LINKAGE_STRUCTURE) showed AMR determinants co-occur in cassette blocks (integron sul1->aadA/dfrA, QRDR clusters) under a determinant-profile-DEDUP clonality proxy. This applies the STRONGER phylogenetic (Mash-clade leave-one-out) control: an integron cassette is itself a MOBILE unit, so co-cassette determinants should predict each other even across held-out clades (real cassette linkage); a collapse would mean the co-occurrence is clade-signature, not cassette-intrinsic.

## Honest caveats
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- Organism = escherichia_coli_shigella (the determinant axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE of the target feature (a common feature present across many clades has cross-clade signal to learn; a clade-restricted accessory feature does not).