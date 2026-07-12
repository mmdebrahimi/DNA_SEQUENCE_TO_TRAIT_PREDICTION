# non-determinant → determinant cross-axis — lineage de-confound (leave-one-clade-out) (2026-07-12)

**Verdict: CROSS_AXIS_GENERALIZES_BEYOND_LINEAGE** — median clade-grouped AUC = 0.913 (vs naive 0.973); 46/60 AMR determinants still predicted at AUC >= 0.7 when their clade is held out.

Mash: 307 E. coli genomes → 118 clades at threshold 0.005 (largest clade 0.254).

Does the E. coli non-determinant -> determinant cross-axis signal survive leave-one-clade-out CV (real beyond lineage), or collapse (lineage-mediated)?

| AMR determinant | n | AUC naive | **AUC clade-grouped** | drop | generalizes |
|---|---|---|---|---|---|
| estX-3 | 26 | 0.998 | **1.0** | -0.001 | YES |
| blaOXA-1 | 48 | 1.0 | **0.998** | 0.001 | YES |
| cmlA1 | 24 | 0.996 | **0.998** | -0.002 | YES |
| aph(3'')-Ib | 76 | 0.996 | **0.994** | 0.002 | YES |
| sul3 | 28 | 1.0 | **0.994** | 0.006 | YES |
| gyrA_S83T | 38 | 0.997 | **0.993** | 0.004 | YES |
| mrx(A) | 76 | 0.997 | **0.993** | 0.004 | YES |
| aac(6')-Ib-cr5 | 48 | 0.989 | **0.992** | -0.003 | YES |
| mph(A) | 77 | 0.999 | **0.991** | 0.008 | YES |
| catB3 | 53 | 0.989 | **0.978** | 0.011 | YES |
| sul2 | 82 | 0.987 | **0.972** | 0.015 | YES |
| aph(6)-Id | 78 | 0.999 | **0.97** | 0.029 | YES |
| dfrA12 | 61 | 0.996 | **0.967** | 0.029 | YES |
| ampC-Kaer | 30 | 0.991 | **0.962** | 0.03 | YES |
| dfrA50 | 21 | 0.97 | **0.961** | 0.009 | YES |
| dfrA14 | 64 | 0.977 | **0.959** | 0.018 | YES |
| sul1 | 108 | 0.976 | **0.952** | 0.023 | YES |
| parC_S80I | 113 | 0.997 | **0.949** | 0.048 | YES |
| fosA5 | 12 | 0.999 | **0.947** | 0.052 | YES |
| aadA2 | 75 | 0.969 | **0.944** | 0.025 | YES |
| oqxB5 | 14 | 0.976 | **0.944** | 0.031 | YES |
| tet(A) | 52 | 0.957 | **0.944** | 0.013 | YES |
| oqxB | 184 | 0.977 | **0.943** | 0.033 | YES |
| qnrB1 | 34 | 0.934 | **0.937** | -0.003 | YES |
| gyrA_S83I | 97 | 0.993 | **0.933** | 0.059 | YES |
| blaCTX-M-15 | 87 | 0.977 | **0.919** | 0.058 | YES |
| oqxA11 | 14 | 0.993 | **0.918** | 0.076 | YES |
| oqxB25 | 23 | 0.967 | **0.915** | 0.053 | YES |
| aac(6')-Ib | 58 | 0.97 | **0.913** | 0.056 | YES |
| aadA1 | 44 | 0.938 | **0.904** | 0.034 | YES |
| oqxB19 | 55 | 0.98 | **0.902** | 0.078 | YES |
| arr-2 | 10 | 0.958 | **0.897** | 0.061 | YES |
| blaKPC-3 | 29 | 0.923 | **0.864** | 0.06 | YES |
| oqxA | 259 | 0.958 | **0.862** | 0.097 | YES |
| blaKPC-2 | 46 | 0.916 | **0.861** | 0.055 | YES |
| aac(3)-IIe | 31 | 0.955 | **0.835** | 0.121 | YES |
| qnrS1 | 32 | 0.884 | **0.835** | 0.049 | YES |
| blaTEM-1 | 94 | 0.917 | **0.815** | 0.102 | YES |
| blaTEM | 26 | 0.79 | **0.806** | -0.015 | YES |
| blaOXA | 46 | 0.891 | **0.79** | 0.102 | YES |
| aph(3')-Ia | 32 | 0.85 | **0.747** | 0.103 | YES |
| blaSHV-1 | 46 | 0.934 | **0.743** | 0.191 | YES |
| gyrA_S83F | 20 | 0.996 | **0.719** | 0.277 | YES |
| aac(3)-IId | 18 | 0.876 | **0.714** | 0.162 | YES |
| arr-3 | 10 | 0.939 | **0.714** | 0.226 | YES |
| dfrA1 | 12 | 0.852 | **0.702** | 0.15 | YES |
| fosA | 287 | 0.91 | **0.674** | 0.235 | no (lineage) |
| ompK36_D135DGD | 32 | 0.962 | **0.672** | 0.289 | no (lineage) |
| gyrA_D87N | 32 | 0.983 | **0.58** | 0.403 | no (lineage) |
| blaSHV-11 | 91 | 0.842 | **0.579** | 0.264 | no (lineage) |
| tet(D) | 24 | 0.853 | **0.558** | 0.296 | no (lineage) |
| blaSHV | 26 | 0.664 | **0.542** | 0.122 | no (lineage) |
| gyrA_S83Y | 10 | 0.842 | **0.529** | 0.312 | no (lineage) |
| blaOXA-9 | 28 | 0.967 | **0.523** | 0.444 | no (lineage) |
| catA1 | 31 | 0.879 | **0.458** | 0.421 | no (lineage) |
| blaSHV-28 | 16 | 0.978 | **0.438** | 0.54 | no (lineage) |
| blaSHV-12 | 71 | 0.961 | **0.296** | 0.666 | no (lineage) |
| aac(3)-IVa | 24 | 0.999 | **clade-concentrated** | None | no (clade-concentrated) |
| aph(4)-Ia | 24 | 0.999 | **clade-concentrated** | None | no (clade-concentrated) |
| ompK35_E42RfsTer47 | 72 | 0.999 | **clade-concentrated** | None | no (clade-concentrated) |

## Literature / mechanism
- The C-core finding (PASS_LINKAGE_STRUCTURE) showed AMR determinants co-occur in cassette blocks (integron sul1->aadA/dfrA, QRDR clusters) under a determinant-profile-DEDUP clonality proxy. This applies the STRONGER phylogenetic (Mash-clade leave-one-out) control: an integron cassette is itself a MOBILE unit, so co-cassette determinants should predict each other even across held-out clades (real cassette linkage); a collapse would mean the co-occurrence is clade-signature, not cassette-intrinsic.

## Honest caveats
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- Organism = klebsiella (the determinant axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE of the target feature (a common feature present across many clades has cross-clade signal to learn; a clade-restricted accessory feature does not).