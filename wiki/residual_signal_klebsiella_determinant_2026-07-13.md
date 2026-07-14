# Residual-signal report — klebsiella / AMR determinant

generated 2026-07-13 · dna_decode.eval.residual_detector (g8) · signal-provenance, NOT a phenotype prediction

**What this is.** For each feature, whether its genotype→axis association SURVIVES leave-one-clade-out (Mash-clade GroupKFold) de-confounding. `generalizes` = real signal beyond lineage; `lineage_mediated` = clade-concentrated (may be clonal structure); `untested` = no cross-axis entry.

- source verdict: `CROSS_AXIS_GENERALIZES_BEYOND_LINEAGE`
- median AUC naive 0.973 → clade-grouped 0.913
- Mash partition: threshold 0.005, 118 clades, largest-clade fraction 0.254
- **tiers:** generalizes=46 · lineage_mediated=14 · untested=0 (60 features)

## Per-feature (residual signal — GENERALIZES first, strongest de-confounded AUC on top)

| feature | tier | family | n | AUC naive | AUC clade-grouped | drop |
|---|---|---|---:|---:|---:|---:|
| `estX-3` | generalizes | other | 26 | 0.998 | 1.0 | -0.001 |
| `blaOXA-1` | generalizes | beta-lactam | 48 | 1.0 | 0.998 | 0.001 |
| `cmlA1` | generalizes | other | 24 | 0.996 | 0.998 | -0.002 |
| `aph(3'')-Ib` | generalizes | aminoglycoside | 76 | 0.996 | 0.994 | 0.002 |
| `sul3` | generalizes | sulfa/trimethoprim | 28 | 1.0 | 0.994 | 0.006 |
| `gyrA_S83T` | generalizes | QRDR | 38 | 0.997 | 0.993 | 0.004 |
| `mrx(A)` | generalizes | other | 76 | 0.997 | 0.993 | 0.004 |
| `aac(6')-Ib-cr5` | generalizes | aminoglycoside | 48 | 0.989 | 0.992 | -0.003 |
| `mph(A)` | generalizes | tet/macrolide | 77 | 0.999 | 0.991 | 0.008 |
| `catB3` | generalizes | phenicol | 53 | 0.989 | 0.978 | 0.011 |
| `sul2` | generalizes | sulfa/trimethoprim | 82 | 0.987 | 0.972 | 0.015 |
| `aph(6)-Id` | generalizes | aminoglycoside | 78 | 0.999 | 0.97 | 0.029 |
| `dfrA12` | generalizes | sulfa/trimethoprim | 61 | 0.996 | 0.967 | 0.029 |
| `ampC-Kaer` | generalizes | other | 30 | 0.991 | 0.962 | 0.03 |
| `dfrA50` | generalizes | sulfa/trimethoprim | 21 | 0.97 | 0.961 | 0.009 |
| `dfrA14` | generalizes | sulfa/trimethoprim | 64 | 0.977 | 0.959 | 0.018 |
| `sul1` | generalizes | sulfa/trimethoprim | 108 | 0.976 | 0.952 | 0.023 |
| `parC_S80I` | generalizes | QRDR | 113 | 0.997 | 0.949 | 0.048 |
| `fosA5` | generalizes | fosfomycin | 12 | 0.999 | 0.947 | 0.052 |
| `aadA2` | generalizes | aminoglycoside | 75 | 0.969 | 0.944 | 0.025 |
| `oqxB5` | generalizes | quinolone-plasmid | 14 | 0.976 | 0.944 | 0.031 |
| `tet(A)` | generalizes | tet/macrolide | 52 | 0.957 | 0.944 | 0.013 |
| `oqxB` | generalizes | quinolone-plasmid | 184 | 0.977 | 0.943 | 0.033 |
| `qnrB1` | generalizes | quinolone-plasmid | 34 | 0.934 | 0.937 | -0.003 |
| `gyrA_S83I` | generalizes | QRDR | 97 | 0.993 | 0.933 | 0.059 |
| `blaCTX-M-15` | generalizes | beta-lactam | 87 | 0.977 | 0.919 | 0.058 |
| `oqxA11` | generalizes | quinolone-plasmid | 14 | 0.993 | 0.918 | 0.076 |
| `oqxB25` | generalizes | quinolone-plasmid | 23 | 0.967 | 0.915 | 0.053 |
| `aac(6')-Ib` | generalizes | aminoglycoside | 58 | 0.97 | 0.913 | 0.056 |
| `aadA1` | generalizes | aminoglycoside | 44 | 0.938 | 0.904 | 0.034 |
| `oqxB19` | generalizes | quinolone-plasmid | 55 | 0.98 | 0.902 | 0.078 |
| `arr-2` | generalizes | other | 10 | 0.958 | 0.897 | 0.061 |
| `blaKPC-3` | generalizes | beta-lactam | 29 | 0.923 | 0.864 | 0.06 |
| `oqxA` | generalizes | quinolone-plasmid | 259 | 0.958 | 0.862 | 0.097 |
| `blaKPC-2` | generalizes | beta-lactam | 46 | 0.916 | 0.861 | 0.055 |
| `aac(3)-IIe` | generalizes | aminoglycoside | 31 | 0.955 | 0.835 | 0.121 |
| `qnrS1` | generalizes | quinolone-plasmid | 32 | 0.884 | 0.835 | 0.049 |
| `blaTEM-1` | generalizes | beta-lactam | 94 | 0.917 | 0.815 | 0.102 |
| `blaTEM` | generalizes | beta-lactam | 26 | 0.79 | 0.806 | -0.015 |
| `blaOXA` | generalizes | beta-lactam | 46 | 0.891 | 0.79 | 0.102 |
| `aph(3')-Ia` | generalizes | aminoglycoside | 32 | 0.85 | 0.747 | 0.103 |
| `blaSHV-1` | generalizes | beta-lactam | 46 | 0.934 | 0.743 | 0.191 |
| `gyrA_S83F` | generalizes | QRDR | 20 | 0.996 | 0.719 | 0.277 |
| `aac(3)-IId` | generalizes | aminoglycoside | 18 | 0.876 | 0.714 | 0.162 |
| `arr-3` | generalizes | other | 10 | 0.939 | 0.714 | 0.226 |
| `dfrA1` | generalizes | sulfa/trimethoprim | 12 | 0.852 | 0.702 | 0.15 |
| `fosA` | lineage_mediated | fosfomycin | 287 | 0.91 | 0.674 | 0.235 |
| `ompK36_D135DGD` | lineage_mediated | other | 32 | 0.962 | 0.672 | 0.289 |
| `gyrA_D87N` | lineage_mediated | QRDR | 32 | 0.983 | 0.58 | 0.403 |
| `blaSHV-11` | lineage_mediated | beta-lactam | 91 | 0.842 | 0.579 | 0.264 |
| `tet(D)` | lineage_mediated | tet/macrolide | 24 | 0.853 | 0.558 | 0.296 |
| `blaSHV` | lineage_mediated | beta-lactam | 26 | 0.664 | 0.542 | 0.122 |
| `gyrA_S83Y` | lineage_mediated | QRDR | 10 | 0.842 | 0.529 | 0.312 |
| `blaOXA-9` | lineage_mediated | beta-lactam | 28 | 0.967 | 0.523 | 0.444 |
| `catA1` | lineage_mediated | phenicol | 31 | 0.879 | 0.458 | 0.421 |
| `blaSHV-28` | lineage_mediated | beta-lactam | 16 | 0.978 | 0.438 | 0.54 |
| `blaSHV-12` | lineage_mediated | beta-lactam | 71 | 0.961 | 0.296 | 0.666 |
| `aac(3)-IVa` | lineage_mediated | aminoglycoside | 24 | 0.999 | None | None |
| `aph(4)-Ia` | lineage_mediated | aminoglycoside | 24 | 0.999 | None | None |
| `ompK35_E42RfsTer47` | lineage_mediated | other | 72 | 0.999 | None | None |

## Gene-family rollup (residual signal by family)

| family | generalizes | lineage-mediated | untested | total |
|---|---:|---:|---:|---:|
| aminoglycoside | 9 | 2 | 0 | 11 |
| beta-lactam | 8 | 5 | 0 | 13 |
| quinolone-plasmid | 8 | 0 | 0 | 8 |
| sulfa/trimethoprim | 7 | 0 | 0 | 7 |
| other | 6 | 2 | 0 | 8 |
| QRDR | 4 | 2 | 0 | 6 |
| tet/macrolide | 2 | 1 | 0 | 3 |
| phenicol | 1 | 1 | 0 | 2 |
| fosfomycin | 1 | 1 | 0 | 2 |

## Honest caveats

- A residual-signal tier is a SIGNAL-PROVENANCE statement, NOT a phenotype prediction — the report says WHERE de-confounded mechanism signal lives, never R/S.
- GENERALIZES = association survives leave-one-clade-out CV (real signal beyond lineage). LINEAGE_MEDIATED = clade-concentrated: MAY be clonal population structure, not a real linkage.
- De-confounding quality is bounded by the Mash-clade partition (threshold + largest-clade fraction in meta); a coarse partition under-removes lineage.
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- Organism = klebsiella (the determinant axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE of the target feature (a common feature present across many clades has cross-clade signal to learn; a clade-restricted accessory feature does not).
