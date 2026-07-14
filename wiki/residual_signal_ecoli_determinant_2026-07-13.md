# Residual-signal report — escherichia_coli_shigella / AMR determinant

generated 2026-07-13 · dna_decode.eval.residual_detector (g8) · signal-provenance, NOT a phenotype prediction

**What this is.** For each feature, whether its genotype→axis association SURVIVES leave-one-clade-out (Mash-clade GroupKFold) de-confounding. `generalizes` = real signal beyond lineage; `lineage_mediated` = clade-concentrated (may be clonal structure); `untested` = no cross-axis entry.

- source verdict: `CROSS_AXIS_GENERALIZES_BEYOND_LINEAGE`
- median AUC naive 0.975 → clade-grouped 0.908
- Mash partition: threshold 0.005, 80 clades, largest-clade fraction 0.353
- **tiers:** generalizes=39 · lineage_mediated=6 · untested=0 (45 features)

## Per-feature (residual signal — GENERALIZES first, strongest de-confounded AUC on top)

| feature | tier | family | n | AUC naive | AUC clade-grouped | drop |
|---|---|---|---:|---:|---:|---:|
| `mrx(A)` | generalizes | other | 98 | 0.997 | 0.997 | 0.0 |
| `cmlA1` | generalizes | other | 13 | 0.996 | 0.993 | 0.003 |
| `mph(A)` | generalizes | tet/macrolide | 98 | 0.993 | 0.992 | 0.001 |
| `aac(6')-Ib-cr5` | generalizes | aminoglycoside | 69 | 0.994 | 0.99 | 0.005 |
| `aadA5` | generalizes | aminoglycoside | 90 | 0.996 | 0.988 | 0.008 |
| `blaOXA-1` | generalizes | beta-lactam | 69 | 0.987 | 0.98 | 0.007 |
| `aph(3'')-Ib` | generalizes | aminoglycoside | 100 | 0.999 | 0.978 | 0.02 |
| `dfrA17` | generalizes | sulfa/trimethoprim | 94 | 0.989 | 0.977 | 0.012 |
| `sat2` | generalizes | other | 10 | 0.981 | 0.977 | 0.004 |
| `gyrA_D87N` | generalizes | QRDR | 140 | 0.997 | 0.976 | 0.02 |
| `parC_S80I` | generalizes | QRDR | 141 | 0.997 | 0.976 | 0.021 |
| `tet(M)` | generalizes | tet/macrolide | 11 | 0.984 | 0.972 | 0.012 |
| `dfrA12` | generalizes | sulfa/trimethoprim | 27 | 0.993 | 0.97 | 0.022 |
| `gyrA_S83L` | generalizes | QRDR | 159 | 0.968 | 0.962 | 0.006 |
| `sul2` | generalizes | sulfa/trimethoprim | 102 | 0.978 | 0.957 | 0.021 |
| `aadA2` | generalizes | aminoglycoside | 31 | 0.964 | 0.955 | 0.009 |
| `catB3` | generalizes | phenicol | 69 | 0.984 | 0.953 | 0.031 |
| `ptsI_V25I` | generalizes | other | 94 | 0.997 | 0.94 | 0.058 |
| `sul1` | generalizes | sulfa/trimethoprim | 107 | 0.957 | 0.931 | 0.025 |
| `uhpT_E350Q` | generalizes | other | 111 | 0.968 | 0.915 | 0.053 |
| `estX-3` | generalizes | other | 12 | 0.982 | 0.91 | 0.072 |
| `aadA1` | generalizes | aminoglycoside | 40 | 0.952 | 0.909 | 0.043 |
| `aph(6)-Id` | generalizes | aminoglycoside | 86 | 0.981 | 0.907 | 0.074 |
| `sul3` | generalizes | sulfa/trimethoprim | 25 | 0.997 | 0.907 | 0.09 |
| `blaCTX-M-15` | generalizes | beta-lactam | 98 | 0.963 | 0.899 | 0.064 |
| `blaCTX-M-55` | generalizes | beta-lactam | 15 | 0.854 | 0.898 | -0.044 |
| `blaCTX-M-27` | generalizes | beta-lactam | 22 | 0.947 | 0.896 | 0.051 |
| `aac(3)-IIe` | generalizes | aminoglycoside | 54 | 0.975 | 0.868 | 0.107 |
| `parE_I529L` | generalizes | QRDR | 88 | 0.998 | 0.852 | 0.146 |
| `qnrS1` | generalizes | quinolone-plasmid | 15 | 0.851 | 0.824 | 0.027 |
| `blaCMY-2` | generalizes | beta-lactam | 15 | 0.81 | 0.819 | -0.009 |
| `catA1` | generalizes | phenicol | 23 | 0.892 | 0.817 | 0.076 |
| `dfrA14` | generalizes | sulfa/trimethoprim | 27 | 0.933 | 0.811 | 0.121 |
| `aph(3')-Ia` | generalizes | aminoglycoside | 23 | 0.922 | 0.809 | 0.113 |
| `blaTEM-1` | generalizes | beta-lactam | 80 | 0.855 | 0.744 | 0.111 |
| `floR` | generalizes | other | 19 | 0.84 | 0.736 | 0.105 |
| `aac(3)-IId` | generalizes | aminoglycoside | 17 | 0.896 | 0.72 | 0.176 |
| `tet(B)` | generalizes | tet/macrolide | 57 | 0.883 | 0.715 | 0.168 |
| `dfrA1` | generalizes | sulfa/trimethoprim | 15 | 0.901 | 0.708 | 0.193 |
| `parE_S458A` | lineage_mediated | QRDR | 33 | 0.968 | 0.682 | 0.285 |
| `blaTEM` | lineage_mediated | beta-lactam | 10 | 0.745 | 0.605 | 0.14 |
| `tet(A)` | lineage_mediated | tet/macrolide | 109 | 0.844 | 0.579 | 0.265 |
| `cyaA_S352T` | lineage_mediated | other | 31 | 0.877 | 0.56 | 0.317 |
| `parE_L416F` | lineage_mediated | QRDR | 11 | 0.987 | 0.245 | 0.742 |
| `parC_E84V` | lineage_mediated | QRDR | 77 | 1.0 | None | None |

## Gene-family rollup (residual signal by family)

| family | generalizes | lineage-mediated | untested | total |
|---|---:|---:|---:|---:|
| aminoglycoside | 9 | 0 | 0 | 9 |
| other | 7 | 1 | 0 | 8 |
| sulfa/trimethoprim | 7 | 0 | 0 | 7 |
| beta-lactam | 6 | 1 | 0 | 7 |
| QRDR | 4 | 3 | 0 | 7 |
| tet/macrolide | 3 | 1 | 0 | 4 |
| phenicol | 2 | 0 | 0 | 2 |
| quinolone-plasmid | 1 | 0 | 0 | 1 |

## Honest caveats

- A residual-signal tier is a SIGNAL-PROVENANCE statement, NOT a phenotype prediction — the report says WHERE de-confounded mechanism signal lives, never R/S.
- GENERALIZES = association survives leave-one-clade-out CV (real signal beyond lineage). LINEAGE_MEDIATED = clade-concentrated: MAY be clonal population structure, not a real linkage.
- De-confounding quality is bounded by the Mash-clade partition (threshold + largest-clade fraction in meta); a coarse partition under-removes lineage.
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- Organism = escherichia_coli_shigella (the determinant axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE of the target feature (a common feature present across many clades has cross-clade signal to learn; a clade-restricted accessory feature does not).
