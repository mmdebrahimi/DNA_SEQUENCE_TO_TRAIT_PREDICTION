# Forward variant-effect method leaderboard (2026-07-14)

Per-protein Spearman(prediction, measured DMS) across the three forward-cell methods. AlphaMissense is human-only (bacterial cells show `—`); ESM2 runs on both but is only populated where a table was built. Higher = better; the learned methods (ESM2 / AlphaMissense) beat deterministic BLOSUM where run.

| protein (assay) | organism | n | BLOSUM62 | ESM2-650M | AlphaMissense | ESM-IF |
|---|---|---:|---:|---:|---:|---:|
| BLAT (BLAT_ECOLX_Stiffler_2015) | E. coli | 4996 | 0.346 | 0.732 | — | — |
| CCDB (CCDB_ECOLI_Tripathi_2016) | E. coli | 1663 | 0.248 | 0.511 | — | — |
| BLAT (BLAT_ECOLX_Firnberg_2014) | E. coli | 4783 | 0.339 | — | — | — |
| BLAT (BLAT_ECOLX_Deng_2012) | E. coli | 4996 | 0.305 | — | — | — |
| IF1 (IF1_ECOLI_Kelsic_2016) | E. coli | 1367 | 0.182 | — | — | — |
| MLAC (MLAC_ECOLI_MacRae_2023) | E. coli | 4007 | 0.182 | — | — | — |
| DYR (DYR_ECOLI_Nguyen_2023) | E. coli | 2916 | 0.152 | — | — | — |
| DYR (DYR_ECOLI_Thompson_2019) | E. coli | 2363 | 0.151 | — | — | — |
| ENVZ (ENVZ_ECOLI_Ghose_2023) | E. coli | 1121 | 0.101 | — | — | — |
| SR43C (SR43C_ARATH_Tsuboyama_2023_2N88) | Arabidopsis | 889 | 0.213 | 0.589 | — | — |
| MBD11 (MBD11_ARATH_Tsuboyama_2023_6ACV) | Arabidopsis | 1155 | 0.196 | — | — | — |
| CP2C9 (CP2C9_HUMAN_Amorosi_2021_abundance) | human | 6370 | 0.333 | — | 0.598 | — |
| TPMT (TPMT_HUMAN_Matreyek_2018) | human | 3648 | 0.240 | — | 0.558 | — |
| PTEN (PTEN_HUMAN_Mighell_2018) | human | 7260 | 0.182 | 0.518 | 0.539 | — |
| MSH2 (MSH2_HUMAN_Jia_2020) | human | 16749 | 0.164 | — | 0.416 | — |
| RL40A (RL40A_YEAST_Mavor_2016) | yeast | 1253 | 0.237 | 0.518 | — | — |
| BBC1 (BBC1_YEAST_Tsuboyama_2023_1TG0) | yeast | 1084 | 0.073 | — | — | — |
| THO1 (THO1_YEAST_Tsuboyama_2023_2WQG) | yeast | 656 | 0.057 | — | — | — |

Learned-vs-deterministic lift (where both present):
- BLAT_ECOLX_Stiffler_2015: ESM2 0.732 − BLOSUM 0.346 = **+0.385**
- CCDB_ECOLI_Tripathi_2016: ESM2 0.511 − BLOSUM 0.248 = **+0.263**
- SR43C_ARATH_Tsuboyama_2023_2N88: ESM2 0.589 − BLOSUM 0.213 = **+0.375**
- CP2C9_HUMAN_Amorosi_2021_abundance: AM 0.598 − BLOSUM 0.333 = **+0.265**
- TPMT_HUMAN_Matreyek_2018: AM 0.558 − BLOSUM 0.240 = **+0.318**
- PTEN_HUMAN_Mighell_2018: AM 0.539 − BLOSUM 0.182 = **+0.357**
- MSH2_HUMAN_Jia_2020: AM 0.416 − BLOSUM 0.164 = **+0.252**
- RL40A_YEAST_Mavor_2016: ESM2 0.518 − BLOSUM 0.237 = **+0.281**
