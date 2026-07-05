# PGx decoder report card (2026-07-05)

_Standing PGx trust surface -- a roll-up, NOT a gate (exit 0 always). No aggregate headline; each cell's honest tier stands alone. CALLING is independently validatable vs GeT-RM (free consensus panel); PHENOTYPE is faithful-to-CPIC (assigned, not measured)._

| gene | trait | GeT-RM core | PharmCAT | func-evidence (A/D/F) | trio Mendelian | residual |
|---|---|---|---|---|---|---|
| CYP2C19 | metabolizer phenotype (PM/IM/NM/RM/UM) | 72/72 (1.0) | 6/6 | A2/D0/F1 | 602/602 | non-core *4/*35 withheld (sentinel v0.1) |
| CYP2C9 | metabolizer phenotype (activity-score) | 73/73 (1.0) | — | A1/D1/F0 | 602/602 | non-core *5/*8/*9/*11 withheld (sentinel v0.1); *6-indel/*61 residual |
| CYP2C8 | star-allele diplotype (*2/*3/*4) — CALLING only (no CPIC phenotype) | 82/82 (1.0) | — | — | — | rare non-core allele mis-called *1 (no sentinel layer v0); no phenotype layer by design |
| CYP3A5 | expressor/non-expressor phenotype (tacrolimus) | 8/8 (1.0) | — | — | — | UNDERPOWERED n=8 (only ~8 GeT-RM CYP3A5 samples overlap 1000G); rare non-core alleles mis-called *1 |
| VKORC1 | warfarin sensitivity (rs9923231) | — | — | A1/D0/F0 | — | — |
| SLCO1B1 | statin myopathy (rs4149056 / *5 521T>C) | — | — | — | — | single-SNP proxy for *5/*15/*17; full SLCO1B1 star typing needs more variants |

## Honest tier per cell

- **CYP2C19:** GeT-RM consensus (independent of consensus tools) + PharmCAT fixtures; phenotype faithful-to-CPIC
- **CYP2C9:** GeT-RM consensus; phenotype faithful-to-CPIC (activity-score)
- **CYP2C8:** GeT-RM consensus (independent of consensus tools); CALLING validated 82/82. NO CPIC metabolizer phenotype — CYP2C8 function is substrate-dependent, so this is a CALLING-only cell (never a PM/IM/NM). Region VCF fetched Docker-free (tabix-over-HTTP).
- **CYP3A5:** REAL GeT-RM CDC multi-lab consensus (independent of the labs); 8/8 core-diplotype incl. *1/*3/*6/*7 (the *7 insertion + *6/*7 non-expressor cases). UNDERPOWERED (n=8). Phenotype faithful-to-CPIC (expressor/non-expressor). First gene outside the CYP2C cluster.
- **VKORC1:** single-SNP genotype->sensitivity (minus-strand encoded); not a star/diplotype system
- **SLCO1B1:** single-SNP genotype->function readout (plus-strand); KNOWLEDGE_BASELINE like VKORC1. NOT an independent star number (rs4149056 IS the truth for a 521 call). CPIC-aligned (simvastatin function is assigned largely from 521T>C).

_Validation axes: GeT-RM = consensus concordance on real 1000G (independent of the consensus tools); PharmCAT = reference-tool fixtures; func-evidence = non-CPIC cross-check of the function assignment (AGREE/DISAGREE/FLAG); trio = Mendelian calling-consistency on 1000G trios. NOT a clinical tool._
