# PGx decoder report card (2026-07-06)

_Standing PGx trust surface -- a roll-up, NOT a gate (exit 0 always). No aggregate headline; each cell's honest tier stands alone. CALLING is independently validatable vs GeT-RM (free consensus panel); PHENOTYPE is faithful-to-CPIC (assigned, not measured)._

| gene | trait | GeT-RM core | PharmCAT | func-evidence (A/D/F) | trio Mendelian | residual |
|---|---|---|---|---|---|---|
| CYP2C19 | metabolizer phenotype (PM/IM/NM/RM/UM) | 72/72 (1.0) | 6/6 | A2/D0/F1 | 602/602 | non-core *4/*35 withheld (sentinel v0.1) |
| CYP2C9 | metabolizer phenotype (activity-score) | 73/73 (1.0) | — | A1/D1/F0 | 602/602 | non-core *5/*8/*9/*11 withheld (sentinel v0.1); *6-indel/*61 residual |
| CYP2C8 | star-allele diplotype (*2/*3/*4) — CALLING only (no CPIC phenotype) | 82/82 (1.0) | — | — | 602/602 | rare non-core allele mis-called *1 (no sentinel layer v0); no phenotype layer by design |
| CYP3A5 | expressor/non-expressor phenotype (tacrolimus) | 8/8 (1.0) | — | — | 602/602 | UNDERPOWERED n=8 (only ~8 GeT-RM CYP3A5 samples overlap 1000G); rare non-core alleles mis-called *1 |
| TPMT | thiopurine phenotype (COMPOUND *3A=*3B+*3C) | 85/85 (1.0) | — | — | 602/602 | rare non-core (*2/*8/*16...) mis-called *1 (no sentinel layer v0) |
| CYP2B6 | efavirenz phenotype (*6-proxy, 516G>T) | 62/62 (1.0) | — | — | 602/602 | single-SNP proxy; *9/*4/other non-core mis-called (785A>G absent from callset) |
| CYP2D6 | metabolizer phenotype (activity-score) — SNP surface only | 46/47 (0.9787) | — | — | 592/602 | hybrid IDENTITY (which of *13/*36/*68) unresolved (PSV analysis, Cyrius-class); hybrid detection sens partial (subtle *36 + opposite-signature *13 missed); non-core SNP alleles (*14/*15/*21/*40/*46) mis-called (no sentinel v0). Structural surface needs a BAM/CRAM |
| VKORC1 | warfarin sensitivity (rs9923231) | — | — | A1/D0/F0 | — | — |
| SLCO1B1 | statin myopathy (rs4149056 / *5 521T>C) | — | — | — | — | single-SNP proxy for *5/*15/*17; full SLCO1B1 star typing needs more variants |

## Honest tier per cell

- **CYP2C19:** GeT-RM consensus (independent of consensus tools) + PharmCAT fixtures; phenotype faithful-to-CPIC
- **CYP2C9:** GeT-RM consensus; phenotype faithful-to-CPIC (activity-score)
- **CYP2C8:** GeT-RM consensus (independent of consensus tools); CALLING validated 82/82. NO CPIC metabolizer phenotype — CYP2C8 function is substrate-dependent, so this is a CALLING-only cell (never a PM/IM/NM). Region VCF fetched Docker-free (tabix-over-HTTP).
- **CYP3A5:** REAL GeT-RM CDC multi-lab consensus (independent of the labs); 8/8 core-diplotype incl. *1/*3/*6/*7 (the *7 insertion + *6/*7 non-expressor cases). UNDERPOWERED (n=8). Phenotype faithful-to-CPIC (expressor/non-expressor). First gene outside the CYP2C cluster.
- **TPMT:** REAL GeT-RM CDC consolidated consensus; 85/85 core-comparable. FIRST compound-allele cell — *3A resolved from two SNPs in cis (*3B+*3C), each alone = *3B/*3C. Phenotype faithful-to-CPIC (thiopurine).
- **CYP2B6:** REAL GeT-RM CDC consolidated consensus; 62/62 on clean *1/*6. SINGLE-SNP *6-proxy (516G>T) — rs2279343 (785A>G) is absent from the 1000G 30x panel so *6 can't be split from *9. Phenotype faithful-to-CPIC.
- **CYP2D6:** The last major pharmacogene. GeT-RM consensus (independent of the consensus tools); 46/47 core-comparable on the SNP-DECODABLE subset (the single miss is a diagnosed structural confound). PRIORITY-ordered per-haplotype resolver (shared-background-aware). Phenotype faithful-to-CPIC (activity-score). Trio-Mendelian 592/602 — the ~2% residual is the structural-confound signature (all homozygous-child). STRUCTURAL SURFACE (read-depth off a real CRAM, dna_decode.pgx.cyp2d6_structural): *5 deletion + *xN duplication CN validated 26/26 on 1000G CRAMs (wiki/cyp2d6_structural_2026-07-06); HYBRID PRESENCE (*13/*36/*68) DETECTED via elevated CYP2D7 depth — sens 0.62 / spec 1.0 / AUROC 0.83 (wiki/cyp2d6_hybrid_2026-07-06; the *68 family detected cleanly). Hybrid IDENTITY (which of *13/*36/*68) still needs PSV analysis (Cyrius-class).
- **VKORC1:** single-SNP genotype->sensitivity (minus-strand encoded); not a star/diplotype system
- **SLCO1B1:** single-SNP genotype->function readout (plus-strand); KNOWLEDGE_BASELINE like VKORC1. NOT an independent star number (rs4149056 IS the truth for a 521 call). CPIC-aligned (simvastatin function is assigned largely from 521T>C).

_Validation axes: GeT-RM = consensus concordance on real 1000G (independent of the consensus tools); PharmCAT = reference-tool fixtures; func-evidence = non-CPIC cross-check of the function assignment (AGREE/DISAGREE/FLAG); trio = Mendelian calling-consistency on 1000G trios. NOT a clinical tool._
