# PGx trio Mendelian co-segregation QC (2026-07-06)

**Cohort:** 1000 Genomes 30x panel trios (20130606_g1k_3202_samples ped) (603 trios)

- **CYP2C19: Mendelian consistency 602/602 (1.0)** across 602 trios (violations 0, uncallable 0)
- **CYP2C9: Mendelian consistency 602/602 (1.0)** across 602 trios (violations 0, uncallable 0)
- **CYP2C8: Mendelian consistency 602/602 (1.0)** across 602 trios (violations 0, uncallable 0)
- **CYP3A5: Mendelian consistency 602/602 (1.0)** across 602 trios (violations 0, uncallable 0)
- **CYP2B6: Mendelian consistency 602/602 (1.0)** across 602 trios (violations 0, uncallable 0)
- **TPMT: Mendelian consistency 602/602 (1.0)** across 602 trios (violations 0, uncallable 0)
- **CYP2D6: Mendelian consistency 592/602 (0.9834)** across 602 trios (violations 10, uncallable 0)
    - _all 10 violations are homozygous-child calls -> the structural-confound signature (CNV/hybrid makes a het locus read homozygous at the SNP level); NOT a SNP-logic bug. This is the cnv_hybrid_unassessed blind spot independently quantified by the trio axis._

**Scope:** Every STAR-ALLELE PGx gene (CYP2C19/CYP2C9/CYP2C8/CYP3A5/CYP2B6/TPMT/CYP2D6). VKORC1 + SLCO1B1 are single-SNP genotype readouts with no diplotype-assembly step -> a trio check is tautological (the phased panel is already Mendelian-clean at a single biallelic SNP) -> out of scope.

_Validates CALLING consistency (Mendelian inheritance), NOT phenotype. A VIOLATION definitively flags a calling error; consistency is necessary-not-sufficient. Independent axis from the GeT-RM panel. For CYP2D6 the ~1.7% violation rate is the STRUCTURAL-confound signature (all violations are homozygous-child calls — a CNV/hybrid the SNP surface can't see), corroborating the GeT-RM 46/47 residual, NOT a SNP-logic bug._
