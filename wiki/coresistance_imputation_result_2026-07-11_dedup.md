# Co-resistance imputation — can unobserved drug-class resistance be imputed? (2026-07-11)

**Verdict: PASS_CORESISTANCE_IMPUTABLE** — 18/21 (organism, drug-class) cells are IMPUTABLE from the other classes (fraction 0.857; PASS bar 0.5). (profile-deduped)

Can a genome's resistance to a drug CLASS be imputed from its OTHER-class determinants, within organism? (the 'predict the unobserved resistance profile' world-model payload)

IMPUTABLE = a genome's held-out drug-class presence is predicted from its OTHER-class determinants, within organism, OOF-AUC 95% CI lower > 0.5. This is the 'impute the unobserved resistance' payload.

| organism | genomes | classes | testable | **imputable** | frac | note |
|---|---|---|---|---|---|---|
| klebsiella | 59 | 17 | 8 | **5** | 0.625 | profile-deduped 307->59 genomes (clonality proxy) |
| escherichia_coli_shigella | 106 | 18 | 12 | **12** | 1.0 | profile-deduped 240->106 genomes (clonality proxy) |
| campylobacter | 20 | 7 | 0 | **0** | 0.0 | profile-deduped 100->20 genomes (clonality proxy) |
| acinetobacter | 16 | 12 | 0 | **0** | 0.0 | profile-deduped 60->16 genomes (clonality proxy) |
| salmonella | 22 | 12 | 1 | **1** | 1.0 | profile-deduped 60->22 genomes (clonality proxy) |

## Per-class imputation AUC (held-out; how well each drug class is imputed from the others)
- **klebsiella**: TETRACYCLINE(0.9151*), SULFONAMIDE(0.8625*), AMINOGLYCOSIDE(0.8428*), NITROFURAN(0.8014*), MACROLIDE(0.7907*), TRIMETHOPRIM(0.6585), COLISTIN(0.5729), RIFAMYCIN(0.5661)
- **escherichia_coli_shigella**: TRIMETHOPRIM(0.9147*), SULFONAMIDE(0.903*), AMINOGLYCOSIDE(0.8947*), QUINOLONE(0.7931*), BETA-LACTAM(0.768*), MACROLIDE(0.7454*), MULTIDRUG(0.7421*), PHENICOL(0.742*)
- **campylobacter**: 
- **acinetobacter**: 
- **salmonella**: SULFONAMIDE(0.7708*)

## Co-resistance NETWORK (drug-class A → classes it predicts, by lift)
### klebsiella
- **AMINOGLYCOSIDE** → BLEOMYCIN(lift 1.28), COLISTIN(lift 1.28), MACROLIDE(lift 1.28), SULFONAMIDE(lift 1.22)
- **COLISTIN** → AMINOGLYCOSIDE(lift 1.28), SULFONAMIDE(lift 1.07), FOSFOMYCIN(lift 1.04), BETA-LACTAM(lift 1.02)
- **MACROLIDE** → SULFONAMIDE(lift 1.4), AMINOGLYCOSIDE(lift 1.28), TRIMETHOPRIM(lift 1.26), NITROFURAN(lift 1.18)
- **NITROFURAN** → TETRACYCLINE(lift 1.28), RIFAMYCIN(lift 1.19), MACROLIDE(lift 1.18), PHENICOL(lift 1.18)
- **RIFAMYCIN** → SULFONAMIDE(lift 1.36), TETRACYCLINE(lift 1.21), NITROFURAN(lift 1.19), TRIMETHOPRIM(lift 1.19)
- **SULFONAMIDE** → BLEOMYCIN(lift 1.48), MACROLIDE(lift 1.4), RIFAMYCIN(lift 1.36), AMINOGLYCOSIDE(lift 1.22)
- **TETRACYCLINE** → BLEOMYCIN(lift 1.31), NITROFURAN(lift 1.28), RIFAMYCIN(lift 1.21), PHENICOL(lift 1.15)
- **TRIMETHOPRIM** → MACROLIDE(lift 1.26), SULFONAMIDE(lift 1.2), RIFAMYCIN(lift 1.19), AMINOGLYCOSIDE(lift 1.11)
### escherichia_coli_shigella
- **AMINOGLYCOSIDE** → BLEOMYCIN(lift 1.33), LINCOSAMIDE(lift 1.33), STREPTOTHRICIN(lift 1.33), SULFONAMIDE(lift 1.28)
- **BETA-LACTAM** → TRIMETHOPRIM(lift 1.1), SULFONAMIDE(lift 1.09), MACROLIDE(lift 1.08), MULTIDRUG(lift 1.08)
- **FOSFOMYCIN** → MULTIDRUG(lift 1.66), QUINOLONE(lift 1.16), MACROLIDE(lift 1.14), PHENICOL(lift 1.1)
- **FOSMIDOMYCIN** → TETRACYCLINE(lift 1.24), MACROLIDE(lift 1.21), PHENICOL(lift 1.18), QUINOLONE(lift 1.08)
- **MACROLIDE** → MULTIDRUG(lift 1.55), TRIMETHOPRIM(lift 1.49), SULFONAMIDE(lift 1.46), PHENICOL(lift 1.31)
- **MULTIDRUG** → FOSFOMYCIN(lift 1.66), MACROLIDE(lift 1.55), PHENICOL(lift 1.5), QUINOLONE(lift 1.38)
- **NITROFURAN** → TETRACYCLINE(lift 1.54), TRIMETHOPRIM(lift 1.47), AMINOGLYCOSIDE(lift 1.21), SULFONAMIDE(lift 1.09)
- **PHENICOL** → BLEOMYCIN(lift 2.26), MULTIDRUG(lift 1.5), TRIMETHOPRIM(lift 1.43), MACROLIDE(lift 1.31)
- **QUINOLONE** → MULTIDRUG(lift 1.38), MACROLIDE(lift 1.23), PHENICOL(lift 1.23), BLEOMYCIN(lift 1.22)
- **SULFONAMIDE** → BLEOMYCIN(lift 1.63), MACROLIDE(lift 1.46), TRIMETHOPRIM(lift 1.44), AMINOGLYCOSIDE(lift 1.28)
- **TETRACYCLINE** → STREPTOTHRICIN(lift 1.68), NITROFURAN(lift 1.54), FOSMIDOMYCIN(lift 1.24), MACROLIDE(lift 1.24)
- **TRIMETHOPRIM** → BLEOMYCIN(lift 1.77), MACROLIDE(lift 1.49), NITROFURAN(lift 1.47), SULFONAMIDE(lift 1.44)
### campylobacter
### acinetobacter
### salmonella
- **SULFONAMIDE** → AMINOGLYCOSIDE(lift 1.29), QUINOLONE(lift 1.29), BETA-LACTAM(lift 1.27), TETRACYCLINE(lift 1.26)

## Honest caveats
- Cross-axis (virulence/stress) NOT available from cached AMR-only runs; AMR drug-class only.
- Cohorts are drug-R/S-SELECTED -> the co-resistance structure reflects the curated cohorts.
- Within-organism de-confounds species; dedup is a crude clonality proxy (Mash-collapse is the follow-on).
- Class presence is DETERMINANT-implied (AMRFinder Class), not measured MIC per drug — a genotype axis.
- Self-distillation from our own caller: associational, NOT causal.