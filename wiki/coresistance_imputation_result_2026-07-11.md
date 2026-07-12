# Co-resistance imputation — can unobserved drug-class resistance be imputed? (2026-07-11)

**Verdict: PASS_CORESISTANCE_IMPUTABLE** — 35/41 (organism, drug-class) cells are IMPUTABLE from the other classes (fraction 0.854; PASS bar 0.5). (raw)

Can a genome's resistance to a drug CLASS be imputed from its OTHER-class determinants, within organism? (the 'predict the unobserved resistance profile' world-model payload)

IMPUTABLE = a genome's held-out drug-class presence is predicted from its OTHER-class determinants, within organism, OOF-AUC 95% CI lower > 0.5. This is the 'impute the unobserved resistance' payload.

| organism | genomes | classes | testable | **imputable** | frac | note |
|---|---|---|---|---|---|---|
| klebsiella | 307 | 17 | 11 | **10** | 0.909 |  |
| escherichia_coli_shigella | 240 | 18 | 14 | **14** | 1.0 |  |
| campylobacter | 100 | 7 | 4 | **1** | 0.25 |  |
| acinetobacter | 60 | 12 | 4 | **3** | 0.75 |  |
| salmonella | 60 | 12 | 8 | **7** | 0.875 |  |

## Per-class imputation AUC (held-out; how well each drug class is imputed from the others)
- **klebsiella**: PHENICOL(0.9782*), SULFONAMIDE(0.972*), TETRACYCLINE(0.9676*), AMINOGLYCOSIDE(0.927*), NITROFURAN(0.8847*), BLEOMYCIN(0.8756*), TRIMETHOPRIM(0.8561*), MACROLIDE(0.8462*)
- **escherichia_coli_shigella**: AMINOGLYCOSIDE(0.9695*), SULFONAMIDE(0.9665*), TRIMETHOPRIM(0.9617*), BETA-LACTAM(0.9271*), QUINOLONE(0.9191*), BLEOMYCIN(0.9138*), MACROLIDE(0.8973*), STREPTOTHRICIN(0.8667*)
- **campylobacter**: QUINOLONE(0.6204*), TETRACYCLINE(0.6038), AMINOGLYCOSIDE(0.5566), BETA-LACTAM(0.2502)
- **acinetobacter**: MACROLIDE(1.0*), STREPTOGRAMIN(1.0*), TRIMETHOPRIM(0.8988*), TETRACYCLINE(0.5404)
- **salmonella**: AMINOGLYCOSIDE(0.9764*), QUINOLONE(0.9667*), PHENICOL(0.9656*), SULFONAMIDE(0.9589*), TRIMETHOPRIM(0.9264*), BETA-LACTAM(0.9207*), TETRACYCLINE(0.8566*), FOSFOMYCIN(0.467)

## Co-resistance NETWORK (drug-class A → classes it predicts, by lift)
### klebsiella
- **AMINOGLYCOSIDE** → BLEOMYCIN(lift 1.57), COLISTIN(lift 1.57), MACROLIDE(lift 1.57), SULFONAMIDE(lift 1.55)
- **BLEOMYCIN** → MACROLIDE(lift 3.04), SULFONAMIDE(lift 1.78), TRIMETHOPRIM(lift 1.76), AMINOGLYCOSIDE(lift 1.57)
- **COLISTIN** → MACROLIDE(lift 2.25), AMINOGLYCOSIDE(lift 1.57), TRIMETHOPRIM(lift 1.52), SULFONAMIDE(lift 1.49)
- **MACROLIDE** → BLEOMYCIN(lift 3.04), COLISTIN(lift 2.25), RIFAMYCIN(lift 1.78), SULFONAMIDE(lift 1.77)
- **MULTIDRUG** → MACROLIDE(lift 1.38), TRIMETHOPRIM(lift 1.16), NITROFURAN(lift 1.11), TETRACYCLINE(lift 1.06)
- **NITROFURAN** → MULTIDRUG(lift 1.11), TETRACYCLINE(lift 1.1), MACROLIDE(lift 1.08), BLEOMYCIN(lift 1.05)
- **PHENICOL** → BLEOMYCIN(lift 1.04), NITROFURAN(lift 1.04), RIFAMYCIN(lift 1.04), TETRACYCLINE(lift 1.04)
- **RIFAMYCIN** → MACROLIDE(lift 1.78), SULFONAMIDE(lift 1.69), TRIMETHOPRIM(lift 1.65), AMINOGLYCOSIDE(lift 1.49)
- **SULFONAMIDE** → BLEOMYCIN(lift 1.78), MACROLIDE(lift 1.77), RIFAMYCIN(lift 1.69), TRIMETHOPRIM(lift 1.63)
- **TETRACYCLINE** → BLEOMYCIN(lift 1.11), NITROFURAN(lift 1.1), MULTIDRUG(lift 1.06), MACROLIDE(lift 1.05)
- **TRIMETHOPRIM** → BLEOMYCIN(lift 1.76), MACROLIDE(lift 1.7), RIFAMYCIN(lift 1.65), SULFONAMIDE(lift 1.63)
### escherichia_coli_shigella
- **AMINOGLYCOSIDE** → BLEOMYCIN(lift 1.34), LINCOSAMIDE(lift 1.34), STREPTOGRAMIN(lift 1.34), STREPTOTHRICIN(lift 1.34)
- **BETA-LACTAM** → QUINOLONE(lift 1.14), TRIMETHOPRIM(lift 1.13), AMINOGLYCOSIDE(lift 1.12), MULTIDRUG(lift 1.12)
- **BLEOMYCIN** → PHENICOL(lift 2.22), MACROLIDE(lift 1.71), TRIMETHOPRIM(lift 1.64), SULFONAMIDE(lift 1.59)
- **FOSFOMYCIN** → MACROLIDE(lift 1.26), QUINOLONE(lift 1.25), MULTIDRUG(lift 1.16), BETA-LACTAM(lift 1.1)
- **FOSMIDOMYCIN** → TETRACYCLINE(lift 1.09), PHENICOL(lift 1.08), BETA-LACTAM(lift 1.06), QUINOLONE(lift 1.06)
- **MACROLIDE** → STREPTOGRAMIN(lift 2.22), BLEOMYCIN(lift 1.71), TRIMETHOPRIM(lift 1.55), SULFONAMIDE(lift 1.53)
- **MULTIDRUG** → PHENICOL(lift 1.52), QUINOLONE(lift 1.37), FOSFOMYCIN(lift 1.16), BETA-LACTAM(lift 1.12)
- **NITROFURAN** → TETRACYCLINE(lift 1.43), TRIMETHOPRIM(lift 1.41), PHENICOL(lift 1.27), AMINOGLYCOSIDE(lift 1.25)
- **PHENICOL** → BLEOMYCIN(lift 2.22), MULTIDRUG(lift 1.52), TRIMETHOPRIM(lift 1.37), MACROLIDE(lift 1.34)
- **QUINOLONE** → MULTIDRUG(lift 1.37), MACROLIDE(lift 1.27), FOSFOMYCIN(lift 1.25), PHENICOL(lift 1.23)
- **STREPTOTHRICIN** → TETRACYCLINE(lift 1.54), SULFONAMIDE(lift 1.43), AMINOGLYCOSIDE(lift 1.34), TRIMETHOPRIM(lift 1.32)
- **SULFONAMIDE** → BLEOMYCIN(lift 1.59), STREPTOGRAMIN(lift 1.59), MACROLIDE(lift 1.53), TRIMETHOPRIM(lift 1.5)
- **TETRACYCLINE** → STREPTOTHRICIN(lift 1.54), NITROFURAN(lift 1.43), MACROLIDE(lift 1.31), BLEOMYCIN(lift 1.3)
- **TRIMETHOPRIM** → BLEOMYCIN(lift 1.64), STREPTOGRAMIN(lift 1.64), MACROLIDE(lift 1.55), SULFONAMIDE(lift 1.5)
### campylobacter
- **AMINOGLYCOSIDE** → TETRACYCLINE(lift 1.45), QUINOLONE(lift 1.38), BETA-LACTAM(lift 0.97)
- **BETA-LACTAM** → TETRACYCLINE(lift 1.0), QUINOLONE(lift 0.99), AMINOGLYCOSIDE(lift 0.97)
- **QUINOLONE** → AMINOGLYCOSIDE(lift 1.38), TETRACYCLINE(lift 1.25), BETA-LACTAM(lift 0.99)
- **TETRACYCLINE** → AMINOGLYCOSIDE(lift 1.45), QUINOLONE(lift 1.25), BETA-LACTAM(lift 1.0)
### acinetobacter
- **MACROLIDE** → STREPTOGRAMIN(lift 3.16), TRIMETHOPRIM(lift 2.63), QUINOLONE(lift 1.07), TETRACYCLINE(lift 1.05)
- **STREPTOGRAMIN** → MACROLIDE(lift 3.16), TRIMETHOPRIM(lift 2.63), QUINOLONE(lift 1.07), TETRACYCLINE(lift 1.05)
- **TETRACYCLINE** → MACROLIDE(lift 1.05), STREPTOGRAMIN(lift 1.05), SULFONAMIDE(lift 1.03), QUINOLONE(lift 1.01)
- **TRIMETHOPRIM** → MACROLIDE(lift 2.63), STREPTOGRAMIN(lift 2.63), QUINOLONE(lift 1.07), SULFONAMIDE(lift 1.03)
### salmonella
- **AMINOGLYCOSIDE** → PHENICOL(lift 1.82), SULFONAMIDE(lift 1.82), TRIMETHOPRIM(lift 1.82), BETA-LACTAM(lift 1.74)
- **BETA-LACTAM** → PHENICOL(lift 2.31), QUINOLONE(lift 2.0), SULFONAMIDE(lift 2.0), TRIMETHOPRIM(lift 1.92)
- **PHENICOL** → TRIMETHOPRIM(lift 3.55), BETA-LACTAM(lift 2.31), SULFONAMIDE(lift 2.22), QUINOLONE(lift 2.0)
- **QUINOLONE** → BETA-LACTAM(lift 2.0), PHENICOL(lift 2.0), SULFONAMIDE(lift 2.0), TRIMETHOPRIM(lift 2.0)
- **SULFONAMIDE** → PHENICOL(lift 2.22), TRIMETHOPRIM(lift 2.22), BETA-LACTAM(lift 2.0), QUINOLONE(lift 2.0)
- **TETRACYCLINE** → SULFONAMIDE(lift 1.65), AMINOGLYCOSIDE(lift 1.61), PHENICOL(lift 1.58), TRIMETHOPRIM(lift 1.58)
- **TRIMETHOPRIM** → PHENICOL(lift 3.55), SULFONAMIDE(lift 2.22), QUINOLONE(lift 2.0), BETA-LACTAM(lift 1.92)

## Honest caveats
- Cross-axis (virulence/stress) NOT available from cached AMR-only runs; AMR drug-class only.
- Cohorts are drug-R/S-SELECTED -> the co-resistance structure reflects the curated cohorts.
- Within-organism de-confounds species; dedup is a crude clonality proxy (Mash-collapse is the follow-on).
- Class presence is DETERMINANT-implied (AMRFinder Class), not measured MIC per drug — a genotype axis.
- Self-distillation from our own caller: associational, NOT causal.