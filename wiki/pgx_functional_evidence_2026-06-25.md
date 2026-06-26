# PGx independent functional-evidence cross-check (2026-06-25)

_Independent (non-CPIC) functional-evidence cross-check of each PGx allele's CPIC function assignment -- the circularity-break. AGREE raises confidence; DISAGREE/FLAG surface where 'faithful-to-CPIC' rests on clinical evidence sequence-signals miss._

**Verdicts:** AGREE 4 / DISAGREE 1 / FLAG 1 / NO_SIGNAL 0  (n=6)

| gene | allele | rsid | CPIC function | variant class | independent signal | verdict |
|---|---|---|---|---|---|---|
| CYP2C19 | *2 | rs4244285 | no function | synonymous_cryptic_splice | Ensembl VEP most-severe = synonymous_variant; the no-function mechanism is a cryptic splice site (aberrant splicing, de Morais 1994) NOT captured by the consequence class | **FLAG** |
| CYP2C19 | *3 | rs4986893 | no function | stop_gained | Ensembl VEP most-severe = stop_gained (p.W212X) -> unambiguous loss of function | **AGREE** |
| CYP2C19 | *17 | rs12248560 | increased function | regulatory_promoter | Documented cis-regulatory promoter variant INCREASING CYP2C19 expression via a new GATA site (Sim 2006); GTEx-liver-significant-eQTL not resolved via the v2 API this run | **AGREE** |
| CYP2C9 | *2 | rs1799853 | decreased function | missense | Ensembl VEP: PolyPhen probably/possibly_damaging, SIFT deleterious(low-conf)/tolerated mixed (p.R144C) -> damaging-leaning | **AGREE** |
| CYP2C9 | *3 | rs1057910 | no function | missense | Ensembl VEP: PolyPhen BENIGN, SIFT mixed (p.I359L) -> predictors UNDER-call | **DISAGREE** |
| VKORC1 | -1639A | rs9923231 | decreased expression (warfarin sensitivity) | regulatory_promoter | Documented cis-regulatory promoter variant LOWERING VKORC1 expression (Rieder 2005); GTEx-liver-significant-eQTL not resolved via the v2 API this run | **AGREE** |

## Notes (the informative cases)

- **CYP2C19 *2 (FLAG):** A consequence-only predictor would UNDER-call this (synonymous); CPIC no-function rests on the documented splice defect -- a real flag, not a defect of the cell.
- **CYP2C9 *3 (DISAGREE):** In-silico predictors call this conservative I->L substitution benign, but CPIC assigns NO function on clinical/functional evidence -- a case where faithful-to-CPIC rests on MORE than sequence prediction. The honest value of this layer.

_Independent signals are ORTHOGONAL to CPIC curation, NOT ground truth: missense -> Ensembl VEP predictors (ML); stop/splice -> consequence class (fact); regulatory -> documented expression effect (measured, primary literature). Small-N per-allele annotation, NOT a concordance-%; GTEx-eQTL confirmation deferred (not asserted)._
