# Determinant-blindness atlas

_Generated 2026-07-21 from the NCBI-PD external-validation cohorts (`scripts/score_ncbipd_extval.ORGANISMS`). Per (organism, drug): of the measured-R isolates, the fraction the deployed cell calls **non-R** because no rule-firing catalog determinant is present._

**Why this exists.** The `/innovate` 2026-07-21 run KILLED the tempting move to *rescue* this false-negative ceiling by scoring the filtered determinants — it dies for expression-driven resistance (azithromycin mtr-efflux, burden gap −0.6) and is a pure **clonal confound** for cumulative-chromosomal resistance (tetracycline pooled burden +6.2 collapses to −1.0 within SNP clusters). The surviving move is to **disclose** the blindness, not hide it. This table is DESCRIPTIVE (not a predictor) → immune to that clonal confound.

**Read it as:** a high `invisible fraction` means the determinant catalog structurally cannot see much of this cell's resistance — a known mechanism gap (efflux / regulatory / porin-loss / multi-locus-cumulative), not a rule bug. `truly-invisible` = the genome carries **zero** determinant token; `rule-limited` = a determinant is present but not one the rule counts. **No aggregate headline** — each cell stands alone.

| organism | drug | n R | **invisible fraction** | invisible (truly / rule-limited) |
|---|---|---|---|---|
| Neisseria gonorrhoeae | azithromycin | 110 | **1.0** | 110 (0 / 110) |
| Neisseria gonorrhoeae | ceftriaxone | 2 | **1.0** | 2 (0 / 2) |
| Neisseria gonorrhoeae | tetracycline | 34 | **0.676** | 23 (0 / 23) |
| Neisseria gonorrhoeae | cefixime | 19 | **0.211** | 4 (0 / 4) |
| Staphylococcus aureus | rifampin | 8 | **0.125** | 1 (0 / 1) |
| Neisseria gonorrhoeae | penicillin | 14 | **0.071** | 1 (0 / 1) |
| Campylobacter | gentamicin | 31 | **0.032** | 1 (0 / 1) |
| Streptococcus pneumoniae | erythromycin | 43 | **0.023** | 1 (0 / 1) |
| Neisseria gonorrhoeae | ciprofloxacin | 94 | **0.011** | 1 (0 / 1) |
| Campylobacter | tetracycline | 66 | **0.0** | 0 (0 / 0) |
| Staphylococcus aureus | ciprofloxacin | 22 | **0.0** | 0 (0 / 0) |
| Streptococcus pneumoniae | tetracycline | 10 | **0.0** | 0 (0 / 0) |

## Notes
- **Descriptive honesty surface, not a predictor** — it reports where the catalog is blind; it does NOT attempt to call those isolates (the `/innovate` burden-rescue that would have is a closed negative: clonal / no-signal).
- A high invisible fraction is expected + honest for known determinant-invisible mechanisms: gono azithromycin (mtr-efflux, 23S-independent — 100% invisible on this cohort), gono tetracycline (chromosomal cumulative), Klebsiella meropenem (porin-loss).
- Reuses the NCBI-PD substrate (provenance-disjoint, NOT methodology-independent — same AMRFinderPlus + same cell). NON-FROZEN cells; the frozen decoder surface is byte-unchanged.
