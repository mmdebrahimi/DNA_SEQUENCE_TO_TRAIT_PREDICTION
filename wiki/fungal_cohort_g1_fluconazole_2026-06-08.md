# Fungal AMR cohort validation (Gate G1) — C. auris fluconazole

> Generated 2026-06-08. Label table: `cauris_g1_subset.assembled.tsv`. Breakpoint: MIC >= 32.0 ug/mL = R (CDC tentative; no formal CLSI/EUCAST C. auris breakpoint exists).
> Caller: deterministic ERG11/FKS1 target-site scan (blastn vs the committed real C. auris reference). Clade = the de-confounding variable (C. auris analogue of bacterial lineage).

## Verdict: **LABEL_LIMITED_FAILURE**

| metric | value |
|---|---|
| isolates (total) | 24 |
| scored (excl. INDETERMINATE/unlabelable) | 24 |
| accuracy | 0.792 |
| sensitivity (recall on MIC-R) | 1.000 |
| specificity | 0.167 |
| TP / TN / FP / FN | 18 / 1 / 5 / 0 |
| INDETERMINATE (no BLAST) | 0 |

## Within-clade de-confounding (per-clade confusion)

| clade | TP | TN | FP | FN | n |
|---|---|---|---|---|---|
| I | 4 | 0 | 0 | 0 | 4 |
| III | 14 | 1 | 5 | 0 | 20 |

## Efflux/aneuploidy discordance (the documented blind spot)

Isolates MIC-R but ERG11-S = candidate non-target resistance (CDR1_efflux_overexpression, ERG11_copy_number_aneuploidy, ERG3_loss_of_function, MDR1_efflux_overexpression, TAC1b_efflux_regulator): (none)

A deterministic target-site scan CANNOT detect these by design. A sub-0.80 sensitivity driven entirely by this set is the documented failure mode (falsifier), not a caller defect.

## Per-isolate

| isolate | clade | MIC | true | pred | bucket | determinants |
|---|---|---|---|---|---|---|
| 3758 | III | 4 | S | S | TN | - |
| 3561_77 | III | 16 | S | R | FP | ERG11:F126L |
| 6057 | III | 16 | S | R | FP | ERG11:F126L |
| 2566 | III | 16 | S | R | FP | ERG11:F126L |
| 4591 | III | 16 | S | R | FP | ERG11:F126L |
| 5268 | III | 16 | S | R | FP | ERG11:F126L |
| 2547_51 | III | 32 | R | R | TP | ERG11:F126L |
| 2546_62 | III | 32 | R | R | TP | ERG11:F126L |
| 4926_60 | III | 32 | R | R | TP | ERG11:F126L |
| 2859_37 | III | 64 | R | R | TP | ERG11:F126L |
| 4928_61 | III | 64 | R | R | TP | ERG11:F126L |
| 5624_2 | III | 64 | R | R | TP | ERG11:F126L |
| 6241_97 | III | 128 | R | R | TP | ERG11:F126L |
| 3406_32 | III | 128 | R | R | TP | ERG11:F126L |
| 4603_48 | III | 128 | R | R | TP | ERG11:F126L |
| 4845_72 | III | 128 | R | R | TP | ERG11:F126L |
| 5006_96 | III | 128 | R | R | TP | ERG11:F126L |
| 5625_3 | III | 128 | R | R | TP | ERG11:F126L |
| 5714 | III | 256 | R | R | TP | ERG11:F126L |
| 5734_5 | III | 256 | R | R | TP | ERG11:F126L |
| 4934_87 | I | 128 | R | R | TP | ERG11:Y132F |
| 4000 | I | 128 | R | R | TP | ERG11:Y132F |
| 5116_74 | I | 128 | R | R | TP | ERG11:Y132F |
| 5233_75 | I | 256 | R | R | TP | ERG11:Y132F |
