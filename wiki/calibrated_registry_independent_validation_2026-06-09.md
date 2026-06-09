# Calibrated AMR registry — INDEPENDENT-cohort out-of-sample validation — 2026-06-09

> The registry was calibrated IN-SAMPLE (N~30). This applies each in-sample rule to a DISJOINT
> second cohort (cohort-1 accessions excluded) + re-calibrates on it. PROMOTION GATE (opt-in ->
> default): OOS acc + sens + SPEC all >= 0.80 AND >= 10 scored strains per class. Config-MATCH is
> a FLAG, not required — a config that generalizes (non-inferior OOS perf) is eligible even if the
> re-calibration tie-break picks a different equally-good config. NCBI labels; AMRFinder pinned.

| organism | drug | scored R/S | in-sample cfg | OOS acc | OOS sens | OOS spec | re-cal cfg | cfg-match flag | PROMOTION |
|---|---|---|---|---:|---:|---:|---|---|---|
| Campylobacter | ciprofloxacin | 15/15 | qrdr_point@1 | 1.0 | 1.0 | 1.0 | qrdr_point@1 | match | ELIGIBLE |
| Klebsiella | ciprofloxacin | 15/15 | qrdr_point@2 | 0.967 | 0.933 | 1.0 | qrdr_point@1 | diff | ELIGIBLE |
| Salmonella | ciprofloxacin | 15/15 | broad@1 | 1.0 | 1.0 | 1.0 | qrdr_point@1 | diff | ELIGIBLE |

## Reading
- **OOS acc/sens/spec** = the in-sample registry rule applied to strains it was
  NOT calibrated on. All >= 0.80 (+ >=10 scored/class) => PROMOTION-eligible.
- **cfg-match flag** = whether re-calibration picks the SAME (counter,threshold). 'diff' is NOT a
  failure — a non-inferior config that generalizes is still eligible (match is a flag, not a gate).

## Honest scope
Second cohort is disjoint by accession but same label source (NCBI AST).
A held-out NCBI cohort is a stronger test than in-sample but still not a different-lab study.
Promotion of any config from opt-in to default remains a deliberate decision on this evidence.