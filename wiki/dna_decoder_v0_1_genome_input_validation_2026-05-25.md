# DNA Decoder v0.1 genome-input validation - 2026-05-25

## Panel

- `562.17721` / `GCA_002201835.1` / expected `R`
- `1438684.3` / `GCA_000692695.1` / expected `R`
- `562.50295` / `GCA_004568615.1` / expected `S`
- `562.7583` / `GCA_001277755.1` / expected `S`

## Results

| Strain | Expected | Cached | Genome-input | Cached p(R) | Genome p(R) | |?| | Cached audit gate | Genome audit gate |
|---|---|---|---|---:|---:|---:|---|---|
| 562.17721 | R | R | R | 0.862065 | 0.861944 | 0.000121 | SUSPEND_CONDITION_4 | SUSPEND_CONDITION_4 |
| 1438684.3 | R | R | R | 0.800415 | 0.803224 | 0.002809 | SUSPEND_CONDITION_4 | SUSPEND_CONDITION_4 |
| 562.50295 | S | S | S | 0.341645 | 0.353244 | 0.011599 | SUSPEND_CONDITION_4 | SUSPEND_CONDITION_4 |
| 562.7583 | S | S | S | 0.132160 | 0.132160 | 0.000000 | SUSPEND_CONDITION_4 | SUSPEND_CONDITION_4 |

## Summary

- prediction concordance: `4 / 4`
- max absolute probability delta: `0.011599`
- cached path and genome-input path both remained canonical audit-aware when an audit sidecar was supplied
- genome-input path preserves audit gate propagation on the same cohort strains
- validation run used `--no-attribution` to isolate core predict-path consistency