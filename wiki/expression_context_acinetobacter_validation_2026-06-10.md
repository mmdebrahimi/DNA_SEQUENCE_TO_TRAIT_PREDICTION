# expression_context independent validation — Acinetobacter x meropenem — 2026-06-10

PRIMARY detector (frozen falsifier rule) on the INDEPENDENT cohort (disjoint from the cached 30).

## Verdict: HOLD

| metric | value |
|---|---|
| n evaluated | 30 (15R/15S) |
| r_rescues (true-R upgraded ABSTAIN->R) | 0 |
| s_upgrades (false-R; must be 0) | 0 |
| abstain_rescue_rate | 0.0 |
| false-upgrade rate | 0.0 |
| false-upgrade Wilson95 upper | 0.2039 |

Gate: PROMOTE iff s_upgrades==0 AND r_rescues>=1 AND n_S>=15. HOLD: r_rescues=0 (detector inert — no true-R upgrade)

Promotion is opt-in only; default-on deferred until n_S materially larger (the Wilson upper bound is still wide at small n_S). Raw BLAST hits + FASTA/ref SHA256 in the JSON sidecar for reproducibility.
