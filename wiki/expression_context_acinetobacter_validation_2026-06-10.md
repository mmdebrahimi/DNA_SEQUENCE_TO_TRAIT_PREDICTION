# expression_context independent validation — Acinetobacter x meropenem — 2026-06-10

PRIMARY detector (frozen falsifier rule) on the INDEPENDENT cohort (disjoint from the cached 30).
**Mechanism-stratified**: the signal targets INTRINSIC-ONLY-R (R with no strong acquired carbapenemase);
acquired-carbapenemase R are resistant via a gene-presence-visible mechanism the signal never targets.

## Verdict: HOLD

| metric | value |
|---|---|
| n evaluated | 30 (15R/15S) |
| **n_target_R** (intrinsic-only-R = the signal's target) | **1** |
| n_acquired_R (acquired carbapenemase = NON-target) | 14 |
| **target_R_rescues** (intrinsic-only-R upgraded ABSTAIN->R) | **0** |
| r_rescues (ALL true-R upgraded — incl. non-target) | 0 |
| s_upgrades (false-R; must be 0) | 0 |
| target_rescue_rate | 0.0 |
| false-upgrade Wilson95 upper | 0.2039 |

Gate: PROMOTE iff s_upgrades==0 AND target_R_rescues>=1 AND n_target_R>=10 AND n_S>=15. HOLD: n_target_R=1 (<10) — UNDERPOWERED on the intrinsic-only-R target subset; 14/15 R are acquired-carbapenemase (non-target). NOT a falsification of the signal — the cohort cannot test it.

**Honest reading:** this is NOT a falsification of the signal. The independent cohort is CONFOUNDED — 14/15 R carry strong acquired carbapenemases (non-target), leaving only n_target_R=1 intrinsic-only-R, far below the 10 needed to test generalization. Intrinsic-only carbapenem-R Acinetobacter is rare in sequenced collections (acquired OXA-23 dominates), so the signal's value on its target subset is not independently testable at adequate power on free opportunistic NCBI AST cohorts found so far. The override stays opt-in/off-by-default. Raw BLAST hits + nearest-distance + FASTA/ref SHA256 in the JSON sidecar for reproducibility.
