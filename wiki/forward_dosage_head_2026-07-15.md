# Forward-cell DOSAGE head — calibrated magnitude + prediction intervals (2026-07-15)

A new output modality for the forward variant-effect cell: turn a rank-score into a **calibrated MAGNITUDE**
prediction with honest uncertainty. The cell's methods (BLOSUM/ESM2/AlphaMissense/ESM-IF) produce a score
whose RANK correlates with the measured effect; the dosage head maps that score to the measured-effect scale
(monotone isotonic calibrator) and wraps it in a **split-conformal prediction interval** with a
pre-registered coverage target — so the decoder can say "this edit's effect is X ± q" with verified coverage,
not just "ranks low." Closer to the north star (a decoder should state magnitude + uncertainty).

## Real result — PTEN, AlphaMissense score → measured DMS magnitude

`scripts/forward_dosage_cell.py`, 7,260 variants, 20 shuffled 50/25/25 fit/calib/test splits, target 80%:

| metric | value | reading |
|---|---:|---|
| **mean held-out coverage** | **0.7993** | target 0.80 — essentially perfect calibration (\|Δ\|=0.0007) |
| mean interval half-width | 1.28 DMS units | interval = calibrated_point ± 1.28 |
| **interval narrowing vs marginal** | **0.20** | the score shrinks the interval 20% vs predict-the-mean |
| mean point Spearman | 0.534 | consistent with AlphaMissense's rank 0.539 |
| **verdict** | **CALIBRATED_DOSAGE** | calibrated AND informative |

## The load-bearing honesty rail (reused from J2's MIC-calibration lesson)

**Split-conformal coverage is guaranteed even for a USELESS model** — the interval simply widens to the
marginal distribution to hit nominal coverage. So coverage alone does NOT prove the score is informative.
The dosage head therefore ALSO reports `interval_narrowing = 1 − q/marginal_q` (how much the score's
conditioning shrinks the interval vs a no-features predict-the-mean baseline), and the verdict requires BOTH
nominal coverage AND narrowing > 0.02. On PTEN the AM score narrows the interval by 0.20 — the calibrated
magnitude reflects real variant-effect signal, not just the marginal spread. (Two offline tests pin this: an
uninformative random x still hits coverage 0.80 but narrows ~0, while an informative x narrows ≫0.)

## Non-duplication (R4)

J2's session already built conformal MIC intervals for the AMR/TB *determinant* decoder
(`scripts/tb_mic_calibration.py` + `hiv_quantitative_calibration._conformal_q`). This dosage head is a
DIFFERENT substrate (forward variant-effect → DMS magnitude, no censoring) but **reuses J2's split-conformal
definition** — `dosage.conformal_q` is byte-equal to `hiv_quantitative_calibration._conformal_q`, asserted
both in the runner and in `tests/test_forward_dosage.py::test_conformal_q_matches_j2_helper`. Same math,
in-package for clean deps; not a re-invention.

## What this adds

The forward cell now has a **magnitude output**, not just a rank: `dna_decode/forward/dosage.py`
(`evaluate_dosage` / `dosage_intervals` / `conformal_q`), method-agnostic (any forward-cell score calibrates
the same way), validated on real DMS with 80% held-out coverage + a genuine informativeness check. 6 offline
tests. Frozen decoder surface (`amr_rules` / `calibrated_amr_rules` / `mic_tiers` / `shipped_decoder_surface`
/ `cohort_manifest`) byte-unchanged (`verify_lock OK`); `dna_decode/forward` NON-frozen. Run:
`uv run python scripts/forward_dosage_cell.py`.
