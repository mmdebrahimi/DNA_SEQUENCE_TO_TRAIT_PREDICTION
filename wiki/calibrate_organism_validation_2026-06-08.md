# calibrate_organism — real-data validation — 2026-06-08

> Validates `dna_decode/eval/calibrate_organism.py` against the cached AMRFinder cohorts from the wider-AMR
> thread. The meta-rule auto-selects (counter, threshold) + intrinsic-family exclusions from a labeled
> cohort by leave-one-out balanced accuracy, and ABSTAINS (verdict EXPRESSION_FLOOR) when no presence-based
> config clears the LOO floor (0.70). Pure logic over cached runs — no Docker, no GPU, no money.
> JSON: `wiki/calibrate_organism_validation_2026-06-08.json`. Unit tests: `tests/test_calibrate_organism.py` (16).

## Result — recovers the right config on every CONTENT/TUNING case; abstains on EXPRESSION

| cohort | drug | verdict | counter@thr | intrinsic excluded | LOO bal-acc | full acc | deployed-rule acc |
|---|---|---|---|---|---:|---:|---:|
| Campylobacter | cipro | **CALIBRATED** | qrdr_point@**1** | — | **1.000** | 1.000 | 0.500 (TUNING miss) |
| Klebsiella | cipro | **CALIBRATED** | qrdr_point@**2** | oqxA, oqxB | **1.000** | 1.000 | 1.000 |
| Salmonella | cipro | **CALIBRATED** | **broad**@1 | — | **1.000** | 1.000 | **0.567** (CONTENT miss) |
| Acinetobacter | meropenem | **EXPRESSION_FLOOR** (abstain) | broad@3 | blaADC, blaOXA-51-family | 0.133 | 0.433 | 0.500 |
| Pseudomonas | meropenem | **EXPRESSION_FLOOR** (abstain) | broad@1 | — | 0.000 | 0.000 | 0.500 |

## What each row demonstrates

- **Campylobacter (TUNING)** — auto-recovers threshold **1** (single gyrA T86I), where the E. coli-tuned
  deployed threshold 2 scored 0.500. ✓
- **Klebsiella (CONTENT, intrinsic efflux)** — auto-recovers qrdr_point@**2** AND auto-excludes the
  intrinsic **oqxA/oqxB** efflux families (≥90% prevalence in both R and S). Matches the hand-built rule
  (1.000) with zero hand-tuning. ✓
- **Salmonella (CONTENT, counter choice)** — auto-switches to the **broad** counter @ threshold 1 to catch
  plasmid qnr, lifting LOO to **1.000** where the deployed QRDR-point-only rule scored **0.567**. This is
  the key proof: the meta-rule selects the COUNTER, not just the threshold — recovering the exact
  organism-specific choice the Klebsiella-vs-Salmonella contrast showed was unavoidable. ✓
- **Acinetobacter + Pseudomonas (EXPRESSION floor)** — the meta-rule still excludes the intrinsics
  (blaOXA-51-family, blaADC) but NO presence config clears the LOO floor (0.70), so it returns verdict
  **EXPRESSION_FLOOR** and recommends ABSTAIN. This is the honest, correct outcome: these organism×drugs are
  expression-driven (ISAba1→OXA-51 overexpression / efflux up-regulation / oprD porin loss), which
  gene-presence fundamentally cannot decode. The tool says "flag, do not predict" instead of shipping a bad
  rule. ✓

## Honest scope / caveats

- LOO balanced accuracy is the honest generalization estimate; intrinsic families are computed per-fold in
  the LOO and once on the full cohort for the returned rule (documented in-sample on N≈30).
- The CALIBRATION_FLOOR (0.70) is a single tunable that decides CALIBRATED vs abstain; it is conservative.
- Acinetobacter's hand-built strength-tier refinement reached acc 0.833 by ALSO excluding the conditional
  OXA-58-like (carried by 9 S + 2 R) — a distinction FINER than the ≥90%-both-classes intrinsic flag can
  auto-discover. So calibrate honestly abstains here rather than matching the hand-curation; closing that
  gap needs strength-tier curation or expression inference, both out of presence-based scope.
- Validates the H3 hypothesis (`calibrate_organism(cohort)` ports the decoder from a ≥15R/15S cohort): on
  3 organisms it matches/beats the hand-curated rule with no tuning; on 2 it correctly abstains. NOT yet
  wired into the deployed `call_resistance` path — this is the validated building block; wiring (let a
  calibrated rule override the default for a recognized organism) is the follow-on.
