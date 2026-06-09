# Soraya --advance — result

- Run: 2026-06-09-eqty-calibrate-brainstorm-fixes
- Family: eukaryotic-trait-decoding-cycle-2026-06-07
- Stop: batch complete (8 actions, under cap 10) · commit 96bffc4 · NOT a plateau (fixed shipped bugs + unblocked promotion)

## Applied the /brainstorm fix batch (2 critical + 3 medium)
- C1 degenerate cohort: INSUFFICIENT_EVIDENCE verdict + MIN_CLASS_COUNT guard + reuse_glob load fix -> registry rebuilt, all 5 entries on valid 15R/15S (Pseudomonas no longer the bogus 7-all-S EXPRESSION_FLOOR).
- C2: loo_balanced_accuracy truly balanced (None if class absent).
- M1: deterministic full-cohort deployed config (no modal/tie-break ambiguity).
- M3: promotion gate = spec floor + min-10/class + config-match-as-flag -> all 3 cipro configs promotion_eligible=True.
- M2: documented + deferred.
- 3 new tests; full suite 947 passed, 0 regressions.

## Next-VOI (not auto-run)
1. Promote the 3 cipro configs opt-in -> default (user sign-off; the eligibility evidence is in).
2. M2: AMRFinder Method-column propagation (latent robustness).
3. Surface organism= + ABSTAIN in the predict CLI.
4. Standing user-gated: Path B G2 GPU run on the workhorse.
