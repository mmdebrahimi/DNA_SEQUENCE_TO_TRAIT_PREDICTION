# Soraya --advance — result

- **Run:** 2026-06-08-gw38-wider-amr-campylobacter
- **Family:** eukaryotic-trait-decoding-cycle-2026-06-07
- **Verdict:** batch complete (4/4 actions) — stop reason: cap reached (planned batch done)

## What ran
Campylobacter × ciprofloxacin wider-AMR validation (`scripts/organism_drug_validate.py`,
`-O Campylobacter`, cohort 30 = 15R/15S, 30 AMRFinder Docker runs). First DIFFERENT-PHYLUM organism
(Campylobacterota vs prior Pseudomonadota).

## Result — FAILS_BAR (TUNING boundary, cleanest finding of the thread)
- Deployed cipro QRDR-POINT **threshold=2** rule: acc 0.500 / **sens 0.000** (calls all S).
- Per-strain perfectly clean: 15/15 R = exactly one `gyrA_T86I`; 15/15 S = zero.
- **threshold=1 → acc 1.000 / sens 1.000 / spec 1.000** (perfect).
- Mechanism caller transfers across the phylum boundary perfectly; only the count threshold is
  organism-specific (E. coli double-mutant tuning vs Campylobacter single gyrA T86I).

## Candidate (documented, NOT wired)
Per-organism cipro threshold (2 Enterobacterales / 1 Campylobacter). Biologically grounded → HIGH
generalization confidence. Wiring is a scoped follow-on (pin E. coli/Klebsiella regressions first).

## Artifacts
- `wiki/campylobacter_ciprofloxacin_validate_2026-06-08.{md,json}`
- Ledger row 30 (`project_state/eukaryotic-trait-decoding-cycle-2026-06-07.md`)
- Commit `b42a683`
- Memory updated: `feedback_intrinsic_genes_break_broad_amr_class_rules.md` (added TUNING boundary type)

## Next-VOI options (not auto-run)
- Wire the per-organism cipro threshold (small scoped code change + regression pins) — turns this
  finding into a product capability.
- More organisms: Pseudomonas × cipro (efflux-heavy, would test threshold + content together);
  Staphylococcus × oxacillin (gram-positive, mecA — but label-unreliability caveat).
- Standing user-gated: Path B G2 GPU run on the workhorse (the cycle's only remaining MVP gate).
