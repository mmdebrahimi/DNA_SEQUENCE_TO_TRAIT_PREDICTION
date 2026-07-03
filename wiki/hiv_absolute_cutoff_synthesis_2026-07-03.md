# HIV v0.2 absolute-cutoff calibration — PI calibrated, NNRTI confirmed, INSTI walled (2026-07-03)

**The deferred "v0.2 absolute calibration" item, done by SOURCING (not fabricating) per-drug clinical
cutoffs.** The PI/NNRTI validators reported an *illustrative* uniform fold≥3 sens/spec; the CLAUDE.md flagged
that per-drug clinical cutoffs "would upgrade to absolute → a v0.2 item, **not fabricated**." This sources
them from Stanford HIVDB's own `DRMcv.R` (`cutoffmat`) — the SAME authoritative script the NRTI cell used —
and re-scores the frozen catalog at each drug's real lower cutoff. Script:
`scripts/hiv_absolute_cutoff_validate.py`.

## Feasibility gate (verified against the fetched DRMcv.R)

| Class | cutoffs in DRMcv.R | v0.2 status |
|---|---|---|
| **PI** | all 8 (FPV/ATV/IDV/SQV 3, **LPV 9**, NFV 3, **TPV 2**, **DRV 10**) | **CALIBRATED** — the genuine new content |
| **NNRTI** | EFV/NVP/ETR/RPV = 3; DOR absent | **CONFIRMED** (cutoffs = the prior illustrative 3); **DOR walled** |
| **INSTI** | none — integrase inhibitors postdate the script | **CUTOFF_UNAVAILABLE** (external wall, not guessed) |

## Result (frozen catalog, PhenoSense fold label, absolute sens/spec/balacc at the DRMcv.R lower cutoff)

**PI (8/8, position-based v0):** FPV 0.786 · ATV 0.841 · IDV 0.843 · LPV 0.789 · NFV **0.880** · SQV 0.789 ·
TPV 0.743 · DRV 0.746 (balacc). Sensitivity is near-perfect (0.986–1.0); **specificity is the story** — it
ranges 0.49 (DRV) → 0.77 (NFV). The low-spec tail is the **deliberate position-based over-call** made
visible at the *real* cutoff (any major-position residue → R, including polymorphisms/revertants). This is
precisely the gap the mutant-specific PI v0.1 catalog (`hiv_pi_mutant_catalog`, 2026-06-23) was built to
close — v0.2 now quantifies it against clinical cutoffs instead of an illustrative one.

**NNRTI (4/5):** EFV 0.908 · NVP 0.944 · ETR 0.728 · RPV 0.724. Because every DRMcv.R NNRTI lower cutoff is
3, these equal the prior illustrative numbers — a **confirmation** that the illustrative fold≥3 was already
the clinical cutoff, not a change. DOR (doravirine) postdates DRMcv.R → walled.

**INSTI (0/5):** every drug CUTOFF_UNAVAILABLE. A free Monogram INSTI clinical cutoff exists in the
literature but is NOT in this canonical free source → external sourcing, deferred. Reported as a wall,
**never guessed**.

Within-B ≈ all for every calibrated drug (consistent with the within-subtype de-confounding finding of the
same day — the calibration is not subtype-driven).

## What it settles

- The PI v0 surface now has **clinically-interpretable absolute sens/spec** at authoritative cutoffs, and
  the position-based over-call is quantified (low spec) — the concrete case for the deployed v0.1
  mutant-specific catalog.
- The honesty rail held: **DOR + all INSTI were not assigned a fabricated cutoff** — the wall is reported.
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9); READ-only over the frozen
  `dna_decode.data.hiv_amr` dispatch. HIV's own trust surface, namespace-separate.

## Reproduce
```bash
uv run python scripts/hiv_absolute_cutoff_validate.py --class all   # PI + NNRTI + INSTI(walled)
uv run pytest tests/test_hiv_absolute_cutoff.py -q                  # 6 offline tests (no network)
# cutoffs: Stanford HIVDB DRMcv.R cutoffmat (GenoPhenoDatasets/DRMcv.R); .Full datasets gitignored
```
