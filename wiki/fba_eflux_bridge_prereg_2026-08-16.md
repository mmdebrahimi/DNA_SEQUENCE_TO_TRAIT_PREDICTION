# E-Flux bridge — PRE-REGISTRATION (written BEFORE any result was computed)

**Date:** 2026-08-16 · **Status:** frozen before execution

This epoch has already produced two self-inflicted verdict reversals from post-hoc choices (a
`min_abs_t` filter that was anti-selective for switchers; a verdict function that compared observed
against nulls instead of against the baseline arm). So the bar is frozen here, in writing, before the
first solve.

## The question

Does **expression-constrained FBA (E-Flux)** improve conditional gene-essentiality prediction over
plain FBA on the same carbon panel?

This is the "bridge experiment" — the structural gap named in
`wiki/genotype_to_phenotype_map_2026-08-13.md` (prediction works at short causal distance, fails
through the network, nothing connects the layers).

## Substrate — and an honest downgrade from this morning

`wiki/fba_bridge_precise_overlap_2026-08-16.json` recorded a **16/25** overlap with PRECISE-NP881's
sole-carbon arm and concluded the blocker was "fully a CODE wall (~1 day of E-Flux work)."

**That conclusion was wrong and is corrected here.** The 346 new NP881 profiles are **not** published
as a processed expression matrix: absent from the `SBRG/precise1k` repo (only `master`, no NP881
directory), absent from GEO (0 hits for either BioProject), no Zenodo deposition. Only **raw SRA
reads** exist. Using NP881 therefore requires downloading and quantifying ~32–82 RNA-seq runs
(~50–250 GB + an alignment pipeline) — not a day of E-Flux.

The substrate-overlap count (16/25) **stands**; it came from SRA run metadata. What was wrong was the
cost/availability inference layered on top of it.

**So this test runs on the panel that is actually available today:** the **11 of our 25 conditions
covered by PRECISE-1K** (`log_tpm_qc.csv`, 60.9 MB, already cached to `D:`), all with ≥2 samples.

## Arms — both on the SAME 11 conditions

The shipped **0.7368 per-cell figure is NOT a valid comparator** here: it was computed on 25
conditions. Comparing an 11-condition E-Flux number against a 25-condition baseline would be the same
class of error as the regulatory-arm verdict bug. Therefore:

| arm | definition |
|---|---|
| **A — baseline** | plain FBA, `apply_carbon_condition` only, recomputed on these 11 conditions |
| **B — E-Flux** | identical, plus expression-derived reaction bounds from PRECISE-1K |

Everything else is held fixed: iML1515, Fitness Browser RB-TnSeq labels (orgId=Keio), essentiality
threshold `fit < -2.0`, `FRAC = 0.01`, the same conditionally-essential two-sided gene subset.

## Frozen decision rule

- **Primary metric:** per-cell agreement (the shipped metric).
- **GO signal:** `per_cell(B) − per_cell(A) ≥ +0.02`.
- **Machado prior confirmed:** `Δ ≤ 0` — i.e. expression constraints do not help, or actively hurt,
  conditional essentiality. Machado & Herrgård 2014 found plain FBA + parsimony as good or better,
  but for **flux**; this would extend it to essentiality.
- **Ambiguous:** `0 < Δ < +0.02`.
- **Committed:** the result is reported in whichever direction it lands, including a null or a
  regression. Secondary readouts (exact-set match, MCC, per-condition) are descriptive only and
  **cannot** be promoted to the headline after the fact.

## E-Flux implementation (frozen before running)

Colijn et al. 2009. Per condition: mean log-TPM across that condition's PRECISE-1K samples → linear
TPM → GPR evaluation (**AND = min**, **OR = sum**) → per-reaction expression score → normalized by the
99th percentile → scales each reaction's bound. Reactions with no GPR, and genes absent from the
expression matrix, are left **unconstrained** (never silently zeroed). Carbon-source exchange bounds
are **exempt** — constraining uptake by expression would confound the condition itself.

## Named caveats (recorded now, not discovered later)

1. **Strain mismatch.** PRECISE-1K is K-12 MG1655; the Fitness Browser labels are orgId=Keio
   (BW25113 parent). Both K-12 derivatives; standard practice, but it is a real seam.
2. **Sample imbalance.** Glucose (621 samples) and glycerol (111) dominate; the other 9 conditions
   have 2–8. The mean expression profile is far better determined for two conditions than for nine.
3. **11 < 25.** This panel is under half the shipped one; a null here does not by itself close the
   bridge on the full panel.
4. **A negative result is a real answer.** Given caveat 2 and the Machado prior, `Δ ≤ 0` is a
   genuinely likely outcome and is pre-committed as informative, not as a failed run.
