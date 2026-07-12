# Genome world-model findings index (2026-07-12)

A single navigational index of the non-frozen world-model research lane — the structure the additive AMR
catalog discards (interactions / quantitative / joint / self-awareness), all extracted from data in hand
(no new labels, GPU, or embeddings; the frozen decoder surface is byte-unchanged across every finding).
Companion to `quantitative_decoder_capability_map_2026-07-12.md` (which indexes calibrated CELLS); this
indexes FINDINGS.

## The findings

| # | finding | verdict | headline number | artifact |
|---|---|---|---|---|
| A | HIV drug-resistance **epistasis** (do mutation interactions beat additive?) | **FAIL_ADDITIVE_SUFFICES** | interactions beat additive on only 2/24 drug-cells; top interaction coefs recover textbook synergy (G140S+Q148H, TAMs, PI 82A+84V) → epistasis is real but rank-redundant | `hiv_epistasis_result_2026-07-11` |
| B | **Quantitative calibration** — coverage-valid prediction INTERVALS (not just R/S) | **CALIBRATED** | 30/34 powered cells INFORMATIVELY calibrated (24 HIV PhenoSense fold-change + 6 TB MIC drugs); coverage 0.893–0.928 | `quantitative_decoder_capability_map_2026-07-12`, `hiv_quantitative_calibration_2026-07-11`, `tb_mic_calibration_panel_2026-07-12` |
| C-core | AMR **determinant co-occurrence** (cassette/integron linkage) | **PASS_LINKAGE_STRUCTURE** | recovers integron sul1→aadA/dfrA cassettes + QRDR clusters + Salmonella floR→aac(3)-IVa; phenotype→genotype inversion works | `determinant_cooccurrence` (2026-07-11) |
| C-deep | class-level **co-resistance imputation** | **PASS_CORESISTANCE_IMPUTABLE** | AUC 0.9–0.97 | `coresistance_imputation` (2026-07-11) |
| C-multiaxis | AMR × **plasmid** × **virulence** joint network | **PASS_MULTIAXIS_LINKAGE** | AMR+plasmid predicts virulence genes AUC 0.79–0.95 (raw); see lineage de-confound below for the honest half | `coresistance_multiaxis`, `virulence_axis_sweep`, `plasmid_axis_sweep` |
| D | deterministic **blind-spot self-awareness** flag | **FLAG_RECOVERS_BLINDSPOT** | a position-novelty flag catches 60% of the HIV EFV catalog blind spot the learned rescue can't (mutant-level-catalog-specific) | `hiv_blindspot_position_novelty` (2026-07-11) |
| L | **lineage de-confound** of the joint network (Mash-clade leave-one-clade-out) | **SPLIT / LINEAGE / GENERALIZES** | determinant↔determinant GENERALIZES 0.908 (E.coli)/0.913 (Kleb) > virulence SPLIT 0.676 > plasmid LINEAGE 0.615; **prediction "mobile plasmids generalize" FALSIFIED** | `genome_world_model_lineage_structure_synthesis_2026-07-12`, `crossaxis_lineage_deconfound*_2026-07-12` |

## The two organizing insights

1. **The additive curated catalog wins even in the continuous regime.** (A) Interactions are rank-redundant
   with an additive model for HIV DR; (B) a determinant-presence model produces coverage-valid MIC/fold
   intervals. The structure the catalog "discards" (epistasis) does not add predictive rank — it confirms the
   catalog is near-complete for the well-studied drugs.

2. **Resistance gene CONTENT is transferable across lineages; the VEHICLE and CONTEXT are lineage-locked.**
   (L) Determinant↔determinant co-occurrence generalizes across held-out Mash clades (0.9+, in E. coli AND
   Klebsiella — specifically the ACQUIRED content; intrinsic/chromosomal genes like fosA collapse), but the
   plasmid backbone and accessory virulence context are clade-fixed. Converges with the ST131 literature.

## Honesty rails that shaped the findings (reusable)

- **Coverage-valid ≠ informative.** Conformal coverage holds even for a useless model (interval → marginal
  spread); the `informative` flag (R²>0.05) separates the TB drugs whose determinants actually predict MIC
  (RIF/INH/EMB/LEV/MXF/ETH) from the coverage-valid-only ones (AMI/KAN/BDQ/LZD).
- **Cross-axis co-occurrence must be lineage-controlled.** The dedup proxy under-controls; leave-one-Mash-
  clade-out is the honest test, and each axis has its OWN verdict (not one global lineage-mediation call).
- **Pre-register the mechanistic prediction, then let the de-confound falsify it** ("mobile ⇒ generalizes"
  was false — plasmid backbones are clade-fixed by stable co-inheritance).
- **The de-confound recovers intrinsic-vs-acquired with no annotation** — cross-lineage generalization AUC is
  a label-free proxy for "is this determinant mobile/acquired vs intrinsic/chromosomal."

## Scope

E. coli/Shigella + Klebsiella (joint axes) + HIV-1 + M. tuberculosis (quantitative). Associational, in-distribution
vs the relevant knowledge base (not independent validation). All local/CPU; frozen decoder surface untouched.
