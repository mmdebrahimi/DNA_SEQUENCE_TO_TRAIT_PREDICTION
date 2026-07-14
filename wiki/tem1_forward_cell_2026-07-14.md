# TEM-1 forward variant-effect cell — the first "edit E. coli → predict phenotype" cell, DMS-validated (2026-07-14)

The decoder's FORWARD direction: make a minor edit (any point mutation in an E. coli protein) → predict the
change in phenotype → validate that exact variant against a FREE, independent, per-variant wet-lab label.
Flagship = **TEM-1 β-lactamase** (`BLAT_ECOLX`), phenotype = **ampicillin growth fitness** — E. coli + an
AMR phenotype (the project's validated home turf) + the fitness-aligned molecular regime (Regime B, where
substitution-severity + learned variant-effect predictors WORK).

## Why this is the right first forward cell (the label wall does NOT bind here)

The project's binding constraint everywhere else is LABELS not models. Deep-mutational-scan (DMS) data breaks
that for the forward/molecular direction: every single point mutation's measured effect on fitness is a free,
independent label. ProteinGym ships 4 TEM-1 β-lactamase DMS assays + a broad E. coli panel, already cached at
`D:/dna_decode_cache/proteingym`. So "edit → predict → validate" is a complete, free loop on E. coli.

## Result — deterministic BLOSUM62 baseline, real per-variant validation

`scripts/tem1_forward_cell.py` scores every single-point variant with the Regime-B predictor
(`dna_decode/forward/variant_effect.py`, BLOSUM62) and correlates against the measured DMS fitness. Signed
Spearman (positive = conservative substitution → preserved fitness), WT residue verified against the
reference (a mismatch is a coordinate/frame error and is COUNTED, never silently scored):

| Assay | phenotype | n single variants | WT-mismatch | Spearman(BLOSUM, DMS) | polarity |
|---|---|---:|---:|---:|---|
| **BLAT_ECOLX_Stiffler_2015** | ampicillin growth (10–2500 µg/mL) | 4,996 | 0 | **+0.347** | ✓ |
| BLAT_ECOLX_Firnberg_2014 | ampicillin fitness | 4,783 | 0 | +0.339 | ✓ |
| BLAT_ECOLX_Deng_2012 | ampicillin fitness | 4,996 | 0 | +0.306 | ✓ |
| DYR_ECOLI_Thompson_2019 (folA/DHFR) | growth (turbidostat) | 2,363 | 0 | +0.151 | ✓ |

**Reproducible ~0.31–0.35 across 3 independent TEM-1 assays**, positive direction confirmed, zero
coordinate-frame errors. folA/DHFR is lower (0.151) — a harder landscape for pure substitution severity.

## UPDATE 2026-07-14 — both named next-lifts LANDED (ESM2 method + genome-level path)

**Lift 1 — ESM2-650M learned upgrade (DONE).** `dna_decode/forward/esm_scorer.py::esm2_logp_table` runs the
650M model once per protein (masked-marginals, CPU, cached to D:) → `predict_effect(method="esm2", esm_table=)`.
Real result on the same 4,996-variant TEM-1 Stiffler DMS:

| method | Spearman(pred, DMS) | n | WT-mismatch | "preserved" precision |
|---|---:|---:|---:|---:|
| BLOSUM62 (deterministic) | +0.347 | 4,996 | 0 | 0.74 |
| **ESM2-650M zero-shot** | **+0.732** | 4,996 | 0 | **0.90** (546 TP / 62 FP) |

**+0.385 absolute lift — more than double the deterministic baseline.** (0.73 is high even for ESM because
TEM-1 β-lactamase is a favorable, deep, well-studied assay; the published ProteinGym *median* is ~0.49.)

**Lift 2 — genome-level edit (DONE).** `dna_decode/forward/genome_edit.py::predict_genome_edit(cds, nt_pos,
ref, alt)` lifts the INPUT from a protein mutation (`M69L`) to a real nucleotide edit in a CDS → codon →
AA change → the Regime-B predictor, classified **silent / missense / nonsense**. Standard genetic code;
REF-base coordinate gate + a DOUBLE check (translated WT AA vs the protein reference) both fail loudly.
Validated by the genetic-code unit tests (real translation logic, not mocks).

## Honest scope (the project's hard-won rails, applied)

1. **Two methods, both DMS-validated.** BLOSUM62 (deterministic, no GPU/network — strong baseline +0.347)
   and ESM2-650M (learned, +0.732). `predict_effect(..., method=)` selects; AlphaMissense would be a further
   drop-in. The learned model earns a large real lift here — the fitness-aligned molecular regime is exactly
   where it should (contrast: on the antagonistic clinical-resistance direction ESM fails, per the
   resistance-conservativeness finding — that direction stays with the Regime-A catalogue).
2. **Regime B only.** This predicts MOLECULAR fitness (enzyme activity/stability). It is NOT an
   organism-level polygenic predictor — `regime="C_organismal"` → `abstain=True` by construction (the
   closed-negative regime; embeddings 0-for-5 de-confounded).
3. **NOT the clinical-resistance direction.** For "does this edit confer clinical resistance" (antagonistic
   selection), raw likelihood/severity scorers FAIL (the resistance-conservativeness finding) — use the
   Regime-A determinant catalogue (`amr_rules` / `tb_amr`) there. The regime tag picks the tool per edit.
4. **The validated quantity is the rank correlation**, reported with every batch; the per-variant
   preserved/damaging tier is a coarse read of the continuous score (illustrative, not the headline).

## What this establishes

The forward "edit → phenotype" capability is now REAL on E. coli, not a moonshot — built by assembling
existing parts (the cached DMS + BLOSUM baseline + a coordinate-checked forward wrapper) in one session, with
a reproducible honest number and free per-variant labels. Two of the three G2P regimes now have a shipped
forward path: Regime A (determinant catalogue, run forward) + Regime B (this cell). Regime C stays an honest
abstention.

**Next lifts (updated 2026-07-14):** ~~ESM2-650M~~ **DONE (+0.732)**; ~~genome-level path~~ **DONE**
(`predict_genome_edit`). Remaining: extend the E. coli panel (CcdB, IF1, folA already cached); wire a REAL
blaTEM nucleotide CDS end-to-end through `predict_genome_edit` (the codon engine is unit-test-validated; a
real-CDS demo is the honest next real-surface touch); optional AlphaMissense method. Frozen decoder surface (`amr_rules` / `calibrated_amr_rules` / `mic_tiers` /
`shipped_decoder_surface` / `cohort_manifest`) byte-unchanged (`verify_lock OK`); `dna_decode/forward` is a
NEW non-frozen package.
