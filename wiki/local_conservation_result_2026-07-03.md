# Frozen-spec local conservation score — the project computes its OWN Site-Independent (purist path)

**Date:** 2026-07-03
**Script:** `scripts/local_conservation_score.py` (`wiki/local_conservation_scores.json`)
**Data (free):** ProteinGym per-protein MSAs (`DMS_msa_files.zip`, 1.4 GB) + its redundancy-based sequence
weights (`DMS_msa_weights.zip`, 43 MB) + per-assay DMS/zero-shot CSVs (Zenodo 15293562). Human Activity
(function) assays. On D:.

## Why (the purist path)

Row 313 (`pg_native_threeway`) still *cited* ProteinGym's published Site-Independent per-mutant score. This
turns that into **"the project computes its own"**: a from-scratch independent-sites conservation score with a
PINNED, hashed spec — reproducing the baseline rather than citing it. Concrete target (from the three-way):
~0.43 median Spearman on function.

## Frozen spec (declared before outcomes; MSA + weights sha256'd per assay)

- **MSA:** ProteinGym a2m. Match columns = per-sequence residues after stripping lowercase + `.` (a2m insert
  states). Protein-position → match-column map built by **walking the focus (first) sequence** (uppercase =
  residue + column; lowercase = residue only; `-` = column only) — NOT `pos − MSA_start` (which is wrong when
  the focus carries inserts/trims; see the bug below).
- **Weights:** ProteinGym's own precomputed per-sequence weights `<UniProt>_theta_<θ>.npy` (no additional
  redundancy filter — avoids double-correction). **Unweighted fallback** (uniform weights) for assays where
  ProteinGym publishes no weights, flagged separately.
- **Frequency:** `f_j(a) = (Σ w·[res==a] + λ) / (Σ w·[res is std AA] + 20λ)` over the 20 AAs, gaps excluded,
  **λ = 0.5**.
- **Score:** `log₂ f_col(mut) − log₂ f_col(wt)` (delta-log-odds; cancels family composition). Positions outside
  the MSA coverage **abstain**.

## Result — reproduces the number AND hits the target

| | n | local conservation (median Spearman) | ProteinGym Site-Independent | reproduction Δ |
|---|---|---|---|---|
| **weighted (ProteinGym weights)** | 7 | **0.451** | 0.459 | **0.0145** |
| unweighted fallback | 12 | 0.480 | — | — |
| **overall function** | 19 | **0.476** | — | — |

- **The implementation is validated:** on the 7 assays where ProteinGym publishes weights, this project's
  independent score **reproduces ProteinGym's Site-Independent to within 0.0145 Spearman per assay** — i.e. the
  same baseline, computed here from the MSA + weights + a pinned formula. The number is now *owned*, not cited.
- **The target is met:** function median **0.476** (19 assays) ≥ the ~0.43 target from the three-way. Even the
  unweighted variant reaches 0.48 on function.
- This closes the review's "purist path": deterministic conservation as a **project-computed**, frozen,
  reproducible score — not a citation.

## Honest scope

- **Weighted coverage is limited by ProteinGym's published weights** (7 of 19 Activity assays have them). For the
  rest, the unweighted fallback is a valid frozen variant (redundancy-uncorrected — it can over-weight deep
  clades — but here it reaches a comparable ~0.48 median; flagged `weighted=false` per assay in the JSON).
- **verify-in-batch caught two real bugs:** (1) the naive `pos − MSA_start` position map gave 7% wt-vs-consensus
  agreement (ncol 294 ≠ protein-length 413 for ADRB2) → fixed to a focus-sequence coordinate walk → 0.047 →
  0.451; (2) `run()` skipped weight-missing assays before the fallback could fire → fixed.
- Assays whose MSA covers only a narrow sub-region (e.g. KCNH2 535–565) or lack an a2m abstain (reported, not
  imputed). Single substitutions only.
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9).

## Reproduce

```bash
uv run python scripts/local_conservation_score.py --selections Activity
uv run pytest tests/test_local_conservation_score.py -q   # 2 offline synthetic tests
# MSAs (DMS_msa_files.zip) + weights (DMS_msa_weights.zip), Zenodo 15293562, extracted on D:
```
