# DepMap fusion + methylation modalities — feature-match law extended to 5 feature types

**Date:** 2026-07-02
**Script:** `scripts/depmap_fusion_methylation.py` (`wiki/depmap_fusion_methylation_scores.json`)
**Data:** DepMap/PRISM 19Q4 (569-line subset) + CCLE_fusions.csv (figshare 20967591) + CCLE RRBS TSS-1kb
methylation (`data.broadinstitute.org/ccle/CCLE_RRBS_TSS_1kb_20180614.txt`, 137 MB). Gitignored on D:.
**De-confounding:** promoted `dna_decode.deconfound` primitives — binary fusion → within-lineage biomarker t;
continuous methylation → clade-centered (primary-tissue-residualized) Spearman.

## What this adds

The multimodal DepMap decoder proved **mutation / copy-number / expression** each carry the de-confounded
drug-response signal for the *matching* mechanism (BRAF-mutation→vemurafenib; EGFR-amplification/expression→
erlotinib). This run adds two feature types that are **invisible to those three**:

- **Gene fusion** (binary presence): a fusion (EML4-ALK) is not an SNV, not a copy gain, and does not raise
  single-gene expression of the partner — so only a fusion-presence feature can capture it.
- **Promoter methylation** (continuous RRBS beta): methylation *silences* a gene with no coding change — so
  only a methylation feature captures the silencing directly.

For each case the matched modality is scored de-confounded, and the "wrong" modalities (mutation / CN /
expression of the **same** biomarker gene) are scored alongside, to show the matched feature type wins.

## Result 1 — FUSION (ALK), directionally unanimous but UNDERPOWERED

PRISM 19Q4 is a 569-line subset. **BCR-ABL1 overlap = 0** (CML lines absent) → imatinib/dasatinib/nilotinib
are **not testable here** (named data wall, below). ALK-fusion overlap = 9 lines.

| Drug | Biomarker | Fusion within-lineage t | ALK mutation | ALK copy-number | ALK expression | Verdict |
|---|---|---|---|---|---|---|
| crizotinib | ALK fusion (n+=9) | **−1.44** | 0.00 | +0.01 | +0.02 | fusion-only, correct dir |
| alectinib | ALK fusion (n+=9) | **−1.25** | 0.00 | +0.02 | −0.01 | fusion-only, correct dir |
| lorlatinib | ALK fusion (n+=8) | **−1.42** | 0.00 | +0.01 | −0.02 | fusion-only, correct dir |

- **Qualitatively textbook:** across all 3 ALK-TKIs, fusion presence is the *only* feature type carrying a
  correct-direction (fusion+ → more sensitive → negative t) de-confounded signal. Mutation, copy-number, and
  expression of ALK are all ~null — an ALK fusion is genuinely invisible to them.
- **Honestly underpowered:** n+≈9 across all lineages, so |t|≈1.3–1.4 (not significant) and no single lineage
  has ≥4 fusion+ cells for a per-lineage delta. This is a **directional** result, not a powered one. It is
  reported as `underpowered: true` in the JSON and not overstated.

## Result 2 — METHYLATION (MGMT, SLFN11), well-powered

~464–474 bridged lines (CCLE cell-name → DepMap_ID via `prism_cells.ccle_name`).

| Drug | Biomarker | Methylation within-lineage ρ | mutation | copy-number | expression | Verdict |
|---|---|---|---|---|---|---|
| temozolomide | MGMT silencing (sensitize) | **−0.171** | 0.00 | +0.04 | +0.169 | MATCH |
| topotecan | SLFN11 silencing (resist) | **+0.103** | 0.00 | −0.00 | −0.221 | MATCH |
| olaparib | SLFN11 silencing (resist) | +0.034 | 0.00 | +0.08 | −0.048 | null (weak) |

- **MGMT → temozolomide** (canonical): high promoter methylation → silenced MGMT → cell can't repair TMZ
  damage → **sensitive** → negative ρ (−0.171), de-confounded by lineage. Mutation/CN are null (MGMT is not
  mutated or copy-altered — it is *epigenetically* silenced).
- **SLFN11 → topotecan**: high methylation → SLFN11 silenced → **resistant** to DNA-damage → positive ρ
  (+0.103), de-confounded. olaparib shows the same sign but weak (+0.034) — the PARP-inhibitor effect is
  genuinely small in this PRISM subset; reported as null, not tuned up.
- **Methylation agrees with expression (opposite sign, same silencing axis):** MGMT methylation −0.171 vs MGMT
  expression **+0.169**; SLFN11 methylation +0.103 vs expression **−0.221**. Methylation and its downstream
  expression readout are consistent — exactly what silencing biology predicts. Methylation is a *valid* feature
  type; it need not beat expression, but it captures the mechanism where mutation and copy-number cannot.

## The law, extended

**Feature type must match mechanism type** now holds across **five** feature types on the same DepMap substrate:

| Feature type | Mechanism it captures | Exemplar (this project) |
|---|---|---|
| point mutation | activating/LoF coding change | BRAF→vemurafenib, TP53→nutlin |
| copy number | amplification/deletion dosage | EGFR/ERBB2→lapatinib |
| expression | transcript-level dependency | EGFR→erlotinib |
| **gene fusion** | oncogenic fusion driver | **ALK→crizotinib/alectinib/lorlatinib** |
| **promoter methylation** | epigenetic silencing | **MGMT→TMZ, SLFN11→topotecan** |

A feature set that omits a mechanism's matching feature type will miss that mechanism no matter how large the
model — the negative result that motivated this whole modality sweep.

## Honest scope / walls

- **Fusion is underpowered in PRISM 19Q4** (ALK n≈9; BCR-ABL n=0). **External data wall** (code-closable): a
  fuller drug-response screen with the CML/heme lines (CTRP or GDSC, which include K562 etc.) would power the
  BCR-ABL→imatinib/dasatinib/nilotinib cases and the ALK cases. The fusion + methylation loaders already
  support this — only the response matrix needs swapping. No new code, just a broader-cohort fetch.
- **In-distribution, not independent.** This is a de-confounded association test on public data, not a held-out
  external validation. It demonstrates *which feature type carries a mechanism*, not a deployable predictor.
- Both new modalities reuse the frozen `dna_decode.deconfound` primitives — no new statistics were invented.
- Frozen bacterial/viral/fungal AMR surface untouched (this is a research-arm cell-line analysis).

## Reproduce

```bash
uv run python scripts/depmap_fusion_methylation.py --data D:/dna_decode_cache/depmap_pilot
uv run pytest tests/test_depmap_fusion_methylation.py -q   # offline synthetic pin (4 tests)
```
