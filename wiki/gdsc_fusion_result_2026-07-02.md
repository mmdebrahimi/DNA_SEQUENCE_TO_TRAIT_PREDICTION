# GDSC fusion decode — the BCR-ABL case powered, the ALK case sharpened

**Date:** 2026-07-02
**Script:** `scripts/gdsc_fusion_decoder.py` (`wiki/gdsc_fusion_scores.json`)
**Data (free, on D:):** GDSC1 + GDSC2 fitted dose-response (release 8.5, `cog.sanger.ac.uk`, 333k + 242k rows) +
DepMap `sample_info.csv` (figshare 20274744, the COSMIC↔DepMap_ID↔lineage bridge) + the same
`CCLE_fusions.csv` the DepMap decoder uses (figshare 20967591).

## Why

The fusion result in `wiki/depmap_fusion_methylation_result_2026-07-02.md` was **directional-only**: the PRISM
19Q4 subset (569 lines) has **zero BCR-ABL1 cell lines** (CML/heme lines absent) and only ~9 ALK-fusion lines.
GDSC screens ~950 lines and **includes the CML lines**, so this run turns the untestable BCR-ABL case into real
numbers and gives the ALK case power.

Join chain: GDSC (keyed by COSMIC_ID) → `sample_info` (COSMIC_ID → DepMap_ID + lineage) → `CCLE_fusions`
(keyed by DepMap_ID). Response metric = **LN_IC50** (lower = more sensitive). De-confounding reuses the promoted
`dna_decode.deconfound.group_centered_biomarker_t` (binary fusion → within-lineage t): does the fusion separate
sensitivity **within a lineage**, not just because BCR-ABL concentrates in the drug-sensitive blood lineage.

## Result — all 5 cases POWERED and de-confounded-significant

| Drug | Fusion | Dataset | n+ | global ΔLN_IC50 | Mann-Whitney p | within-lineage t | within-lineage Δ |
|---|---|---|---|---|---|---|---|
| **imatinib** | BCR-ABL1 | GDSC1 | 9 | −4.02 | ~0 | **−6.5** | −3.53 (blood) |
| **dasatinib** | BCR-ABL1 | GDSC2 | 10 | −8.34 | ~0 | **−14.95** | −7.61 (blood) |
| **nilotinib** | BCR-ABL1 | GDSC2 | 9 | −7.58 | ~0 | **−9.09** | −5.88 (blood) |
| **ponatinib** | BCR-ABL1 | GDSC1 | 8 | −6.14 | ~0 | **−4.83** | −4.53 (blood) |
| **crizotinib** | ALK | GDSC2 | 15 | −2.19 | 2e−05 | **−4.09** | −3.25 (lymphocyte) |

(ΔLN_IC50 = fusion+ mean − fusion− mean; negative = fusion+ more sensitive. within-lineage Δ = the same
difference computed inside the single lineage that carries the fusion.)

### What changed vs the PRISM-only result

- **BCR-ABL → ABL-TKIs went from UNTESTABLE (n=0) to the strongest fusion signal in the project.** Four ABL
  inhibitors (imatinib / dasatinib / nilotinib / ponatinib), each with n=8–10 BCR-ABL+ lines, each with a
  4–8 log-unit global sensitization and a within-lineage t of **−4.8 to −15**. The signal is not the lineage
  confound: comparing BCR-ABL+ (CML) against BCR-ABL− (AML/ALL) lines **within the blood lineage**, the BCR-ABL+
  lines are still 3.5–7.6 LN_IC50 units more sensitive.
- **ALK → crizotinib sharpened from directional to powered:** n=9 → 15, within-lineage t −1.44 → **−4.09**
  (de-confounded within the lymphocyte/ALCL lineage, Δ −3.25). The directional PRISM result is confirmed.

## The feature-match law, now with a powered fusion cell

A gene fusion is invisible to point-mutation, copy-number, and single-gene expression (shown null for ALK in the
PRISM cross-check). GDSC now makes the *positive* half of that claim strong: where the mechanism **is** a fusion
(BCR-ABL, ALK), fusion presence is a decisive, de-confounded predictor of the matching TKI's response. The
"feature type must match mechanism type" law holds across all five feature types, and the fusion cell is no
longer the weak one.

## Honest scope

- **GDSC LN_IC50 is not comparable across GDSC1 and GDSC2** (different concentration ranges) — each case names
  its dataset; no cross-assay pooling. imatinib/ponatinib are GDSC1-only; dasatinib/nilotinib/crizotinib GDSC2.
- **In-distribution association, not a held-out external validation.** This proves which feature type carries a
  mechanism on a broad public screen; it is not a deployed predictor.
- **n+ is still modest in absolute terms** (8–10 BCR-ABL lines exist in CCLE at all), but the effect is so large
  and the Mann-Whitney + within-lineage t so significant that the case is powered, not directional.
- Reuses the frozen `dna_decode.deconfound` primitives — no new statistics. Frozen bacterial/viral/fungal AMR
  surface untouched (leak guard 9/9).

## Reproduce

```bash
# (one-time) GDSC1/2 xlsx -> parquet cache already on D:; sample_info + CCLE_fusions fetched.
uv run python scripts/gdsc_fusion_decoder.py
uv run pytest tests/test_gdsc_fusion_decoder.py -q   # 5 offline synthetic tests
```
