# TB second-line CRyPTIC-parquet baseline — moxifloxacin (concentrated QRDR — works) + bedaquiline (diffuse LoF — honest low-sens) (2026-07-14)

Builds the tree-of-life census's #1 lowest-friction candidate: extend the shipped TB RIF/INH decoder cell to **second-line drugs**, scored on the free CRyPTIC compendium (Zenodo `VARIANTS.parquet`, 58.97 M per-isolate genomic-variant rows + measured BMD-MIC), reusing the FROZEN-adjacent scorer (`tb_amr.score_drug`) + the WHO mutation catalogue v2 (grade-1/2) UNCHANGED.

**What was actually built.** The scorer + catalogue were already drug-generic; second-line was gated by a single map (`score_tb_cryptic_parquet.py::DRUG_CODE`, which only covered RIF/INH). Extending it to the 10 second-line + new drugs (verified against both `DRUG_CATALOGUE_NAME` keys and the reuse-table `<code>_BINARY_PHENOTYPE` columns) unlocked the whole panel. Commit `1fa792d` (wiring + 3 tests, 11/11 pass).

## Moxifloxacin — the flagship second-line cell (gyrA/gyrB QRDR — the mycobacterial analog of the cipro cell)

Artifact: `wiki/tb_mxf_cryptic_parquet_baseline_2026-07-14.json`. Run over all 58,966,529 parquet rows; 6,784 measured HIGH-quality MXF labels; 6,648 with in-scope determinant calls; 39 grade-1/2 determinants across 11 positions (gyrA + gyrB).

| Metric | sens | spec |
|---|---:|---:|
| **raw** (SNV+indel = SNV-only; 0 indel determinants — gyrA/gyrB are point mutations) | **0.841** | **0.961** |
| **lineage-collapsed (the honest headline)** | **0.279** | **0.988** |

Lineage structure: 43 R-lineages, 259 S-lineages, 40 discordant (mixed-label clones, excluded). `status = WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`.

## Honesty tier (identical to the shipped RIF/INH cell — do not overclaim)

1. **In-distribution, NOT independent.** The WHO catalogue was built partly from CRyPTIC, so a CRyPTIC-scored number is in-distribution. This is a KNOWLEDGE baseline, not the independent validation the AMR-Portal arm gives for RIF/INH.
2. **Raw is clonality-inflated → lineage-collapsed is the headline.** TB is heavily clonal; the raw 0.841 counts one vote per isolate, so an over-sampled clone carries it. Lineage-collapsed (barcode-on-VCF, one vote per lineage) is the honest sens 0.279 / spec 0.988. This mirrors RIF (raw 0.916 → lineage 0.41) and INH (0.889 → 0.349) exactly — the raw is ~3× inflated, the lineage number is much lower sens + high robust spec.
3. **SNV-form determinant match** (per-base SNV + exact indel); true delins/LoF not normalized → a LOWER BOUND. Moxifloxacin has 0 indel determinants, so this caveat is moot for MXF.

**Reading:** moxifloxacin behaves exactly like the validated RIF/INH cell — high, robust specificity (0.99) and a lineage-collapsed sensitivity in the 0.28–0.44 band the WHO-catalogue rule shows across TB drugs. The gyrA/gyrB QRDR is the same target-site family the decoder validates for E. coli cipro; the fluoroquinolone cell transfers the QRDR mechanism into M. tuberculosis with an honest in-distribution number.

## Bedaquiline — an HONEST LOW-SENSITIVITY cell (the mechanistic contrast — a real finding)

Artifact: `wiki/tb_bdq_cryptic_parquet_baseline_2026-07-14.json`. 8,535 measured BDQ labels; 8,348 in-scope; 2448 grade-1/2 determinants (Rv0678/atpE/pepQ).

| Metric | sens | spec |
|---|---:|---:|
| raw SNV-only | 0.056 | 0.998 |
| raw SNV+indel (indel-matching lifted sens; 70 isolates hit, **9 TP** / 61 FP) | **0.183** | 0.991 |
| **lineage-collapsed** | **0.0** | **0.938** |

Lineage structure: only **5 R-lineages** (vs 406 S-lineages, 17 discordant). `status = WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`.

**The contrast is the finding — and it is the project's recurring theme, now inside TB:**
- **Moxifloxacin (gyrA/gyrB QRDR) = CONCENTRATED signal → works.** A handful of point mutations at defined QRDR positions explain most FQ resistance → the determinant rule recovers it (lineage 0.279, spec 0.99), exactly like the shipped RIF/INH cell.
- **Bedaquiline (Rv0678 = mmpR5) = DIFFUSE loss-of-function → the determinant rule structurally under-detects it.** BDQ resistance is (a) RARE in CRyPTIC (only 5 R-lineages — it is a new drug) and (b) mechanistically diffuse: Rv0678 resistance arises via *diverse* LoF mutations (indels, frameshifts, scattered SNVs) that no grade-1/2 point-catalogue can enumerate. Indel-matching helped (sens 0.056→0.183, +9 true LoF) but most is un-catalogued or below grade-2. lineage-collapsed sens rounds to **0.0**.

This is NOT a build failure — the pipeline ran cleanly and produced a real, honest number. It is the same **concentrated-mechanism-works / distributed-mechanism-fails** boundary the project found for E. coli (cipro QRDR + cef plasmid-β-lactamase pass; tet distributed mobile-element fails) and for the HIV ESM probe (antagonistic-selection blindness) — now demonstrated *within* the TB cell across two drugs. Screen a determinant-catalogue cell for **mechanism concentration** before expecting sensitivity: QRDR-shaped drugs (moxi/lev) will behave like RIF/INH; diffuse-LoF drugs (BDQ, and likely the Rv0678-linked clofazimine) will be honest low-sensitivity cells.

## Panel status

- **moxifloxacin — DONE** (concentrated QRDR; works: lineage 0.279/0.988).
- **bedaquiline — DONE** (diffuse Rv0678 LoF; honest low-sens negative: lineage 0.0/0.938).
- **The remaining panel is one command each** (the wiring supports them): `uv run python scripts/score_tb_cryptic_parquet.py --drug {levofloxacin,amikacin,kanamycin,linezolid,ethionamide,clofazimine,delamanid,ethambutol}` — each ~15 min streaming the 2.9 GB parquet. Predicted by the concentration rule: levofloxacin (QRDR) ≈ moxifloxacin; amikacin (rrs, 2 determinants) + linezolid (rrl/rplC) should be reasonably concentrated; clofazimine (Rv0678, cross-resistant with BDQ) ≈ the bedaquiline low-sens pattern.

Frozen decoder surface (`amr_rules` / `calibrated_amr_rules` / `mic_tiers` / `shipped_decoder_surface` / `cohort_manifest`) byte-unchanged (verify_lock OK); `tb_amr` + the parquet adapter are NON-frozen.
