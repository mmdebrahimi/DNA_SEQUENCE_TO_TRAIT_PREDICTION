<!-- memo-schema: 0.4 -->
# ΔΔG → drug-resistance fold-change: predictor accuracy + calibration (supported memo, 2026-07-09)

> `/research` run (executed by hand under Soraya, L1). Slug: `ddg-to-fold-change-calibration-drug-resistance-2026-07-09`. Intake-validated. Search partially policy-blocked on the antimicrobial-resistance+MIC angle; the thermodynamic identity + oncology-kinase analog carried the result. Rejected/verify-needed rows: `..._unsupported.md`.

## Research Context (problem anchor)

Verbatim topic: *"How well do computed protein ΔΔG values (FoldX/Rosetta/AlphaFold-based + stability ΔΔG) correlate with and predict measured antimicrobial/antiviral drug-resistance fold-change (IC50/MIC), and what calibration methods map computed ΔΔG → measured fold-change?"*

This tests **`/hypothesise` #2** for the dna_decode genome world model: *manufacture physics labels (ΔΔG) to move drug resistance out of the antagonistic-selection regime (where learned likelihood scorers provably lose) into the molecular-property regime (where learning wins ~2.1×), sidestepping the LABELS wall.*

## Supported claims

**HIGH confidence (physical law + WebFetched primary):**

1. **The ΔΔG→fold-change map is EXACT, not a fit — no calibration is needed for competitive binding.** Fold-change in Ki/Kd = `exp(ΔΔG / RT)`; `RT·ln(10) = 1.373 kcal/mol` per 10× at 300 K (1.42 at 310 K). Self-computed and cross-checked against a fetched primary: *"within 1.4 kcal/mol of experiment (i.e., within a 10-fold error in Kd change at 300 K)"* — [flex_ddg ligand paper, PMC6311686](https://pmc.ncbi.nlm.nih.gov/articles/PMC6311686/). The Abl-kinase study independently uses *"a ten-fold reduction ... (1.36 kcal mol⁻¹)."*

2. **Therefore the bottleneck is COMPUTING ΔΔG, and the error exponentiates.** A predictor's kcal/mol RMSE maps to a *multiplicative* fold-change error: **1.0 → 5.4×, 1.46 → 11.6×, 1.9 → 24×** (300 K). This is the single most important number for the project.

3. **Rosetta flex_ddg predicting a mutation's effect on LIGAND binding: RMSE 0.99 kcal/mol (Pearson 0.39, n=86) / RMSE 1.46 (Pearson 0.25, n=134)** — [PMC6311686](https://pmc.ncbi.nlm.nih.gov/articles/PMC6311686/), WebFetched. Note the *low Pearson (0.25–0.39)*: even at a "good" RMSE, rank-ordering mutations is weak.

**MEDIUM confidence (real primary, abstract-via-search — full text paywalled/redirected):**

4. **FEP+ (expensive alchemical) classified 144 clinical Abl-kinase point mutations as resistant vs susceptible at 88% accuracy**, using the 10-fold (1.36 kcal/mol) threshold — [Hauser et al., Comms Biol 2018](https://www.nature.com/articles/s42003-018-0075-x). This is the **direct proof-of-concept** that ΔΔG→resistance calibration works — for a *target-site competitive* kinase inhibitor.

5. **Alchemical FEP kinase:inhibitor MUE ≈ 1.0 kcal/mol** (CDK2/JNK1/p38/Tyk2) — the accuracy ceiling of the expensive method; still ~5× fold-change uncertainty per claim 2.

6. **Best cheap deep-learning stability surrogate RaSP: PCC 0.62, MAE 0.94 kcal/mol** (zero-shot) — [eLife 82593](https://elifesciences.org/articles/82593). **flex ddG protein–protein interface: R 0.58, reliable >1 kcal/mol hotspot ID** — [Barlow et al., JPCB 2018](https://pubs.acs.org/doi/10.1021/acs.jpcb.7b11367). **Cheap ddG from homology/predicted structure collapses toward a null model below ~45% template identity** — [CSBJ 2022](https://www.csbj.org/article/S2001-0370(22)00542-6/pdf).

## Bottom line for the hypothesis

- **The calibration question is basically solved by physics** (fold = exp(ΔΔG/RT)); the real question is *predictor accuracy*, and it is **regime-split**:
  - **Cheap ddG (FoldX/Rosetta/RaSP, RMSE ~1–1.9 kcal/mol)** → good enough to *classify* resistance hotspots (>1 kcal/mol), **not** to predict quantitative MIC fold-change (5–24× error).
  - **FEP-grade (MUE ~1.0 kcal/mol)** → ~88% resistant/susceptible classification on clinical mutations, but expensive and demonstrated on *competitive kinase binding*.
- **Scope wall (the load-bearing caveat):** a binding-ΔΔG model is blind to the *dominant* non-binding resistance mechanisms — **efflux, β-lactamase hydrolysis, target amplification, promoter/regulatory changes, kcat/Km shifts**. It fits **target-site-competitive** cells (β-lactam↔PBP, protease/kinase/RT inhibitors) and NOT the project's efflux/regulatory cells.
- **Net:** #2 is *viable but narrowed* — it does NOT give a general MIC predictor, but it could give a **physics-labelled classifier for target-site-competitive resistance** that (a) needs no wet-lab labels, (b) sidesteps the likelihood-model anti-alignment, (c) is a *third* signal distinct from both the catalog and the closed embedding track. The falsifiable next step is cheap (FoldX/Rosetta cartesian_ddg on a handful of known target-site DRMs vs their measured fold-change; check if it clears the >1 kcal/mol hotspot bar).

## Decisions for Human Confirmation

| # | Decision / candidate use | Verification needed | Confidence |
|---|---|---|---|
| 1 | Scope hypothesis #2 to **target-site-competitive** resistance cells only (kinase/protease/RT/PBP), not efflux/regulatory. | Confirm which of the project's cells are competitive-binding (HIV PR/RT + INSTI = yes; efflux/β-lactamase = no). | high |
| 2 | Use **cheap ddG (Rosetta cartesian_ddg / RaSP) as a hotspot CLASSIFIER**, reserve FEP+ for a small validation set — do NOT expect quantitative MIC from cheap ddG. | Pilot: FoldX/Rosetta on ~10 known HIV RT/PR DRMs vs measured PhenoSense fold-change; does it clear the >1 kcal/mol bar + beat the catalog's blind spot? | high |
| 3 | Treat the ΔΔG→fold-change map as **exact physics (no learned calibration layer)** — model only the ΔΔG computation. | None — it's `exp(ΔΔG/RT)`; verified 1.373 kcal/mol/decade at 300 K. | high |
| 4 | The antimicrobial-MIC calibration literature is **policy-search-walled here** — retrieve it on an unfiltered surface before committing. | Re-run the ΔΔG→MIC antimicrobial search on a machine without the usage-policy filter; find the antimicrobial (not just oncology) analog. | medium |
