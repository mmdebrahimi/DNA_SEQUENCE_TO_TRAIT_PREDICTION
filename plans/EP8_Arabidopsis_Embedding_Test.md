# EP-8 — Arabidopsis flowering-time embedding test (QUEUED — gated on compute)

> **Status:** QUEUED 2026-06-07 (Path B of the ratified eukaryotic cycle). **BLOCKED on a money/hardware
> decision** — needs a ≥24GB GPU or a cloud budget. Fires the moment the user confirms compute.
> **Anchors on:** `research_outputs/eukaryotic-multimodal-substrate-feasibility-2026-06-07.md`.
> **Why it matters:** this is the **first true YES/YES/YES embedding-niche substrate** the project has
> found — sampling-independent quantitative label (flowering time) + NO curated mechanism catalog +
> organism-specific depth (~1,003 accessions). The frozen-foundation-model thesis went 0-for-3 on AMR /
> pathotype / carbon-util; Arabidopsis flowering-time is its **fairest possible shot**. A clean PASS would
> be the first evidence embeddings earn their keep; a clean FAIL would make the thesis 0-for-4 on a
> substrate that met every precondition — decisive either way.

## The compute gate (the blocker — user's decision)
- Plant DNA foundation model (PlantCaduceus / AgroNT / NT-multispecies-v2) needs a **≥24GB GPU** (the
  Caduceus paper used RTX 3090 / A100). The project's GTX 860M (4GB) **cannot** run it.
- Precision 7780's RTX 3500 Ada (~12GB) is borderline — PlantCaduceus's 512bp windows MIGHT fit; needs a
  VRAM check. Else paid cloud (A100 spot ~$1-2/hr) = **MONEY GATE → user must approve a budget.**
- **Soraya will not provision paid compute autonomously.** This EP stays parked until the user names a
  ≥24GB machine OR approves a cloud budget.

## Terminal claim (when fired)
On the AraGWAS flowering-time cohort (1,003 accessions, 10°C; 10.7M SNPs), a plant-DNA-FM embedding +
regression head predicts flowering time with R² beating (a) a SNP-PRS baseline and (b) a population-
structure/kinship-only baseline by a meaningful margin under clade-stratified CV — the embedding-niche
PASS bar. Within-lineage diagnostic (reuse `within_lineage_diagnostic.py` logic) confirms mechanism, not
population structure (the same trap that sank cipro NT).

## Substrate (verified)
- Data: AraPheno/AraGWAS — flowering time 1,003 (10°C) / 971 (16°C) accessions, public; 10.7M-SNP matrix
  over 2,029 accessions. Fully downloadable (no access gate, unlike UKB/AoU).
- Label: quantitative phenotype (growth-chamber measurement) — sampling-independent. ✅
- No curated mechanism catalog (polygenic) — the embedding niche. ✅

## Build steps (fire on compute confirmation)
1. Confirm GPU (≥24GB) or cloud budget — **the gate**.
2. Download AraGWAS genotype matrix + flowering-time phenotype (AraPheno).
3. Embed accessions with a plant DNA-FM (PlantCaduceus first) on the GPU; OR a SNP-PRS / GBLUP baseline.
4. Clade-stratified CV (kinship-aware folds) + within-lineage diagnostic + the two baselines (PRS, kinship-only).
5. Verdict: embedding R² vs baselines → the embedding-niche PASS/FAIL. Capstone artifact either way.

## Falsifier
Embedding does NOT beat the SNP-PRS / kinship baseline, OR within-lineage = population-structure → the
embedding thesis is 0-for-4 on a substrate that met every precondition; conclude embeddings don't earn
their keep for solo-scale G2P and close that frontier permanently.
