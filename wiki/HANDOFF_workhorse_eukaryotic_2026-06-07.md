# Workhorse Handoff — Eukaryotic Cycle, Path B (Arabidopsis embedding test)

> For the **Precision 7780 (RTX 3500 Ada)** — the GPU machine. Sync = **git pull** (NOT Gmail/zip).
> Your scope is **ONLY Path B (Phase 2)**: the GPU-bound Arabidopsis flowering-time embedding test. The
> laptop owns ALL of Path A (fungal AMR) + all data-prep — do NOT touch fungal AMR or duplicate it.
> Full plan: `plans/Eukaryotic_DualMachine_Coordination.md` + `plans/EP8_Arabidopsis_Embedding_Test.md`.

## DO NOT (avoid duplication)
- Do NOT work on fungal AMR / C. auris / `dna_decode/data/fungal_amr.py` / `scripts/fungal_erg11_caller.py` — that's the laptop's (Path A).
- Do NOT edit Path-A ledger rows or fungal artifacts. Append only Path-B rows + Path-B result packets.
- Do NOT start until **Gate G1 is shared** (laptop posts the fungal-AMR result packet) AND the user confirms compute — see "Preconditions".

## ▶ GO — G1 SHARED + Path B PRE-STAGED (updated 2026-06-08, laptop side)
- ✅ **Gate G1 (Path A) is DONE + shared** — C. auris fluconazole decoder validated (method transfers,
  sens 1.0; LABEL_LIMITED_FAILURE acc 0.79 — documented label-limitation, not a defect). The phase gate
  that held Path B is now RELEASED. See `wiki/fungal_ep7_g1_closeout_2026-06-08.md`. Path B is GO.
- ✅ **Path B is PRE-STAGED — execute, don't hunt.** Full executable spec (every URL pinned + verified,
  baselines + CV + G2 PASS/FAIL frozen): **`plans/EP8_PathB_PreStage_Manifest.md`** (REVISED 2026-06-08
  post-brainstorm — read the header + §8 open decisions first). Workhorse order:
  1. `uv run python scripts/fetch_arabidopsis_pathb.py` (confirms committed phenotype labels + HEAD-checks
     genotype URLs), then `--download` to pull the 19.2 GB genotype VCF to D:.
  2. **Run the CPU-only §0.5 G2 dry-manifest FIRST — it gates all GPU work.** Do NOT embed until every
     dry-manifest check is green (accession join, window/coord table, N-fraction QC, matched variant matrix,
     group labels). Decide the §8 open design choices (estimand / window-selection rule / GPU budget) at this
     point.
  3. Only then run the §1a phenotype-agnostic embedding → baselines → leave-one-group-out CV → G2 verdict.
  - NOTE the framing fix: the gate is the **phenotype-agnostic** embedding (not the curated FLC/FRI/FT panel —
    that's now a secondary diagnostic only). PASS needs a paired-bootstrap CI excluding 0, not a point margin.
- ✅ Phenotype labels already COMMITTED at `data/arabidopsis/` (FT10 n=1162, FT16 n=1122; 1122 with both).
- Pre-commit (ratified): a clean **G2-FAIL does NOT auto-close** the embedding frontier (KEEP-OPEN).

## Status (updated 2026-06-07 — laptop side)
- ✅ **Workhorse identity CONFIRMED:** personal **Precision 7780 (RTX 3500 Ada ~12GB)**, NOT Bombardier/DLP.
  Path B safety gate cleared; personal code may run there.
- ✅ **VRAM-fit gate PRE-RESOLVED by the laptop (skip task 1):** PlantCAD inference fits 12GB decisively.
  4 sizes (l20=20M / l24=40M / l28=112M / l32=225M); authors target the RTX 3090; 225M weights ~0.5GB fp16
  @ 512bp ctx. The "24-80GB" figure was TRAINING cost (8×H100), NOT inference. **Use `PlantCaduceus_l32`
  (225M, richest) — it fits the RTX 3500 Ada 12GB for frozen-embedding extraction with no quantization;**
  drop to l28/l24 only if throughput is the bottleneck. HF: `kuleshov-group/PlantCaduceus_l32`.
- ✅ **Compute greenlit:** local 12GB GPU approved; paid cloud A100 DEFERRED (Databricks high-GPU exists but
  is money-gated — do NOT fire it without explicit user OK). Since 12GB suffices, no cloud is needed for v0.

## Preconditions (check before any GPU work)
1. `git pull` — confirm you have `plans/EP8_Arabidopsis_Embedding_Test.md` + this handoff + the laptop's
   pre-staged AraGWAS manifest + baseline spec (laptop produces these in Phase 1 pre-stage).
2. ✅ VRAM-fit — DONE (see Status above); start at task 2.
3. ✅ Compute greenlit — local 12GB only (see Status above).
4. **Still required:** Gate G1 shared (laptop posts the fungal-AMR result packet) before you START — per the
   share-and-decide phase-gate design. You may `git pull` + read now, but hold GPU work until G1 lands.

## Your tasks (Path B, Phase 2 — in order)
1. ✅ **VRAM-fit check — DONE by laptop** (PlantCaduceus_l32 225M fits 12GB; see Status). Skip.
2. **Download** AraGWAS genotype matrix (10.7M SNPs / 2029 acc) + AraPheno flowering-time phenotype
   (1,003 acc @ 10°C). Public, no access gate. (Laptop pre-stages the exact URLs/manifest.)
3. **Baselines (no GPU):** SNP-PRS + kinship/population-structure-only regression on flowering time
   (clade-stratified CV). These are the bars the embedding must beat.
4. **Embedding (GPU):** embed accessions with the fitting plant DNA-FM → regression head → same
   clade-stratified CV folds.
5. **Within-lineage diagnostic:** confirm the embedding predicts mechanism, not population structure
   (reuse the logic of `scripts/within_lineage_diagnostic.py` — the trap that sank cipro NT).
6. **Result packet:** `wiki/phase2_arabidopsis_result_<date>.md` — embedding R² vs both baselines +
   within-lineage verdict = the embedding-niche **PASS/FAIL** (Gate G2). Push to origin.

## Verdict contract (Gate G2)
- **PASS:** embedding R² beats BOTH baselines by a meaningful margin AND within-lineage > population-structure
  → first evidence the frozen-FM thesis earns its keep (after 0-for-3 on AMR/pathotype/carbon-util).
- **FAIL:** embedding ≤ baselines OR within-lineage = structure → thesis 0-for-4 on a substrate that met
  every precondition → close the embedding frontier permanently. Document honestly either way.

## Handoff back
Push the result packet + a one-line ledger row (`project-state-row`, Path-B only) → notify the laptop side
via the next session. The laptop + you then DECIDE at Gate G2 whether to iterate (different FM / cloud A100)
or conclude.
