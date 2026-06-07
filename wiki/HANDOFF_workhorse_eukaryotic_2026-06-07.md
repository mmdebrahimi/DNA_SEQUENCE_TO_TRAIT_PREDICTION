# Workhorse Handoff — Eukaryotic Cycle, Path B (Arabidopsis embedding test)

> For the **Precision 7780 (RTX 3500 Ada)** — the GPU machine. Sync = **git pull** (NOT Gmail/zip).
> Your scope is **ONLY Path B (Phase 2)**: the GPU-bound Arabidopsis flowering-time embedding test. The
> laptop owns ALL of Path A (fungal AMR) + all data-prep — do NOT touch fungal AMR or duplicate it.
> Full plan: `plans/Eukaryotic_DualMachine_Coordination.md` + `plans/EP8_Arabidopsis_Embedding_Test.md`.

## DO NOT (avoid duplication)
- Do NOT work on fungal AMR / C. auris / `dna_decode/data/fungal_amr.py` / `scripts/fungal_erg11_caller.py` — that's the laptop's (Path A).
- Do NOT edit Path-A ledger rows or fungal artifacts. Append only Path-B rows + Path-B result packets.
- Do NOT start until **Gate G1 is shared** (laptop posts the fungal-AMR result packet) AND the user confirms compute — see "Preconditions".

## Preconditions (check before any GPU work)
1. `git pull` — confirm you have `plans/EP8_Arabidopsis_Embedding_Test.md` + this handoff + the laptop's
   pre-staged AraGWAS manifest + baseline spec (laptop produces these in Phase 1 pre-stage).
2. **VRAM-fit gate (first task):** can a plant DNA-FM (PlantCaduceus first; 512bp windows) run inference on
   the RTX 3500 Ada (~12GB)? If YES → proceed. If NO → STOP + report; the user decides cloud budget (money gate).
3. User has explicitly greenlit Phase 2 compute (the eukaryotic cycle's money/hardware gate).

## Your tasks (Path B, Phase 2 — in order)
1. **VRAM-fit check** (above). Report pass/fail + the model that fits.
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
