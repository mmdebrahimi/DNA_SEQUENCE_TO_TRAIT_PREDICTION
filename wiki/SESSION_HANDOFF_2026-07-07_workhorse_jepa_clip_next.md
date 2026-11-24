# Workhorse handoff — J-track next steps (JEPA / CLIP "world-model sidecar") · 2026-07-07

**From:** DNA-11 (dna_decode laptop, GTX 860M — no viable ML GPU)
**To:** Workhorse (RTX 3500 Ada Laptop GPU)
**Branch:** `mosafer` (this is the J-track / workhorse-GPU scratch line, where the J2 ESM-2 bank lives)
**Non-dup checked:** yes — this ADDS the ranked next-step menu; it does not redo the banked J2 work.

---

## TL;DR

- **J2 protein-representation Phase 1 is DONE and banked here** (`wiki/j2_esm2_result_2026-07-05.md`): ESM2-650M zero-shot on ProteinGym DMS → **median |Spearman| 0.491** (shuffled 0.011), 40 assays, on your RTX 3500 Ada. Clears the 0.45 pass + 0.48 stretch bar. *The learned-representation thesis provably works at the molecular (protein-variant-effect) layer.*
- **Everything past that is GPU-gated and pre-registered to be FALSIFIED**, not confirmed. The honest prior is strong skepticism: learned-embedding → **organism phenotype** is a **0-for-5 de-confounded negative** across this project. The J-track's only honest home is **representation / imputation / molecular-effect**, NOT organism-trait prediction.
- **Recommended next GPU move: J3 falsifier first (cheap, decisive), then J4 CLIP head only if J3 passes.** Do NOT jump to J4/J5 scale — the critical path is `J3 → (pass) → J4`.
- **Money gate:** any *cloud* compute or paid API = **user approval required** (inviolable). Local RTX 3500 Ada runs are free → proceed. Training (fine-tune / real-encoder JEPA) is free on-GPU but gate the *decision to spend the time* on "does beating the bar power a decoder feature."

---

## Context — where the whole project is (so you have the frame)

The validated product is a **deterministic multi-kingdom genotype→trait decoder** (bacterial/viral/fungal AMR + human PGx + genome-map). The **full current decoder suite lives on `main`** (mirrored to `mosfaer`) — including **CYP2D6 now decoder-complete for short-read WGS** and HLA B*57:01. This `mosafer` branch is the **older J-track scratch line** and is intentionally behind on decoder code; work the J-track here, pull decoder context from `main`/`mosfaer` if needed.

The binding constraint everywhere is **labels, not models**. The 8-gate negative-results map (`wiki/negative_results_map_2026-06-13.md` on main) screens any new dataset before labor.

---

## What is DONE (do not redo)

| Piece | Status | Evidence |
|---|---|---|
| **J1** genome-JEPA prototype (EMA target-encoder + stop-grad + anti-collapse) | ✅ verified (synthetic, CPU) | `scripts/genome_jepa_prototype.py` + `tests/test_genome_jepa_prototype.py` (6 offline) |
| **J2-protein Phase 1** ESM2-650M zero-shot masked-marginals on ProteinGym DMS | ✅ **PASS 0.491** | `wiki/j2_esm2_result_2026-07-05.md`; runner `scripts/esm_zeroshot_dms.py` |
| DMS falsifier harness (proxy) | ✅ 0.417 (AlphaMissense) | `scripts/dms_learned_model_falsifier.py` (+5 offline tests) |

---

## Ranked next steps (GPU-gated)

### ★ #1 — J3 falsifier (do this FIRST; it gates everything above it)
**Question:** does a learned JEPA/ESM representation + a small head **beat BOTH** (a) frozen-FM embeddings **AND** (b) an LD/statistical baseline, on a **de-confounded** downstream — masked-SNP imputation OR variant-effect ranking?

- **Comparators (E19 gauntlet):** frozen-FM embedding (`nt.h5`/`dnabert2.h5`) + linear head · LD/statistical baseline (the imputer) · k-mer/CNN · shuffled-label + PCA-only + held-out-clade negative controls.
- **Split:** lineage/ancestry-blocked, clade-stratified CV. **Reuse `scripts/within_lineage_diagnostic.py`** to rule out population structure.
- **PASS bar (pre-registered — freeze BEFORE the run):** JEPA/ESM rep beats frozen-FM **and** LD on the **full metric surface** (not one cherry-picked metric) by a meaningful margin, AND the within-lineage diagnostic rules out structure.
- **FAIL → honest close:** the learned-representation frontier for solo-scale G2P is closed (consistent with 0-for-5); the deterministic scan + LD imputer stay the product. **This is a legitimate, bankable outcome — a clean FAIL is a win.**
- **Cost:** mostly laptop-feasible once J2 embeddings are cached; the embedding-generation pass is the GPU part.

### #2 — J4 CLIP head (the actual "world-model sidecar" — ONLY if J3 passes)
Contrastive (InfoNCE) alignment of the seq/protein representation ↔ a **molecular/expression phenotype modality**, in a shared space.
- **Gate:** (a) your GPU (free) + (b) **paired de-confounded labels** — a seq-rep ↔ molecular-phenotype pairing that survives lineage residualization. Sourcing (b) is the real blocker; propose candidate modalities in the reply (e.g. ProteinGym-adjacent function assays, expression panels) and we de-confound-screen them together.
- **Honest scope:** CLIP earns its keep on the **molecular** axis where J2 already works — NOT on organism traits.

### #3 — genome-JEPA J2 "real encoder" (parallel GPU fork)
Swap the J1 synthetic windows → **real sequence windows**, init from an open FM (Caduceus / NT), on a de-confounded organism slice (Arabidopsis 1001G / openSNP EUR / a bacterial panel). Feeds J3 with a genome (not protein) representation. Same pre-registered falsifier as #1.

---

## Reproduce the banked J2 (workhorse env — exact)

Your machine specifics that already worked (from the banked result):
- venv: `.venv` in the repo; runner `scripts/esm_zeroshot_dms.py`.
- **HuggingFace weight payloads 403 on your network** → the runner falls back to direct `fair-esm` checkpoints from `dl.fbaipublicfiles.com`; `fair-esm` is installed into the venv.
- Scratch on `I:\scratch_fea\b0652085_tmp\` (keep heavy caches OFF the repo drive).

```powershell
$env:HF_HOME='I:\scratch_fea\b0652085_tmp\j2_esm2\hf_cache'
$env:TORCH_HOME='I:\scratch_fea\b0652085_tmp\j2_esm2\torch_home'
$env:HF_HUB_DISABLE_XET='1'
$env:HF_HUB_DISABLE_SYMLINKS_WARNING='1'
.\.venv\Scripts\python.exe scripts\esm_zeroshot_dms.py `
  --data-dir I:\scratch_fea\b0652085_tmp\j2_esm2\proteingym `
  --model facebook/esm2_t33_650M_UR50D --max-assays 40
```

Data: ProteinGym v1.3 substitutions benchmark (official Zenodo record) + the leaderboard `pg_spearman_dms.csv`. Field context to beat/match: ESM2-650M 0.484 · GEMME 0.484 · TranceptEVE-L 0.475.

---

## Discipline rails (do not skip — this project's anti-confound spine)

1. **Freeze the falsifier BEFORE the GPU run.** No post-hoc bar-moving.
2. **De-confound every claim.** Within-lineage / clade-blocked split + the negative-control gauntlet. An overall AUROC that conflates lineage + mechanism is not a result.
3. **Beat the DOMAIN-KNOWLEDGE baseline, not just k-mer.** Beating k-mer ≠ working.
4. **A clean FAIL is a shippable result.** The honest-close branch is pre-authorized.
5. **Frozen AMR surface is byte-locked** (`amr_rules.py` + `calibrated_amr_rules.json`) — the J-track must not touch it.

---

## Branch note (so nothing gets lost)

- `mosafer` (this branch) = J-track + older PGP-UK/ClinVar scratch. **18 commits unique to here** (J2 bank + PGP-UK cohort rows 342 + ClinVar Mendelian) that are NOT on `main`/`mosfaer`.
- `mosfaer` = up-to-date mirror of `main`'s full decoder suite (merged 2026-07-07, `f75bfdf`).
- **Consolidating `mosafer` ↔ `mosfaer` is an open user/authority decision** — flagged, not executed. If you bank J3/J4 results here, they stay on `mosafer` until that consolidation call is made.

## Reply requested

When you pick up: (1) which #1/#2/#3 you're running, (2) the pre-registered falsifier you froze, (3) the candidate paired-label modality for J4 (if you reach it) so we de-confound-screen it before you spend GPU time.
