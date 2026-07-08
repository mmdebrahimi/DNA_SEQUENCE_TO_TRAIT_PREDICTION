# J2 Phase 2b — LoRA fine-tune ESM-2 (genuine training), pre-registered 2026-07-08

The deeper world-model bet. Phase 2 beats 0.48 by scaling/combining **frozen** open weights (inference).
Phase 2b asks: **does LoRA-adapting the ESM-2 backbone teach it transferable variant-effect biology that
generalizes to proteins it never trained on?** That is the genuine "improve the world model" claim.

## The honest bar — cross-protein, leakage-controlled (NOT within-assay CV)

- **Split PROTEINS into folds** so no protein's assays appear in both train and test
  (`split_proteins` — verified leakage=0 on the real 217-assay reference: fold 0/5 = 175 train / 42 test
  assays, 148 / 38 proteins, protein-leak **0**).
- **LoRA-adapt** ESM-2 on the TRAIN proteins' variants (pairwise margin ranking loss vs the wet-lab DMS).
- **Evaluate masked-marginals on the HELD-OUT proteins** with BOTH the base model and the adapted model,
  reusing the Phase-1 drift-guarded `score_assay` (eval is byte-identical to the zero-shot path).
- **HEADLINE = `finetuned_median − zeroshot_median` on the SAME held-out fold.**
  - `> +0.005` (a real margin, ideally across folds) → fine-tuning improved the world model.
  - `≤ 0` → it only memorized train proteins; the frozen zero-shot path (Phase 2) is the better product.

**Why not within-assay CV:** training a head on 80% of an assay and predicting its own 20% mostly measures
the readout regressor, not the backbone's knowledge (R2 framing check). Comparing to the global 0.48 is also
dishonest here — the fold is a subset, so the base-vs-adapted delta on the *identical* held-out set is the
fair comparator.

## Anti-overfit / honesty rails

- Fold assignment is **frozen by seed** before any training; report every fold, no cherry-picking.
- The comparator (base model) is evaluated on the **exact same** held-out assays as the adapted model.
- Peft/torch are lazy-imported; the leakage-safe split + pair sampler are pure + unit-tested (5 tests in
  `tests/test_j2_phase2b_lora.py`). The training/eval loop is real but **GPU-run** — the number is the
  user's kernel run, not claimed here.

## Turnkey (free Kaggle/Colab T4/P100; data attached per the runbook Step 0)

```bash
# 0. no-GPU sanity: confirm the leakage-safe split for this fold
python j2_phase2b_lora_finetune.py --dry-run --data-dir /kaggle/input/<slug> --nfolds 5 --fold 0

# 1. one fold (GPU): LoRA-adapt on train proteins, report zero-shot vs finetuned on held-out proteins
pip install peft
python j2_phase2b_lora_finetune.py --data-dir /kaggle/input/<slug> --nfolds 5 --fold 0 \
    --model facebook/esm2_t33_650M_UR50D --dtype float16 --rank 8 --epochs 2 --out j2b_fold0.json
# repeat --fold 1..4 (or shard across kernels) -> average the per-fold deltas for the headline.
```

Upload BOTH `j2_phase2b_lora_finetune.py` and `j2_phase1_esm2_proteingym.py` (the scaffold reuses the
Phase-1 scoring core). ESM2-650M + LoRA rank-8 fits a free T4 comfortably; 3B is the stretch.

## Scope / provenance

Self-contained runner: `notebooks/j2_phase2b_lora_finetune.py`. Reuses the drift-guarded Phase-1 core.
Committed to side branch `soraya-j2-phase2` (main = DNA-11's frontier; merge is user authority). Frozen AMR
surface byte-unchanged. Genome-embedding scaling stays CLOSED. Gated value: Phase-2b matters most if
Phase-2 (3B/ensemble) does NOT already clear the decoder's needs — run Phase-2 first, then decide.
