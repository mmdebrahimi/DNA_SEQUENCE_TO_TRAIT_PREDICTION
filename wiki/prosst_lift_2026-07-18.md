# ProSST structure hybrid — the lift is REAL and POWERED (Step-6 run, LOCAL)

**Date:** 2026-07-18 · **Model:** `AI4Protein/ProSST-2048` (transformers, CPU, masked-LM log-ratio) ·
**Structures:** ProteinGym's own pre-quantized structure sequences (vocab 2048) · **Script:**
`scripts/prosst_lift.py` · **Data:** `data/processed/prosst_lift_checkpoint.jsonl`

## Question

The modality-hybrid sweep found `ESM2 + ProSST` +0.05 paired vs ESM2-650M (win 87%) — the biggest single
modality lever — but only from ProteinGym's PRECOMPUTED `ProSST-2048` column. This Step-6 run computes ProSST
scores **ourselves** and checks (1) our scorer reproduces the column and (2) the `ESM2 + our-ProSST` hybrid
lift reproduces. **Ran entirely LOCALLY on CPU — no Kaggle** (see "Local, not Kaggle" below).

## Result — reproduction EXACT, lift REAL and POWERED

**(1) Our scorer is exact.** Over **56 proteins**, our ProSST reproduces ProteinGym's `ProSST-2048` column at
**Spearman 1.0000 (min 1.0000)** — given the same pre-quantized structure tokens + the same model, our wiring
(structure-token +3 shift + `[1,…,2]` wrap, CLS/EOS strip, masked-LM log-ratio) + the force-tie fix (below)
reproduce the reference byte-for-byte. The scorer is provably correct.

**(2) The hybrid lifts — powered, and on every phenotype.** Paired median Δ (ESM2+our-ProSST − ESM2):

| phenotype | n | median Δ | win-rate |
|---|---:|---:|---:|
| **Expression** | 9 | **+0.1028** | 9/9 |
| **Stability** | 22 | **+0.0893** | 21/22 |
| Binding | 4 | +0.0847 | 4/4 |
| OrganismalFitness | 12 | +0.0201 | 11/12 |
| Activity | 9 | +0.0135 | 7/9 |
| **overall** | **56** | **+0.0668** | **52/56 (93%), sign-p 1.1e-11** |

Median hybrid |Spearman| **0.594** vs ESM2 0.515 vs ProSST-alone 0.585. **Every category is positive**, the
lift is largest on the structure-dominated phenotypes (Expression +0.103, Stability +0.089) exactly as the
modality-per-category finding predicted, and the overall sign test is **1.1e-11** (highly significant, n=56).

**This is the inverse of the MSA-Transformer evolution run** (`wiki/msa_transformer_lift_2026-07-17.md`),
which reproduced its column (0.84) but did NOT lift (+0.0008, 52%). Here the structure hybrid reproduces
EXACTLY (1.0) AND lifts strongly (+0.067, 93%, p=1e-11). **Structure is the modality lever; the fast-evolution
path was not.** The forward-cell's biggest measured lever is now realized end-to-end in our own pipeline.

## Honest scope

- **Reproduction is exact BY CONSTRUCTION** — we score with ProteinGym's own pre-quantized structure tokens
  (not an independently quantized structure), so 1.0 proves our SCORER wiring is correct, not that an
  independent structure would give the same tokens. The deployable finding is the **LIFT** (+0.067, powered),
  which is a genuine paired comparison against ESM2 on real DMS.
- For a NOVEL protein (no pre-quantized tokens), `prosst_variant_table(pdb_path=…)` needs `torch_geometric`
  for `PdbQuantizer` (absent on this host — the same wall as ESM-IF); run the quantizer on a Linux/GPU host
  or reuse ProteinGym's tokens. The transformer forward itself runs on CPU.
- N=56 (pre-quantized ∩ ESM2-table ∩ seq_len≤400). Powered for the aggregate + Stability; the smaller
  categories (Binding n=4, Activity/Expression n=9) are directional but consistently positive.

## Two real bugs fixed by verify-in-batch

1. **Missing MLM decoder weight → force-tie.** The ProSST checkpoint OMITS `cls.predictions.decoder.weight`,
   expecting it tied to the input embeddings. Newer `transformers` no longer auto-ties (the config uses a
   legacy tie key) → it loads a RANDOM head → garbage logits → **reproduction 0.013**. Force-tying
   (`model.cls.predictions.decoder.weight = model.get_input_embeddings().weight`) → **reproduction 1.0**.
2. **Canonical ss wiring.** The structure tokens shift **+3** (past `<pad>/<cls>/<eos>`) and wrap `[1,…,2]`;
   logits strip CLS/EOS (`[:, 1:-1]`) so residue at 1-based pos → 0-based `pos-1` (from the repo's
   `zero_shot/proteingym_benchmark.py`). The draft's `[0,…,0]` wiring gave ~0 reproduction.

## Local, not Kaggle

The plan tagged Step 6 as an attended Kaggle run. It ran LOCALLY instead: `transformers` is installed and the
ProSST model caches to D:, and the ProSST **transformer forward needs only `transformers` + torch — NOT
`torch_geometric`** (that is only for the `PdbQuantizer`). ProteinGym's pre-quantized structure tokens
(Google-Drive `proteingym_benchmark.zip`, `structure_sequence/2048/{DMS_id}.fasta`) skip the quantizer
entirely, so the whole reproduce-the-column-then-lift validation runs on CPU. No Kaggle, no data-upload bug.

## Shipped

- `dna_decode/forward/prosst_scorer.py` — canonical ss wiring + the force-tie fix (the seam is now a WORKING
  local scorer, not just a mock-tested stub).
- `scripts/prosst_lift.py` — reads ProteinGym pre-quantized `.fasta` structures; the validation harness.
- Frozen decoder surface byte-unchanged (`verify_lock OK`).
