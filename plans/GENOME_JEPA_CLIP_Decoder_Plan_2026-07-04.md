# Genome JEPA / CLIP Embedding-Predictor — Decompose + Plan + Pre-registered Falsifier (2026-07-04)

> **Session split:** the LD **imputer** (statistical masked-SNP prediction) is being built in the *DNA AI
> Decode 11* session. **This** track builds the **neural / world-model counterpart**: a genome-**JEPA**
> (predict masked-window *embeddings* from context) — and its CLIP variant. Same task as the imputer, one
> level up (latent, not genotype). Companion: `plans/GENOME_WORLD_MODEL_Scoping_Brief_2026-07-04.md`.
> **Status:** J1 prototype BUILT + VERIFIED (this session); J2+ workhorse-GPU-gated.

## 0. The three paradigms, precisely (answers "is it JEPA?")
- **JEPA (LeCun's world model) = YES, this is it.** Predict the **embedding** of a masked/target region from
  the **context embedding**, in latent space, with an **EMA target-encoder + stop-gradient** to avoid
  representation collapse. It does NOT reconstruct the input. I-JEPA (images), V-JEPA (video) → **genome-JEPA**
  (masked genomic window). "Predicts the embeddings" = JEPA. ✅
- **MAE / masked-LM** (DNABERT/NT) — predict the masked *tokens* (input space). Different; what the frozen FMs do.
- **CLIP** — contrastive *alignment of two modalities* (InfoNCE), e.g. sequence↔phenotype. Not masked prediction.
- **Honest link:** JEPA is the **neural masked-predictor**; the LD imputer is the **statistical masked-predictor**.
  Both exploit haplotype/local structure — the exact signal that is a *feature* for imputation and a *confound*
  for phenotype-prediction (the 0-for-5 wall). So JEPA's honest home is representation/imputation, NOT phenotype.

## 1. Honest prior (integrity rail — do not bury)
dna_decode learned-embedding→PHENOTYPE is **0-for-5** de-confounded negative. JEPA is **representation
learning**, not phenotype prediction — so it is not directly refuted by that. BUT its *value* is unproven:
a JEPA rep only earns its keep if it beats (a) the **frozen FM embeddings** (nt.h5/dnabert2.h5) AND (b) an
**LD/statistical baseline** on a concrete downstream (masked-SNP imputation or variant-effect) under a
de-confounded split. Default prior after 0-for-5: skepticism. The plan is built to **falsify**, not confirm.

## 2. J1 — the prototype (BUILT + VERIFIED this session, CPU)
`scripts/genome_jepa_prototype.py` + `tests/test_genome_jepa_prototype.py`. Implements the real JEPA
mechanism: `context_encoder` (online) + `target_encoder` (EMA, stop-grad) + `predictor`; mask a contiguous
window; predict its target-embeddings from the pooled context; latent MSE; EMA update; anti-collapse monitor.
**Verified on synthetic Markov-DNA (learnable local structure):**
- loss 0.87 → 0.18 (trains), target-variance 0.176 > 0 (**no collapse** — the load-bearing JEPA property),
- **shuffled negative control** stays higher (0.33) → it captures **real local dependency**, not a constant.
This proves the MECHANISM is correct. It is NOT a genome model (synthetic data, CPU, tiny).

## 3. Decompose — build families
| Fam | Deliverable | Depends | Gate |
|---|---|---|---|
| **J1 prototype** | correct JEPA mechanism (EMA/stop-grad/anti-collapse), verified | — | ✅ DONE (laptop) |
| **J2 real encoder** | swap synthetic → real sequence windows; init from an open FM (Caduceus/NT) or train small | J1 | **compute (GPU)** + data (D:) |
| **J3 downstream falsifier** | JEPA rep vs frozen-FM + LD baseline on masked-SNP imputation / variant-effect | J2 | de-confounded slice |
| **J4 CLIP head** | contrastive align JEPA seq-rep ↔ a phenotype/expression modality | J2 | paired de-confounded labels |
| **J5 scale** | longer context, bigger corpus, in-silico perturbation ("world model" proper) | J2-J4 | **compute/$** |

**Critical path:** J1 ✅ → J2 (GPU) → J3 (falsifier). J4/J5 only if J3 shows the rep has value.

## 4. --plan — phased (each gated + falsifiable)
- **Phase 1 (DONE, laptop, $0):** J1 prototype — mechanism verified. *This session.*
- **Phase 2 (workhorse GPU / cloud):** J2 — real-sequence JEPA (init from an open FM), on a de-confounded
  organism slice (Arabidopsis 1001G / openSNP EUR / a bacterial panel). *Gate:* GPU (RTX 3500 Ada VRAM check
  or cloud budget = **money gate → user approval**).
- **Phase 3 (laptop-feasible once J2 embeddings are cached):** J3 falsifier — does the JEPA rep + a small head
  beat frozen-FM + LD on the pre-registered downstream? Reuse the de-confound gate + CI-aware falsifier.
- **Phase 4 (compute-gated):** J4 CLIP head + J5 scale/perturbation — only if Phase 3 passes.

## 5. Pre-registered falsifier (freeze BEFORE the GPU run — the project's anti-confound discipline)
- **Task:** masked-SNP imputation accuracy (LOO) OR ProteinGym/AlphaMissense variant-effect ranking (Spearman).
- **Comparators (E19 gauntlet):** frozen-FM embedding (nt.h5/dnabert2.h5) + linear head · LD/statistical
  baseline (the imputer) · k-mer/CNN · shuffled-label + PCA-only + held-out-clade negative controls.
- **Split:** de-confounded, lineage/ancestry-blocked, clade-stratified CV (reuse `within_lineage_diagnostic`).
- **PASS bar:** JEPA rep beats BOTH frozen-FM AND the LD baseline on the full metric surface (not one metric)
  by a meaningful margin, AND the within-lineage diagnostic rules out population structure.
- **FAIL → honest close:** if JEPA does not beat frozen-FM + LD, the learned-representation frontier for solo-
  scale G2P is closed (consistent with the 0-for-5 phenotype evidence); the deterministic scan + LD imputer
  remain the product.

## 6. Reproduce (J1)
```bash
uv run python scripts/genome_jepa_prototype.py            # trains + asserts mechanism (structured vs shuffled)
uv run pytest tests/test_genome_jepa_prototype.py -q      # 6 offline mechanism tests
```

## 7. One-line
The JEPA decoder is real and its mechanism is now built + verified on the laptop (J1). Everything past the
mechanism (real sequence, value on a downstream) is compute-gated and **pre-registered to be falsified** —
because after 0-for-5, a learned representation must *prove* it beats the frozen FM + the LD imputer, not be
assumed to.
