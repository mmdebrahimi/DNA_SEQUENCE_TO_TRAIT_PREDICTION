# The modality-hybrid test — is there a move past the ESM2-650M molecular ceiling?

**Date:** 2026-07-17 · **Substrate:** ProteinGym `pg_zeroshot` per-variant precomputed scores (canonical
reference implementations, 99 model columns) · **Script:** `scripts/forward_modality_hybrid_sweep.py` ·
**Data:** `wiki/forward_modality_hybrid_2026-07-17.json`

## Question

The learned world model works in exactly ONE regime — **molecular fitness (DMS)**. At scale its sequence
baseline is **ESM2-650M**, and two facts were already established:

- **Scale is a dead end** — ESM2 peaks at 650M; 3B/15B regress (`feedback_g2p_decoder_regime_boundary`).
- **The headroom is MODALITY, not parameters** — +evolution (MSA) and +structure.

So the highest-VOI world-model move is not a bigger model — it is a **modality upgrade**, and specifically a
**HYBRID**: does combining *orthogonal* signals (learned sequence + explicit evolution + structure) beat the
best single modality? That is the project's hybrid thesis (deterministic ⊕ learned), one level down.

## Method

For each of the 95 ProteinGym assays where every modality column is present (structure + MSA available),
compute abs-Spearman(model, `DMS_score`) with **mid-rank ties** over ≥20 shared non-NaN variants, then a
**paired** per-protein comparison against ESM2-650M. Naive hybrids = **rank-average** of ≥2 component score
vectors, each oriented by ProteinGym's fixed *higher=fitter* convention (verified in-run: 94/95 baseline
assays positive) — **no label is fit, only the standard orientation**.

**Pre-registered bar** (a candidate "beats ESM2-650M"): PAIRED **median Δ > 0** AND **win-rate ≥ 60%** AND a
two-sided **sign-test p < 0.05**. Falsifiable: if nothing clears it, the single-modality ceiling is real.

## Result — 10 candidates clear the bar

Baseline **ESM2-650M median abs-Spearman = 0.4926 (N=95)** — matches the full-benchmark 0.49275 (N=194), so
the modality subset is not cherry-picked.

| model | modality | median | Δ vs ESM2 | win-rate | sign-p | beats |
|---|---|---:|---:|---:|---:|:--:|
| **HYB_ESM2+GEMME+ProSST** | seq⊕evo⊕struct | **0.5468** | **+0.056** | **0.905** | ~0 | ✅ |
| HYB_GEMME+ProSST | evo⊕struct | 0.5330 | +0.052 | 0.821 | ~0 | ✅ |
| VenusREM | canonical hybrid | 0.5309 | +0.039 | 0.695 | 2e-4 | ✅ |
| HYB_ESM2+ProSST | seq⊕struct | 0.5304 | +0.050 | 0.874 | ~0 | ✅ |
| ProSST-2048 | structure | 0.5203 | +0.028 | 0.611 | 0.040 | ✅ |
| S3F | structure | 0.5043 | +0.012 | 0.600 | 0.049 | ✅ |
| **HYB_ESM2+GEMME** | seq⊕evo | **0.5004** | **+0.022** | **0.842** | ~0 | ✅ |
| HYB_ESM2+SaProt | seq⊕struct | 0.4997 | +0.029 | 0.853 | ~0 | ✅ |
| SaProt_650M_AF2 | structure | 0.4988 | +0.015 | 0.632 | 0.013 | ✅ |
| HYB_ESM2+MSA_T | seq⊕evo | 0.4754 | +0.013 | 0.663 | 0.002 | ✅ |
| — ESM2_650M — | sequence | 0.4926 | base | — | — | |
| GEMME | evolution | 0.4783 | −0.010 | 0.421 | 0.15 | ✗ |
| ESM2_3B | sequence (bigger) | 0.4493 | −0.022 | 0.368 | 0.013 | ✗ |
| ESM2_15B | sequence (bigger) | 0.4446 | −0.038 | 0.368 | 0.013 | ✗ |
| EVmutation | evolution | 0.4234 | −0.055 | 0.263 | ~0 | ✗ |
| Site_Independent | evolution | 0.3865 | −0.061 | 0.200 | ~0 | ✗ |
| ESM2_8M | sequence (tiny) | 0.2327 | −0.183 | 0.105 | ~0 | ✗ |

## Three findings

1. **Scale is dead — re-confirmed on paired terms.** ESM2 8M→35M→150M→**650M (peak)**→3B (−0.022,
   p=.013)→15B (−0.038, p=.013). Bigger sequence models *lose*. The `feedback_g2p_decoder_regime_boundary`
   ESM2-peaks-at-650M finding reproduces here, paired, from the per-variant scores.

2. **The modality lever is LIVE.** Structure (ProSST +0.028, SaProt +0.015, S3F +0.012) and the canonical
   retrieval-hybrid VenusREM (+0.039) each beat ESM2-650M paired. The ProteinGym leaderboard ordering
   (structure/retrieval > sequence) holds on this project's own baseline.

3. **The hybrid thesis is validated — and this is the novel result.** A **naive rank-average** of orthogonal
   modalities beats ESM2-650M on **84–90%** of proteins:
   - `ESM2+GEMME` (**seq⊕evo**): +0.022, **win 84%**, p≈0 — *even though GEMME alone LOSES to ESM2* (−0.010).
     Two signals each ≈ or below ESM2 combine to beat it: that is the orthogonality payoff.
   - `ESM2+ProSST` (**seq⊕struct**): +0.050, win 87%.
   - `ESM2+GEMME+ProSST` (**all three**): **+0.056, win 90.5%** — the top of the board, above every single
     method incl. published SOTA (VenusREM 0.531) on this subset.

   The signature is the **win-rate**: hybrids improve *consistently per protein* (84–90%), not via a few big
   wins — the hallmark of a real ensemble.

   **This refutes the prior note** "ESM+AlphaMissense ensembling gives no paired lift" — that was two
   sequence-ish predictors. **ESM (learned sequence) + GEMME (explicit evolution)** are genuinely orthogonal
   and the lift is real, paired, and robust.

## Deployability (checked before claiming a cell upgrade)

The naive rank-average hybrid **RANKS, needs no calibrator, no DMS, no label** — same deployability class as
the shipped inverse; it does **not** hit the calibrator-transfer wall. What it needs at run time:

| hybrid | needs at run time | infra bar |
|---|---|---|
| `ESM2+GEMME` | ESM2 forward pass + an **MSA** (jackhmmer/hhblits) → GEMME | MSA search |
| `ESM2+ProSST` | ESM2 + a **3D structure** (AlphaFold DB) → ProSST | structure lookup |
| 3-way | both of the above | MSA + structure |

Cheapest deployable upgrade: **`ESM2+GEMME`** (sequence+evolution, +0.022, 84% win). Best overall: the 3-way
(+0.056, 90% win). Both are above the wheel-only BLOSUM v0 / precomputed-ESM2 v0.5 infra tier.

## Honest scope

- **This deepens the ONE working regime (molecular fitness).** It does NOT break the organism-level
  polygenic wall (a closed negative) nor the antagonistic clinical-resistance inversion. The world model's
  reach is unchanged; its *ceiling within its one regime* is raised and the lever is named.
- Scores are ProteinGym's canonical reference implementations (standard benchmark use). The contribution is
  the **paired** analysis against this project's own ESM2-650M baseline + the naive-hybrid construction +
  the deployable-rank framing.
- N=95 is the structure+MSA-available subset; the baseline on it (0.4926) ≈ the full 0.49, so it is
  representative, not selected.

## Shipped

- `dna_decode/forward/variant_effect.py::rank_average_hybrid` — pure, label-free rank-average of ≥2
  precomputed score tables (the deployable form of finding 3); `predict_effect(..., method="hybrid",
  hybrid_tables=[...])`.
- `wiki/forward_modality_hybrid_2026-07-17.{md,json}` (this memo + data).
- Frozen decoder surface byte-unchanged (`verify_lock OK`).
