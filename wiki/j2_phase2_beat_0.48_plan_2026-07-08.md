# J2 Phase 2 — beat ESM-2 ProteinGym 0.48 on FREE compute (pre-registered, 2026-07-08)

> ## STATUS 2026-07-09: **RESOLVED — bar already met; Lever 1 FALSIFIED. Do not run this plan as written.**
>
> The pre-registration below is preserved verbatim (it is a pre-registration; it gets annotated, never
> rewritten). What actually happened, from a real Kaggle T4 run over all 217 assays
> (`wiki/proteingym_esm2_650m_full_2026-07-09.md`):
>
> **1. The PASS bar is already cleared by the FROZEN 650M, with zero levers applied.** Our run scored
> **median 0.490** over 217/217 assays (windowing on, no ensemble, no 3B) against the pre-registered bar
> of `median |Spearman| > 0.48`. Only 3/217 assays have negative rho, so the signed median equals the
> |rho| median — apples-to-apples. **Phase 2's goal is met without spending a single lever.**
>
> **2. Lever 1 ("ESM2-3B … 650M ≈0.48 → 3B ≈0.51 (published). HIGH — the certain beat") is WRONG, and
> wrong in direction.** Recomputed from ProteinGym's own per-assay CSV (vendored at
> `wiki/refs/proteingym_v1.3_DMS_level_Spearman.csv`):
>
> | ESM2 | per-assay median | published `Average_Spearman` |
> |---|---|---|
> | **650M** | **0.484** | 0.414 |
> | 3B | 0.467 | 0.406 |
> | 15B | 0.438 | 0.400 |
>
> **650M is the PEAK of the ESM2 family.** Running Lever 1 would burn 1–3 GPU-h to score *below* the bar,
> and would read as a project failure when it is a known, published regression. Independently corroborated
> on humsavar: 650M → 3B is worse on 4/7 proteins, mean paired delta −0.021
> (`wiki/esm_am_ensemble_paired_2026-07-09.md`).
>
> **3. Lever 2 (long-protein windowing) is already done** and is included in the 0.490 —
> `scripts/kaggle_proteingym_sweep.py` windows proteins >1022 aa; 217/217 scored, 0 dropped, 0 errors.
>
> **4. Lever 3 (650M + 3B ensemble) now has a weak prior**, since it would ensemble the peak checkpoint
> with a strictly worse one. Not formally falsified — cross-*checkpoint* rank-averaging is a different
> pairing from the cross-*method* ESM+AlphaMissense ensemble that showed no paired lift — but it is no
> longer a "med-HIGH" bet. If run, judge it on the **paired** per-assay delta vs best-single, never on a
> difference of medians (that error is what produced the false "ENSEMBLE helps" verdict; see `b82aa79`).
>
> **Where the real headroom is** (official leaderboard; ESM2-650M is rank 45 / 97):
>
> | direction | ceiling | delta vs 650M |
> |---|---|---|
> | same modality (single sequence) | 0.458 (VespaG, rank 16) | +0.044 |
> | **+ structure** (ProSST, rank 3) | **0.507** | **+0.093** |
> | + MSA / retrieval (rank 1–2) | 0.518 | +0.104 |
>
> Scale is dead; naive score-averaging is dead. The remaining levers are **modality** (structure, MSA) and
> **training** (Phase 2b LoRA) — not more parameters.

**Goal:** improve our genomics world model in the ONE regime where a learned model is *supposed* to win —
molecular-property (protein variant-effect) — by beating the published ESM2-650M ProteinGym zero-shot
number (**median |Spearman| ≈ 0.48**) on a **free** T4/P100 kernel. No money, no training in this phase
(inference + ensembling only).

## Why this track (honest regime placement)

Per the project's regime-boundary finding ([[feedback_g2p_decoder_regime_boundary]]): curated-catalog →
deterministic wins; organism-polygenic → neither (genome embeddings are **0-for-4 de-confounded**, a
CLOSED negative — do NOT scale on a bigger GPU); **molecular-property → learned wins (~2.1× BLOSUM).**
J2/ESM-2 on ProteinGym is the molecular-property cell. This is the world-model direction where free
compute genuinely helps.

## Three levers (implemented 2026-07-08; CPU-tested; GPU run is the user's)

| # | Lever | Flag | Expected effect | Certainty |
|---|---|---|---|---|
| 1 | **ESM2-3B in fp16** | `--model facebook/esm2_t36_3B_UR50D --dtype float16` | 650M ≈0.48 → 3B ≈**0.51** (published). Fits a free T4 16 GB in fp16 (~6 GB weights). Pure inference. | **HIGH** — the certain beat |
| 2 | **Long-protein windowing** | `--long-mode window` | Recovers assays with proteins >1022 aa (currently DROPPED). **Verified on the real reference: 16/217 assays (7.4%) exceed 1022 aa** — the exact recovery scope. **Completeness/honesty, NOT a guaranteed median lift** — long proteins may score lower; the number becomes more complete, not necessarily higher. | med (direction unknown by design) |
| 3 | **Model ensemble** (650M + 3B) | two runs `--keep-scores` → `--ensemble-merge` | z-score rank-average across checkpoints — the classic ProteinGym top approach; reliably ≥ best single model. Free (just more inference). | med-HIGH |

## Pre-registered PASS bar (derived, not asserted — R2)

- **PASS = median |Spearman| > 0.48** (the `STRETCH` constant) on the **FULL joinable cohort** (all ~217
  assays in one shard or merged shards — NOT a cherry-picked subset), with the shuffled control < 0.05.
- Report the number **honestly whatever it is.** Lever 1 alone should clear it per published 3B numbers; if
  it does not on our joinable subset, that itself is the finding (our cohort ≠ the published cohort).
- **Anti-overfit rail:** we do NOT tune `--long-mode`/ensemble membership against the score. Pick the config
  a priori (3B + window + {650M,3B} ensemble), run once, report.

## Pre-run guard (no GPU — run FIRST, added 2026-07-08)

Confirm the ProteinGym data is correctly attached BEFORE spending GPU minutes:

```bash
python j2_phase1_esm2_proteingym.py --self-test --data-dir /kaggle/input/<slug>
# -> "{have_csv}/{total} assays have their per-assay CSV attached; N exceed maxlen ..."
#    exit 0 = READY; exit 1 = data not attached (fix the dataset mount before the GPU cell).
```
The data-loading path is R3-verified against the REAL official ProteinGym reference (217 assays parse;
16 exceed 1022 aa). A green self-test means the GPU run won't fail at load time.

## Turnkey commands (free Kaggle kernel — GPU T4×2 or P100, Internet ON for weights)

```bash
# Lever 1 (the certain beat) — one cell:
python j2_phase1_esm2_proteingym.py --data-dir /kaggle/input/<proteingym-slug> \
    --model facebook/esm2_t36_3B_UR50D --dtype float16 --long-mode window \
    --out j2_3b.json
# -> prints median |Spearman|; PASS if > 0.48.

# Lever 3 (ensemble 650M + 3B) — three cells (last needs no GPU):
python j2_phase1_esm2_proteingym.py --data-dir /kaggle/input/<slug> --model facebook/esm2_t33_650M_UR50D \
    --dtype float16 --long-mode window --keep-scores --out j2_650m.json
python j2_phase1_esm2_proteingym.py --data-dir /kaggle/input/<slug> --model facebook/esm2_t36_3B_UR50D \
    --dtype float16 --long-mode window --keep-scores --out j2_3b.json
python j2_phase1_esm2_proteingym.py --ensemble-merge j2_650m.json j2_3b.json
```

Data: upload `D:/dna_decode_cache/proteingym` (217 assays) once as a private Kaggle Dataset (see
`notebooks/J2_PHASE1_RUNBOOK.md` Step 0). **D: was disconnected on the build host 2026-07-08** — reconnect
it (or confirm the dataset is already on Kaggle) before the run.

## Honest scope / what this is NOT

- This is **stronger + more complete + ensembled use of free open weights**, NOT a novel or fine-tuned
  model. The claim is "our harness beats 0.48 on free compute", not "we built a better protein LM".
- **Phase 2b (deferred, gated):** a real supervised probe / LoRA fine-tune with per-assay CV — genuine
  training, worth doing ONLY if beating 0.48 powers a concrete decoder feature (unchanged gate from the
  runbook). Not in this phase.
- **Not verified locally:** the actual number needs a GPU + the ProteinGym data — the user's kernel run.
  Everything CPU-testable (window math, ensemble combine, merge wiring) is unit-tested (14 tests) +
  CLI-smoked. Genome-embedding scaling stays CLOSED (do not reopen on GPU).

## Provenance

Runner: `notebooks/j2_phase1_esm2_proteingym.py` (+ canonical `scripts/esm_zeroshot_dms.py`, kept in sync;
drift guard `tests/test_j2_phase1_notebook.py` now pins `window_for_position` + `combine_variant_scores`
too). Committed non-collidingly to side branch **`soraya-j2-phase2`** (main = DNA-11's frontier); merge to
main is the owner's/user's call.
