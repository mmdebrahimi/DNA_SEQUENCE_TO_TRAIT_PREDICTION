# ESM2-650M at scale: the learned molecular ceiling vs the shipped wheel-only default (2026-07-17)

**The learned method (ESM2-650M) transforms both molecular cells from weak-to-modest to strong across
~190 ProteinGym proteins. Tonight's tiny-N ESM headlines generalize; the shipped `blosum62` default does
not.** Computed on a free Kaggle T4 (~1.5 h, 201 masked-marginal tables, 0 skipped).

| capability | shipped default (`blosum62`, wheel-only) | **learned (ESM2-650M, GPU)** | ratio |
|---|---:|---:|---:|
| **INVERSE** — rank inverse materially beats a random pick | 13.5% (N=200) | **72.9%** (N=188) | **5.4×** |
| **INVERSE** — any positive edge | 59.0% | **93.6%** | — |
| **FORWARD** — median \|Spearman\| vs measured DMS | 0.201 (N=209) | **0.493** (N=194) | **2.4×** |

## What this settles

Tonight (2026-07-17) two things were established on FOUR proteins: the ESM rank inverse beat the null 4/4,
and the ESM forward Spearman was ~0.52. The overnight N=200 blosum sweeps then showed the **shipped default**
is weak (inverse 13.5% material; forward median 0.20) — an honest bound that corrected an overclaim. The
open question was whether the *learned* headline was itself a 4-protein fluke.

**It is not.** At N≈190:
- the ESM inverse materially beats a random pick on **73%** of proteins
  (median margin +36.5%), vs 14% for blosum;
- the ESM forward reaches a median \|Spearman\| of **0.49** (matching the
  known ProteinGym ESM2-650M benchmark ~0.49), vs 0.20 for blosum.

So the molecular cells' real capability is REAL and general — it just lives in the learned method, which
needs a GPU. The wheel-only default a GPU-less `pip install` user gets is the weak floor; ESM is the strong
ceiling, now measured at scale rather than on a hand-picked handful.

## What ships / what doesn't

- The 201 ESM2-650M masked-marginal tables (17.9 MB) are NOT committed (too large for git). Regenerate them
  with the Kaggle kernel `emanueleebrahimi/proteingym-esm2-650m-tables` (fair-esm masked-marginal,
  ~1.5 h on a free T4) or locally with a GPU via `scripts/tem1_forward_cell.py --method esm2`.
- The shipped CLI default STAYS `blosum62` (no GPU assumption). The Python API path
  `predict_effect(..., method="esm2", esm_table=...)` / `propose_edits(..., method="esm", esm_table=...)`
  is the strong path for anyone with a GPU or a precomputed table.
- Honest scope: ProteinGym single-substitution assays; censored assays excluded; the inverse RANKS, never
  doses (the magnitude form remains non-deployable). Blosum still WINS on convenience (instant, offline) and
  is not useless — it just is not a reliable design tool on its own.

## Provenance
`wiki/proteingym_inverse_sweep{,_esm}_2026-07-17.json` · `wiki/forward_{blosum,esm}_proteingym_2026-07-17.json`
· scripts `forward_inverse_proteingym_sweep.py --method {blosum62,esm}` +
`forward_blosum_proteingym_sweep.py --method {blosum62,esm}` · Kaggle kernel + `scripts/` (this session).
