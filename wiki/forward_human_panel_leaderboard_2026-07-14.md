# Forward cell — human DMS panel + 3-method leaderboard + AlphaMissense routing (2026-07-14)

Completes the Phase-6 human forward cell: a 4-protein human DMS panel, AlphaMissense wired into the router
demo, and a per-protein BLOSUM / ESM2 / AlphaMissense leaderboard across both organism classes.

## Human DMS panel — AlphaMissense beats BLOSUM on every protein

| protein | UniProt | assay | n | BLOSUM62 | AlphaMissense | lift |
|---|---|---|---:|---:|---:|---:|
| CYP2C9 | P11712 | Amorosi_2021 (abundance) | 6,370 | 0.333 | **0.598** | +0.265 |
| TPMT | P51580 | Matreyek_2018 | 3,648 | 0.240 | **0.558** | +0.318 |
| PTEN | P60484 | Mighell_2018 (function) | 7,260 | 0.182 | **0.539** | +0.357 |
| MSH2 | P43246 | Jia_2020 | 16,749 | 0.164 | **0.416** | +0.252 |

All join AlphaMissense at offset 0 with 100% coverage + 0 WT-mismatch. AM roughly 2–3× the deterministic
baseline across all four.

## 3-method leaderboard (`scripts/forward_leaderboard.py`) — 13 proteins, both kingdoms

| protein (assay) | organism | n | BLOSUM62 | ESM2-650M | AlphaMissense |
|---|---|---:|---:|---:|---:|
| BLAT (Stiffler_2015) | E. coli | 4996 | 0.346 | **0.732** | — |
| CcdB (Tripathi_2016) | E. coli | 1663 | 0.248 | **0.511** | — |
| CYP2C9 (Amorosi_2021) | human | 6370 | 0.333 | — | **0.598** |
| TPMT (Matreyek_2018) | human | 3648 | 0.240 | — | **0.558** |
| **PTEN (Mighell_2018)** | human | 7260 | 0.182 | **0.518** | **0.539** |
| MSH2 (Jia_2020) | human | 16749 | 0.164 | — | **0.416** |

(+ 7 more E. coli BLOSUM-only rows.) **PTEN is the full 3-method comparison:** BLOSUM 0.182 « ESM2 0.518 ≈
AlphaMissense 0.539 — both learned methods ~0.52–0.54, AM edges ESM2 slightly. The pattern:

- **Bacteria:** ESM2 ≫ BLOSUM (AlphaMissense N/A — human-only).
- **Human:** AlphaMissense ≳ ESM2 ≫ BLOSUM.
- Learned methods beat deterministic BLOSUM by **+0.25 to +0.39** everywhere they run.

## AlphaMissense in the router (`scripts/forward_router_demo.py`)

The router now demonstrates the human learned path: **PTEN K13F → damaging (AM 1.0)**, **PTEN T2A → preserved
(AM 0.05)**, alongside the bacterial ESM2 (TEM-1) + WHO-catalogue (TB) + abstain (organismal) regimes. The
router's Regime-B picks BLOSUM / ESM2 / AlphaMissense per protein by organism.

## Status

The forward "edit → predict molecular effect" cell now spans **two organism classes with three methods**,
per-protein-validated on measured DMS: E. coli (BLOSUM + ESM2) and human (BLOSUM + ESM2 + AlphaMissense),
plus the regime router that picks the right predictor. 24 forward tests. Frozen decoder surface byte-unchanged
(`verify_lock OK`); `dna_decode/forward` NON-frozen. Regenerate: `uv run python scripts/forward_leaderboard.py`.
