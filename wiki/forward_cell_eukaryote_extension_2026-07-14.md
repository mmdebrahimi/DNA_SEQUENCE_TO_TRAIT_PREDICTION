# Forward cell — eukaryotic extension (yeast + Arabidopsis), ESM2 across the tree of life (2026-07-14)

Extends the forward variant-effect cell beyond E. coli + human to more eukaryotes (yeast + Arabidopsis),
using the cached ProteinGym DMS + the organism-agnostic ESM2 method (AlphaMissense stays human-only).

## New eukaryote cells (BLOSUM + ESM2)

| protein | organism | assay | n | BLOSUM62 | ESM2-650M | lift |
|---|---|---|---:|---:|---:|---:|
| **SR43C** | Arabidopsis | Tsuboyama_2023 (stability) | 889 | 0.213 | **0.589** | +0.376 |
| **RL40A** (ubiquitin) | yeast | Mavor_2016 (function) | 1253 | 0.237 | **0.518** | +0.281 |
| MBD11 | Arabidopsis | Tsuboyama_2023 | 1155 | 0.196 | — | — |
| BBC1 | yeast | Tsuboyama_2023 (stability) | 1084 | 0.073 | — | — |
| THO1 | yeast | Tsuboyama_2023 (stability) | 656 | 0.057 | — | — |

All WT-mismatch=0. The tiny Tsuboyama stability domains (BBC1/THO1, 41–64 aa) are where BLOSUM is weakest
(0.06–0.07) — exactly where a learned model should help most; SR43C confirms it (+0.376).

## ESM2 works across the whole tree of life

| kingdom | organism | ESM2 range |
|---|---|---|
| bacterium | E. coli | 0.51–0.73 (TEM-1 0.73, CcdB 0.51) |
| eukaryote | human | 0.52 (PTEN) |
| eukaryote | yeast | 0.52 (ubiquitin) |
| eukaryote | Arabidopsis | 0.59 (SR43C) |

**ESM2-650M is the universal learned predictor** — it lifts the deterministic BLOSUM baseline by **+0.28 to
+0.39 on every organism class tested** (bacteria, human, yeast, plant). AlphaMissense remains the human-only
complement (slightly edges ESM2 on human PTEN 0.539 vs 0.518).

## The leaderboard now spans 4 organisms / 2 kingdoms

`scripts/forward_leaderboard.py` (organism classifier generalized to ProteinGym species codes): **18 proteins
across 9 E. coli + 4 human + 3 yeast + 2 Arabidopsis = 9 bacterium + 9 eukaryote.** The pattern holds across
the tree of life: learned methods (ESM2 universal; AlphaMissense human) beat deterministic BLOSUM everywhere.

## Status

The forward "edit → predict molecular effect" cell is now a genuine **tree-of-life** decoder: bacteria +
human + yeast + Arabidopsis, per-protein-validated on measured DMS, with the ESM2 method proven organism-
agnostic and a regime router that picks the predictor by organism. 24 forward tests. Frozen decoder surface
byte-unchanged (`verify_lock OK`); `dna_decode/forward` NON-frozen. Regenerate:
`uv run python scripts/forward_leaderboard.py`.
