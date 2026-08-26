# The ParD epistasis anomaly is 5x smaller than published — a mutation-order pooling artifact (2026-08-25)

**Correction to `wiki/forward_epistasis_sweep_2026-07-27.md`.** Its headline anomaly — joint-ESM2 scoring
being catastrophically worse than the additive null on ParD, **Δ = −0.283** — is inflated by pooling across
mutation orders. The honest **within-order Δ is −0.053**.

The sweep's *conclusion* is unchanged and in fact strengthened: the additive null remains the robust
deployed default. What changes is the magnitude of the one exception, and the story told about it.

No new compute: this is re-analysis of the committed sweep JSON plus the cached ProteinGym assays.

## The arithmetic (established — not a hypothesis)

Per-order Spearman, n-weighted, vs the pooled figure the sweep reported:

| protein | metric | within-order | pooled | pooling gain |
|---|---|---:|---:|---:|
| GFP | additive | 0.018 | 0.101 | +0.082 |
| GFP | joint | 0.029 | 0.087 | +0.058 |
| HIS7 | additive | 0.273 | 0.342 | +0.069 |
| HIS7 | joint | 0.282 | 0.347 | +0.065 |
| **ParD** | **additive** | **0.259** | **0.543** | **+0.284** |
| **ParD** | joint | 0.206 | 0.260 | +0.054 |

| protein | within-order Δ | pooled Δ (published) |
|---|---:|---:|
| GFP | +0.010 | −0.014 |
| HIS7 | +0.009 | +0.005 |
| **ParD** | **−0.053** | **−0.283** |

**Two facts fall straight out, both pure arithmetic:**

1. **The published −0.283 overstates the within-order effect ~5x.** ParD's pooled additive ρ (0.543) is
   higher than its ρ at *any* individual order (0.301 / 0.356 / 0.120) — the signature of a confounder.
   Mutation order is that confounder.
2. **Joint's pooling gain is near-constant across all three proteins** (+0.054, +0.058, +0.065), while
   **additive's swings 4x** (+0.069 → +0.284). So the anomaly is not "joint collapses on ParD"; it is
   "*additive* gets an unusually large pooling bonus on ParD."

This is the same shape as the AMR side's clonality inflation: a grouping variable the predictor tracks
inflates the pooled metric. There the group was the clone; here it is the mutation order.

## What explains additive's outsized bonus — two hypotheses, one FALSIFIED

**H1 (mine, FALSIFIED).** The additive score is a SUM of k terms, so it scales with k by construction; it
should harvest cross-order signal wherever fitness declines steeply with k. Measured Spearman(k, fitness):

| protein | ρ(k, fitness) | additive pooling gain |
|---|---:|---:|
| GFP | **−0.476** | +0.082 |
| ParD | −0.410 | **+0.284** |
| HIS7 | −0.122 | +0.069 |

**GFP has the STEEPER slope and the far SMALLER gain.** H1 is dead as stated — the rank slope is not the
governing quantity.

**H2 — SINCE CONFIRMED (2026-08-25, `wiki/forward_epistasis_h2_confirmed_2026-08-25.md`):** tested WITHIN protein across order-subsets (which removes the protein-identity confound below) it holds at Spearman +0.973 / +0.995 on GFP / HIS7 over 26 subsets each, and two negative controls make it causal — shuffling ORDER LABELS drives the gain to exactly 0.000 on all three proteins, while shuffling FITNESS WITHIN order leaves it intact. The pooling gain is **entirely** between-order structure. Original underpowered reading follows:

**H2 (CONSISTENT, UNDERPOWERED — not established).** What should govern pooling inflation is not the slope
but how much of the total fitness variance sits *between* orders — η²(k):

| protein | η²(k) | additive pooling gain |
|---|---:|---:|
| ParD | 0.432 | +0.284 |
| GFP | 0.254 | +0.082 |
| HIS7 | 0.016 | +0.069 |

Monotone — but on **n=3 proteins**, where a monotone ordering arises by chance ~17% of the time. This is a
hypothesis consistent with three points, **not** a mechanism. Testing it needs more proteins with
multi-order coverage.

**H3 (CONFOUNDED, cannot be tested here).** Mutation density: every Δ < −0.05 occurs above 3% k/L, and
ParD (L=93) is the **only** protein in the panel above 3%. Density is therefore perfectly confounded with
protein identity — n=1 in the danger zone. The one-instance-is-not-a-mechanism trap; recorded, not claimed.

## Consequences

- **`predict_multi_effect`'s additive null stays the deployed default** — unchanged, and better supported:
  the worst measured within-order penalty for using additive is −0.053, not −0.283.
- **Report the WITHIN-ORDER Δ**, not the pooled one, whenever a sweep spans mutation orders. The pooled
  figure answers a different question ("can the score rank across orders?") and flatters any predictor
  that scales with k.
- **The "degrades with order" observation survives** the correction: ParD's per-order Δ still runs
  +0.054 → −0.086 → −0.126 across k=2,3,4. That trend is within-order and is not a pooling artifact.

## Reproduce

```bash
uv run python scripts/epistasis_pooling_check.py     # pure re-analysis; no GPU, no network
```
Sources: `wiki/forward_epistasis_sweep_2026-07-27.json` + `D:/dna_decode_cache/epistasis/*.csv`.
Frozen AMR surface untouched (this is the non-frozen forward cell).
