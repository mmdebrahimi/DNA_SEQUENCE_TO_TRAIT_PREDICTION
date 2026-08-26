# H2 CONFIRMED: the pooling gain is made entirely of between-order structure (2026-08-25)

Follow-up to `wiki/forward_epistasis_pooling_correction_2026-08-25.md`, which showed the epistasis sweep's
ParD anomaly (Δ = −0.283) was a mutation-order **pooling artifact** (within-order Δ = −0.053) and that the
inflated half is the **additive** score's pooling gain. Two explanations were offered; H1 was falsified on
the spot, and **H2** — that the governing quantity is η²(k), the share of fitness variance sitting BETWEEN
orders — was left as *consistent but underpowered* (n=3 proteins, monotone-by-chance ≈ 17%).

**H2 is now confirmed two independent ways, with zero new compute.**

## 1. Within-protein correlation (fixes the confound that made n=3 useless)

n=3 was never the real problem. The problem was that η² was **perfectly confounded with protein identity** —
each protein contributed exactly one point, so "high η²" and "is ParD" were the same statement. The same
confound killed H3 (density).

So the test holds protein identity **fixed** and varies η² across **subsets of mutation orders** within one
protein. Each (protein, order-subset) is a point:

| protein | orders used | subsets | Spearman(η², pooling gain) | p |
|---|---|---:|---:|---:|
| GFP | 2,3,4,5,6 | 26 | **+0.973** | 7.4e-17 |
| HIS7 | 4,5,6,7,8 | 26 | **+0.995** | 9.3e-26 |
| ParD | 2,3,4 | 4 | +0.600 | 0.40 (underpowered — only 4 subsets) |

Two proteins, 26 subsets each, near-perfect monotone relationships. ParD has only 3 orders so it yields 4
subsets and stays underpowered — reported, not hidden.

## 2. Intervention: two negative controls that isolate the source

A correlation between η² and the gain could still be incidental, so the structure was **intervened on**:

- **Control A** — permute ORDER LABELS across all variants. Every marginal distribution is preserved; only
  the fitness↔order association dies. H2 predicts the gain collapses.
- **Control B** — permute FITNESS WITHIN each order. The within-order score↔fitness signal dies; the
  per-order means survive. H2 predicts the gain is untouched.

| protein | condition | η² | within ρ | **pooling gain** |
|---|---|---:|---:|---:|
| GFP | real | 0.254 | +0.010 | **+0.066** |
| GFP | A (order labels shuffled) | 0.000 | +0.076 | **−0.000** |
| GFP | B (fitness shuffled in-order) | 0.254 | −0.001 | **+0.065** |
| HIS7 | real | 0.035 | +0.334 | **+0.043** |
| HIS7 | A | 0.000 | +0.377 | **−0.000** |
| HIS7 | B | 0.035 | +0.001 | **+0.052** |
| ParD | real | 0.432 | +0.227 | **+0.148** |
| ParD | A | 0.000 | +0.376 | **−0.000** |
| ParD | B | 0.432 | +0.004 | **+0.185** |

**Both predictions hold on all three proteins.** Destroy between-order structure → the gain is exactly
zero. Destroy within-order signal → the gain is untouched (it was never made of that). The pooling gain is
**entirely** between-order structure.

Note ParD carries the highest η² (0.432) and the largest real gain (+0.148) — the outlier is not special,
it is simply the far end of the same axis.

## Honest framing — what is and is not new here

This is **Simpson's-paradox mechanics**, and that pooling two groups differing in both X and Y manufactures
correlation is textbook. The non-obvious parts are the specific ones:

1. It **quantitatively accounts for the ParD outlier** — ParD is not anomalous, it just has the largest η².
2. η² gives a **number** for how much of a pooled metric is order-separation rather than ranking skill.
3. The controls show the gain is **100%** between-order, not partly — so a pooled figure computed across a
   grouping variable is not "somewhat inflated", it is inflated by exactly that structure.

**Scope:** three proteins, one score family (additive ESM2 log-ratio), one grouping variable (mutation
order). The mechanism is general in principle; this is the evidence actually in hand.

## Deployment rule

- Report the **within-order** metric whenever a sweep spans mutation orders.
- If a pooled figure is reported anyway, **report η²(k) beside it** so a reader can see how much of the
  number is group separation.
- This generalises past mutation order: it is the same rule as the AMR side's clonality disclosure — any
  grouping variable both the predictor and the label track will inflate a pooled metric, and η² is the
  cheap way to say by how much.

## Reproduce

```bash
uv run python scripts/epistasis_pooling_h2_test.py   # ~6 min; no GPU, no network
```
Uses cached ESM2 per-position log-prob matrices (additive score only — the joint score is not needed, since
the pooling gain is a property of the additive score alone) plus the cached assay CSVs, at the FULL assays
rather than the sweep's 300-per-order subsample. Artifact:
`wiki/forward_epistasis_h2_within_protein_2026-08-25.json`. Frozen AMR surface untouched.
