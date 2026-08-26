# H2: the pooling gain is made of between-order structure (2026-08-25)

> **AMENDED 2026-08-25 — this memo OVERCLAIMED; read
> `wiki/forward_epistasis_eta2_identifiability_2026-08-25.md` first.** Four corrections, all load-bearing:
> **(1)** The governing quantity is **NOT** η² of fitness. A pooling gain needs BOTH the label and the
> predictor to separate by group, and the additive score's own between-order η² was never measured. Measured,
> the two are near-collinear (ρ +0.982 / +0.999) and the two well-powered proteins point in OPPOSITE
> directions when partialled (GFP score-side +0.815 / fitness +0.020; HIS7 fitness +0.559 / score **−0.135**).
> The supportable claim is the JOINT one — aligned between-order separation in label AND predictor — and which
> side carries it is **not separably identifiable** here. **(2)** The p-values below are **not inferential**:
> the 26 subsets are nested and share most of their variants, so n_eff ≈ the number of ORDERS (5). Treat them
> as descriptive. The honest resampling unit is leave-one-order-out, and the relationship IS robust to it
> (ρ stays +0.955..+1.000). **(3)** "confirmed", "causal", "entirely", "100%" are **overstated** — Control A
> re-deals the same group sizes from one shuffled pool, so gain → 0 *by construction*; it checks the
> arithmetic, not a mechanism. **(4)** Control B did **not** leave the gain "untouched" — it rose on 2 of 3
> (HIS7 +0.043→+0.052, ParD +0.148→+0.185), which is explainable but was asserted past.
>
> The original text is preserved below unedited except for these markers. What still stands: **the pooling
> gain is between-order structure, and the ParD anomaly is a pooling artifact** (`wiki/forward_epistasis_pooling_correction_2026-08-25.md`).

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
| GFP | 2,3,4,5,6 | 26 | **+0.973** | 7.4e-17 *(NOT inferential — see banner)* |
| HIS7 | 4,5,6,7,8 | 26 | **+0.995** | 9.3e-26 *(NOT inferential — see banner)* |
| ParD | 2,3,4 | 4 | +0.600 | 0.40 (underpowered — only 4 subsets) |

> These p-values assume the 26 subsets are independent. They are not — they are nested and share most of
> their variants. n_eff ≈ 5 orders. Use the leave-one-order-out sweep in
> `wiki/forward_epistasis_eta2_identifiability_2026-08-25.md` instead.

Two proteins, 26 subsets each, near-perfect monotone relationships. ParD has only 3 orders so it yields 4
subsets and stays underpowered — reported, not hidden.

## 2. Intervention: two negative controls that isolate the source

A correlation between η² and the gain could still be incidental, so the structure was **intervened on**:

- **Control A** — permute ORDER LABELS across all variants. Every marginal distribution is preserved; only
  the fitness↔order association dies. H2 predicts the gain collapses.
- **Control B** — permute FITNESS WITHIN each order. The within-order score↔fitness signal dies; the
  per-order means survive. H2 predicts the gain is untouched. *(AMENDED: it was not — the gain ROSE on 2 of
  3, HIS7 +21% and ParD +25%. Explainable, since gain = pooled − within and B zeroes a positive within term
  that was diluting the difference; but the prediction as stated did not hold exactly.)*

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
zero. Destroy within-order signal → the gain survives. The pooling gain is dominated by between-order
structure.

> **AMENDED.** Two walk-backs. **Control A is near-definitional** — it re-deals the SAME group sizes from one
> shuffled pool, so every group becomes an i.i.d. draw from the same population, pooled and within estimate
> the same quantity, and the gain → 0 *by construction*. What A establishes is that the arithmetic is
> implemented correctly: a pipeline check, not a causal intervention. **Control B did NOT leave the gain
> "untouched"** — it rose on 2 of 3 (HIS7 +0.043→+0.052 = +21%, ParD +0.148→+0.185 = +25%). The direction is
> explainable and should have been explained: B zeroes the within-order ρ (HIS7 +0.334→+0.001, ParD
> +0.227→+0.004), and since gain = pooled − within, removing a positive within term that was *diluting* the
> difference raises the gain. "entirely" and "100%" below are withdrawn.

Note ParD carries the highest η² (0.432) and the largest real gain (+0.148) — the outlier is not special,
it is simply the far end of the same axis.

## Honest framing — what is and is not new here

This is **Simpson's-paradox mechanics**, and that pooling two groups differing in both X and Y manufactures
correlation is textbook. The non-obvious parts are the specific ones:

1. It **quantitatively accounts for the ParD outlier** — ParD is not anomalous, it just has the largest η².
2. η² gives a **number** for how much of a pooled metric is order-separation rather than ranking skill.
   *(AMENDED: η²(label) ALONE is under-specified — the predictor's own η² matters equally and is not
   separable from it here. Disclose both.)*
3. ~~The controls show the gain is **100%** between-order~~ — **WITHDRAWN**, see the control amendment above.
   The controls show the gain is dominated by between-order structure; Control A cannot establish the
   proportion because its collapse is definitional.

**Scope:** three proteins, one score family (additive ESM2 log-ratio), one grouping variable (mutation
order). The mechanism is general in principle; this is the evidence actually in hand.

## Deployment rule

- Report the **within-order** metric whenever a sweep spans mutation orders. *(AMENDED: state the estimand
  first — this is right for "which scorer ranks at fixed mutation burden?", which IS the deployed question,
  but not for "can this score rank a realistic mixed-order pool?", where cross-order ranking is the target.)*
- If a pooled figure is reported anyway, **report η²(k) beside it** so a reader can see how much of the
  number is group separation. *(AMENDED: disclose η² of the **label AND the predictor** — a pooled metric is
  inflated only where the two are aligned, and one alone cannot tell you.)*
- This generalises past mutation order: it is the same rule as the AMR side's clonality disclosure — any
  grouping variable both the predictor and the label track will inflate a pooled metric. *(AMENDED: η² is
  the cheap way to FLAG exposure, not to "say by how much". `dna_decode/eval/clonality.py` computes a
  corrected estimator with uncertainty — one lineage vote, `wilson_ci`, `effective_lineage_n`, DISCORDANT
  clusters excluded rather than majority-voted. η² is a scalar flag with no CI and no corrected estimate;
  the corrected number is the within-group metric itself.)*

## Reproduce

```bash
uv run python scripts/epistasis_pooling_h2_test.py   # ~6 min; no GPU, no network
```
Uses cached ESM2 per-position log-prob matrices (additive score only — the joint score is not needed, since
the pooling gain is a property of the additive score alone) plus the cached assay CSVs, at the FULL assays
rather than the sweep's 300-per-order subsample. Artifact:
`wiki/forward_epistasis_h2_within_protein_2026-08-25.json`. Frozen AMR surface untouched.
