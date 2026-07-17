# Does the molecular inverse generalize? — cross-protein boundary map (2026-07-17)

**Inverse design WORKS on 4/5 proteins. The LEARNED oracle earns its keep on 3/5. Magnitude is certifiable on 2/5.**

The blaTEM falsifier passed at +53%. That was **one protein**. This runs the identical falsifier
across every protein with a cached ESM2 table — E. coli x2, human, yeast, Arabidopsis.

## Result

| protein | organism | fwd rank | esm err/span | vs null | win | vs BLOSUM | win | interval | verdict |
|---|---|---:|---:|---:|:-:|---:|:-:|:-:|---|
| TEM-1 beta-lactamase | E. coli | 0.73 | 0.0307 | +58.9% | 6/6 | +52.9% | 6/6 | 0/6 | **oracle earns keep** |
| CcdB toxin | E. coli | 0.51 | 0.0128 | +77.1% | 6/6 | +12.5% | 2/6 | 0/6 | works, *BLOSUM suffices* |
| PTEN | human | 0.52 | 0.0262 | +49.3% | 6/6 | +33.1% | 6/6 | 6/6 | **oracle earns keep** |
| RL40A (ubiquitin) | yeast | — | 0.0408 | +41.1% | 5/6 | +26.9% | 3/6 | 0/6 | **FAIL vs null** |
| SR43C | Arabidopsis | — | 0.0355 | +50.0% | 6/6 | +57.4% | 6/6 | 4/6 | **oracle earns keep** |

*err/span normalizes by each protein's measured-effect span — spans differ >10x (RL40A 0.75 vs
PTEN 8.57), so absolute errors are not comparable across rows. `win` = paired per-split wins;
all three methods face the identical partition on each split.*

## Three findings

### 1. Inverse SELECTION generalizes; the LEARNED oracle does NOT

4/5 proteins beat the no-oracle null — the inverse lands ~2x
closer to a target than guessing (err/span 1.3–4.1% vs the null's 5.2–7.5%). But only
3/5 beat **BLOSUM62**. On CcdB the 1992 substitution
matrix ties ESM2-650M (2/6 paired wins, +12.5%); on RL40A it is a coin flip (3/6).
**The blaTEM PASS was not representative.** Where BLOSUM suffices that is an *engineering win*, not
a failure: no GPU, no model, no precomputed table — just run the matrix.

### 2. Inverse utility does NOT track forward rank quality (the pre-registered question, answered)

The natural assumption is that a better forward ranker inverts better. It is **false**:

| protein | forward Spearman | earns its keep vs BLOSUM? |
|---|---:|---|
| PTEN | **0.518** | **yes — 6/6 paired wins, +33.1%** |
| CcdB | **0.5115** | **no — 2/6 paired wins, +12.5%** |

Near-identical rank quality, opposite verdicts. So you cannot read a leaderboard Spearman and
conclude the oracle will be useful for design on that protein — it must be measured per protein.
*(n=5: this falsifies the assumption; it does not establish what the real predictor is.)*

### 3. Selection quality and magnitude certifiability are ORTHOGONAL — measured in both directions

| protein | selection vs BLOSUM | informative interval |
|---|---|---|
| blaTEM | **best in the sweep** (+52.9%) | **0/6 — certifies nothing** |
| PTEN | +33.1% | **6/6 — fully certified** |

This project already knew *ranks well ≠ pins the dose* (CcdB). The sweep shows the **converse** too:
*pins the dose ≠ selects well*. They are independent axes, so a design tool must report both — and
an interval that merely **brackets** the target proves nothing, because split-conformal coverage
holds even for a useless model.

## Scope

only blaTEM has a committed CDS -> only it is restricted to the single-nt-accessible (real genome-edit) space; the rest use the full DMS variant set (protein-level). Within a protein every method faces the SAME pool, so each verdict is valid; cross-protein margins must carry `candidate_space`.

Test: PAIRED per-split comparison (esm2/blosum62/null face the IDENTICAL partition on each split): a win must sweep every paired split AND clear the 25% material margin. An unpaired range test would discard the pairing and is far too conservative.

## What this licenses

A **per-protein-gated** selection inverse: run this falsifier on the target protein FIRST, and let
it choose the scorer (BLOSUM where BLOSUM suffices — most proteins — and ESM only where it earns
its keep). Do **not** ship a blanket 'ESM2 inverse' on the strength of blaTEM; on 2/5 proteins here
that would burn a GPU to match a substitution matrix, and on RL40A it would not beat guessing.
