# Is the inverse deployable on a protein with NO DMS? (2026-07-17)

**Yes — but only the RANK version. It works on 4/4 usable proteins; the learned oracle earns its keep on 3/4.**

## The problem with the passing magnitude result

the calibrator maps score -> measured effect, so it needs the TARGET protein's DMS -- and if you have that you already know every variant's effect. Cross-protein transfer is impossible by construction, not merely inaccurate: the assays do not share a scale (CcdB's entire range [-9.00,-2.00] lies below TEM-1's minimum -3.56, so a TEM-1 calibrator cannot express a CcdB value at all).

So the magnitude inverse — the one that PASSED at +53% on blaTEM — is confined to a narrow niche:
*I scanned half the positions, now design at the unscanned ones*. It cannot serve a novel protein.

## The deployable alternative

Ask for a **rank**, not a dose: *propose an edit at the p-th percentile of damage among the
reachable edits*. That needs **no calibrator, no DMS, no measured label** — only the oracle's own
score ordering. Graded by the percentile the proposal ACTUALLY lands at in the measured
distribution. Errors are in percentile points, so they are comparable across proteins with
incommensurate assay scales — the property the magnitude version structurally lacks.

Random-pick null (exact, no RNG): top-1 **0.3250**, best-of-5 **0.0922** percentile points.

| protein | organism | esm top-1 | esm best-of-5 | BLOSUM best-of-5 | vs null | vs BLOSUM | verdict |
|---|---|---:|---:|---:|---:|---:|---|
| TEM-1 beta-lactamase | E. coli | 0.1344 | 0.0375 | 0.0489 | +59.3% | +23.3% | works, *BLOSUM suffices* |
| PTEN | human | 0.1270 | 0.0405 | 0.0705 | +56.1% | +42.5% | **oracle earns keep** |
| RL40A (ubiquitin) | yeast | 0.2238 | 0.0600 | 0.1309 | +34.9% | +54.2% | **oracle earns keep** |
| SR43C | Arabidopsis | 0.1964 | 0.0520 | 0.1395 | +43.6% | +62.7% | **oracle earns keep** |
| CcdB toxin | E. coli | — | — | — | — | — | **EXCLUDED — censored** (79% tied at -2.0) |

## What this means for what ships

**Ship the rank inverse.** It needs nothing but the protein sequence, works on every usable assay
here, and its error is ~2-5 percentile points with 5 proposals. **Do not ship the magnitude
inverse**: it requires the very data that would make it unnecessary.

The two versions do NOT agree on which proteins they suit, which is why both were measured:

| | magnitude inverse | rank inverse |
|---|---|---|
| blaTEM | star (+52.9% vs BLOSUM) | *BLOSUM suffices* (+23.3%) |
| RL40A | **fails vs null** | works (+34.9%) |

A per-protein gate is therefore not optional — the right scorer, and even whether the learned
oracle is worth loading at all, changes protein by protein and question by question.

## Honest limits

- **It ranks, it does not dose.** *Near the top of the damaging tail* — never *fold-change 4.2*.
- **A censored assay is excluded, not scored.** CcdB (79.3% of variants tied at its ceiling) has no
  well-defined percentile; ungated it posted −159% vs null, which reads as an oracle failure and is
  not one — the metric is simply undefined there.
- **Regime B (molecular fitness) only.** Not clinical resistance, where this scorer class is below
  chance.
- **n=4 proteins.** Enough to falsify *utility tracks rank*; not enough to model what does predict it.
