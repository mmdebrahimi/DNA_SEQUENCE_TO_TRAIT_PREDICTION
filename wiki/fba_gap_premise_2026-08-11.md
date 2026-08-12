# Track C's premise, tested before building it (2026-08-11)

Track C of the design epoch proposes protein-function prediction over dark-matter genes → candidate
missing reactions → better FBA accuracy. That is a large build, so its premise was tested first, with a
falsifier written down before the run.

> **The premise.** FBA's **false negatives** — genes essential in vivo that the model calls dispensable —
> are missed because the reactions they control are **blocked or touch a dead-end metabolite**, i.e. the
> model is incomplete. The alternative is that they are missed for reasons gap-filling cannot reach:
> regulation, kinetics, moonlighting, or medium.

**Pre-registered falsifier:** FN gap-adjacency must sit materially above the TN rate among genes the model
calls dispensable, Fisher two-sided p < 0.05. Equal rates falsify it.

## Result — SUPPORTED

Yeast / **iMM904**, the only cross-organism essentiality cell currently SCORED. Gold standard: SGD
(1,215 essential genes). 905 model genes scored.

| cell | n | gap-adjacent | rate |
|---|---|---|---|
| TP *(essential, called essential)* | 43 | 6 | 14.0% |
| FP | 67 | 6 | 9.0% |
| **FN** *(essential, called dispensable)* | **92** | **65** | **70.7%** |
| TN *(dispensable, called dispensable)* | 703 | 307 | 43.7% |

**FN 70.7% vs TN 43.7%, Fisher two-sided p = 1.2 × 10⁻⁶.** Among the genes FBA calls dispensable, the ones
it is *wrong* about are far more likely to sit next to a hole in the model. Gap structure predicts **which**
essential genes the model misses.

## The cheaper lever was checked first, and it does not dissolve the result

43.9% of iMM904's reactions are blocked — but blocked *under the model's default medium*. Medium
mis-specification is a config fix; a structural gap needs the whole Track C build. These are separable by
opening every exchange:

| | reactions | share |
|---|---|---|
| blocked, default medium | 692 / 1577 | 43.9% |
| **blocked with every exchange open** — structural | **553 / 1577** | **35.1%** |
| unblocked purely by opening the medium | 139 | **20.1% of the blocked set** |

**~80% of the blocked set is structural.** The cheap lever is real but minor — worth doing first (setting
the organism's standard medium is nearly free and would recover ~20%), and it does not explain the
enrichment away.

## What this does and does not establish

**Does:** the premise Track C rests on is not wishful. The errors are concentrated exactly where the model
is incomplete, and that incompleteness is mostly structural rather than a medium artifact.

**Does not:**

- **This is correlation, not repair.** It shows gaps *mark* the false negatives; it does not show that
  filling a specific gap flips a specific gene. The decisive test remains a measured **MCC delta after
  gap-filling** — that is the actual Track C deliverable, and it is untouched by this probe.
- **The mechanism is near-tautological in one direction.** A gene whose reactions are all blocked carries
  zero flux, so knocking it out cannot reduce growth and it *must* be called dispensable. What is
  informative is that this holds **more** for truly-essential genes (70.7%) than for truly-dispensable ones
  (43.7%) — the tautology alone would not produce that gap.
- **Gap-adjacency is coarse** — "controls ≥1 blocked or dead-end-touching reaction". A sharper version
  would weight by how much of the gene's flux capacity is affected.
- **One organism.** P. aeruginosa is `MODEL_WALLED` (no GEM in BiGG) and S. aureus is `LABEL_WALLED`
  (iYS854 uses `USA300HOU_####` ids, NTML is JE2 `SAUSA300_####` — a crosswalk away), so yeast is the only
  cell that can currently carry this test.

## Recommended sequencing for Track C

1. **Set yeast's standard medium and re-score** — nearly free, recovers ~20% of the blocked set, and
   sharpens the baseline every later number is measured against.
2. **Then** the function-prediction build, measured by MCC delta on this same gold standard.

Doing (2) before (1) would attribute to gap-filling whatever (1) would have fixed for nothing.

## Reproduce

```bash
uv run python scripts/fba_gap_premise_check.py --organism yeast
# exit 0 = premise supported, 1 = falsified, 2 = not scorable
```

Sidecar: `wiki/fba_gap_premise_yeast_2026-08-11.json`. Tests: `tests/test_fba_gap_premise.py`.
