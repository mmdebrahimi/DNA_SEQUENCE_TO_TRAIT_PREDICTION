# FBA gets the gene right and the *condition* wrong (2026-08-12)

A new validation cell for the FBA decoder, and the first one that measures the property strain design
actually depends on.

Every essentiality number in this repo so far has been **single-condition** — one medium, one E/N call per
gene. That measures whether a model knows a gene matters *at all*. It cannot measure the thing a designer
relies on: that a gene is dispensable on glucose and required on succinate, so deleting it is free in one
process and lethal in another.

## The substrate — and the wall that wasn't there

`scripts/fba_carbon_growth_validate.py` recorded a 2026-08-03 probe concluding a two-sided E. coli
phenotype set was **"SI-locked"**. That was re-checked and it is **wrong**. The Europe PMC
`supplementaryFiles` REST endpoint serves both candidate papers' supplements freely and scriptably:

```
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3261703/supplementaryFiles
```

That is the route the earlier probe missed — it had tried publisher pages and the PMC web UI, which are
JS/Cloudflare-gated, exactly as the Kosuri fetch was. **Europe PMC's REST API is not.**

The gold standard is now committed at
`data/raw/ecoli_conditional_essentiality/orth2011_table_s1_conditional_essentiality.tsv`
(Orth et al. 2011, *Mol Syst Biol* 7:535, Supplementary Table 1 — CC-BY, redistributed with attribution):

**1,075 E. coli K-12 genes × 4 minimal media** — glucose aerobic, glucose anaerobic, lactate aerobic,
succinate aerobic — each an experimental E/N, with **68 genes conditionally essential** (essential in ≥1
medium, dispensable in ≥1 other). Two properties make it the right substrate:

1. **Two-sided by construction.** A conditionally-essential gene is its own control — same gene, same GPR,
   only the medium moves. A model cannot score well here by calling everything dispensable.
2. **It carries its own reproduction gate.** The supplement also ships the paper's own iJO1366 FBA calls,
   so the pipeline is checked against a published result before any new number is trusted.

## Result

**Reproduction gate passed** — the paper's own iJO1366 calls reproduce sane per-condition numbers, and the
recomputed conditionally-essential subset (68) matches the supplement's own flag exactly.

| | glucose aer | glucose anaer | lactate aer | succinate aer |
|---|---|---|---|---|
| **iJO1366** *(paper's own, the gate)* MCC | 0.6981 | 0.6692 | 0.6759 | 0.6648 |
| **iML1515** *(successor, ours)* MCC | **0.7428** | **0.7244** | **0.7013** | **0.7010** |

Per condition, iML1515 improves on its predecessor by ~+0.04 MCC across all four media. Then the
conditional metric:

| | exact-set match | per-cell agreement |
|---|---|---|
| null — always dispensable | **0/68** | **0.5588** |
| null — always essential | 0/68 | 0.4412 |
| iJO1366 *(paper's own)* | 4/68 = **5.9%** | 0.5735 |
| **iML1515** *(ours)* | 3/67 = **4.5%** | 0.5709 → **lift +0.0121** |

**Both models reproduce the conditional switch for about 1 gene in 20, and beat a constant predictor by
about one percentage point.**

### What that means

The per-condition MCC of 0.70–0.74 is real, but it is carried almost entirely by genes that are *always*
essential or *never* essential. On the 68 genes where the answer actually depends on the medium — the only
genes where a condition-aware model earns its keep — there is close to nothing.

**And eight years of model development did not move it.** iML1515 is the maintained successor to iJO1366;
it is better on every single-condition metric and *slightly worse* on the switch (4.5% vs 5.9%). Whatever
improved between the two reconstructions, it was not conditional resolution.

The null control is what makes this legible. A per-cell agreement of 0.57 reads like signal until you
notice that predicting "dispensable everywhere" scores 0.5588 on the same subset — because most
conditionally-essential genes are essential in only one or two of the four media. The module refuses to
report the metric without its null (`constant_baselines`), and two tests pin the values.

## Honest limits

- **In-distribution**: a published knowledge baseline, not an independent-lab claim.
- **Media are approximations** — the reconstruction's own M9 mineral background with the carbon source
  swapped and the oxygen bound set. They are not calibrated to the assay's exact medium.
- **Lactate is scored as L-lactate** (`EX_lac__L_e`); the assay may have used D,L-lactate.
- iML1515 is the *successor* of the model the paper scored, so the gap between the gate row and our row is
  a **model** difference, not a reproduction failure.
- 1,064 of 1,075 gold-standard genes are present in iML1515; 11 are not scored.

## Why this matters for the design epoch

Track A ships growth-coupled strain design — it proposes knockouts that make a product obligatory. Those
proposals are only as trustworthy as the model's ability to say **"this gene is required *here* and not
*there*."** This cell measures exactly that capability and finds it is ~5%.

That is not a reason to distrust the design cell's outputs wholesale — the designs it found were validated
against a published OptKnock-lineage result, and the growth-coupling test is a flux argument rather than a
per-gene essentiality call. It *is* the right place to point model improvement, and it is a far better
target than the single-condition metric that Track C's premise check showed is at its honest ceiling.

## Reproduce

```bash
uv run python scripts/fba_conditional_essentiality_validate.py
```

Sidecar: `wiki/fba_conditional_essentiality_ecoli_2026-08-12.json`.
Tests: `tests/test_fba_conditional_essentiality.py` (15, offline — the gold standard is committed).
