# Growth-coupled strain design — Succinate exchange (Escherichia coli K-12 MG1655, 2026-08-07)

Model **iML1515** (Escherichia coli K-12 MG1655) · target **`EX_succ_e`** · growth floor **10% of wild type** (0.0877 /h of 0.876997 /h).

**Growth coupling** means every flux distribution that grows at the floor also secretes the product — production is obligatory, not merely allowed. That is what makes an engineered strain stable: selection for growth becomes selection for production.

## Baseline (wild type)

- min product flux **0.0** · max **15.429274** → **POSSIBLE**

## Search

- exhaustive single knockouts over 1516 genes; pairs are a BOUNDED heuristic over the top 40 singles (780 pairs), NOT an exhaustive double-knockout search
- genes scanned **1516**; non-viable at the growth floor **196**; singles evaluated **1320**; pairs evaluated **780**

## Result — 0 growth-coupled design(s)

**No growth-coupled design found** under this search. The closest candidates (still uncoupled — the cell can avoid producing) were:

| knockouts | min product | max product | growth (/h) | coupling |
|---|---|---|---|---|
| `b4132` | 0.0 | 15.429274 | 0.876997 | POSSIBLE |
| `b3617` | 0.0 | 15.429274 | 0.876997 | POSSIBLE |
| `b0383` | 0.0 | 15.429274 | 0.876997 | POSSIBLE |
| `b1465` | 0.0 | 15.429274 | 0.876997 | POSSIBLE |
| `b2799` | 0.0 | 15.429274 | 0.876997 | POSSIBLE |

A negative result here is informative: it bounds what single/paired knockouts can do for this target on this medium, and points at medium or pathway changes instead.

## Scope

STOICHIOMETRIC prediction. A coupled design is a HYPOTHESIS FOR THE BENCH, not a validated strain: FBA does not model regulation, enzyme kinetics, toxicity, metabolic burden, or whether the knockout strain is constructible.
