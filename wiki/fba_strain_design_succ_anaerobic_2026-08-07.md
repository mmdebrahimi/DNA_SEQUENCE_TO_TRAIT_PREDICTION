# Growth-coupled strain design — Succinate exchange (Escherichia coli K-12 MG1655, 2026-08-07)

Condition: **ANAEROBIC (O2 uptake closed)**.

Model **iML1515** (Escherichia coli K-12 MG1655) · target **`EX_succ_e`** · growth floor **10% of wild type** (0.015754 /h of 0.15754 /h).

**Growth coupling** means every flux distribution that grows at the floor also secretes the product — production is obligatory, not merely allowed. That is what makes an engineered strain stable: selection for growth becomes selection for production.

## Baseline (wild type)

- min product flux **0.005256** · max **16.26593** → **OBLIGATORY**

> **The wild type is ALREADY growth-coupled here** (guaranteed 0.005256), so `OBLIGATORY` alone is not evidence of a design. Every design below must **beat** that floor; **2096** knockout sets were coupled only by inheriting it and are NOT counted.

## Search

- exhaustive single knockouts over 1516 genes; pairs are a BOUNDED heuristic over the top 40 singles (780 pairs), NOT an exhaustive double-knockout search
- genes scanned **1516**; non-viable at the growth floor **200**; singles evaluated **1316**; pairs evaluated **780**

## Result — 0 growth-coupled design(s)

**No growth-coupled design found** under this search. The closest candidates (still uncoupled — the cell can avoid producing) were:

| knockouts | min product | max product | growth (/h) | coupling |
|---|---|---|---|---|
| `b3290` | 0.005256 | 16.26593 | 0.15754 | OBLIGATORY |
| `b2094` | 0.005256 | 16.26593 | 0.15754 | OBLIGATORY |
| `b2497` | 0.005256 | 16.26593 | 0.15754 | OBLIGATORY |
| `b3093` | 0.005256 | 16.26593 | 0.15754 | OBLIGATORY |
| `b1488` | 0.005256 | 16.26593 | 0.15754 | OBLIGATORY |

A negative result here is informative: it bounds what single/paired knockouts can do for this target on this medium, and points at medium or pathway changes instead.

## Scope

STOICHIOMETRIC prediction. A coupled design is a HYPOTHESIS FOR THE BENCH, not a validated strain: FBA does not model regulation, enzyme kinetics, toxicity, metabolic burden, or whether the knockout strain is constructible.
