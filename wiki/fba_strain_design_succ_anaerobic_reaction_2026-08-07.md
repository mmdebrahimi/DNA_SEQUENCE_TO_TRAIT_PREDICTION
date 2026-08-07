# Growth-coupled strain design — Succinate exchange (Escherichia coli K-12 MG1655, 2026-08-07)

Condition: **ANAEROBIC (O2 uptake closed)**.

Model **iML1515** (Escherichia coli K-12 MG1655) · target **`EX_succ_e`** · growth floor **90% of each strain's OWN max growth** (wild type max 0.15754 /h) · knockout level **reaction**.

**Growth coupling** means every flux distribution that grows at the floor also secretes the product — production is obligatory, not merely allowed. That is what makes an engineered strain stable: selection for growth becomes selection for production.

## Baseline (wild type)

- min product flux **0.047306** · max **3.675306** → **OBLIGATORY**

> **The wild type is ALREADY growth-coupled here** (guaranteed 0.047306), so `OBLIGATORY` alone is not evidence of a design. Every design below must **beat** that floor; **3604** knockout sets were coupled only by inheriting it and are NOT counted.

## Search

- exhaustive single reaction knockouts over 2266 candidates; pairs are a BOUNDED heuristic over the top 40 singles: 780 pairs, 816 triples. NOT exhaustive at depth>1
- reaction candidates scanned **2266**; non-viable at the growth floor **258**; singles evaluated **2008**; pairs **780**; triples **816**

## Result — 0 growth-coupled design(s)

**No growth-coupled design found** under this search. The closest candidates (still uncoupled — the cell can avoid producing) were:

| knockouts | min product | max product | growth (/h) | coupling |
|---|---|---|---|---|
| `ETHAtex` | 0.047306 | 3.675306 | 0.15754 | OBLIGATORY |
| `R15BPK` | 0.047306 | 3.675306 | 0.15754 | OBLIGATORY |
| `TRDR` | 0.047306 | 3.675306 | 0.15754 | OBLIGATORY |
| `NNAM` | 0.047306 | 3.675306 | 0.15754 | OBLIGATORY |
| `PGMT` | 0.047306 | 3.675306 | 0.15754 | OBLIGATORY |

A negative result here is informative: it bounds what single/paired knockouts can do for this target on this medium, and points at medium or pathway changes instead.

## Scope

STOICHIOMETRIC prediction. A coupled design is a HYPOTHESIS FOR THE BENCH, not a validated strain: FBA does not model regulation, enzyme kinetics, toxicity, metabolic burden, or whether the knockout strain is constructible.
