# Growth-coupled strain design — Succinate exchange (Escherichia coli K-12 MG1655, 2026-08-07)

Condition: **ANAEROBIC (O2 uptake closed)**.

Model **iML1515** (Escherichia coli K-12 MG1655) · target **`EX_succ_e`** · growth floor **90% of each strain's OWN max growth** (wild type max 0.15754 /h) · knockout level **reaction**.

**Growth coupling** means every flux distribution that grows at the floor also secretes the product — production is obligatory, not merely allowed. That is what makes an engineered strain stable: selection for growth becomes selection for production.

## Baseline (wild type)

- min product flux **0.047306** · max **3.675306** → **OBLIGATORY**

> **The wild type is ALREADY growth-coupled here** (guaranteed 0.047306), so `OBLIGATORY` alone is not evidence of a design. Every design below must **beat** that floor; **24** knockout sets were coupled only by inheriting it and are NOT counted.

## Search

- exhaustive single reaction knockouts over 5 candidates; pairs are a BOUNDED heuristic over the top 5 singles: 10 pairs, 10 triples. NOT exhaustive at depth>1
- reaction candidates scanned **5**; non-viable at the growth floor **0**; singles evaluated **5**; pairs **10**; triples **10**

## Result — 1 growth-coupled design(s)

| knockouts | guaranteed product flux | gain over wild type | max product flux | growth (/h) |
|---|---|---|---|---|
| `LDH_D`, `ALCD2x`, `PFL` | **9.263835** | +9.216529 | 13.158454 | 0.081588 |

## Scope

STOICHIOMETRIC prediction. A coupled design is a HYPOTHESIS FOR THE BENCH, not a validated strain: FBA does not model regulation, enzyme kinetics, toxicity, metabolic burden, or whether the knockout strain is constructible.
