# CYP2C19 caller validation -- pharmcat (2026-06-25)

**Honesty tier:** FAITHFUL_TO_PHARMCAT (reference tool's own test fixtures; in-distribution, NOT independent)

- Core diplotype concordance: **6/6** (1.0)
- Core phenotype concordance: **6/6** (1.0)
- Non-core blind-spot cases (excluded from headline): 2

## Core fixtures (alleles in *1/*2/*3/*17)

| fixture | expected | predicted | match | phenotype (pred) | phenotype match |
|---|---|---|---|---|---|
| s1s1 | *1/*1 | *1/*1 | OK | Normal Metabolizer | OK |
| s1s17 | *1/*17 | *1/*17 | OK | Rapid Metabolizer | OK |
| s1s1rs12248560missing | *1/*1 | *1/*1 | OK | Normal Metabolizer | OK |
| s1s2 | *1/*2 | *1/*2 | OK | Intermediate Metabolizer | OK |
| s2s2 | *2/*2 | *2/*2 | OK | Poor Metabolizer | OK |
| s2s3 | *2/*3 | *2/*3 | OK | Poor Metabolizer | OK |

## Non-core cases -- v0.1 sentinel layer WITHHOLDS rather than mis-calls (2/2 correctly withheld)

| fixture | expected (PharmCAT) | core-proxy | phenotype_status | sentinel |
|---|---|---|---|---|
| s1s35 | *1/*35 | *1/*1 | phenotype_withheld | *35 |
| s1s4b | *1/*4b | *1/*17 | phenotype_withheld | *4 |

_Core concordance is on the v0 SNP-defined set (*1/*2/*3/*17). Non-core fixtures are now WITHHELD by the v0.1 sentinel layer (phenotype_status=phenotype_withheld) rather than silently mis-called, and do NOT count toward the headline. The GeT-RM INDEPENDENT consensus number is a DATA-ACCESS step (labels in paper supplements; the 1000G tooling is proven via Docker bcftools -- see wiki/pgx_1000g_population_2026-06-25); this harness consumes it via --source getrm --expected-tsv. NOT a clinical tool._
