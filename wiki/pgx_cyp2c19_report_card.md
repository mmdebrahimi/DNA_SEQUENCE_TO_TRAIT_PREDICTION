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

## Non-core blind-spot cases (caller covers core SNP set only -- mis-call EXPECTED)

| fixture | expected (PharmCAT) | predicted (ours) | flags |
|---|---|---|---|
| s1s35 | *1/*35 | *1/*1 |  |
| s1s4b | *1/*4b | *1/*17 |  |

_Core concordance is on the v0 SNP-defined set (*1/*2/*3/*17). Non-core fixtures are reported as blind-spot cases (the caller is EXPECTED to mis-call a non-core star to a *1-substituted diplotype + flag it) and do NOT count toward the headline. The GeT-RM INDEPENDENT number (1000 Genomes VCFs) needs tabix/bcftools + a fetch on a Linux/Docker host -- the named follow-up; this harness consumes that cohort identically via --expected-tsv. NOT a clinical tool._
