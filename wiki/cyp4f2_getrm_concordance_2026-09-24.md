# CYP4F2 single-SNP surface vs GeT-RM CDC consensus (2026-09-24)

Scored **64** 1000G samples with CDC multi-lab consensus CYP4F2 diplotypes. Bar and predictions frozen beforehand: `wiki/cyp4f2_star_acceptance_bar.json`.

This cell carries the SAME *"never an independent star-concordance number"* claim as SLCO1B1, plus an **AF-corroboration** claim cited as a validation leg. Both are tested.

## 1. Star-NAME axis

- exact star-name agreement **52/64 = 0.8125**
- 12/64 samples carry `*2`, which this caller cannot see (it reads exactly one SNP, rs2108622)
- **verdict: STAR_PROXY_PARTIAL** · registered prediction (~0.81, 12 disagreements) **HELD**

## 2. FUNCTION axis (a CONTROL, not validation)

- `*3` carriership: sens **1.0** · spec **1.0** · PPV **1.0** (n=64), dosage-exact 64/64
- definitional: GeT-RM's `*3` assignment derives from the one SNP we read.

## 3. SAFETY axis

- truth `*3` carriers: **22**; called Normal Function anyway: **0**
- *Scope:* *2's FUNCTIONAL status is NOT settled by this run and must not be asserted from memory. CPIC's warfarin dosing algorithm uses *3; whether *2 carries an independent dose effect is a separate literature question. So a *2 carrier called *1 is a NAMING finding, and is only a SAFETY finding if *2 is reduced-function -- which this measurement does not establish.

## 4. AF-corroboration axis (the second shipped claim)

Computed over the **full 3202-sample panel**, not the 64-sample overlap.

| pop | asserted | measured | delta | within +/-0.10 |
|---|---|---|---|---|
| EUR | 0.29 | 0.2773 | -0.0127 | yes |
| EAS | 0.79 | 0.2128 | -0.5772 | **NO** |

- **verdict: AF_CLAIM_WRONG**

All measured superpopulation frequencies:

| pop | n | ALT alleles | freq |
|---|---|---|---|
| AFR | 893 | 140 | 0.0784 |
| AMR | 490 | 242 | 0.2469 |
| EAS | 585 | 249 | 0.2128 |
| EUR | 633 | 351 | 0.2773 |
| SAS | 601 | 502 | 0.4176 |
| ALL | 3202 | 1484 | 0.2317 |

## Which reference allele drives the naming gap

| allele | samples | name agrees |
|---|---|---|
| `*1` | 54 | 49 |
| `*3` | 22 | 15 |
| `*2` | 12 | 0 |

## Scope limits

- *2's FUNCTIONAL status is NOT settled by this run and must not be asserted from memory. CPIC's warfarin dosing algorithm uses *3; whether *2 carries an independent dose effect is a separate literature question. So a *2 carrier called *1 is a NAMING finding, and is only a SAFETY finding if *2 is reduced-function -- which this measurement does not establish.
- the rsID defining *2 is NOT asserted here; what is load-bearing and verifiable is that this caller genotypes exactly ONE SNP (rs2108622) and *2 is not it
- GeT-RM consensus is caller-derived, not wet-lab -- agreement-with-reference, never correctness
- n=64; *2 appears on 12 samples and *3/*3 on 3 -- per-allele claims below those counts are underpowered
- AF is computed over the FULL 3202-sample panel (not the 64-sample truth overlap), which is the right denominator for a population-frequency claim and a DIFFERENT denominator from every other number here
- truth diplotypes that do not parse as two haplotypes are counted and EXCLUDED, never defaulted
