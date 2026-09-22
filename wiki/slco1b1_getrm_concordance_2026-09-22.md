# SLCO1B1 single-SNP surface vs GeT-RM CDC consensus (2026-09-22)

Scored **88** 1000G samples with CDC multi-lab consensus SLCO1B1 diplotypes. Bar frozen beforehand: `wiki/slco1b1_star_acceptance_bar.json`.

The cell's docstring says validation is *"never an independent star-concordance number"*. That is half right, and the two halves are reported separately below.

## 1. Star-NAME axis (genuinely testable -- the caller cannot see 388A>G)

- exact star-name agreement **35/88 = 0.3977**
- 53/88 samples carry at least one allele the caller structurally cannot identify (`*1A`/`*1B`/`*14`/`*15`/`*17`/`*21`)
- **verdict: STAR_PROXY_NOT_REFERENCE_AGREEING**

## 2. FUNCTION axis (a CONTROL, not validation)

GeT-RM's star assignment itself derives from genotypes including 521T>C, so this is close to self-comparison. A near-perfect number here is **expected and is not evidence**.

- reduced-function carriership: sens **1.0** · spec **1.0** · PPV **1.0** (n=88), dosage-exact 88/88

## 3. SAFETY axis (the clinical one)

- truth reduced-function carriers: **19**
- of those, called Normal Function anyway: **0**

## Which reference allele drives the naming gap

| allele | samples | name agrees |
|---|---|---|
| `*1` | 62 | 35 |
| `*1B` | 23 | 0 |
| `*15` | 15 | 0 |
| `*14` | 8 | 0 |
| `*1A` | 8 | 0 |
| `*17` | 4 | 0 |
| `*21` | 3 | 0 |
| `*5` | 1 | 0 |

## Honest limits

- GeT-RM consensus is caller-derived, not wet-lab -- agreement-with-reference, never correctness
- the FUNCTION axis shares its defining SNP with the reference's own allele assignment and is therefore a control, not validation
- n=88; *5 (1 haplotype), *17 (4) and *21 (3) are individually underpowered and must not carry a per-allele claim
- truth diplotypes that do not parse as two haplotypes are counted and EXCLUDED, never defaulted
