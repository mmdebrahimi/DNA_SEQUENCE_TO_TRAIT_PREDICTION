# BV-BRC MIC 4-drug feasibility census — 2026-05-18

Phase 2 entry: counts per-drug feasibility of building N=150 cohorts from BV-BRC AST data at two label-quality bars.

## Methodology

Strains gated through 6 stages, computed at TWO label-quality bars:

- **Strict-MIC** = HIGH_R + HIGH_S only (median MIC >=4x CLSI-R or <=CLSI-S/4; 4x safety margin)
- **Relaxed-MIC** = strict + DECISIVE_R + DECISIVE_S (median MIC outside borderline gray zone but no 4x margin)

Strains pass through:
1. E. coli broth-microdilution AST rows for the drug
2. Distinct genome IDs
3. With parseable MIC value (numeric, not NA)
4. Strict-MIC OR relaxed-MIC classification
5. With assembly_accession (NCBI-downloadable)
6. Passing assembly QC (contig_count <= 500, n50 >= 50000)

Breakpoints: CLSI 2024 + EUCAST 14.0 (E. coli).

**Why two bars:** the audit framework (mechanism x MIC x opacity merge with SUSPEND gate) is the downstream cohort-quality gate. Strict-MIC pre-filters on a 4x safety margin upstream — biology-aware noise reduction, but doing the audit framework's job twice. Relaxed-MIC accepts any MIC clearly outside the gray zone and lets the audit framework flag noisy strains downstream. Both numbers reported so the user can pick the bar that matches the project's failure-tolerance.

## Per-drug census

### Ciprofloxacin

Breakpoints: CLSI R>=2.0, S<=0.5 / EUCAST R>=1.0, S<=0.25

| Stage | Strict-MIC | Relaxed-MIC |
|---|---:|---:|
| 1. Total E. coli broth-MIC AST rows | 1,793 | 1,793 |
| 2. Distinct genome IDs | 1,750 | 1,750 |
| 3. With parseable MIC value | 1,287 | 1,287 |
| 4. MIC classification pass | 513 | 518 |
| 5. With assembly_accession | 21 | 21 |
| 6. Passing assembly QC | 21 | 21 |

Tier distribution (all distinct genome IDs, stage 2):

| Tier | Count |
|---|---:|
| BORDERLINE | 717 |
| NO_MIC | 463 |
| HIGH_S | 396 |
| HIGH_R | 117 |
| AMBIGUOUS | 52 |
| DECISIVE_S | 5 |

**Strict-MIC final:** R=17, S=4, Total=21
**Relaxed-MIC final:** R=17, S=4, Total=21

**Verdict:** **Strict-MIC: INFEASIBLE at any N>=60** (R=17, S=4)
**Verdict:** **Relaxed-MIC: INFEASIBLE at any N>=60** (R=17, S=4)

### Ceftriaxone

Breakpoints: CLSI R>=4.0, S<=1.0 / EUCAST R>=2.0, S<=1.0

| Stage | Strict-MIC | Relaxed-MIC |
|---|---:|---:|
| 1. Total E. coli broth-MIC AST rows | 383 | 383 |
| 2. Distinct genome IDs | 354 | 354 |
| 3. With parseable MIC value | 264 | 264 |
| 4. MIC classification pass | 164 | 164 |
| 5. With assembly_accession | 69 | 69 |
| 6. Passing assembly QC | 68 | 68 |

Tier distribution (all distinct genome IDs, stage 2):

| Tier | Count |
|---|---:|
| HIGH_R | 161 |
| BORDERLINE | 97 |
| NO_MIC | 90 |
| AMBIGUOUS | 3 |
| HIGH_S | 3 |

**Strict-MIC final:** R=66, S=2, Total=68
**Relaxed-MIC final:** R=66, S=2, Total=68

**Verdict:** **Strict-MIC: INFEASIBLE at any N>=60** (R=66, S=2)
**Verdict:** **Relaxed-MIC: INFEASIBLE at any N>=60** (R=66, S=2)

### Tetracycline

Breakpoints: CLSI R>=16.0, S<=4.0 / EUCAST: no breakpoints (ECOFF only)

| Stage | Strict-MIC | Relaxed-MIC |
|---|---:|---:|
| 1. Total E. coli broth-MIC AST rows | 357 | 357 |
| 2. Distinct genome IDs | 322 | 322 |
| 3. With parseable MIC value | 154 | 154 |
| 4. MIC classification pass | 31 | 31 |
| 5. With assembly_accession | 2 | 2 |
| 6. Passing assembly QC | 1 | 1 |

Tier distribution (all distinct genome IDs, stage 2):

| Tier | Count |
|---|---:|
| NO_MIC | 168 |
| BORDERLINE | 123 |
| HIGH_S | 17 |
| HIGH_R | 14 |

**Strict-MIC final:** R=1, S=0, Total=1
**Relaxed-MIC final:** R=1, S=0, Total=1

**Verdict:** **Strict-MIC: INFEASIBLE at any N>=60** (R=1, S=0)
**Verdict:** **Relaxed-MIC: INFEASIBLE at any N>=60** (R=1, S=0)

### Gentamicin

Breakpoints: CLSI R>=16.0, S<=4.0 / EUCAST R>=4.0, S<=2.0

| Stage | Strict-MIC | Relaxed-MIC |
|---|---:|---:|
| 1. Total E. coli broth-MIC AST rows | 1,394 | 1,394 |
| 2. Distinct genome IDs | 1,363 | 1,363 |
| 3. With parseable MIC value | 806 | 806 |
| 4. MIC classification pass | 647 | 649 |
| 5. With assembly_accession | 136 | 136 |
| 6. Passing assembly QC | 134 | 134 |

Tier distribution (all distinct genome IDs, stage 2):

| Tier | Count |
|---|---:|
| HIGH_S | 629 |
| NO_MIC | 557 |
| BORDERLINE | 146 |
| HIGH_R | 18 |
| AMBIGUOUS | 11 |
| DECISIVE_S | 2 |

**Strict-MIC final:** R=2, S=132, Total=134
**Relaxed-MIC final:** R=2, S=132, Total=134

**Verdict:** **Strict-MIC: INFEASIBLE at any N>=60** (R=2, S=132)
**Verdict:** **Relaxed-MIC: INFEASIBLE at any N>=60** (R=2, S=132)

## Cross-drug summary

### Strict-MIC (HIGH_R + HIGH_S only; 4x safety margin)

| Drug | R | S | Total | N=150 per-class? | N=100 per-class? | N=60 smoke? |
|---|---:|---:|---:|---|---|---|
| Ciprofloxacin | 17 | 4 | 21 | no | no | no |
| Ceftriaxone | 66 | 2 | 68 | no | no | no |
| Tetracycline | 1 | 0 | 1 | no | no | no |
| Gentamicin | 2 | 132 | 134 | no | no | no |

### Relaxed-MIC (HIGH + DECISIVE; lean on audit framework for noise)

| Drug | R | S | Total | N=150 per-class? | N=100 per-class? | N=60 smoke? |
|---|---:|---:|---:|---|---|---|
| Ciprofloxacin | 17 | 4 | 21 | no | no | no |
| Ceftriaxone | 66 | 2 | 68 | no | no | no |
| Tetracycline | 1 | 0 | 1 | no | no | no |
| Gentamicin | 2 | 132 | 134 | no | no | no |

## What this informs

- **Stage 2 cohort building**: which drugs have enough strict-MIC strains for N=150 per-class cohorts.
- **Drug-parameterize the audit infrastructure** ([[Candidate-1-framing]]): which 4th-mechanism-class drug (gentamicin = aminoglycoside) clears the smoke threshold for the architectural-finding falsifier.
- **Compute budget allocation**: only drugs that pass N=150 per-class warrant Databricks NT cache populates.

## What this does NOT do

- Does not count MLST coverage on the feasible subset (separate diagnostic at `scripts/diagnose_bvbrc_mlst_gaps.py`)
- Does not check NCBI Datasets API availability for each accession (downloadable in practice may be lower than `assembly_accession`-present)
- Does not run AMRFinder mechanism audit on the feasible subset (separate step; only worthwhile post-cohort-build)
