# dna-amr cross-source validation — NCBI Pathogen Detection — 2026-06-07

> Recommendation #3 (2026-06-06 brainstorm): validate the deterministic AMR caller against an
> INDEPENDENT label source (NCBI Pathogen Detection), zero accession overlap with BV-BRC cohorts.
> Headline question: does the per-drug rule add value over vanilla AMRFinder on UN-tuned data?

- Source: NCBI Pathogen Detection `PDG000000004.6094` (https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Escherichia_coli_Shigella/latest_snps/Metadata/PDG000000004.6094.metadata.tsv)
- Cohort: 22 E. coli, 22 with AMRFinder runs; non-overlap with BV-BRC enforced by construction
- AMRFinder image (pinned): `ncbi/amr:4.2.7-2026-03-24.1`

## dna-amr per-drug rule vs naive AMRFinder (independent labels)

| Drug | N | dna-amr acc | sens | spec | naive acc | sens | spec | Δacc |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ciprofloxacin | 22 | **0.955** | 1.0 | 0.917 | 0.864 | 1.0 | 0.75 | +0.091 |
| ceftriaxone | 22 | **0.864** | 1.0 | 0.727 | 0.636 | 1.0 | 0.273 | +0.228 |
| gentamicin | 22 | **1.0** | 1.0 | 1.0 | 0.591 | 1.0 | 0.25 | +0.409 |
| tetracycline | 22 | **0.909** | 1.0 | 0.8 | 0.909 | 1.0 | 0.8 | +0.0 |

## Discordance taxonomy (dna-amr failure modes vs independent labels)

| Drug | FN (R missed: efflux/porin/regulatory/low-level) | FP (called R, susceptible: label/expression/borderline) |
|---|---:|---:|
| ciprofloxacin | 0 | 1 |
| ceftriaxone | 0 | 3 |
| gentamicin | 0 | 0 |
| tetracycline | 0 | 2 |

## Interpretation

- Δacc > 0 means the per-drug policy (threshold + Subclass refinement) beats vanilla "any drug-class determinant -> R" on data the rule was NOT tuned on — the product-value headline.
- Honest scope: NCBI Pathogen Detection aggregates public AST submissions (different source/curation than BV-BRC, NOT a controlled different-lab study). Closes the same-database gap, not the same-methodology gap.
- Non-overlap with all BV-BRC cohorts is guaranteed by construction (selection excluded the 176 BV-BRC accessions).