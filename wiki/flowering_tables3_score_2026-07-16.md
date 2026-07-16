# Arabidopsis flowering cell — scored on Zhang 2020 Table S3 (N=854)

*Generated 2026-07-16 by `scripts/flowering_tables3_score.py`. Verdict: **SCORED_BEATS_NULL**.*

## Headline

| metric | value |
|---|---|
| pooled accuracy | **0.733** |
| best constant predictor (null) | 0.502 (always `early`) |
| sensitivity (late) | 0.967 |
| specificity (early) | 0.501 |
| **structure-group-weighted accuracy** | **0.710** (null 0.676) |
| groups beating their own null | 7/9 |

Split: `FT16 > cohort median` = **70.25 days**. Observed 425 late / 429 early.

## Who is missing (non-random dropout)

163 of 1017 accessions have FT16_mean='NA' and cannot be scored. The dropout is NOT random wrt genotype: 9.8% of dropped accessions carry a deleterious FRI vs 24.1% overall -- i.e. FUNCTIONAL-FRI (late-candidate) accessions are preferentially unphenotyped. This shifts the class balance the null baseline is computed against.

## Why sensitivity and specificity diverge (it is not the label)

| direction | n | precision | mean FT16 |
|---|---|---|---|
| `FRI loss-of-function -> early` | 229 | **0.939** | 49.9d |
| `FRI functional -> late` | 625 | **0.658** | 83.5d |

The rule is strong in one direction and weak in the other. **Losing FRI reliably makes a plant early**
(93.9%); **having FRI does not make it late** (65.8%). Functional FRI is *necessary but not sufficient* — which is precisely what the cell's two-locus rule says (you also need a strong FLC), and precisely
what this FRI-only run cannot use.

A near-perfect sensitivity beside a coin-flip specificity is normally the signature of a **wrong
label surrogate** in this project. It is not that here, and the false-positive spread is what rules
it out: the 214 functional-FRI-yet-early accessions span 32.0–70.2 days with a median of 59.2d, far below the 70.2d cut. A thresholding/label artifact would pile them up
*at* the boundary. These are real early plants that really do carry a functional FRI — the FLC route
and the polygenic residue, exactly the mechanisms the cell names as beyond its reach.

## The population-structure correction (read this before the pooled number)

The paper reports non-functional FRI is overrepresented in central/western Europe and rare elsewhere
(their Figure S3, P<0.01). So FRI genotype correlates with ancestry, and the pooled accuracy partly
measures *'can you recognise a central European accession?'* — the plant analogue of the clonality
inflation corrected elsewhere in this project. Per-STRUCTURE-group accuracy, each against **its own**
null:

| STRUCTURE group | n | accuracy | its null | beats null? | observed late/early |
|---|---|---|---|---|---|
| central_europe | 157 | 0.682 | 0.720 | no | 44/113 |
| south_sweden | 132 | 0.848 | 0.826 | **yes** | 109/23 |
| admixed | 112 | 0.696 | 0.616 | **yes** | 43/69 |
| spain | 94 | 0.553 | 0.511 | **yes** | 48/46 |
| western_europe | 90 | 0.789 | 0.767 | **yes** | 21/69 |
| germany | 79 | 0.823 | 0.722 | **yes** | 22/57 |
| italy_balkan_caucasus | 66 | 0.712 | 0.712 | no | 47/19 |
| asia | 55 | 0.564 | 0.527 | **yes** | 29/26 |
| north_sweden | 47 | — | — | *unscorable (one class only)* | 47/0 |
| relict | 22 | 0.727 | 0.682 | **yes** | 15/7 |

## Scope limits (load-bearing)

- FRI-ROUTE ONLY: Table S3 has no FLC status, so the cell's two-locus AND collapses to its FRI route -- the weaker, MEDIUM-confidence one. The FLC route (the cell's distinctive claim, the Da(1)-12 class) is NOT tested by this run.
- FLC assumed functional for every accession. Documented-rare, not absent; every FLC-route accession present is a guaranteed miss.
- HABIT/direction call, not quantitative days-to-flower. The paper reports FRI/FLC explains only part of long-day variation; the residue is polygenic + environmental.
- 163 of 1017 accessions have FT16_mean='NA' and cannot be scored. The dropout is NOT random wrt genotype: 9.8% of dropped accessions carry a deleterious FRI vs 24.1% overall -- i.e. FUNCTIONAL-FRI (late-candidate) accessions are preferentially unphenotyped. This shifts the class balance the null baseline is computed against.
- IN-DISTRIBUTION, NOT INDEPENDENT: the cell's catalogue and this label both trace to the same literature. This measures faithfulness of the rule to the paper's own data -- it is NOT an out-of-distribution validation.

## Threshold sensitivity

Re-scored at the upper tertile (86.25 days): accuracy 0.580 vs null 0.670 → **FAILS_NULL_BASELINE**. The median cut is used for the headline because it makes the constant-predictor baseline maximally hard (~50/50); the tertile row shows whether the verdict rides on that choice.

## Provenance

- Substrate: `data/arabidopsis/zhang2020/tpj14716-sup-0012-TableS3.tsv` — Zhang L & Jimenez-Gomez JM (2020) Functional analysis of FRIGIDA using naturally occurring variation in Arabidopsis thaliana. The Plant Journal 103:154-165. doi:10.1111/tpj.14716. Table S3, CC-BY 4.0.
- Cell reference-integrity guard: `True`
- Phenotype: FT16_mean (days to first flower, long days 16C) -- the paper's own phenotype
