# DGRP learned-path test — the 5th de-confounded negative (with a new substrate lesson)

**Date:** 2026-07-02
**Script:** `scripts/dgrp_learned_decoder.py` (`wiki/dgrp_learned_scores.json`)
**Data (free):** DGRP2 genotypes (NCSU dm6 SNPs, via the reliable Aertslab mirror `resources.aertslab.org/DGRP2/`
— the canonical NC State server was down) + DGRPool phenotype `StarvationRes` (id 2798, JSON API). On D:.

## The test

Does a **multivariate genotype→trait model beat population structure** on a *Drosophila* quantitative-trait
panel, once de-confounded? This is the animal-panel instance of the project's recurring embedding/learned-path
test — negative on yeast growth, Arabidopsis flowering, cipro within-lineage, and pathotype. DGRP was the one
untouched high-VOI model-organism candidate. Method mirrors `yeast_growth_decoder.py`, reusing the promoted
`dna_decode.deconfound` primitives: build a line×SNP matrix, derive structure, compare a naive whole-genome
ridge (`cv_r2`) against the de-confounded `within_group_r2` (5-fold inside each structure group, scored on
group-centered residuals), with a within-group permutation null.

## Result — negative, and DGRP is a *degenerate substrate* for the structure test

| metric | value |
|---|---|
| N lines (genotype × StarvationRes) | 197 |
| SNPs (MAF≥0.05, thinned) | 12,365 |
| naive whole-genome ridge `cv_r2` | **0.0154** |
| de-confounded `within_clade_r2` (PCA-ward, 3 strata) | **0.0349** |
| permutation-null p95 (abs corr) | 0.14 (→ r² floor ~0.02) |
| **verdict** | **LEARNED_PATH_NEGATIVE_UNDER_DECONFOUNDING** |

Both the naive (0.015) and de-confounded (0.035) r² sit at the permutation-noise floor. The learned path shows
**no usable genotype→trait signal** on DGRP StarvationRes.

### The new lesson: DGRP has no discrete population structure

Euclidean average-linkage clustering (the method that gave yeast/Arabidopsis their discriminating clades)
collapses DGRP into **one clade of 191/197 lines** (`discrete_structure_degenerate: true`). DGRP is a single
Raleigh-derived population — it lacks the discrete structure that makes the clade-de-confounding test sharp.
So the test had to use **PCA-continuous structure** (top-5 genotype PCs → ward, 3 usable strata) instead. The
mechanism of the negative is therefore *distinct* from the earlier cases:

- **Yeast / Arabidopsis:** naive r² was HIGH, then **collapsed** under de-confounding (the model had learned
  structure). Classic signal-vs-structure negative.
- **DGRP:** there is **no naive signal to inflate** (0.015) **and no discrete structure to de-confound
  against** — a "no-signal + wrong-substrate-for-the-test" negative.

## Honest caveats (do not overclaim)

- **Not a claim that DGRP StarvationRes is unpredictable in principle.** StarvationRes is heritable; proper
  DGRP genomic prediction uses **all ~2–4M SNPs + a GRM/GBLUP**, not a 12k-thinned ridge at α=10, N=197. The
  naive 0.015 is a statement about *this thinned-SNP multivariate model*, whose job was to establish a signal
  the de-confounding could then attack. It didn't clear the bar — which is the point of the test — but a
  full-genome GBLUP would post a higher naive r².
- **Single trait, single thinning.** A fuller sweep (multiple traits, all SNPs, tuned α) is deferred; it would
  sharpen the *naive* number but is very unlikely to flip the *de-confounded* verdict, given the 0-for-4 prior
  and DGRP's lack of discrete structure.

## Where this leaves the learned track

**The learned/embedding path is now 0-for-5 under de-confounding** (yeast, Arabidopsis, cipro within-lineage,
pathotype, DGRP). This was the expected, informative outcome — it re-confirms that the project's validated
product is the **deterministic determinant-scan decoder** (AMR/TB/fungal/HIV/SARS + the human PGx lane + the
new C. elegans ben-1 cell), not a learned genotype→trait model. DGRP adds the substrate lesson: a panel without
discrete population structure is degenerate for the clade-de-confounding test — reach for PCA-continuous
structure there. See `wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md` (this is the 5th data point).

## Reproduce

```bash
uv run python scripts/dgrp_learned_decoder.py --build   # streams the VCF -> SNP matrix cache, then runs
uv run pytest tests/test_dgrp_learned_decoder.py -q      # 3 offline synthetic tests
# genotype mirror: resources.aertslab.org/DGRP2/NCSU/final/dm6/DGRP2.source_NCSU.dm6.final.SNPs_only.vcf.gz
# phenotype: https://dgrpool.epfl.ch/phenotypes/2798.json  (StarvationRes)
```
