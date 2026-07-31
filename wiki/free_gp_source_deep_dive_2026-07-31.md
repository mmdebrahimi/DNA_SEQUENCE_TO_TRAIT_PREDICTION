# Deep-dive: a FREE paired genotype+phenotype source for a real decoding-validation arm (2026-07-31)

**Goal:** find a free, paired (genotype+phenotype) source that survives this project's screening rails,
then build a real decoding-validation arm on it. The binding constraint has always been LABELS not
models: every prior candidate hit one of the 8 rejection gates (`wiki/negative_results_map_2026-06-13.md`)
or failed the embedding-niche three-part test (sampling-independent lab label + no curated catalog +
organism depth ≥100), and the population-structure confound killed 4 learned-decoder attempts
(Arabidopsis flowering-time, cipro within-lineage, pathotype, HIV-ESM).

## Candidates surveyed

| source | depth | phenotype | verdict |
|---|---|---|---|
| **Bloom et al. 2013 yeast cross** | 1,008 segregants | 46 quantitative growth traits (colony size) | **RECOMMENDED — confound-free** |
| 1011 natural yeast isolates (Peter 2018) | 1,011 strains × 223 traits | growth/life-history | viable but HARD (population structure = Arabidopsis-like confound) |
| Yeast proteome/transcriptome QTL (~1000 isolates) | ~889 isolates | protein/mRNA abundance | natural isolates → same structure confound |
| Bacterial AMR / TB (BacGWASim etc.) | thousands | binary resistance | SATURATED / done (10 provdisjoint cells + TB) |
| ProteinGym DMS | 217 assays | protein fitness | molecular (the `forward` cell already uses it — not organism phenotype) |

## Why Bloom 2013 is the breakthrough substrate

The 1,008 segregants come from a SINGLE genetic cross (lab strain BY × wine strain RM → recombinant
offspring). This design is the one thing every prior candidate lacked: **it removes population structure
by construction.** Every segregant is a ~50/50 recombinant mosaic of the same two parental genomes, and
meiotic recombination randomizes which parent contributes each locus — so there is NO clade / ancestry /
sampling stratification for the genotype to secretly track. Any cross-validated genotype→phenotype signal
is therefore a GENUINE causal link, not the confound that inflated (and then sank) the Arabidopsis and
cipro-within-lineage attempts.

### Screen against the 8 rejection gates — Bloom PASSES all
1. **circular-label** — NO: colony size is a direct lab measurement, not a tool output the model competes with.
2. **study==class** — NO: one study; the phenotype varies WITHIN the panel by genotype (that is the point).
3. **sampling-defined** — NO: growth is measured, not a sampling attribute.
4. **surveillance-domination** — NO: 1,008 unique recombinants, no dominant clone.
5. **assembly-attrition** — N/A: genotype is 11,623 markers (parental origin), no assembly needed.
6. **MIC-censoring** — NO: continuous colony size, not censored MIC.
7. **provenance-not-separable** — N/A: provenance is the two known parents; cross-validation is over segregants.
8. **dedup-collapses-balance** — NO: unique segregants, balanced by design.

### Three-part embedding-niche test — Bloom PASSES all
- **sampling-independent lab label:** ✅ colony size under controlled media.
- **no curated catalog:** ✅ quantitative yeast growth QTL has no simple determinant catalog (unlike AMR/PGx).
- **organism depth ≥100:** ✅ 1,008 same-species recombinants.

**Plus the population-structure confound — SIDESTEPPED** (the cross design), which no prior candidate did.

## What decoder this substrate admits (honest framing)

Bloom's genotype is MARKER-based (parental origin, 1/−1), and there is NO curated determinant catalog for
these growth traits. So:
- The **deterministic curated-catalog decoder** (AMR/PGx style) does NOT apply here (no catalog).
- A **genomic-prediction / learned decoder** DOES apply and is the natural fit — and, uniquely, this is a
  FAIR test of it because the confounds are removed.

Bloom 2013 itself showed the detected loci explain "nearly the entire additive contribution to heritable
variation" → genomic prediction demonstrably WORKS on this data. So the pilot arm is expected to produce a
genuine POSITIVE (cross-validated r² > 0) — the project's FIRST clean genotype→phenotype positive after a
string of confounded negatives. That is the value: it flips the project from "0-for-N confounded" to "here
is a working decoder + a baseline on a confound-free substrate."

## The decoding-validation arm (design)

- **Layer 1 (pilot, this build):** cross-validated genomic prediction — ridge/BLUP on the marker genotype
  → colony-size phenotype, per trait, report held-out predictive r² (the standard genomic-prediction
  baseline). Establishes the arm + a baseline. Highly-heritable traits (e.g. metal/oxidative stress) should
  give clear r².
- **Layer 2 (research payoff, follow-on):** does a DNA foundation-model embedding of each segregant's actual
  SEQUENCE (parental-allele-substituted) BEAT the marker-ridge baseline on this clean substrate? This is the
  genuinely novel test — whether the "AI" adds value once confounds are gone — and needs reconstructing
  per-segregant sequence from parental alleles (heavier; deferred).

## Recommendation
Build the Bloom-2013 genomic-prediction arm (Layer 1) as the project's first confound-free
decoding-validation arm. 1011 natural isolates is the harder natural-population follow-on (needs the same
within-clade de-confounding that Arabidopsis required). Data is free (no DUA); download + format confirmed
before the build.

Sources: Bloom et al. 2013 Nature 494:234-7 (doi:10.1038/nature11867); Peter et al. 2018 Nature 556:339
(1011 genomes); consolidated yeast trait resource Mol Syst Biol 2025 (doi:10.1038/s44320-025-00136-y).
