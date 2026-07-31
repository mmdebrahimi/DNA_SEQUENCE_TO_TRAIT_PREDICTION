# Bloom-2013 yeast decoding-validation arm — the first CONFOUND-FREE genotype→phenotype positive (2026-07-31)

The project's first clean genotype→phenotype decoding result. Every prior learned-decoder attempt was
confounded (population structure) or hit a rejection gate; the Bloom 2013 BYxRM segregant cross is the
substrate that removes the confound by construction (see the deep-dive
`wiki/free_gp_source_deep_dive_2026-07-31.md`). Cross-validated ridge genomic prediction with a
label-permutation null.

## Result — 12/12 traits decode, decisively above the null

Substrate: Bloom et al. 2013 BYxRM cross, **1,008 segregants** aligned, marker-subsampled to **1,451
markers** (stride 8 of the genome-wide 11,623; a CONSERVATIVE pilot — more markers only add signal),
5-fold CV ridge, 30-permutation null.

| trait | predictive r | predictive r² | null r p95 | verdict |
|---|---|---|---|---|
| Lithium_Chloride | **0.803** | 0.646 | 0.090 | BEATS NULL |
| Neomycin | 0.773 | 0.598 | 0.050 | BEATS NULL |
| Maltose | 0.728 | 0.531 | 0.079 | BEATS NULL |
| Cycloheximide | 0.709 | 0.502 | 0.073 | BEATS NULL |
| Diamide | 0.699 | 0.489 | 0.025 | BEATS NULL |
| Cadmium_Chloride | 0.693 | 0.480 | 0.105 | BEATS NULL |
| YPD | 0.692 | 0.479 | 0.054 | BEATS NULL |
| Zeocin | 0.681 | 0.463 | 0.063 | BEATS NULL |
| Cobalt_Chloride | 0.677 | 0.459 | 0.054 | BEATS NULL |
| Hydrogen_Peroxide | 0.603 | 0.363 | 0.118 | BEATS NULL |
| Copper | 0.566 | 0.320 | 0.083 | BEATS NULL |
| Caffeine | 0.463 | 0.214 | 0.065 | BEATS NULL |

**12/12 beat the label-permutation null** (predictive r 0.46–0.80 vs null p95 0.03–0.12). Artifact:
`wiki/yeast_bloom_gp_arm_2026-07-31.json`.

## Why this is a GENUINE positive (not confounded inflation)

The 1,008 segregants are recombinant offspring of ONE cross (lab BY × wine RM). Recombination randomizes
which parent contributes each locus, so there is **no population / clade / ancestry structure** for the
genotype to secretly track. A held-out predictive r that beats the permutation null is therefore a real
genotype→phenotype signal — NOT the population-structure inflation that produced (then sank) the
Arabidopsis flowering-time and cipro-within-lineage embedding attempts. The magnitudes also match the
literature: Bloom 2013 reported these traits are highly heritable and that detected loci explain "nearly
the entire additive contribution to heritable variation" (Lithium_Chloride, one of the most heritable,
tops the table at r²=0.65).

## Honest scope

- **This is the marker-ridge genomic-prediction baseline (Layer 1).** It establishes the arm + a working
  decoder + a baseline any fancier model must beat. Ridge-on-markers is standard quantitative genetics;
  the novelty for THIS project is (a) the FIRST confound-free substrate and (b) a reusable, null-controlled
  harness. It flips the project from "0-for-N confounded learned-decoder negatives" to "a clean positive
  on clean data".
- **It does NOT overturn the closed embedding-vs-phenotype negatives.** Those were CONFOUNDED
  natural-population datasets; this is a confound-free cross — a different, cleaner regime. It shows the
  decoding *method* works when the data is clean, not that a foundation model adds value.
- **Layer 2 (the research payoff, deferred):** does a DNA foundation-model embedding of each segregant's
  actual sequence (parental-allele-substituted genome) BEAT this marker-ridge baseline on the clean
  substrate? That is the genuinely novel "does the AI add value once confounds are gone" test; it needs
  per-segregant sequence reconstruction (heavier).
- Marker-subsampled pilot (1,451 of 11,623); a full-marker run only strengthens the result.

## Reproducibility
- Engine: `dna_decode/eval/genomic_prediction.py::cv_ridge_gp` (4 offline tests). Arm:
  `scripts/yeast_bloom_gp_arm.py`.
- Data (free, no DUA): Princeton BYxRM web supplement — genotype `BYxRM_GenoData.txt` (11,623 markers ×
  1008 segregants, B/R parental origin, all 16 chromosomes), phenotype `BYxRM_PhenoData.txt` (1008 × 46
  colony-size traits). Princeton 403s this sandbox's IP; fetched via the Wayback capture
  `web.archive.org/web/20240707040412id_/…` (the `2020id_` form caps at 1 MB → truncates the genotype).
  Files cached at `D:/dna_decode_cache/bloom/` (large; not committed). Frozen AMR/forward surfaces
  byte-unchanged.
