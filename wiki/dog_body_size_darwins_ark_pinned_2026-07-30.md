# Dog body-size catalog — pinned + functionally validated on Darwin's Ark (2026-07-30)

**Verdict: the dog body-SIZE cell is genuinely buildable on the free imputed-SNV substrate — the exact
opposite of the coat-colour outcome.** All four dominant body-size loci (IGF1/HMGA2/STC2/GHR) were pinned
via the OMIA→canFam4-liftover→.bim pipeline, found in-panel at their exact lifted positions, and
**functionally validated** against owner-reported height (Q121 z-score, N=3276). The combined 4-locus
additive polygenic score correlates with height at **r=+0.619 (R²=0.383)** — ~38% of cross-breed height
variance from 4 SNPs, squarely in the literature's "~half from a handful of loci" range.

## Why this substrate works where coat colour failed

Coat colour failed because its causal variants are **indels/SVs/imputation gaps** (CBD103 3bp-del, ASIP
SINE, MLPH frameshift) — structurally absent from a biallelic-SNV panel. Body-size causal variants are
**SNPs**, present in the 29M-SNV imputed panel. Same free Dryad cohort, opposite feasibility — because the
determinant *variant type* differs, not the trait's tractability. (Reusable lesson: pick the substrate to
the variant TYPE, not just the trait.)

## Pinned catalog (verified against the REAL .bim, not from memory)

Pipeline per locus: OMIA/literature canFam3.1 coord → UCSC canFam3ToCanFam4 liftover (pyliftover, 0-based)
→ exact `.bim` match (canFam4, IDs `chr:pos:ref:alt`) → **functional validation**: big-allele dosage vs
Q121 height z, N=3276, sign as the literature predicts.

| locus | gene | canFam3.1 source | canFam4 panel SNP | big allele | single-SNP r(dose,height) |
|---|---|---|---|---|---|
| **IGF1** | insulin-like growth factor 1 | chr15:41,221,438 (Sutter07/Rimbault13) | `chr15:41513523:G:A` | G | **+0.505** |
| **HMGA2** | high-mobility-group AT-hook 2 | chr10:8,348,804 (Rimbault13) | `chr10:8703415:G:A` | G | **+0.542** |
| **STC2** | stanniocalcin 2 | gene chr4:39,151,951–39,165,514 (NCBI 489112) | `chr4:40070215:T:A` | T | +0.369 |
| **GHR** | growth-hormone receptor | gene chr4:66,705,544–66,845,096 (NCBI 403721) | `chr4:67710295:C:T` | C | +0.299 |

Notes on the two CFA4 loci: the CanMap *tag-SNP* positions (39.2M / 67.0M canFam3.1) sit **outside** the
STC2/GHR gene bodies, so the first peak-centered scan under-read them; re-scanning the actual **gene
windows** (lifted from the NCBI gene spans) recovered the real signal. The STC2 hit sits ~17 kb downstream
of the gene body — matching the literature "SNP ~20 kb downstream of STC2" (Rimbault 2013). The IGF1 SINE
insertion itself is absent (structural, as expected); its intron-2 SNP is in complete LD and tags it.

Every dosage→height relationship is a clean monotonic dose-response, e.g. HMGA2: AA −1.19 → AG −0.65 →
GG +0.32; IGF1: AA −0.57 → AG +0.11 → GG +0.67.

## Combined polygenic score (the capstone)

Equal-weight sum of big-allele dosages across the 4 independent loci (0–8) vs Q121 height, N=3276:

- **r = +0.619, R² = 0.383.** Monotonic across the whole range:

| score | mean height-z | n |
|---|---|---|
| 0 | −1.15 | 32 |
| 2 | −1.13 | 192 |
| 4 | −0.62 | 317 |
| 5 | −0.09 | 462 |
| 6 | +0.14 | 721 |
| 8 | +0.77 | 500 |

## Honest scope

- This is a **RELATIVE size rank** (a polygenic score), NOT a calibrated absolute-height predictor. Q121 is
  a covariate-adjusted quantile-normalised z-score, so the validated claim is "more big-alleles → taller
  RANK", not "20 inches". Absolute height needs a raw-inches label + breed covariate (a v0.1 item).
- v0 = the 4 dominant loci. IGF1R/SMAD2 (secondary) and FGF4 (leg-length chondrodysplasia retrogene, a
  structural insertion) are OUT; the module ABSTAINS on FGF4-affected leg length rather than guessing.
- Faithful-to-literature: applies published loci + measured directions; not a new GWAS.

## Shipped artifact

- `dna_decode/pigment/dog_body_size.py` — the pinned+validated catalog (`SIZE_LOCI`), the additive
  `polygenic_size_score(dosages)` scorer, the relative-rank map, and `reference_integrity_ok()` (offline
  catalog+rule guard). `RULES_VERSION = dog-body-size-v0.1.0`.
- `tests/test_dog_body_size.py` — 12 offline tests (synthetic dosages; no genotype data / no network).
- Reproduction: liftover + `.bim` grep + functional scan use `dna_decode/pigment/plink_io.py` on the
  Darwin's Ark `DarwinsDogs_2024_N-3277_canfam4_gp-0.70_biallelic` set + the Q121 height TSV.

Frozen AMR/forward surfaces byte-unchanged (this cell imports nothing from them). This completes step 1
(pin the catalog coords) + most of step 2 (validate vs Q121 height). Step 3 — shipping the full
`typing:dog:morphology` cell (a `predict()` / CLI / cell-registry contract, plus the 4 binary/ordinal
morphology traits Q124/Q127/Q128/Q245) — is the natural next increment.
