# Single-locus falsification cells — photic sneeze + asparagus anosmia (2026-07-08)

Two new deterministic single-SNP cells built to TEST single-locus sufficiency against free openSNP
self-report labels. Both are moderate-effect GWAS hits (not Mendelian), so the honest expectation is
near-chance — the same calibration contrast the existing **cilantro** cell provides. Directions are
**sourced, not fabricated** (WebSearch, 2026-07-08).

## Cells

| trait | SNP | gene | risk/direction (sourced) | source |
|---|---|---|---|---|
| photic sneeze | rs10427255 | ZEB2-region | **C** allele → sneezer (OR~1.32) | [Eriksson 2010 PLoS Genet 6:e1000993](https://pmc.ncbi.nlm.nih.gov/articles/PMC2891811/) |
| asparagus anosmia | rs4481887 | OR2M7-region | **A** allele → smeller, **GG** → anosmic | [Eriksson 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2891811/) + [Pelchat 2011 Chem Senses 36:9](https://academic.oup.com/chemse/article/36/1/9/442551) |

## Result — vs the majority-class null (the load-bearing comparison)

| trait | n | accuracy | majority-null | verdict | sens / spec |
|---|---|---|---|---|---|
| photic | 108 | 0.602 | 0.574 | **NEAR_CHANCE** | 0.85 / 0.42 |
| asparagus | 61 | 0.574 | 0.852 | **BELOW_NULL** | 0.67 / 0.56 |

**Both FALSIFY single-locus sufficiency for these traits** (the point of the exercise):
- **Photic** clears the null by only +2.8 pp — the OR-1.32 C allele is common, so the dominant-risk rule
  over-calls "sneezer" (sens 0.85 / spec 0.42) and barely beats "predict everyone non-sneezer."
- **Asparagus** lands *below* the null: the anosmia allele (GG) is ~42% frequent, but only ~15% self-report
  anosmia — so the GG→anosmic rule massively over-predicts anosmia. Flipping the direction is worse (~0.43),
  confirming the sourced direction is correct but the single locus simply doesn't determine the trait here.

This is the intended contrast: **earwax (rs17822931) is strong-Mendelian and clears the null; photic +
asparagus + cilantro are weak GWAS hits that do not.** A single associated SNP ≠ a deterministic decoder.

## Integrity note (circular-label guard)

The openSNP self-reports for both traits contained **genotype-referencing values** ("photic sneezer with the
snp", "gg - but... i can smell it", "rs10427255", bare "cc") — these are contaminated by the very genotype
they would validate. `single_snp_traits._mentions_genotype` **drops** them (including leading bare-genotype
tokens) so the labels stay independent. Without this guard the concordance would be spuriously inflated.

## Honest tier

- Self-reported labels (near-independent, non-circular after the guard, noisy) — PILOT/DEMO, NOT a lab assay.
- Moderate-effect GWAS SNPs — **near-chance / below-null is the CORRECT outcome**, not a bug.
- NOT a clinical tool.

Reproduce: `uv run python scripts/single_snp_opensnp_validate.py --trait photic` (or `asparagus`). Rules +
binners: `dna_decode/data/single_snp_traits.py`. Tests: `tests/test_single_snp_photic_asparagus.py` (6).
