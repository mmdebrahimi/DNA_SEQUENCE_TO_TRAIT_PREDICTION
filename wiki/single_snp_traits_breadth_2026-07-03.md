# Visible/organismal-trait breadth — 3 new single-locus openSNP cells (2026-07-03)

Extends the deterministic gene→trait decoder (eye colour rs12913832, ABO rs8176719) with **three new
single-SNP human-trait cells**, reusing the archived openSNP dump + the eye-colour cell's zip-streaming
machinery. The outcome is a **calibrated spread, not three clean wins** — which is the honest, informative
result: the deterministic decoder wins where the biology + label support it, is label-limited where they
don't, and correctly FAILS on a weak-effect association. Modules: `dna_decode/data/single_snp_traits.py`
(sourced rules) + `scripts/single_snp_opensnp_validate.py` (generalized scorer).

## Result (openSNP self-report label, PILOT/DEMO tier, deterministic strand-agnostic rules)

| Cell | rsid / gene | tier | n scored | acc | pos sens | neg spec | verdict |
|---|---|---|---|---|---|---|---|
| **Lactase persistence** | rs4988235 / MCM6-LCT | strong Mendelian | 339 | **0.779** | 0.527 (intolerant) | 0.874 (tolerant) | **PILOT WIN** — rule tracks the self-report |
| **Earwax type** | rs17822931 / ABCC11 | strong Mendelian | 247 | 0.765 | 0.136 (dry) | **0.994** (wet) | **LABEL-LIMITED** — rule correct, label noisy |
| **Cilantro soap** | rs72921001 / OR6A2 | weak-association *contrast* | 51 | 0.235 | 0.75 (soapy) | 0.191 | **CORRECTLY FAILS** (by design) |

## The calibrated reading (maps onto the decoder-regime boundary)

- **Lactase = a real deterministic pilot win.** The −13910*T persistence allele is common in the
  European-dominated openSNP cohort, so BOTH classes are powered; the recessive-non-persistence rule tracks
  self-reported lactose intolerance at acc 0.78. Textbook rule, non-circular self-report — a clean regime-1
  demonstration (deployed-rule integration, as always flagged; not new biology).
- **Earwax = the rule is right, the LABEL is the weak link.** Specificity 0.994 (the rule correctly
  identifies the rare truly-dry `TT`/`AA` homozygotes; the ~4% genetic-dry rate matches the European
  ABCC11 dry-allele frequency). But 27% of users *self-report* dry — far above the ~1-4% genetic rate — so
  the excess is non-genetic dryness (cleaning/age/subjective) → sensitivity collapses to 0.136. This is the
  project's **"high-spec/low-sens → suspect the label"** pattern at the finest grain: the GENOTYPE is the
  trustworthy output, the self-report can't score it in this ancestry.
- **Cilantro = the deterministic decoder CORRECTLY fails.** rs72921001 is a *weak GWAS association*, not a
  Mendelian rule; the risk allele is common, so a deterministic "carrier → soapy" call over-calls massively
  (FP 38 vs TP 3) → acc 0.235, below chance. Included ON PURPOSE as a calibration contrast: it demonstrates
  the decoder distinguishes a strong single-locus trait from a weak association rather than "winning"
  everywhere. A deterministic rule is the wrong tool for a weak-effect locus — the honest boundary.

## What it settles

- The deterministic single-locus decoder **generalizes to more organismal traits** — but the result is
  honestly REGIME-DEPENDENT: it wins on a strong-single-locus trait with a clean, common-allele,
  well-powered label (lactase); it is label-limited when the causal allele is rare in the cohort's ancestry
  (earwax); it correctly fails on a weak association (cilantro). Not "regime-1 wins everywhere."
- **Across-species (cat/dog/cattle coat colour) remains data-walled** — these are HUMAN cells (openSNP);
  the animal visible-trait cells still need a free, non-circular, joinable genotype×phenotype table (the
  horse-coat lesson; horse was only unblocked by the Sarcidano PMC mirror). Not attempted here (external
  wall, unchanged).
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9). openSNP `.zip` on D:
  (gitignored). Label tier = self-report (PILOT/DEMO), like the eye-colour cell.

## Reproduce
```bash
uv run python scripts/single_snp_opensnp_validate.py --trait all   # earwax / lactase / cilantro
uv run pytest tests/test_single_snp_traits.py -q                   # 6 offline tests (synthetic zip, no D:)
# data: OpenSNP archive dump (archive.org/details/opensnp_data_dumps, 2017-12-08) on D:
```
