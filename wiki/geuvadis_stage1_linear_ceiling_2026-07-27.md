# GEUVADIS Stage-1: the linear cis-eQTL ceiling (project-owned organism-multimodal number, 2026-07-27)

**The number the organism-multimodal DNA-encoder arm has to beat — measured on our own data, free.**

## Result

Single best cis-eQTL SNP → gene expression, cross-individual, on **652 genes** (chr1–3) × **462 GEUVADIS
individuals** across 5 populations (CEU/FIN/GBR/TSI/YRI), evaluated **pooled vs WITHIN population**:

| | mean \|Spearman ρ\| |
|---|---|
| **Pooled** (all individuals) | **0.286** |
| **Within-population** (each pop vs its own null, sample-weighted) | **0.291** |
| **Inflation** (pooled − within) | **−0.005** (essentially zero) |
| genes with within-pop ρ ≥ 0.1 | 99.7% |

Per-chromosome (chr1 296 / chr2 207 / chr3 149 genes) the numbers are flat to ±0.006 — a stable estimate.

## Interpretation (the honest, slightly surprising finding)

**The linear cis-eQTL ceiling is ρ ≈ 0.29, and it is genuinely DE-CONFOUNDED** — a single validated causal
cis-eQTL SNP predicts expression *within* each population as well as pooled (within is even fractionally
higher). This is the biologically correct result: cis-eQTLs are largely population-shared causal variants, so
they carry real cross-individual signal that does not depend on population structure.

This **refines** the project's R3 population-structure lesson rather than contradicting it: the pooled-inflation
confound (Arabidopsis +23 pp pooled → +3.4 pp within) is a property of **polygenic / embedding** predictors
that can *learn* population structure as a shortcut — **not** of a single validated causal SNP. So the correct
"linear ceiling" for the multimodal comparison is this de-confounded ρ ≈ 0.29, not an inflated pooled number.

**Consequence for the organism-multimodal question:** the DNA-encoder arm has to beat ρ ≈ 0.29 cross-individual
to justify adding it. The field's SOTA on this exact task (Nat Genet 2023; Variformer 2026) shows sequence
models **tie** the linear cis-eQTL / elastic-net ceiling cross-individual and do not generalize to unseen loci.
So the DNA arm is expected to **match, not beat**, this ceiling — the project-owned confirmation of row 572's
resolution that organism-level multimodal is closed, now with *our own* number instead of a cited one.

## Honesty / scope

- **Single-best-SNP is a CONSERVATIVE FLOOR of the linear ceiling.** A full elastic-net over the cis-window
  (PrediXcan-style, multiple SNPs) would be somewhat higher than ρ ≈ 0.29. The de-confounding conclusion
  (within ≈ pooled) is what matters and is robust to this.
- **In-sample eQTL selection.** The best SNP per gene was chosen on the EUR373 GEUVADIS analysis, so the
  absolute pooled/within values are in-sample optimistic. The **pooled-vs-within CONTRAST** (the de-confounding
  signal) is robust to this — both arms share the same optimism.
- **This is Stage-1 (the ceiling).** Stage-2 — actually running a DNA-encoder arm (Enformer/Borzoi-style ref→alt
  delta) over GRCh37 cis-window sequences and comparing to this ceiling within-population — is the GPU/Kaggle
  step that would *directly* confirm the tie. Data for it (GRCh37 genome + GENCODE v19) is already downloaded.

## Reproduce

```
uv run python scripts/geuvadis_stage1_linear_ceiling.py --chroms 1,2,3
```
Inputs (free, on D:): `GD462.GeneQuantRPKM` + `E-GEUV-1.sdrf` + `EUR373.gene.cis.best` + matched genotype VCFs.
Evaluator: `dna_decode/organism_multimodal/deconfound_eval.py` (pooled vs within-population Spearman vs null).
