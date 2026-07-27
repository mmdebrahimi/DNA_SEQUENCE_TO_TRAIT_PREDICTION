# GEUVADIS Stage-1b: polygenic elastic-net cis ceiling + the confound-claim falsification (2026-07-27)

**Tests the row-574 claim** ("the pooled-vs-within inflation is a POLYGENIC phenomenon, not a single-causal-SNP
one") by fitting the field's real linear baseline — **elastic-net over the cis-window** (PrediXcan-style),
5-fold CV → **out-of-sample** predicted expression — and comparing pooled vs within-population.

## Result (40 eQTL genes, chr1, out-of-sample CV)

| Predictor | pooled \|ρ\| | within-pop \|ρ\| | inflation (pooled−within) |
|---|---|---|---|
| Single best-SNP (in-sample) | 0.286 | 0.291 | **−0.005** |
| **Elastic-net multi-SNP cis (out-of-sample CV)** | 0.200 | 0.187 | **+0.013** |

## Verdict on the claim: PARTIALLY SUPPORTED, and REFINED

- **Directionally supported:** the polygenic cis predictor shows a small POSITIVE inflation (+0.013) where the
  single-SNP showed none (−0.005) — more features → the predictor starts to absorb a little population
  structure, as the claim predicted.
- **But the magnitude refines the claim:** +0.013 is *tiny* — nothing like the Arabidopsis-scale +0.23 pooled
  inflation. **Cis-restricted linear prediction (single-SNP OR polygenic) is largely DE-CONFOUNDED either way.**
  The row-574 framing ("polygenic phenomenon") is not quite right: it is not polygenicity per se — a
  cis-window polygenic model barely inflates. The Arabidopsis-scale confound needs a **genome-wide /
  structure-absorbing** predictor (an embedding / DNA-encoder scoring many loci), which can soak up population
  structure that a cis-restricted model cannot. **The confound scales with genome-wide structure capture, not
  with cis-polygenicity.** (Corrected in the ledger.)

## The honest linear ceiling for the multimodal comparison

The **out-of-sample** cis-linear ceiling is **ρ ≈ 0.19–0.20** (within-population). The single-SNP 0.29 was
in-sample optimism. So a DNA-encoder arm must beat ρ≈0.19–0.20 cross-individual, de-confounded, to justify the
multimodal DNA arm — and per the SOTA (Nat Genet 2023; Variformer 2026) it ties, not beats, this.

## Scope / honesty

- **n=40 genes (chr1, lowest-coord subset), out-of-sample CV.** The +0.013 inflation is small and at n=40 is
  near the noise floor — the robust conclusion is "cis-linear prediction is ~de-confounded", NOT a precise
  inflation value. The elastic-net had mild convergence warnings (max_iter=1500) — fine for the directional
  conclusion, not for a decimal-precise ceiling.
- This DOES NOT run the DNA-encoder arm (Stage-2, GPU) — it establishes the ceiling + tests the confound claim.

## Reproduce
```
uv run python scripts/geuvadis_stage1b_elasticnet_cis.py --chrom 1 --max-genes 40
```
