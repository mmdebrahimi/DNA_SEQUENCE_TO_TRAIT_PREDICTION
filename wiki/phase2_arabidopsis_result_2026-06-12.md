# Gate G2 — Arabidopsis flowering-time embedding result (canonical, 2026-06-12)

**Verdict: RESOLVED — stable mixed / ambiguous, embedding NEGATIVE under de-confounding.**
The embedding-niche hypothesis (H2) is **falsified** in its best-designed test.

Source: `wiki/pathb_databricks_handoff_2026-06-12.md` (workhorse: Databricks `Eng-GPU-Cluster-01`,
`PlantCaduceus_l32`, FT10 cohort n=1003 / 9 PC groups). Ported to this machine via Downloads handoff;
the raw run JSON/MD live on the Databricks volume (`.../RCA Analysis/outputs/runs/`), not in this repo.

## What ran

- Primary model `PlantCaduceus_l32` on a real Linux GPU (T4) — the env boundary that blocked canonical G2
  on the Windows workhorse is **cleared**. Probe green (mamba_ssm, CUDA, tokenizer, forward pass all OK).
- Three real bounded packets at `window_budget=64`, seeds 42 / 7 / 99 (FT10 full1003).
- `window_budget=128` attempted, failed on T4 CUDA OOM (hardware ceiling, not a scientific result).

## Results (3 seeds, all `window_budget=64`)

| metric | embedding (s42/s7/s99) | structure-only | winner |
|---|---|---|---|
| global r² | −0.046 / −0.036 / −0.031 | −0.449 | embedding (less-bad) |
| spearman | 0.213 / 0.221 / 0.226 | 0.484 | **structure-only** |
| **within-group r²** (de-confounded) | **−0.173 / −0.110 / −0.128** | **+0.038** | **structure-only** |

Stable across all three seeds — same direction every time.

## Interpretation (the honest read)

- **Both models have negative global r²** — neither is a usable flowering-time predictor.
- The embedding only "wins" **global r²**, the metric most polluted by population structure. It captures
  *which sub-population* an accession belongs to, which correlates with flowering time but is not causal.
- On the two metrics that test for real signal — **rank correlation (spearman)** and especially
  **within-group r²** (predict flowering time *within* a kinship/PC group, i.e. with population structure
  removed, which the EP8 manifest designated as the PRIMARY de-confounder) — the embedding **loses**, and
  its within-group r² is **negative across all three seeds**.
- Conclusion: the PlantCaduceus embedding learned **population structure, not the flowering-time causal
  signal.** This is the same failure mode seen in AMR (NT learned lineage, not mechanism; within-lineage =
  chance). See `wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md`.

## Why this is the load-bearing test (not just one more negative)

Per the project's three-part embedding-niche criterion, Arabidopsis flowering time was the **best-designed
shot** for embeddings — it satisfies all three conditions AMR fails:
1. sampling-independent label ✓ (common-garden flowering-time measurement, not a clinical/sampling artifact)
2. no curated catalog ✓ (no AMRFinder-equivalent gene list to beat)
3. organism-specific depth ✓ (1003 Arabidopsis accessions)

All three met, de-confounded result negative. The embedding hypothesis took its fair best shot and did not
clear the bar.

## Do NOT

- **Do not** retry `window_budget=128` on a bigger/paid GPU to "rescue" this. A negative *within-group* r²
  is a signal-vs-structure problem; more window will not convert it positive. Scaling up is a **money gate**
  spent on a question already qualitatively answered at the de-confounding level.

## Gate status

- **G2 = RESOLVED (FAIL for the embedding hypothesis; clean negative).**
- With G1 already MET (fungal determinant-scan transfers, sens 1.0; label-limited spec), the eukaryotic
  cycle terminal (G1 ∧ G2 both resolved) is **REACHED** — cycle complete.
- Strategic consequence: the validated, shippable artifact is the **deterministic decoder suite** (10
  lineage-disclosed SCORED AMR cells), not the foundation-model embedding bet.
