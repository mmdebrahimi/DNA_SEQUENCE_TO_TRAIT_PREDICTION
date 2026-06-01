# Pathotype Learned-Model — State + Bake-off Plan — 2026-05-31

> Synthesis after the LT/ST pipeline-validation smoke. Decides the next real step: a small ExPEC-vs-EPEC representation bake-off, NOT scaling the per-CDS pipeline.

## Smoke result (pipeline validation)

`scripts/pathotype_laptop_pipeline.py` ran end-to-end on 8 strains (4 LT / 4 ST von Mentzer ETEC): download → Bakta annotate → `.ffn` CDS → NT v2 100M embed (200 genes/strain, CUDA) → mean-pool → LOSO. **Pipeline PASS** (all 8 embedded, ~5,100 CDS each).

**LOSO AUROC = 0.0000 — a degenerate small-N artifact, NOT a model verdict.** Per-strain probas were two-valued (all neg → 0.5606, all pos → 0.4184), perfectly rank-inverted: the documented N=8 balanced-LOSO mechanism where XGBoost predicts the training-prior and it inverts against the held-out label (CLAUDE.md gene-presence + calibration lessons). The LT/ST label is also circular (toxin genes the representation sees). So the number is uninformative by construction; the smoke's only goal — prove the chain runs on the laptop — was met.

## Decision: ExPEC-vs-EPEC bake-off at N≈20–30, do NOT scale yet

Scaling the per-CDS Bakta+NT pipeline to ~250 strains ≈ 40 laptop-hours of annotation before answering the actual question. Instead, one small experiment on **non-circular** labels with the project's **mandated classical control**:

| Representation | Cost | Role |
|---|---|---|
| **k-mer baseline** (`loso_kmer.run_kmer_xgboost_loso`) | cheap (no Bakta, no GPU) | the **mandated control** (CLAUDE.md:211 — FM must beat classical by ≥3pp) — run FIRST |
| **whole-genome windowed NT** mean-pool | medium (GPU, no Bakta) | drops the 40hr annotation bottleneck |
| **per-CDS Bakta + NT** | expensive (Bakta ~10min/genome) | only if it earns ≥3pp over the cheaper two |

**Contrast:** ExPEC (isolation-site labels = independent) vs EPEC (DECA-curated). Both are non-circular, unlike LT/ST.

**Pre-committed stopping rule:** if k-mer matches or beats NT, stop treating embeddings as default (the project's own Phase-2 redesign trigger). Cache all NT embeddings so pooling variants don't re-pay GPU.

## Known risks to control in the bake-off
- **tet precedent:** NT mean-pooling failed (0.40) on distributed mobile-element signal; pathotype virulence is partly distributed → real risk NT dilutes. The k-mer/gene-presence control is the honest check.
- **Confounding:** ExPEC (Salipante/Hazen US clinical) vs EPEC (Hazen DECA) differ by study/geography/year → plain LOSO may reward batch/population structure. Flag + prefer lineage/ST-aware splits; report the confound.
- **EPEC label independence** needs an operational definition before trusting it.
- **Small-N degeneracy:** check the score distribution (two-valued = degenerate) before reading any AUROC; N≈24 is better than 8 but still small.

## ETL note
ExPEC/EPEC clean rows in Horesh F1 are WGS-master accessions (Salipante JSIS…, Hazen AIEY…), not GCA — need either GCA resolution or direct ENA WGS-set FASTA fetch. The bake-off representations all need only assembly FASTA, so the fetch path is the only ETL gate.
