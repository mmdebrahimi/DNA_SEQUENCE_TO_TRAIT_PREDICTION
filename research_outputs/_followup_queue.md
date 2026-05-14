# Research Follow-up Queue (V1.5 invocation)
<!-- queue-schema: 0.1 -->

> Generated 2026-05-13. Stale-days threshold: 30. Source memos scanned: 1.
> This queue surfaces "Decisions for Human Confirmation" rows from supported memos —
> NOT a promotion list. Human review + per-memo Promotion Gates remain required before
> any number lifts into rules / wiki / code.

## Summary

- Source memos scanned: 1
- Schema-drift memos (skipped): 0
- Total raw candidate rows extracted: 5
- Unique candidates after dedup: 5
- Active candidates (≤ 30 days old): 5
- Stale candidates (> 30 days old): 0
- Cross-flavor candidates (public + internal): 0

## Active Candidates

| Claim | Numeric value | Units | Locator(s) | Candidate use / Verification needed | Confidence | Sources | Problem-anchor(s) | Notes |
|---|---:|---|---|---|---|---|---|---|
| Random Forest selected as standard downstream classifier on frozen DNA-FM embeddings | qualitative | — | https://pmc.ncbi.nlm.nih.gov/articles/PMC12663285/ | **Candidate use:** at the dna_decode decision gate, evaluate Random Forest alongside XGBoost as the head on frozen NT embeddings; field consensus is RF for DNA-FM downstream tasks. **Verification needed:** confirm the author/title of the Nature Communications source (author-identity-uncertain flag); confirm RF outperforms XGBoost on dna_decode's specific 12-strain cipro mini-cohort empirically before adopting. | medium | sota-bacterial-amr-prediction-small-cohorts-2026-05-13 | SOTA architectures for predicting bacterial AMR from genome sequence on small cohorts (2024-2026), with attention to: classifier head choice on frozen foundation-model embeddings, alternatives to XGBoost (deep tabular, MLP, attention pooling), and fine-tuning vs frozen-embedding trade-offs on small (12-150 strain) cohorts | author-identity-uncertain (confidence downgraded high→medium per source-text identity advisory) |
| XGBoost ranks ~8.30 vs TabPFN ~4.88 on small datasets (≤1250 samples, 57-dataset benchmark) | 8.30 | rank (lower=better) | https://arxiv.org/abs/2305.02997 | **Candidate use:** add TabPFN to the dna_decode decision-gate comparison (alongside AMRFinder/k-mer/gene-presence/clade-only/NT-XGBoost). TabPFN is specifically designed for small-sample regimes; XGBoost has known disadvantage there. **Verification needed:** check TabPFN's input-dimension limit (it has well-known constraints around ~100 features in early versions); confirm it can ingest 512-dim NT embeddings either via PCA or via TabPFN v2's larger limits. | high | sota-bacterial-amr-prediction-small-cohorts-2026-05-13 | SOTA architectures for predicting bacterial AMR from genome sequence on small cohorts (2024-2026), with attention to: classifier head choice on frozen foundation-model embeddings, alternatives to XGBoost (deep tabular, MLP, attention pooling), and fine-tuning vs frozen-embedding trade-offs on small (12-150 strain) cohorts | — |
| GBDT performance positively correlated with samples-to-features ratio (NN-class opposite) | qualitative | — | https://arxiv.org/abs/2305.02997 | **Candidate use:** flag in dna_decode CLAUDE.md as a known limitation — current architecture (frozen 512-dim NT embeddings + XGBoost) is in the GBDT-disfavored regime (N=12-150, D=512, samples << features). Use this as motivation for either (a) feature reduction (mean-pool to per-gene aggregate) or (b) classifier swap (RF / TabPFN / MLP). **Verification needed:** measure actual samples-to-features ratio at decision-gate time; check whether dimension reduction (PCA, learned linear projection) improves XGBoost vs leaving raw embeddings. | high | sota-bacterial-amr-prediction-small-cohorts-2026-05-13 | SOTA architectures for predicting bacterial AMR from genome sequence on small cohorts (2024-2026), with attention to: classifier head choice on frozen foundation-model embeddings, alternatives to XGBoost (deep tabular, MLP, attention pooling), and fine-tuning vs frozen-embedding trade-offs on small (12-150 strain) cohorts | — |
| E. coli cipro + AMRFinderPlus + k-mer + XGBoost = >90% accuracy on 256-genome cohort (5-fold CV) | >90 | % accuracy | https://pubmed.ncbi.nlm.nih.gov/39320197/ | **Candidate use:** sets the classical-baseline floor that dna_decode's NT-XGBoost approach must beat at the decision gate. Anything ≤90% on a similar-sized cohort would indicate NT is NOT adding signal beyond classical features. **Verification needed:** confirm Talamantes-Becerra cohort size + drug + CV method match dna_decode's intended comparison; the 256-genome cohort is larger than the planned 12-strain mini, so direct comparison requires running on a 150-250 strain cohort to be apples-to-apples. | high | sota-bacterial-amr-prediction-small-cohorts-2026-05-13 | SOTA architectures for predicting bacterial AMR from genome sequence on small cohorts (2024-2026), with attention to: classifier head choice on frozen foundation-model embeddings, alternatives to XGBoost (deep tabular, MLP, attention pooling), and fine-tuning vs frozen-embedding trade-offs on small (12-150 strain) cohorts | — |
| For ciprofloxacin in E. coli, SNP tables may outperform gene presence-absence (resistance is mutation-driven) | qualitative | — | https://pmc.ncbi.nlm.nih.gov/articles/PMC11684616/ | **Candidate use:** when dna_decode runs the classical baselines at the decision gate, prioritize a SNP-table feature variant (gyrA / parC / parE specific positions) over generic gene-presence for cipro specifically. Different drugs may want different feature types. **Verification needed:** confirm dna_decode's current `classical_baselines.py` supports SNP-table features OR add a SNP-feature variant; cross-check against AMRFinderPlus mutation-aware output schema. | high | sota-bacterial-amr-prediction-small-cohorts-2026-05-13 | SOTA architectures for predicting bacterial AMR from genome sequence on small cohorts (2024-2026), with attention to: classifier head choice on frozen foundation-model embeddings, alternatives to XGBoost (deep tabular, MLP, attention pooling), and fine-tuning vs frozen-embedding trade-offs on small (12-150 strain) cohorts | — |

## Stale Candidates

(None — only one source memo, captured today.)

## Schema-drift Memos

(None — all source memos carry a recognized schema marker.)

## Footer notes

- Verbatim values are preserved from source memos. This skill does NOT edit, rewrite, or interpret claim text.
- Public-source rows (research-intake) and internal-source rows (athena-intake) preserve their original locator semantics. Internal doc IDs do not resolve at a public URL.
- Promotion Gate (per memo type) remains the human-confirmation step that lifts candidates into rules/wiki. This queue does NOT bypass it.
