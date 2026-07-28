# Essentiality E3 — human learned complement (BAGEL, 2026-07-28)

**Rounds out "any organism":** the same cheap-feature classifier (aa-composition + length + conserved-core
score), 5-fold CV, now trained on the **human** BAGEL CEG2/NEG gold-standard. Sequences from UniProt human
reviewed proteome (20,416; symbol→seq by GN=). $0, no GPU.

## Result (human, n=1539: 659 essential / 880 non-essential)

| model | CV AUROC |
|---|---|
| conserved-core (E.coli-tuned, human transfer) | 0.572 |
| logistic regression | 0.825 |
| **HistGradientBoosting** | **0.911** |

**AUROC lift +0.339** (0.572 → 0.911). Tail-recovery: on the 1445 genes the conserved-core misses (567
essential), the learned model scores **AUROC 0.898** — it recovers essentiality signal the E.coli-tuned core
is blind to in human.

## Reading it — verdict LEARNED_COMPLEMENT_EARNS_KEEP (decisively)

- The lift is far bigger than E. coli (+0.339 vs +0.098) for a clear reason: the conserved-core catalogue is
  **E. coli-tuned**, so it transfers poorly to human (0.572), leaving huge headroom — while the learned model
  is *trained on human* and captures the human-specific essential core (proteasome/spliceosome/etc.) the
  catalogue can't see. This is exactly the "per-organism learned complement" the plan predicted.

## Honesty (H8) — the absolute numbers are NOT cross-organism comparable
- **BAGEL CEG2 vs NEG are curated to be cleanly SEPARABLE** reference sets (designed for CRISPR-screen
  normalization: "always-essential" vs "never-essential" genes), so 0.911 is on an EASIER task than the
  genome-wide E. coli Goodall benchmark (all genes, learned 0.795). Do NOT read 0.911 as "genome-wide human
  essentiality accuracy" — it is CEG2-vs-NEG reference-set discrimination. The valid cross-organism claim is
  the DIRECTION (learned ≫ core in both), not the absolute magnitudes.
- 5-fold CV, per-protein-independent features; tail-recovery (0.898 on core-missed genes) confirms the signal
  is sequence-intrinsic, not core-score reuse.

## The essentiality cell now spans both organisms with a learned complement
| organism | conserved-core | learned (E3) | benchmark |
|---|---|---|---|
| E. coli | 0.697 | **0.795** (+0.098) | Goodall genome-wide (harder) |
| human | 0.572 | **0.911** (+0.339) | BAGEL CEG2/NEG (separable reference) |

## Reproduce
Inline pipeline (UniProt human FASTA + gene_info + BAGEL on D:); same features as `essentiality_e3_learned.py`.
