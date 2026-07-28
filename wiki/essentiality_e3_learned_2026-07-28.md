# Essentiality E3 — the learned complement earns its keep (E. coli, 2026-07-28)

**Does a supervised classifier beat the conserved-core decoder + recover the tail it misses?** YES — and
**cheaply**: aa-composition + protein length + the conserved-core score, 5-fold stratified CV on the Goodall
gold-standard. **No GPU, no Kaggle, $0** (the R2 pre-bar "try cheap features first" call — ESM2 turned out
unnecessary for the E3 verdict).

## Result (E. coli, n=3776: 354 essential / 3422 non; base rate 9.4%)

| model | features | CV AUROC |
|---|---|---|
| conserved-core (deterministic baseline) | function catalogue | 0.697 |
| logistic reg | aa-composition only | 0.721 |
| logistic reg | aa + length + core-score | 0.783 |
| **HistGradientBoosting** | aa + length + core-score | **0.795** |

**AUROC lift: +0.098 (0.697 → 0.795).** The learned complement is a strictly better ranker.

## The "lift past the ceiling" test, done correctly

The naive recall-at-matched-specificity comparison is **INVALID here** — the conserved-core score is
degenerate (most genes score exactly 0), so specificity-quantile thresholds collapse (verify-in-batch caught
a spurious "core recall = 1.0"). The valid metric is **tail recovery**: on the **3586 genes the conserved-core
MISSES** (score < 2, its non-essential prediction), 220 are truly essential — and the learned model scores
**AUROC 0.685 on that tail** (null 0.5). **So the learned complement recovers essentiality signal in exactly
the region the deterministic core is blind to** — the honest form of "lift recall past the conserved-core
ceiling."

## Verdict: LEARNED_COMPLEMENT_EARNS_KEEP
+0.098 AUROC overall AND 0.685 AUROC on the core-missed tail, from cheap CPU features. This validates the
plan's E3 hypothesis (a learned complement beats the deterministic core) — and shows it doesn't even need
ESM2/GPU to do so for E. coli. (ESM2 embeddings on Kaggle remain an optional polish that MIGHT lift further;
not required for the verdict.)

## Honesty (H8)
- 5-fold CV (no leakage); features are per-protein independent. AUROC is the valid comparison; recall@spec is
  not reported (degenerate baseline score).
- The core-score is included as a feature, so the +0.098 is the ADDITIONAL signal aa-composition+length bring
  ON TOP of the catalogue (aa-composition alone already beats the core at 0.721).
- E. coli only; the human arm (BAGEL labels) is the natural next run (human protein sequences needed).

## Reproduce
`scripts/essentiality_e3_learned.py` (needs D: Goodall XLSX + E. coli feature table + protein FASTA). $0, CPU.
