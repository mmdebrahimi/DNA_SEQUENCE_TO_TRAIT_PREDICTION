# Forward cell Phase 6 — human proteins + AlphaMissense method (PTEN cell, 2026-07-14)

The organism-class jump (long-horizon roadmap Phase 6): the forward variant-effect cell now handles **human
proteins**, which unlocks **AlphaMissense** as a third predictor method — exactly the regime where the
bacterial cell couldn't use it (AM is human-proteome-only by design).

## What was built

- `dna_decode/forward/am_scorer.py` — `load_am_for_uniprot` + `am_table_for_mutants` (join the cached
  `am_filtered.tsv` by UniProt accession, applying the DMS→UniProt position **offset**) + `am_tier` (AM's
  own published class thresholds: benign ≤ 0.34, pathogenic ≥ 0.564).
- `predict_effect(method="alphamissense", am_table=…)` — `raw_score = 1 − AM_pathogenicity` (flipped so
  higher = benign = preserved, consistent with the BLOSUM/ESM sign); a variant not AM-covered raises loudly.
- The forward cell (`--method alphamissense --uniprot P60484 --offset 0`) + the router (`am_table`) thread it.

## Flagship: PTEN (P60484), ProteinGym Mighell_2018 — 7,260 variants, 100% AM-covered, WT-mismatch 0

| method | Spearman(pred, DMS) | note |
|---|---:|---|
| BLOSUM62 (deterministic) | **0.182** | substitution severity |
| **AlphaMissense** (learned, human-only) | **0.540** | **+0.357 over BLOSUM** |
| ESM2-650M (learned) | *background run in flight* | appended on completion |

**AlphaMissense triples the deterministic baseline on PTEN** — the learned human variant-effect signal, on
the canonical human tumor-suppressor. The join is clean (offset 0, every mutant AM-covered, WT verified).

## Why this matters

- **New organism class:** the forward "edit → predict molecular effect" cell is no longer E. coli-only — it
  works on human proteins with the same interface (`predict_effect` / the cell / the router).
- **AlphaMissense unlocked:** the method the bacterial cell *couldn't* use (human-only) is now a first-class
  predictor, complementing ESM2 (which works on both bacteria and human). The router's Regime-B can pick
  BLOSUM / ESM2 / AlphaMissense per protein.
- **Roadmap Phase 6 opened:** eukaryotic proteins are now in scope for the forward cell.

## Honesty rails

- AlphaMissense is human-proteome-only (Cheng et al. 2023) — this is complementary coverage, not a universal
  method; the cell picks the method by organism.
- The validated quantity is the rank correlation; the per-variant tier uses AM's published class thresholds.
- Offset + WT-integrity gates: the DMS→UniProt offset is explicit and the WT residue is verified (0 mismatch).

21 → 21+ forward tests (4 new AM). Frozen decoder surface (`amr_rules` / `calibrated_amr_rules` / `mic_tiers`
/ `shipped_decoder_surface` / `cohort_manifest`) byte-unchanged (`verify_lock OK`); `dna_decode/forward` is a
NEW non-frozen package. Run: `uv run python scripts/tem1_forward_cell.py --dms-id PTEN_HUMAN_Mighell_2018
--method alphamissense --uniprot P60484 --offset 0`.
