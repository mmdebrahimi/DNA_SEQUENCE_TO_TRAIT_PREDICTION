# Protein variant-effect cell — genotype→MOLECULAR phenotype (mapping the deterministic/learned boundary)

**Date:** 2026-07-02
**Script:** `scripts/proteingym_variant_effect.py` (`wiki/proteingym_variant_effect_scores.json`)
**Data (free, unit-level):** PTEN VAMP-seq abundance DMS (Matreyek et al. 2018), fetched from **MAVEDB**
`urn:mavedb:00000102-0-1` (5,477 single-amino-acid variants × measured protein abundance). On D:.

## Why this cell exists

The session's earlier answer to "more genotype→phenotype in humans/higher animals" identified the tractable
frontier: when the **organism-level** trait is polygenic (the learned path is 0-for-5 under de-confounding),
shift to a **molecular** phenotype where the genotype→phenotype map is short and the label is free + unit-level.
This is that cell: **protein variant → measured function**, on the canonical human PTEN abundance DMS.

The deterministic feature is the **BLOSUM62 substitution score** of each missense change (authoritative, from
Biopython — no hand-transcribed matrix). Conservative substitution (high BLOSUM) → abundance preserved;
disruptive (low BLOSUM) → abundance reduced.

## Result

| metric | value |
|---|---|
| missense variants scored | 5,083 |
| **Spearman(BLOSUM62, PTEN abundance)** | **0.217** |
| nonsense mean abundance | 0.177 |
| missense mean abundance | 0.717 |
| direction sanity (nonsense ≪ missense) | ✅ |

- The deterministic substitution-severity score has **real but modest** signal (ρ = 0.22): conservative
  substitutions preserve PTEN abundance, disruptive ones destabilize it.
- **Direction sanity is perfect**: 201 nonsense (truncating) variants average abundance 0.177 vs 0.717 for
  missense — truncation destabilizes the protein, exactly as expected. This validates the assay + the encoding.

## The point — this maps the deterministic/learned boundary

This is the **first modality in the project where the LEARNED path clearly beats the deterministic rule.** The
published ProteinGym leaderboard puts zero-shot models (ESM-1v / EVE / TranceptEVE) in the **~0.4–0.5 Spearman**
range on PTEN abundance — roughly double the BLOSUM62 baseline. (Cited, **not** run here — those need a GPU.)

That is the **opposite** of the project's AMR / PGx / ben-1 cells, and it sharpens the law the project has been
building:

| regime | winning decoder | why | project cells |
|---|---|---|---|
| curated high-effect catalog exists | **deterministic** determinant scan | a few variants dominate; the rule is the biology | AMR, TB, fungal, HIV, SARS, PGx, ben-1 |
| organism-level polygenic trait | *neither* (learned 0-for-5 de-confounded) | thousands of tiny effects + population-structure confound | yeast, Arabidopsis, cipro, pathotype, DGRP |
| **molecular property, no catalog** | **learned** (evolutionary/structural context) | distributed signal that self-supervised models capture; BLOSUM is a weak floor | **this cell (PTEN variant effect)** |

So "more human/higher-animal G2P" is real at the **molecular** level — and it's the one place a learned model is
the right tool, not the deterministic scan. The deterministic BLOSUM baseline (ρ=0.22) is the honest floor that
quantifies how much the learned models actually add.

## Honest scope

- **Single assay, single deterministic feature.** BLOSUM62 is a deliberately simple baseline; Grantham distance
  or a physicochemical model would likely edge it up but stay well below the learned models. The cell's job is to
  establish the deterministic floor + the modality, not to compete with ESM.
- **The learned numbers are cited from ProteinGym, not reproduced here** (GPU). No fabricated figure — the
  contrast is stated as the published ~0.4–0.5 range.
- MAVEDB also hosts a **CYP2C9** activity DMS (`urn:mavedb:00000095-a-1`) — a natural future bridge to the
  project's existing CYP2C9 PGx cell (star-allele metabolizer call ↔ per-variant molecular effect).
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9); this is a new research cell.

## Reproduce

```bash
uv run python scripts/proteingym_variant_effect.py    # reads D:/dna_decode_cache/proteingym/pten_scores.csv
uv run pytest tests/test_proteingym_variant_effect.py -q   # 2 offline synthetic tests
# data: https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000102-0-1/scores
```
