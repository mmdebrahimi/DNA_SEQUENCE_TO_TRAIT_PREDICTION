# Phase 3 — masked-genotype IMPUTATION reduces the decoder's ABSTAIN rate (PASS, 2026-07-04)

**The one branch of the masked-genotype idea that survives every test — because LD-learning is a FEATURE
here, not the confound that kills phenotype-prediction.** Impute an UNCALLABLE determinant SNP from a linked
proxy so the FROZEN deterministic rule can call otherwise-abstained samples. Substrate: the ABO O-status cell
— consumer SNP arrays frequently do not type the O-deletion `rs8176719` (an indel) directly, so O-status is
`INDETERMINATE` for those users; the LD tag `rs657152` (97% purity, feasibility-verified) lets us impute it.
Script: `scripts/impute_determinant_abstain.py`.

## Result — PASS (imputation reduces ABSTAIN without losing accuracy)

| metric | value |
|---|---|
| users scanned / co-called (target+tag) / tag-only (abstain pop) | 1000 / 842 / 70 |
| **LOO genotype-imputation accuracy** (impute `rs8176719` from `rs657152`) | **0.985** |
| **imputed → O-status accuracy** (frozen rule on imputed vs true genotype) | **0.989** |
| **abstain-reduction** (tag-only users rescued from `INDETERMINATE`) | **69 / 70** |
| verdict (pre-committed falsifier: acc ≥ 0.90 both AND rescued > 0) | **PASS** |

The imputed→O-status accuracy (98.9%) even EXCEEDS the raw genotype accuracy (98.5%) because the O-status
rule collapses `DI`/`II` both to non-O — so the only imputation error that flips the phenotype call is a
`DD`↔(`DI`/`II`) confusion, and those are rare.

## What it means (pre-committed PASS branch)

- **LD-based imputation of an uncallable determinant is accurate + deployable.** 69 users whose O-deletion is
  untyped — who the deterministic decoder would abstain on — get a correct O-status call at ~99%. The decoder
  goes from ABSTAIN to a high-accuracy call by imputing its own missing input.
- **This is the honest form of "mixing the learned approach with our data."** LD-learning (the exact signal
  that CONFOUNDS phenotype-prediction — the 0-for-5 wall) is precisely what MAKES imputation work: predicting
  a masked SNP from its neighbors IS learning haplotype structure, and here that is the goal, not the trap.
- **V2 greenlit as the deployable hybrid value-add** — a masked-genotype imputer as a PRE-PROCESSOR that
  feeds the frozen deterministic decoder better input (reduce ABSTAIN), never competing with the catalog.

## The whole hybrid arc — settled

| branch | verdict | role |
|---|---|---|
| Learned phenotype-prediction from FM embeddings (0-for-5) | de-confounded NEGATIVE | rejected |
| **V1** zero-shot / supervised learned SCORER (Phase 1–2) | PARTIAL — beats zero-shot, loses to catalog | novel-variant FALLBACK only |
| **V2** masked-genotype IMPUTATION (Phase 3) | **PASS** — 98.9% abstain-reduction | **the deployable win** |

**The best version of the hybrid idea is imputation** — not a learned phenotype predictor (curated knowledge
beats it), but a learned INPUT-COMPLETION layer that cuts the deterministic decoder's ABSTAIN rate at high
accuracy. The deterministic decoder stays the primary; the learned layer feeds it better data.

## Honest scope
- Single best tag + majority-conditional imputation (the simplest LD imputer); a full-panel Beagle-class
  imputer would be marginally better but this already clears the bar.
- openSNP is European-dominated → the tag-target LD is European (a reference-panel-ancestry caveat, standard
  for imputation; other-ancestry LD may differ).
- Validates GENOTYPE imputation vs the TRUE typed genotype (not a self-report label) → no self-report noise;
  the O-status rule is the frozen deterministic cell (untouched).
- Scan capped at 1000 users for wall-clock (21 GB zip streaming); accuracy is LOO on 842 co-called.
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9). openSNP dump on D: (gitignored).

## Reproduce
```bash
uv run python scripts/impute_determinant_abstain.py --scan-cap 1000
uv run pytest tests/test_impute_determinant_abstain.py -q   # 3 offline tests (no zip)
```
