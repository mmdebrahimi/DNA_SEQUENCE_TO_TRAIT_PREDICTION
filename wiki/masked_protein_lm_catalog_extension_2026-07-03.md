# Masked-protein-LM "extend the catalog" test — a decisive NEGATIVE with a clear mechanism (2026-07-03)

> **AMENDED 2026-07-03 (adversarial review): this result is UNDERPOWERED / INCONCLUSIVE, not a decisive
> NO_WIN, and the mechanistic "resistance ⊥ conservation" claim is a plausible CONCERN, not an established
> law.** Two flaws: (1) the blindspot AUC (0.293) compared R variants vs only n=2 benign variants at
> DIFFERENT positions — ESM LLR is position-conservation-dominated, so that AUC largely measures position
> conservation, not R-vs-S; (2) the Spearman (0.073) is WITHIN the all-R set (magnitude ordering among
> resistant variants), not the R-vs-S separation the hypothesis needs — no S set was scored. The proper
> test (HIV RT/PR/IN single-mutant isolates, rich R+S at attributable per-variant folds) is run separately
> (`wiki/esm_hiv_resistance_matched_test_2026-07-03.md`) and supersedes this verdict. The genotype-level
> LD-confound argument (below) stands for PURE self-supervised models; a supervised head is a separate case.

**The user's masked-prediction idea ("mask parts, predict the missing part — CLIP/world-model style"),
tested at the granularity where it is CAUSAL (protein, not genotype), against the sharpest target (does it
fill the deterministic AMR catalog's documented gaps?). Phase-1 result: `NO_WIN` for AMR — and the reason is
a real regime-boundary refinement, not a dead end.** Script: `scripts/esm_catalog_extension_test.py`.

## What was tested (the good, cheap, decisive version of the idea)

Not "replace the deterministic decoder with a learned embedding" (that lost 0-for-5 de-confounded). Instead:
does a masked protein-LM (**ESM-2 650M**, zero-shot masked-marginals, Meier 2021 — literally "predict the
masked residue") rank the SARS-CoV-2 Mpro catalog's OWN MISSED true-resistant variants (the documented FN
blindspot) as more deleterious than benign polymorphisms? If yes → a learned layer that plugs the catalog's
holes = the hybrid the user envisioned. Real ESM inference on the committed NC_045512 Mpro protein
(H41/C145 catalytic + E166 nirmatrelvir verified); label = CoV-RDB measured fold (the cell's own label).

## Result — clean null for AMR

| metric | value | reading |
|---|---|---|
| Spearman(−LLR, log fold) over 29 variants | **0.073** | ESM does NOT order Mpro variants by resistance magnitude (null) |
| blindspot-rescue AUC (missed-R vs benign) | **0.293** | most catalog-missed R variants score LESS deleterious than benign |
| verdict | **`NO_WIN`** | zero-shot masked protein-LM does not extend the AMR catalog |

The scoring is functioning (the large-effect misses F140A fold 21 → LLR −1.71, S144E fold 480 → LLR −1.45
DO register as deleterious); the signal simply does not track resistance across the set.

## Why — the mechanistic finding (resistance ⊥ conservation)

A masked protein-LM measures **evolutionary tolerance / fitness** — "does this substitution break the
protein." **Drug resistance is a DIFFERENT axis**: a resistance mutation changes drug BINDING while
*preserving* enzyme function — so it is, by design, evolutionarily TOLERATED (low ESM deleteriousness). The
resistance signal is therefore **orthogonal** to the fitness signal a masked-LM learns. ESM's ~0 Spearman on
Mpro AMR fold is not a bug or a capacity limit — the SAME model reaches ~0.4–0.5 on FUNCTION/fitness DMS
assays (ProteinGym; the project's own DMS arc). The gap between 0.5-on-function and 0.07-on-resistance IS the
axis mismatch.

**This sharpens the regime map:**

| molecular phenotype | masked protein-LM / conservation |
|---|---|
| function / activity / stability / abundance (DMS) | **WINS** (validated last session: ESM/EVE/GEMME ~0.5; our own conservation 0.476) |
| **drug resistance** (evade binding, keep function) | **BLIND** (this test: ~0 correlation) — it is curated-determinant knowledge, not a fitness signal |

And it retroactively explains WHY the deterministic curated catalog is the right tool for AMR: resistance is
specific, non-fitness, curated knowledge a general masked-LM cannot infer.

## Consequence for the sequenced plan ("both, protein first")

Phase 1 did NOT clear the bar → the **genotype world-model (Phase 2) is NOT greenlit**. Two independent
reasons now converge: (a) at the genotype granularity the masked objective learns LD/population structure
(the confounder — the LD-vs-purifying-selection distinction) so it re-hits the 0-for-5 wall; (b) at the
protein granularity, where the paradigm IS causal, it captures FUNCTION but is BLIND to RESISTANCE, so it
does not add to the AMR product. The masked-prediction paradigm's value for this project is already
realized where it applies — molecular function (the conservation/AlphaMissense arc) — and does not intersect
the curated resistance catalogs.

## Honest scope
- Benign negative set is small (2 sourced polymorphisms; CoV-RDB is R-enriched) → the AUC is a first-pass
  indicator, BUT the null Spearman over 29 variants is the robust, well-powered readout.
- Model is 650M zero-shot; a 3B/15B ESM is the only residual scale-up, but the axis-mismatch is mechanistic
  (fitness ≠ evasion), so a bigger fitness model is not expected to close a resistance gap. HIV RT/PR/IN
  (rich R+S) would be the powered confirmation IF the paradigm had shown any signal here — it did not.
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9); READ-only over the frozen
  catalog. torch-cpu + fair-esm installed (removable; ESM weights on D:).

## Reproduce
```bash
HF_HOME=D:/hf_cache TORCH_HOME=D:/dna_decode_cache/torch \
  uv run python scripts/esm_catalog_extension_test.py --model esm2_t33_650M_UR50D
uv run pytest tests/test_esm_catalog_extension.py -q   # offline helper tests (no model)
```
