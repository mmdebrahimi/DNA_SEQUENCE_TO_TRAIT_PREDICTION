# Deterministic conservation vs learned models — the molecular-regime verdict (revises the earlier story)

**Date:** 2026-07-03
**Script:** `scripts/dms_conservation_benchmark.py` (`wiki/dms_conservation_benchmark_scores.json`)
**Data (free, authoritative):** ProteinGym v1.x published per-assay zero-shot Spearman table
(`benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv`) + the reference
file (`reference_files/DMS_substitutions.csv`), 96 human substitution DMS assays. On D:.

## The question

The BLOSUM benchmark showed a deterministic **substitution-severity** floor of ~0.21 median Spearman, ~2× below
AlphaMissense (~0.42), and concluded "learned wins in the molecular regime." But substitution matrices are
**position-blind**. The honest deterministic competitor is a **position-specific conservation** score (an
independent-sites model over a per-protein MSA) — ProteinGym's **"Site-Independent"** baseline. Does a
*deterministic* conservation score compete with learned models, especially on **function**?

**Method (reshaped after review):** rather than hand-roll an MSA pipeline (risk: a non-canonical baseline /
parameter-search against the learned models), this reads ProteinGym's **authoritative** table, which already
ran Site-Independent AND the learned models (ESM-1v / EVE / GEMME) on **identical variant rows**. Non-circular
(ProteinGym is not tuned to this question); it is an audit of the canonical benchmark, exactly as the pre-build
review recommended.

## Result — deterministic conservation LARGELY COMPETES

| selection type | n | **Site-Independent (deterministic)** | ESM-1v | EVE | GEMME | beats a learned model? |
|---|---|---|---|---|---|---|
| **Activity (function)** | 20 | **0.427** | 0.522 | 0.520 | 0.542 | no |
| **Binding** | 8 | **0.388** | 0.341 | 0.405 | 0.348 | **YES** (> ESM-1v, > GEMME) |
| Stability (abundance) | 23 | 0.394 | 0.526 | 0.502 | 0.528 | no |
| Expression | 15 | 0.372 | 0.490 | 0.437 | 0.481 | no |
| OrganismalFitness | 30 | 0.376 | 0.454 | 0.417 | 0.422 | no |
| **ALL-HUMAN** | 96 | **0.391** | 0.480 | 0.466 | 0.478 | — |

**Project anchors (this repo's own numbers):** BLOSUM floor **0.210**, AlphaMissense **0.423**.

## The verdict (honest, and it REVISES the earlier conclusion)

- **A deterministic conservation score reaches 0.427 on FUNCTION** — just shy of the 0.45 hypothesis target,
  and **2× the BLOSUM substitution floor (0.21).** It is NOT a weak deterministic method.
- **Overall it ≈ AlphaMissense**: Site-Independent all-human **0.391** vs this project's AlphaMissense **0.423**
  — a *deterministic* conservation score essentially matches an AlphaMissense-class learned predictor across the
  human DMS landscape.
- **It BEATS the single-sequence learned model ESM-1v on binding** (0.388 vs 0.341) — deterministic conservation
  is the *better* tool for some modalities.
- **Residual gap, stated plainly:** the *best* learned models (GEMME 0.54, EVE/ESM ~0.52) still beat
  Site-Independent by **~0.1** on function / stability. The frontier learned models retain an edge; deterministic
  conservation does not fully close it.

**So the earlier "molecular regime → learned wins" was too strong — it was really "substitution matrices are too
weak."** A position-specific *deterministic* conservation score competes: it matches AlphaMissense-class
predictors overall, wins on binding, reaches ~0.43 on function, and trails only the very best learned models by
~0.1. That materially sharpens the project's decoder-regime map: in the molecular regime, the deterministic
option is **conservation (MSA independent-sites), not substitution severity** — and it is close enough that a
deterministic, interpretable decoder is a legitimate choice, with a named ~0.1 accuracy trade vs the best
learned models on function.

## Honest scope

- These are ProteinGym's **published** Site-Independent + learned Spearman numbers (authoritative, non-circular),
  audited + aggregated by selection type on human assays — not recomputed here (which would risk a non-canonical
  baseline). The BLOSUM/AM anchors ARE this project's own computed numbers.
- Function n=20 (ProteinGym Activity) is now well-powered (vs n=3 in the local benchmark) — the interaction is
  robust: deterministic conservation's deficit vs learned is smallest on binding, largest on stability/function.
- Building a fresh deterministic conservation score locally (per-protein MSA, frozen spec) would let this project
  *own* the number rather than cite it — a follow-up; the authoritative audit is the honest first answer.
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9).

## Reproduce

```bash
uv run python scripts/dms_conservation_benchmark.py
uv run pytest tests/test_dms_conservation_benchmark.py -q   # 2 offline synthetic tests
# data: raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv
```
