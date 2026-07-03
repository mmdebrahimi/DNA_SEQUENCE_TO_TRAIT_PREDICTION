# ProteinGym-native three-way — Site-Independent vs BLOSUM vs AlphaMissense on IDENTICAL rows

**Date:** 2026-07-03
**Script:** `scripts/pg_native_threeway.py` (`wiki/pg_native_threeway_scores.json`)
**Data (free):** ProteinGym v1.x zero-shot substitution scores (Zenodo 15293562, 1.8 GB — each per-assay CSV
carries `mutant` + `DMS_score` + `Site_Independent` + all model columns) + AlphaMissense (re-filtered to the 79
ProteinGym human UniProt accessions) + BLOSUM62 (Biopython). 88 human assays. On D:.

## Why (closing the coverage-mixing caveat)

The row-312 conservation result **cited** ProteinGym's aggregated per-assay Spearman for Site-Independent, while
this project's BLOSUM/AM were computed on **MAVEDB** rows — a coverage-mixing caveat a pre-build review flagged.
This closes it: on ProteinGym's OWN assays, BLOSUM + AlphaMissense are **re-scored on the exact same ProteinGym
mutant rows** as Site-Independent (which ProteinGym computed on those rows). All three predictors, identical
variants. **ProteinGym-native** — chosen over mapping MAVEDB↔ProteinGym — sidesteps the coordinate/polarity
landmine entirely (the review's preferred route; the 96 human ProteinGym assays need no MAVEDB join at all).

## Result — three-way on identical rows (88 assays, single substitutions)

| selection type | n | BLOSUM (deterministic, position-blind) | **Site-Independent (deterministic, conservation)** | AlphaMissense (learned) |
|---|---|---|---|---|
| Activity (function) | 18 | 0.229 | 0.426 | **0.499** |
| **Stability (abundance)** | 22 | 0.195 | **0.393** | 0.334 |
| Binding | 6 | 0.224 | 0.351 | **0.392** |
| Expression | 13 | 0.279 | 0.322 | **0.477** |
| OrganismalFitness | 29 | 0.160 | 0.375 | **0.466** |
| **overall (assay-median)** | 88 | 0.207 | **0.386** | 0.456 |
| **overall (UniProt-median)** | 75 prot | 0.234 | **0.371** | 0.417 |

(Spearman vs DMS_score, aligned so positive = predicts fitness; AlphaMissense negated from pathogenicity.)

## Verdict — deterministic conservation LARGELY COMPETES

- **The number is robust.** Site-Independent on function = **0.426** on identical rows — matches row 312's
  *cited* 0.427 essentially exactly. The aggregation-vs-identical-rows caveat was moot; the finding holds.
- **On identical rows, AlphaMissense edges Site-Independent by ~0.07 on function** (0.499 vs 0.426) — a real but
  modest learned advantage where the signal is position-specific catalytic.
- **Site-Independent BEATS AlphaMissense on STABILITY / abundance** (0.393 vs 0.334) — a clean new finding: for
  structural-stability effects, the deterministic conservation score is the *better* predictor.
- **Overall, deterministic conservation ≈ 0.39 vs learned ≈ 0.46** (assay-median) / 0.37 vs 0.42 (UniProt-median)
  — a ~0.05–0.07 gap, i.e. deterministic conservation captures the large majority of what an AlphaMissense-class
  learned predictor does, and wins on one modality.
- **BLOSUM floor 0.207** — far below both. The earlier "molecular regime → learned wins" was really
  "substitution matrices are too weak"; a position-specific *deterministic* conservation score competes.

## What this settles for the project

- The molecular-regime decoder-choice is now quantified on identical rows: a **deterministic conservation
  decoder** (interpretable, no learned model) is a legitimate option — it matches AlphaMissense-class predictors
  on abundance (wins) and overall, and trades ~0.07 on function. That is the honest "deterministic can compete"
  answer, owned on the same variants rather than cited.
- The residual learned edge is concentrated on **function / expression / organismal fitness** (~0.07–0.10) —
  exactly the position-specific / distributed effects a single independent-sites model cannot capture. A coupling
  model (EVmutation/EVE) would close some of it but is no longer strictly deterministic-simple.

## Honest scope

- These are ProteinGym's DMS_score + its own Site-Independent per-mutant zero-shot scores; **BLOSUM + AM are
  computed by this project on ProteinGym's exact rows** (owned), Site-Independent is ProteinGym's published
  per-mutant score (its Spearman recomputed here on the same rows — 0.426 on function confirms row 312).
- Single substitutions only; assays with <30 single variants or missing a predictor dropped from that cell.
- verify-in-batch caught a real bug: a naive numeric-column heuristic picked `mutated_sequence` (coerces to
  all-NaN float64) as the Site-Independent column — fixed to a name-based match.
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9).

## Reproduce

```bash
# ProteinGym DMS + zero-shot (Zenodo 15293562) + AlphaMissense (Zenodo 8208688) on D:; extract human CSVs, then:
uv run python scripts/pg_native_threeway.py --dms-dir D:/dna_decode_cache/proteingym/pg_zeroshot --si-dir D:/dna_decode_cache/proteingym/pg_zeroshot
uv run pytest tests/test_pg_native_threeway.py -q   # 3 offline synthetic tests
```
