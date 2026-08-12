# The yeast essentiality number was a medium-mismatch artifact (2026-08-11)

Follow-on from `wiki/fba_gap_premise_2026-08-11.md`, which recommended setting the organism's standard
medium **before** building anything. Doing it changed a published metric.

## The mismatch

Gene essentiality is **medium-dependent**. An amino-acid biosynthesis gene is essential on minimal medium
and dispensable on rich — that is biology, not model error. So the medium the model runs on has to match
the medium the labels were measured in.

It did not:

| | condition |
|---|---|
| **SGD gold standard** (deletion collection, inviable nulls) | **YPD — rich** |
| **iMM904 default medium** | glucose **minimal** — 10 open exchanges, **all 31** amino-acid / nucleobase exchanges closed |

A minimal-medium model *must* call biosynthesis genes essential that are dispensable when the nutrient is
supplied. The prediction from that direction is specific: the errors should be **false positives**, and
they should shrink under a rich medium.

## Result — the prediction holds

Same model, same labels, same gene set. Only the medium differs:

| medium | wt growth | TP | **FP** | FN | TN | **MCC** | precision | recall | accuracy |
|---|---|---|---|---|---|---|---|---|---|
| minimal *(model default)* | 0.288 | 43 | **67** | 92 | 703 | **0.2524** | 0.391 | 0.319 | 0.824 |
| **YPD-like rich** *(label-matched)* | 0.984 | 34 | **13** | 101 | 757 | **0.3773** | **0.723** | 0.252 | 0.874 |

**False positives fall 81% (67 → 13); precision nearly doubles (0.39 → 0.72); MCC rises ~50% relative.**
The discrimination tier moves **WEAK → MODERATE**. No new biochemistry — a config change.

Recall dips slightly (0.319 → 0.252, FN 92 → 101), which is expected: supplying nutrients also rescues a
few genes the model previously got right for the wrong reason.

## Why this matters beyond one number

**The remaining error is now almost entirely false negatives** — 101 of 114. That is exactly the error
class gap-filling targets, so fixing the medium first does not just improve the score, it *cleans the
signal Track C is aimed at*. Had the function-prediction build run first, it would have been credited with
the 54 false positives the medium fix removes for free.

This is the sequencing the premise check recommended, and it is now measured rather than argued.

## What shipped

- `dna_decode/fba/medium.py` — `RICH_MEDIUM_EXCHANGES` (20 amino acids + 6 nucleobases/nucleosides),
  `rich_medium()` / `apply_rich_medium()`. **Additive**: the reconstruction's own carbon source and oxygen
  bounds are kept and only supplements are opened. Exchanges a given model lacks are skipped, not raised,
  so it ports across reconstructions.
- `dna_decode/fba/essentiality_labels.py` — `ESSENTIALITY_LABEL_CONDITION` records the growth condition
  **with the label source**, since the gold standard is what knows how it was measured.
- `scripts/fba_essentiality_validate.py` — `--medium {label_matched,default,rich}`, defaulting to
  `label_matched`. The artifact records `medium_mode`, `medium_condition` and the supplement count, so no
  future reader has to guess which medium produced a number.
- An organism with **no** recorded condition falls back to the model default — it does not guess a medium.

## Limits

- The rich medium is **YPD-*like***, not a calibrated YPD: uniform 10 mmol/gDW/h uptake for every
  supplement, no vitamins/lipids, and no attempt to match yeast-extract composition. A curated YPD would
  likely shift the number again.
- Only yeast has a recorded label condition, because it is the only cross-organism cell that is SCORED.
  E. coli/Keio (`fba_keio_validate.py`) is a separate path and was not touched.
- The superseded 0.252 is kept in `wiki/fba_essentiality_yeast_2026-08-03.json` as history; the current
  number is `wiki/fba_essentiality_yeast_2026-08-11.json`.

## Reproduce

```bash
uv run python scripts/fba_essentiality_validate.py --organism yeast          # label-matched (rich)
uv run python scripts/fba_essentiality_validate.py --organism yeast --medium default   # the old number
```

Tests: `tests/test_fba_medium.py`.
