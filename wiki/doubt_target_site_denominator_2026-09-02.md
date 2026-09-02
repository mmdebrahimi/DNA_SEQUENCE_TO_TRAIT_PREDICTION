# The doubt layer's completeness screen transfers to the target-site arm — one vocabulary, not two

**Open question (F-A candidate 2): "extend the completeness screen to the target-site catalogs — one
vocabulary or two?" Answer: ONE.** The purity signature is well-formed on the target-site arm and it
fires, naming exactly one candidate after a family-wise correction over 638 units.

Measured rather than built, because the same question has twice been settled by a census here: the
per-cell regime field looked like ceremony until the census showed the column was not constant, and NNRTI
curation looked obvious until its recovery was measured against the incumbent.

---

## The shape test

The AMR screen ranks determinant **families the deployed rule cannot represent** by purity: carriers
labelled R and never S. That signature needs three things, or it reports structure that is not there:

| condition | why it matters | HIV NNRTI / EFV |
|---|---|---|
| a negative class | "zero S carriers" is otherwise true by construction | **yes** — 1,170 S of 2,168 |
| units the rule cannot already represent | else the flag can never fire on anything actionable | **yes** — 638 substitutions outside the 8 catalogued positions with ≥5 carriers |
| purity must discriminate | if every unit is pure, purity separates nothing | **yes** — only 12 of 638 are pure |

Run on the **best case** deliberately: HIV NNRTI is the one target-site cell with a free, independent,
isolate-level wet-lab label (Stanford PhenoSense fold-change). A failure there would have generalised
downward and settled the question as "two vocabularies" for a few minutes' work.

## What it found

**`V179F` — 15 carriers, all resistant, p = 8.8 × 10⁻⁶**, the only unit of 638 surviving the family-wise
correction. Position 179 is **not** among the deployed catalog's positions (100, 101, 103, 106, 181, 188,
190, 230).

**Independent corroboration inside the repo:** the 2026-09-01 curation measurement — a completely
different method, multivariate OLS on log fold-change — named `V179D` (12 carriers) and `V179E` (3) among
its blind-spot drivers. Position 179 surfacing twice by unrelated routes is worth more than either alone.

## What this does NOT license

**It is not a curation recommendation.** Data-derived NNRTI curation was measured three ways on
2026-09-01 and **declined**: every variant recovered less of the blind spot than the free position-novelty
flag already does, and the best-scoring one deleted canonical `Y181C`. Nothing here reopens that.

The screen is **L2**, and L2 qualifies a call without competing with it. The actionable form of this
result is that an isolate carrying `V179F` and receiving a susceptible call should carry a doubt flag
saying that call is the least trustworthy — not that the L1 catalog should change. That distinction is
the layer's entire design claim, and it is what makes a positive here safe.

## Honest limits

- **One cell.** Measured on HIV NNRTI only, chosen as the best case. A negative would have generalised
  downward; this positive does not generalise upward, and each other target-site cell needs its own check.
- **Different measurement from the AMR arm.** R/S here is a threshold on continuous fold-change at the
  sourced `DRMcv.R` cutoff (3.0); the AMR arm uses categorical AST. Same vocabulary, not the same label.
- **In-distribution.** Catalog and label both trace to Stanford HIVDB, so this is not an independent test
  of the catalog — it is a screen for what the catalog does not cover.
- **15 carriers is small.** The p-value survives correction over the 638 units actually tested, which is
  the right denominator, but the finding rests on 15 isolates.
- **The screen was not wired into the shipped doubt block for this arm.** This probe answers the design
  question; wiring it is separate work with its own augment-only verification.

## Reproduce

```bash
uv run python scripts/doubt_target_site_denominator_probe.py
```

Reads the gitignored Stanford dataset and skips cleanly when it is absent. Touches no frozen file.
