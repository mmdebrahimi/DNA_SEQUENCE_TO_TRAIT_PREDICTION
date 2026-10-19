# MPRA substrate census: three rungs, three ceilings, one of them forgotten

**2026-10-12.** Every figure below is MEASURED from the raw file on this date, not estimated and not
carried from a paper. Machine-readable twin: `wiki/mpra_substrate_census_2026-10-12.json`.

## Why this artifact exists

The GLM `ceiling:` criterion sat `blocked:untried-arm` on a declared 4-decade-span arm whose stated
blocker was **"no substrate with longer context; unblock = acquiring a substrate, a USER call."** That
blocker was wrong in a specific and instructive way: a **321,123-fragment library at 48-499 bp** was
already cached in the SAME directory as the 150 bp tiles the project trains on, and
`dna_decode/glm/genomewide.py` even mentions it in a comment as a thing a caller might accidentally mix
in. Nobody had asked whether it was the substrate the arm needed.

**The honest scope of that error.** The LITERAL 4-decade arm is still substrate-gated: four decades of
span needs a ~10,000 bp window and the fragment library tops out at 499 bp. So the blocker was NARROWER
than stated, not simply false. What was wrong was the implication that no longer-context data existed.

**A substrate blocker is a claim about what you looked for.** This census is the durable record so the
next session reads measured properties instead of re-deriving them, or re-discovering them in a year.

## The three rungs

All three live in `D:/dna_decode_cache/mpra/`.

| rung | file | n | length (bp) | sequence-model ceiling | single-rep agreement |
|---|---|---|---|---|---|
| peak tiles | `GSE144621_peak_tile_expression_formatted_std.txt.gz` | 46,713 rows; **44,106** after the published filters | all exactly **150** | **0.9410** | 0.7945 |
| fragments, LB | `GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz` | **321,123** | 48-499, median 244, p90 290 | **0.7928** | 0.4583 |
| fragments, M9 | `GSE144621_U00096.2_frag-rLP5_M9_expression.txt.gz` | **318,457** | 48-499, median 244 | **0.8014** | 0.4731 |
| sigma70 designed | `GSE108535_sigma70_variant_data.txt.gz` (+ `..._barcode_mapping.txt.gz`) | **10,898 usable** (join rate 1.00) | all exactly **150** | not computable — the two published tables correlate at exactly 1.0 (one monotone transform of the other, NOT replicates) | — |

**CORRECTED 2026-10-13 — the first version of this table reported the wrong column as the ceiling.**
`noise_ceiling()` returns two different quantities and they were conflated: `single_replicate_agreement`
(how well rep1 tracks rep2) and `ceiling_for_a_sequence_model` (the Spearman-Brown step-up of that,
square-rooted — because the published readout is the MEAN of two replicates, which is more reliable than
either alone). The fragment rung's ceiling is **~0.79-0.80, not ~0.46-0.47**, and the consequence is
strategic, not cosmetic — see rung 2 below.

**Ceilings are per-rung and MUST NOT be inherited.** 0.941 belongs to the peak tiles alone; it is a
property of that library's selection-for-signal, not of "bacterial promoter data". Quoting 0.941 against
fragment data would be the shared-key error one level up.

## Rung 1 — peak tiles (the current substrate)

- 46,713 rows; the two published filters (`category == "tile"` and `len(seq) == 150`) leave **44,106**.
- **The size lever is EXHAUSTED here.** 46,713 is the entire file and 44,106 are already in use, so the
  lever the literature ranks first — training-set size, which improves log-linearly and is reported
  unsaturated even at the replicate ceiling — cannot be pulled on the 0.941 substrate. More data means
  leaving this rung, and that is what makes the fragment rung's lower ceiling a real cost rather than a
  footnote.
- Best measured arm on it: `conv-3width-13k-params-150bp` at Spearman 0.3234, which is **34% of 0.941**.
  A single GC number reaches 0.3112 and every learned arm's paired CI against GC includes zero.

## Rung 2 — the fragment libraries (longer context, much noisier labels)

- 98.84% of LB fragments are longer than 150 bp; 359 distinct lengths.
- **The trade is MILD, not steep — this reverses the first version of this census.** ~7.3x the rows and
  3.3x the median length, against a sequence-model ceiling that falls only from **0.9410 to 0.7928 (LB) /
  0.8014 (M9)**. The first draft read the single-replicate agreement (0.4583 / 0.4731) as the ceiling and
  therefore described a halving of headroom; the true figure is a ~15% reduction. **The direction of the
  recommendation flips with it:** 7.3x the data for 85% of the headroom is an attractive trade, where
  7.3x the data for 49% of it was not. The mechanism still holds — peak tiles were selected at expression
  peaks so their replicates agree more closely, while fragments are random genomic pieces — but averaging
  two replicates recovers most of what single-replicate noise costs.
- Max span at a 3 bp floor is **2.22 decades** here versus 1.70 on the tiles — so the multi-scale
  question CAN be asked wider on this rung, even though the literal 4-decade arm still cannot run.
- Two conditions (rich LB, minimal M9) over near-identical fragment sets make this also a
  condition-contrast substrate, which no current arm uses.

## Rung 3 — the sigma70 designed library (constructed variation)

> **RETRACTED AND REWRITTEN 2026-10-13. Every number in the first version of this section was wrong,
> from one cause.** It reported "only 3,821 (35.1%) join to a sequence" and concluded the 65% gap was
> *absent sequence data*. It is not: the expression table writes UP-element names with HYPHENS
> (`gourse-326fold-up`) and the barcode map writes them with UNDERSCORES (`gourse_136fold_up`).
> Normalising the separator joins **all 10,898 at rate 1.00**. The figure 3,821 is exactly what a naive
> join returns — and `dna_decode/glm/expression.py` already existed, already fixed this, already
> enforced a `MIN_JOIN_RATE = 0.95` floor, and `tests/test_glm_expression.py` already pinned
> `rep.naive_join == 3821` with the comment *"the naive-join counterfactual is the point of the fix"*.
> A second loader was written beside it and reproduced the exact failure the first one documents.
> **The canonical loader is `dna_decode.glm.expression.load_pairs`. There is no other.**

- **10,898 designed promoters, all 10,898 usable** — join rate **1.00**, zero dropped for length, zero
  for non-ACGT, 10,898 distinct sequences. All exactly 150 bp, so the existing encoder runs unchanged.
- Design is a near-complete factorial, **3 x 8 x 8 x 8 x 8 = 12,288 cells** over UP_element x Minus35 x
  Spacer x Minus10 x Background, with **10,898 filled (88.7%)**. `UP_element` has **three** levels
  (`noUP`, `gourse-136fold-up`, `gourse-326fold-up`) — the first version claimed it collapsed to one,
  because the naive join had kept only the `noUP` arm. **It is 5-way crossed, not 4-way.**
- The leakage-proof split is `expression.leave_element_out(pairs, column)`, which holds out entire
  element IDENTITIES on **any of the five axes** — the bacterial analogue of a held-out-chromosome
  split. Disjointness is BY CONSTRUCTION (membership is a pure function of the element label), so no
  sequence can land on both sides.
- Expression is continuous and uncensored: **10,898 distinct values, mode-share 0.0001, ~615-fold
  range** — clears the project's own `assay_degeneracy` bar (>25% mode-share / <20 distinct levels) by
  a wide margin.
- **Measured headline on this rung:** positional one-hot ridge reaches Spearman **~0.59** on held-out
  Minus10 while GC alone reaches **~0.14** (`tests/test_glm_expression.py`, pinned). Unlike the peak-tile
  rung, the crude baseline does NOT win here — that contrast is a finding, and it is why this rung is
  worth more than its row in the first draft suggested.
- **Why the leakage controls are non-optional on this rung:** a random forest on deliberately
  IRRELEVANT features scored 0.75 AUC on a published saturation-mutagenesis benchmark and
  "reproduc[ed] the prediction quality of the best CAGI challenge submissions", via leakage from
  adjacent-segment splits within one element. Specify the split and an irrelevant-feature control BEFORE
  fitting anything.

## Three parsing traps, all of which cost a run, and all three already solved in the codebase

0. **The sigma70 join loses 65% of the data while looking perfectly healthy** — hyphen-vs-underscore in
   the UP-element name. A naive join returns 3,821 of 10,898 and every downstream statistic computes
   cleanly on it. `expression.py::_norm_key` is the fix (and it PROTECTS the mutation arrow `->`, which a
   blanket `replace("-", "_")` would also rewrite); `load_pairs` enforces `MIN_JOIN_RATE = 0.95` so the
   failure cannot go quiet. **This trap was documented, fixed, and test-pinned BEFORE it was hit again** —
   see the retraction banner under rung 3.

1. **These files are SPACE-delimited, not tab-delimited** — and the project ALREADY KNEW. A
   `split("\t")` parse returns **zero rows** and reads as an empty dataset rather than failing; it
   surfaced only because the script crashed downstream on an empty list. **But the fix was already in
   the codebase:** `dna_decode/glm/genomewide.py::sniff_delimiter` decides tab-vs-whitespace per file
   and its own docstring says *"within GSE144621 the frag files are space-delimited and the peak-tile
   file is tab-delimited. A global assumption parses one of them to zero rows."* The real lesson is
   therefore not "these files are space-delimited" but **"a one-off parse re-learns a trap the shared
   loader already handles"** — reuse `sniff_delimiter` + `_split` rather than splitting inline.
2. **`load_peak_tiles` wants `D:/dna_decode_cache/mpra`, not its parent.** Pointed at
   `D:/dna_decode_cache` it raises `GenomeWideDataError: ... missing and downloading is disabled`, which
   reads like absent data rather than a wrong path.

## What this census does NOT establish

- No sigma70 replicate ceiling. GSE108535's two published tables correlate at **exactly 1.0** — one is a
  monotone transform of the other, not a replicate pair — so "what fraction of achievable" is **not
  computable** for rung 3 from these files. A number on this rung therefore has no denominator, which is
  a real limit on comparing it to the 0.9410 peak-tile results.
- ~~Whether the 65% of sigma70 variants lacking sequences are recoverable from another supplementary
  file.~~ **CLOSED 2026-10-13 — the question was an artifact of a broken join.** Nothing is missing;
  all 10,898 join once the hyphen/underscore separator is normalised.
- Whether training on the fragment rung transfers to the tile rung. That is a different and possibly
  stronger design than treating the rungs as separate criteria, and nothing here tests it.
