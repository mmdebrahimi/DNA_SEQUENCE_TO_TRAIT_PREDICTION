# The human lentiMPRA panel is acquired, verified, and screened (2026-10-27)

**What this closes.** `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` was recorded as
blocked on data acquisition. The 2026-10-27 analysis
(`wiki/ceiling_arm_decision_2026-10-27.md`) found that wall was a property of our drives, not of the
world. This memo records the fetch, what it actually contains, and the two findings that change the
arm's shape. **The arm is NOT run; nothing is claimed about its performance.**

Source: **Zenodo 10558183**, Agarwal et al., *Massively parallel characterization of transcriptional
regulatory elements*. **Licence CC-BY-4.0, open access** — read from the API, not assumed, and
`scripts/fetch_lentimpra.py` REFUSES if the licence ever reads as anything else.

## CORRECTION to a figure this project was carrying

Earlier in the session I reported the panel as **">680,000 sequences"**, taken from the paper's
abstract. Measured from the deposit itself:

| quantity | measured |
|---|---|
| designed sequences (Table S3, three sheets) | 483,627 |
| rows carrying measured activity (Table S6) | 182,859 |
| **joinable (sequence + activity + fold)** | **159,795** |

The usable labelled set is **~4.3x smaller** than the figure quoted. 680,000 is **not reproducible
from these tables** and this memo does not claim to know what it counts (plausibly all libraries and
both orientations, unverified). **Quote 159,795.**

Element length is exact rather than approximate. The sequence column reads verbatim:

> `230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)`

So **200 nt of biology** inside a 230 nt oligo. The loader returns the sequence AS PUBLISHED and
records that column name, because a silent adapter trim would change what every downstream number is
about.

## What was fetched, and what was deliberately not

Four supplementary tables, 74.19 MB, **all four MD5-verified** against the publisher's manifest,
landed at `D:/dna_decode_cache/mpra_human/` (the D: rule; C: was at 88%, D: had 908 GB free).
Namespace-separate from `D:/dna_decode_cache/mpra/` so a human table can never be mistaken for a
GSE144621 *E. coli* one.

**The two 3.5 GB / 3.97 GB dumps were NOT downloaded, on measured grounds.** Their zip central
directories were read over HTTP **range requests** (Zenodo returns 206), so their contents were
listed without transferring them:

| file | entries | model checkpoints | data |
|---|---|---|---|
| `final_dump.zip` | 303 | **270 files, 3.83 GB** | ~10 MB |
| `final_joint_dump.zip` | 303 | **270 files, 4.33 GB** | ~10 MB |

Both are ~99.7% trained LegNet weights, and the ~10 MB of TSVs in them are the authors' own **model
predictions** (`scores/final/<cell>/.../testfold_N.tsv`), not the measurement library. So 7.5 GB
would have bought weights. `--select dumps --yes` remains available; the default does not.

## The fetcher's verification is the point of it, and it fired on the first real run

A first hand-rolled `curl` of two tables wrote **a 20,899-byte Zenodo "application-unavailable" page
and a 150-byte nginx "502 Bad Gateway" page, under genuine `.xlsx` filenames**, because Zenodo was
returning 502/503/504 intermittently. Fed to a parser those are indistinguishable from "this table
has no usable rows". This repo has already paid for that exact shape once: the cached GPS assembly
cohort was 100% unusable (HTTP 403 pages and truncated gzips under genome filenames).

Both corrupt files were deleted and `scripts/fetch_lentimpra.py` was built with **whole-file MD5**
verification, retries with backoff, and deletion of anything that fails. On the real run, **Table
S3's first attempt returned the 504 page, was caught by the size check, deleted, retried, and
verified on attempt 2.**

**A HEAD request or an opening-bytes sniff would be the WRONG guard**, and this is measured, not
stylistic: a truncated gzip decompresses its opening blocks perfectly and only raises at the END, so
reading the first 2 KB passes exactly the files that later fail.

## Two findings that change the arm's shape

### 1. WTC11 does not join — the intersection is EMPTY, not small

| cell | designed | measured | joined | rate |
|---|---|---|---|---|
| HepG2 | 164,307 | 61,463 | 61,463 | **1.0000** |
| K562 | 243,778 | 98,332 | 98,332 | **1.0000** |
| WTC11 | 75,542 | 23,064 | **0** | **0.0000** |

Its Table S6 activity rows are keyed by **Ensembl gene id** (`ENSG00000000003`); its Table S3 design
rows use **element names** (`seq1000_F`). Different entity types. Table S7 (the joint library)
rescues only **699 of 23,064 (3.03%)**, because S7's ENSG-style names are its 2,396 `promoter` rows
and nothing else; a `WTC11_` prefix join gives **0**.

`load_cell("WTC11")` **REFUSES** rather than returning the overlap, for the reason
`expression.load_pairs` enforces `MIN_JOIN_RATE` after the sigma70 hyphen-vs-underscore gap returned
3,821 of 10,898: a 3%-complete cohort has the same shape as a complete one. Relaxing the bar is
tested to NOT manufacture rows — the intersection is empty at any threshold.

**Open, not closed:** WTC11's 23,064 measured elements are real and their sequences exist somewhere
(GEO, or a table not in this deposit). Recovering them would add 14% to the labelled set.

### 2. A published incumbent benchmark exists on the IDENTICAL folds

Table S6 ships the authors' own **10-fold assignment** (near-equal folds, carried through by the
loader) plus per-fold scores for five methods. That means a leakage-safe split is handed to us AND
any result is directly comparable without re-deriving a split:

| method | HepG2 (mean of 10) | K562 | WTC11 |
|---|---|---|---|
| Biochemical (cell-type matched) | 0.7320 | 0.7235 | 0.7124 |
| Biochemical (universal) | 0.7365 | 0.7404 | 0.7448 |
| MPRAnn | 0.7665 | 0.8082 | 0.7572 |
| Sei | 0.7914 | 0.8150 | 0.8212 |
| Enformer | 0.8306 | 0.8497 | **0.8335** |
| **MPRALegNet** | **0.8333** | **0.8665** | 0.7946 |

**The metric is NOT stated in the sheet.** `published_incumbents()` returns
`metric="UNSPECIFIED_IN_SOURCE"` and a test pins that, because comparing a Spearman of our own
against an unknown readout is a readout-type mismatch rather than a comparison. Resolving it (from
the paper's methods) is a prerequisite to any numeric claim against this table.

Worth noting for later: **MPRALegNet beats Enformer on HepG2 and K562 but LOSES on WTC11** — the one
cell type we cannot currently load. That asymmetry is pinned by test so it is not smoothed away.

## Screens

**Ten rejection gates: `CLEARS`** — 5 pass, 5 not-applicable, 0 trip
(`wiki/lentimpra_gate_screen_2026-10-27.json`). G1 and G3 are the judgment gates and carry written
evidence; the screen refuses a verdict without it. All four previously-committed hand verdicts still
re-derive, so this packet caused no regression.

**G6 (assay degeneracy) passes decisively**, measured on the real activity values and screened on the
WORSE of the two cells:

| cell | distinct levels | mode share (bar <=0.25) | at min / at max |
|---|---|---|---|
| HepG2 | 7,192 | 0.000586 | 1 / 1 |
| K562 | 5,494 | 0.000814 | 1 / 1 |

Exactly one value at each extreme — the opposite of the CcdB failure mode (79.3% tied at a ceiling,
which posted the forward/inverse sweep's best number for that reason).

**`sd` is exactly 1.0000 in both cells**, so `Observed log2(RNA/DNA)` has been **standardized**, not
left as a raw ratio. Recorded because it decides what an absolute margin would mean and makes
cross-cell comparison of raw values suspect.

## Substrate consequences

- **Uniform 230 bp**, all 483,627 designed sequences. So this panel needs **none** of the masked
  pooling the *E. coli* fragment rung required — `encoder.py`'s uniform-length contract applies
  directly, and `require_uniform_length` will not refuse it.
- **Span:** `log10(200/3) = 1.82` decades on the biological element, or `log10(230/3) = 1.88`
  including adapters. Comparable to the tile rung's 1.70, below the fragment rung's measured 1.91
  median. **This panel is not a span lever**; its value is the SUBSTRATE CHANGE (human, constructed,
  with published incumbents).
- **Regime:** constructed variation, learned representation, molecular-adjacent readout — the regime
  `eval/regime.py` records as WORKING, not the closed natural-population zero-shot negative.

## Limits

Acquisition and screening ONLY. **No arm was run and no bar was frozen**; a gate PASS bounds whether
a usable label exists, never whether the work is worth doing. The 680,000 figure is unexplained
rather than refuted. WTC11 is unresolved. The published benchmark's metric is unknown, so no number
here may be compared against it yet. The two 7.5 GB dumps were characterised by reading their
central directories, not by extracting them, so the claim about their contents is a claim about
their file listing.

## Reproduce

```bash
uv run python scripts/fetch_lentimpra.py --list            # manifest + MD5s, no download
uv run python scripts/fetch_lentimpra.py --select tables   # 74 MB, verified
uv run python scripts/screen_candidate_gates.py --verify   # the 4 committed hand verdicts
uv run pytest tests/test_lentimpra.py -q                   # 22 tests, ~2m12s
```

`openpyxl` is required and was installed **venv-only** via `uv pip install`, leaving the frozen
`uv.lock` byte-identical (md5 `892c08339f82bd1a23f898821adb9735` before and after). This expires
CLAUDE.md's recorded reason for declining it — *"the workbook it would parse is itself absent"* — a
constraint that has now lifted (R9).
