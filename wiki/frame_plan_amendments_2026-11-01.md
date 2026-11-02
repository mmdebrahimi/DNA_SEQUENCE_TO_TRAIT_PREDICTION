# The loader regression is fixed, and the ordering divergence turned out to be real

<!-- doc-status: active -->
<!-- measured-on: 2026-11-01 -->

<!-- Execution record for plans/Clear_The_Three_Blockers_The_Frame_Conditioning_Review_Found_Plan/.
     Artifacts: wiki/frame_plan_window_ordering_2026-11-01.json (new, from
     scripts/frame_window_ordering_probe.py) and wiki/frame_plan_feasibility_2026-11-01.json
     (pre-existing, deliberately untouched). Review and dispositions in that plan folder's
     brainstorm.md. -->

## DECISIONS NEEDED

1. **Should the guard's new exemption row be backed by a test that pins the mention count?** The row
   added to `ALLOWED_MENTIONS["GSE144621"]` carries a comment recording that `tile_controls.py` has
   exactly one prose mention of the accession. A comment cannot fail, so if that count grows nothing
   breaks. A test asserting the exact occurrence count would convert the observation into a guard.
   Not done, because it adds a test the plan did not ask for.
2. **Three questions the frame plan still carries, unchanged by this work.** Which frame granularity
   is primary (`in_window_any` at 66.47% versus `in_window_full` at 19.13%); whether
   `arms.native150 = 0.5076` stays at its measured value if the frame turns out to explain it; and
   whether the 689 sequence-duplicate rows should be excluded from the rung's splits, which would
   make every committed number on that rung incomparable.

## CORRECTIONS

**The ordering divergence is real, and it is now a number instead of an assertion.** The frame plan
claimed to reproduce the upstream `extract_promoter_regions` exactly. That was a contradiction,
because the upstream function threads a running `last_end` through gene features in FILE order and
the cached annotation's coding-sequence rows are not ascending by start on either strand. Measured
by `scripts/frame_window_ordering_probe.py`:

| ordering | windows | bases covered |
|---|---|---|
| file order | 3,490 | 624,812 |
| sorted ascending | 3,496 | 625,255 |

The symmetric difference is **26 windows**, and **146 of 43,061 uniquely-located tiles change their
in-window status**. Small — 0.74% of windows and 0.34% of tiles — but non-zero, so the claim of exact
reproduction could not stand. Row ordering is now an explicit parameter defaulting to sorted, every
pinned figure is labelled with the ordering it came from, and the resulting claim is scoped to *the
corpus rule as specified at 200 bp with a 5 bp intergenic gate* rather than *as the authors ran it*.
A further limit is recorded in the artifact: the authors read **GenBank**, not this GFF3, so file
order here is a different file's order and faithful reproduction is not well-defined from this input.

**A superseded figure replaced in four places.** Stated-coordinate agreement was quoted as
`63/2000 = 3.15%`, a first-2,000-row spot check. Measured across every uniquely-located tile it is
**1193/43061 = 2.77%**. Same conclusion — a coincidence rate — better number. Corrected in the frame
plan's Codebase Context, Pre-Change Baseline and Risk Flags, and in the index bullet that restated it.

**An index bullet that would have outlived the claim it restates.** `wiki/plans-index.md` asserted
that the frame plan's Step 1 "replicates the authors' `extract_promoter_regions` exactly". The index
is what a reader scans to see what a plan decided without opening it, so a false bullet there is
worse than a false sentence in the body. Rewritten in the same step that changed the body, which is
the standing rule that a correction is finished only when every surface restating the claim carries it.

## What actually got fixed

**The red ownership guard is green.** `tests/test_substrate_loader_ownership.py` went from 1 failed
and 7 passed to **8 passed**. `dna_decode/glm/tile_controls.py` no longer opens the substrate; it
calls `genomewide.load_peak_tiles` and filters the result, so it is a consumer of the registered
owner rather than a second reader. The guard exists because the fragment files in that series are
space-delimited while the peak-tile file is tab-delimited, and a second loader inherits none of the
owner's delimiter sniffing.

**Routing changed no data, which is why the fix cost nothing.** Through the owner the counts are
**936 `random` and 470 `neg_control` at 150 bp, total 1,406** — byte-identical to what the second
parser produced.

**The fix needed both acts, in order, and that was verified rather than assumed.** With the parser
deleted the guard **still failed**, on the single surviving docstring mention, because its criterion
is the bare text pattern `\bGSE\d+` over whole files. Only then was the exemption row added, which is
the evidence that the row does real work rather than covering a live reader. An exemption alone would
have been the wrong fix; a deletion alone was insufficient.

**Seven test fixtures were widened rather than seven assertions relaxed.** The owner enforces a
minimum of eight columns and requires `peak_start`, `peak_end` and `strand`, which the five-column
fixtures lacked. The review's estimate of "roughly a ten-line change" was wrong for exactly this
reason, and the fixture rebuild was the bulk of the work. Two tests were added: a syntax-tree
assertion that no import names `gzip` and no attribute targets it, and a pin on the deliberate
narrowing of the signature from "file or directory" to "directory only".

## Verification

- Affected-file subset: **160 passed**. Ownership guard **1 failed → 8 passed**. Control loader
  **9 → 11 tests**, all passing with widened fixtures and unrelaxed assertions. New ordering probe:
  **8 tests**, including a non-vacuity check that the two orderings genuinely differ on a constructed
  case, and pins on both replicated upstream quirks that fail if a future edit tidies them.
- `scripts/contract_number_audit.py` still returns `NOTHING_TO_ADJUDICATE 0/215 across 25 cells`.
- The frozen five, `uv.lock`, `genomewide.py` and `expression.py` all show an empty diff.
- `scripts/promoteratlas_arm_dryrun.py` still passes its embedded-payload byte-identity assertion.
- In-package importers of `genomewide` are exactly `atlas_embed.py`, `encoder.py` and the newly
  routed `tile_controls.py` — three, up from two, and the third is the intended addition.
- The amended frame plan validates clean and its wave grouping is unchanged at
  `[[1,2,4],[3],[5],[6],[7]]`, with all 7 steps and all 7 named outcomes surviving.

## Honest limits

The frame-conditioning experiment itself has **not** run. This work made its plan executable and
cleared a regression that was blocking nothing but trust; the three cells, the acceptance bar and
every measured question remain ahead. The ordering probe bounds the SHAPE of the divergence on this
GFF3, not its exact magnitude in the corpus the authors actually built. And the guard's exemption is
per-module, so substrate-reading code added to `tile_controls.py` later will not be caught — a
residual the guard's own comments already name for its two other exempted modules.
