# Moving substrate does NOT establish relief from the Markov-chain defeat — and my primary geometry measured the wrong thing

**Status:** `DIAGNOSTIC`. Graded by `wiki/human_falsifier_bar_2026-11-21.json`, frozen before any
number below existed. Grid: `wiki/human_falsifier_margin_2026-11-21.json` (8/8 cells, committed
verbatim before this reading). Verdict recomputed by
`uv run python scripts/human_falsifier_margin.py --verdict-only` ->
`wiki/human_falsifier_verdict_2026-11-21.json`.

## DECISIONS NEEDED

**One, and it is about my own bar rather than about the result.** The bar named the FULL 230 nt
sequence as the PRIMARY geometry, for comparability with the measured activity labels and the
2026-11-01 encoder arm. The data shows that choice measures the wrong thing for this question (see
CORRECTIONS). I have **not** swapped primary and secondary, because swapping a frozen choice after
seeing which is more favourable is precisely what a pre-registration exists to prevent. Whether to
re-freeze a successor bar with the 200 nt element as primary is a scope decision and it is yours.

## CORRECTIONS

**The frozen PRIMARY verdict is real and uninformative, and the grid proves it.** On the full 230 nt
sequence the reading is `CHAIN_IS_WEAKER_ON_HUMAN`, unanimous across both cells and both feature
modes. But the publisher's sequence is `15nt 5' adaptor - 200nt element - 15nt 3' adaptor`, and
stripping those constant flanks collapses the null's distinguishability while barely moving its
noise floor:

| cell / mode | null, full 230 | null, element 200 | null drop | ceiling move |
|---|---|---|---|---|
| HepG2 / both | 0.9591 | 0.5295 | **-0.4296** | +0.0279 |
| HepG2 / kmer | 0.8757 | 0.5106 | **-0.3651** | +0.0054 |
| K562 / both | 0.9396 | 0.7495 | -0.1901 | +0.0092 |
| K562 / kmer | 0.8088 | 0.7711 | -0.0377 | +0.0015 |

On HepG2 almost the entire signal was the adaptors: the null falls to 0.5295 — essentially its own
floor of 0.5388 — once 15 bases are removed from each end.

**That drop alone is suggestive but underdetermined, so the mechanism was MEASURED rather than
asserted** (`scripts/human_adaptor_probe.py` -> `wiki/human_adaptor_probe_2026-11-21.json`,
`ADAPTOR_ARTIFACT_MECHANISM_CONFIRMED`). The first draft of this memo claimed the flanks were
constant on the strength of a column description; both halves are now measured:

- **The flanks are perfectly constant.** Across 3,000 natural sequences per cell the 5' 15 nt is
  **1 distinct value** (`AGGACCGGATCAACT`, share 1.0000) and the 3' 15 nt likewise
  (`CATTGCGTGAACCGA`), identical in both cell types. Per-position max-base share is **1.000** on both
  flanks against **0.262** (HepG2) and **0.328** (K562) in the interior.
- **The Markov chain never reproduces them.** Of 120 generated sequences per cell, **0 match either
  flank exactly**; the mean matching 5' prefix is 3.69 of 15 bases (max 8) for HepG2 and 3.62 (max 7)
  for K562.

**This refutes the obvious rival explanation.** Stripping removes 30 of 230 positions (13%), so the
weaker separation could in principle be a shorter-sequence effect rather than anything about
adaptors. It is not: positions that are **identical in every natural sequence** carry zero natural
variability, so they cannot contribute to a discriminator's ability to characterise natural sequence
— which is exactly why the self-control's floor barely moves. They can only supply the
constant-versus-not signal a stochastic generator fails.

So a 3-mer Markov chain cannot reproduce two constant flanks, the discriminator reads that, and
**none of it is about regulatory structure.** The bar itself flagged this risk in advance (*"constant flanks are trivially learnable
and could inflate distinguishability for a reason that is not about regulatory structure"*) and then
made the flanked geometry primary anyway. That is an internal inconsistency in my own
pre-registration, exposed by the measurement rather than by review.

## The result

Each row compares the human **within-substrate margin** (`null - self_control_ceiling`) against the
committed bacterial margin, with the tolerance being the human self-control's own measured seed
spread. Raw levels are never compared across substrates.

| geometry | mode | HepG2 margin | K562 margin | bacterial | verdict |
|---|---|---|---|---|---|
| full 230 (frozen PRIMARY) | kmer | +0.3742 | +0.2633 | +0.0543 | `CHAIN_IS_WEAKER_ON_HUMAN` |
| full 230 (frozen PRIMARY) | both | +0.4482 | +0.4122 | +0.0814 | `CHAIN_IS_WEAKER_ON_HUMAN` |
| element 200 (secondary) | kmer | **+0.0036** | **+0.2241** | +0.0543 | `CHAIN_IS_COMPARABLE` (cells contradict) |
| element 200 (secondary) | both | **-0.0093** | **+0.2129** | +0.0814 | `CHAIN_IS_COMPARABLE` (cells contradict) |

**The geometry that actually bears on the question returns `CHAIN_IS_COMPARABLE`, and it does so
because the two cells of the same panel CONTRADICT each other.** HepG2 says
`CHAIN_IS_STRONGER_ON_HUMAN` (its margin is *below* bacterial, by more than its seed spread); K562
says `CHAIN_IS_WEAKER_ON_HUMAN`. The bar's conservative resolution — fixed in advance, so a split
could not be settled by picking the agreeable cell — reports `CHAIN_IS_COMPARABLE` and records the
contradiction.

## What this does and does not answer

**The original question was: is the Markov-chain defeat a property of the bacterial substrate or of
generators?** On the geometry that isolates regulatory structure, **it is not established either
way.** Two cell types from one panel, sharing an assay and an adaptor design, disagree in opposite
directions by more than their own measured noise.

**One striking row, reported without being leaned on.** On HepG2's 200 nt elements the null sits
essentially AT its own noise floor (margin +0.0036 for kmer, **-0.0093** for both) while the uniform
control still reaches 0.996 — so the discriminator works and a 3-mer Markov chain is nonetheless
near-indistinguishable from real human regulatory elements there. If that generalised it would mean
the bacterial 0.1186 gap is not a generator-universal property. **K562 contradicts it directly**, so
it does not generalise within this panel, let alone beyond it.

**The cell disagreement is UNEXPLAINED and is not explained away.** Both cells share the assay, the
adaptors and the fold scheme. Their libraries differ in size (164,307 designed for HepG2 against
243,778 for K562) and presumably in composition, but nothing here establishes that as the cause.

**This says nothing about whether a generator would succeed on human.** No generator was run, by
design — the current one is `GENERator-v2-prokaryote-1.2b-base` and a human arm would need a new base
model, confounding substrate with model family. This characterised the OPPONENT only.

## Instrument validity

The uniform control reaches **0.978-1.000** in all 8 cells, so the discriminator separates random
sequence from natural everywhere and no row is `INSTRUMENT_INVALID`. Discriminator capacity was held
at the bacterial run's settings (k=4, bins=6, 120 generated, 1200 natural, seed 7, Markov order 3)
because the margin involves the self-control, so changing capacity would move both terms. The
fit/eval split is the publisher's own 10-fold assignment with fold 1 held out, disjointness asserted
on content rather than inferred from fold ids.

## Honest limits

- One human panel (Agarwal 2025 lentiMPRA), two cell types, 230 bp. ONE substrate change, not a
  substrate sweep: the panels differ in assay, length, kingdom and library construction, and nothing
  here says WHICH difference carries any effect.
- No human replicate ceiling exists (`wiki/human_frame_shape_2026-11-16.json`: no per-replicate
  columns, so no Spearman-Brown step-up). The bacterial 0.941 is NOT borrowed, so the two substrates'
  results are not expressible on a common normalised scale even as margins.
- A k-mer-plus-positional discriminator sees LOCAL composition while the published generator failures
  name long-range organisation, so any low score remains necessary-NOT-sufficient.
- The published incumbents are unavailable as a comparator
  (`wiki/lentimpra_metric_resolution_2026-10-28.md`, `UNRESOLVED_TABLE_IDENTITY_UNVERIFIED`).
- One eval fold (fold 1) and one generation seed per cell. The tolerance is the self-control's seed
  spread, which bounds discriminator-split noise but not generation-sampling noise on the null.
- The forecast recorded 0.5 for `CHAIN_IS_WEAKER_ON_HUMAN`; the secondary geometry returned
  `CHAIN_IS_COMPARABLE`, so the forecast was neither confirmed nor refuted.
- The adaptor probe explains the **full-230 rows only**. It says nothing about the element-200 rows,
  where the cells genuinely disagree, and it is not evidence about generation.
