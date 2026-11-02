# Pre-execution review — Frame Conditioning plan

<!-- doc-status: active -->
<!-- measured-on: 2026-11-01 -->

<!-- Pre-execution review of `technical-plan.md` in this folder, run after the plan was saved and
     before any implementation step wrote a file. Evidence artifact:
     `wiki/frame_plan_feasibility_2026-11-01.json` (re-runnable via
     `scripts/frame_plan_feasibility_probe.py`). SINGLE-MODEL REVIEW -- the second-opinion tool was
     unavailable (`codex` exited 127, not on PATH; LM Studio UNREACHABLE and
     `lmstudio-auditor.sh` absent), so this carries anchoring-bias risk: the reviewer wrote the
     plan. Two of the four critical findings were established by measurement rather than reasoning,
     which is the partial mitigation; the rest should be re-reviewed independently when a second
     opinion is available. -->

## DECISIONS NEEDED

1. **Faithful-but-ill-defined replication vs an explicitly chosen rule** (critical issue C2 below).
   My lean: an explicit sorted rule with the divergence documented, and the claim scoped to *"the
   corpus rule as specified at 200 bp / 5 bp"* rather than *"as the authors ran it"*.
2. **Whether C1 is fixed inside this plan or as a separate change.** It is a red guard inherited
   from the preceding `--advance` run, not caused by this plan.
3. **Whether it matters that this study cannot lower the headline number.** Every outcome either
   leaves `arms.native150 = 0.5076` intact or adds a caveat, so it bounds the *interpretation* and
   cannot function as a check on the value. Probably fine for a limit study, but it should be
   decided deliberately rather than discovered later.

## CORRECTIONS

**A regression I introduced in the preceding run is live, and my "142 tests pass" claim did not
cover it.** See C1. The claim was true of the five test files I chose to run and the guard that
fails was outside that set — the standing lesson that *a zero-regressions claim inherits its
suite's scope*, demonstrated on myself.

**Two figures in the plan are superseded by the full-scale probe.** The plan cites stated-coordinate
agreement as `63/2000 = 3.15%`, a first-2,000-row spot check; measured over all 43,061
uniquely-located tiles it is **1193/43061 = 2.77%**. Same conclusion (a coincidence rate), better
number. Trivial-overlap share is **2,359 (8.2% of `in_window_any`)**, not the 2,540/9.0% quoted from
the rounded histogram.

**The plan's sentrux anti-drift check is unenforceable as written** and was corrected during
execution: `.sentrux/baseline.json` holds nine aggregate scalars (`total_import_edges: 2244`,
`cycle_count: 3`, …) with no per-edge record, so a compare cannot distinguish an unintended new
import edge into `genomewide.py` from an intended one — and four new leaf modules move the edge
count regardless. Replacement check: assert the set of files under `dna_decode/` carrying an actual
`genomewide` **import** is exactly `{atlas_embed.py, encoder.py}` (verified: the only two, at
`atlas_embed.py:169` and `encoder.py:333`, both importing `spearman`; two further files match a text
search for the name but only in docstrings).

## 1. Critical issues

### C1 — a pre-existing guard is RED, and it is my defect `[grounded]`

`tests/test_substrate_loader_ownership.py::test_each_substrate_accession_is_read_by_exactly_ONE_module`
fails:

> `GSE144621 is named in ['glm/tile_controls.py'] as well as its owner glm/genomewide.py`

The guard's recorded trap is that the fragment files are SPACE-delimited while the peak-tile file is
TAB-delimited, so a second loader that does not inherit `sniff_delimiter` can parse one of them to
zero rows — which reads as *"no signal"* rather than *"no data"*. `tile_controls._rows` hardcodes
`split("\t")`: correct for this file today, wrong by exactly the named mechanism. This is the
second-loader failure the resource-ownership rail exists to prevent, written while the plan beside
it discussed `genomewide.py`'s read-only status.

**The fix dissolves the violation rather than exempting it.** Measured: `load_peak_tiles` returns all
**46,713** rows including **936 `random` + 470 `neg_control` at 150 bp** — identical counts to what
the second loader produced. So `load_control_tiles` should be a filter over the owner's output and
the parser deleted. An `ALLOWED_MENTIONS` entry would be the WRONG fix: the module genuinely reads
the file, and the exemption list is for prose-only mentions.

Baseline suite: **14 failed / 7227 passed / 180 skipped** (797.9 s). `CLAUDE.md` records the 14 as
*"all absent DATA, NONE a code defect"* — now wrong by one.

### C2 — Step 1 contradicts itself and would fail its own acceptance test `[grounded]`

Step 1 says *"reproduce `extract_promoter_regions` exactly"*. The upstream function threads a running
`last_end` through features in **file order**. Measured: GFF3 CDS rows are **not** ascending by start
within either strand (`+` False, `-` False). The feasibility probe SORTS before threading, so the
pinned 3,496 windows / 625,255 bp / 28,622 in-window come from a sorted replica; a faithful unsorted
replica threads `last_end` differently and will not reproduce them.

Compounding it: the authors parsed **GenBank**, not GFF3 — a different source with a different
order — so "faithful" is not well-defined from this input at all. The probe now records the sorted
choice as an explicit DIVERGENCE rather than a reproduction, and Step 1 must do the same: make
ordering a named parameter, document the divergence, and re-derive the three pinned figures under
whichever order is chosen.

### C3 — both additive cells lack the control the subtractive experiment got right `[inferred]`, key premise `[grounded]`

- **C2 cell (regression adjustment).** Adding an informative feature to the GC baseline strengthens
  it, which *mechanically* shrinks the atlas margin whether or not the atlas advantage is
  positional. Shrinkage is not absorption unless it exceeds what a **matched-informativeness,
  non-positional** feature produces. The subtractive experiment had a random-direction control for
  exactly this; the additive one has no analogue.
- **C3 cell (stratification).** The strata are **not exchangeable** — measured: in-window IQR
  **0.4735** vs out-of-window **0.2710**, medians 0.4594 vs 0.3962, and the sd *inverts*
  (3.380 vs 5.535, outlier-driven), so the shapes differ and not merely the scales. A rank
  correlation on a compressed outcome attenuates for BOTH channels, so a smaller within-stratum
  margin is expected for reasons unrelated to frame absorption. The plan's pre-commitment to
  re-measure the baseline within the stratum handles the LEVEL but not this attenuation of the
  MARGIN.

**Fix for both:** a size-matched RANDOM stratification as the null, and each stratum's expression
spread reported beside its margin. The bar must NAME these before it is frozen, or the experiment
gets graded without them.

## 2. Medium issues

- **M1** `in_window_any` has no stated threshold `[grounded]`. At `>0` it fires on 28,622; a small
  epsilon lowers it. Name the threshold.
- **M2** The binary collapses a graded variable `[grounded]`. Histogram: 14,915 at 0.0, 9,664 at
  1.0, ~18.5k spread between. Consider grading on `overlap_fraction`, or a three-way
  none/partial/full stratification. **The specific worry that `in_window_any` was dominated by
  trivial overlaps is REFUTED** — only 8.2% overlap by ≤0.2.
- **M3** Duplicate sequences bite harder under stratification `[inferred]`. 689 rows share a
  sequence with another row; in a ~14.4k stratum a train/test-straddling duplicate matters more than
  in the 43k cohort.
- **M4** The sentrux anti-drift check — see CORRECTIONS; already replaced.

## 3. Deferred issues

- The 689 duplicate sequences as a substrate-wide question: changing the split would make every
  committed number on this rung incomparable, so it is out of scope here.
- `pos_control` geometry incompatibility (0 of 106 rows at 150 bp) — a recorded limit, unchanged.
- The `doc_validity` grammar mismatch (12 of 724 wiki documents carry an authored `doc-status`
  header; 0 parse) — separate, and the parser lives in the shared cross-project skill tree.

## 4. What's solid

- **The cross-device anchor.** All 10 per-seed values bit-identical (max |Δ| 0.00e+00) and per-tile
  predictions at `r = 1.000000` against the original GPU run. Unusually strong, survived review
  untouched, and it is the right thing to gate every new number on.
- **Locating by sequence rather than stated coordinates.** 2.77% agreement and offsets in families
  make this unambiguous.
- **Grader in Wave 0, with the reachable outcome set required to EQUAL the check order** rather than
  merely be contained in it — the shape that catches a dead validity branch.
- **Bar frozen after feasibility, before results; 0.02 inherited rather than invented.**
- **The pre-commitment that the measured 0.5076 is never retro-edited.**
- **Excluded-tile arithmetic reconciliation** (emitted + excluded = input).

## 5. Open tradeoffs

Carried to DECISIONS NEEDED at the top, where a reader will actually see them.

## AMENDED 2026-11-01 — two of this review's own claims were wrong, and its open question is settled

**"Roughly a ten-line change" was wrong.** The module change is small, but the owner enforces a
minimum of eight columns and requires `peak_start`, `peak_end` and `strand`, which the five-column
synthetic fixtures lacked. Seven of the nine tests used them, so the fixture rebuild was the bulk of
the work. The estimate is left above as written; this is the correction.

**"An ALLOWED_MENTIONS entry would be the WRONG fix" was half wrong.** True only of adding a row
*instead of* deleting the parser. False as a blanket claim: the guard's criterion is the bare text
pattern `\bGSE\d+` over whole files, and the module names the accession in a docstring sentence that
survives the fix. Measured during execution — with the parser deleted the guard **still failed** on
that one mention. Clearing it required both acts in order, and the row was added only after
observing the failure, which is what makes it evidence rather than assumption.

**The ordering divergence this review flagged is real.** Measured at file order 3,490 windows /
624,812 bases against sorted 3,496 / 625,255 — a symmetric difference of 26 windows and 146 of
43,061 uniquely-located tiles changing in-window status. So the "exact reproduction" claim genuinely
could not stand, and row ordering is now an explicit parameter with every pinned figure labelled.

**Two figures in this review are superseded.** Stated-coordinate agreement is 1193/43061 = 2.77%
across every uniquely-located tile, and the trivial-overlap share is 2,359 or 8.2%. Both already
appear corrected in the CORRECTIONS section above; noted here so the chain is complete.

**Outcome.** The red ownership guard is green: 1 failed / 7 passed → 8 passed. Full record:
`wiki/frame_plan_amendments_2026-11-01.md`.

## 6. Recommended next step

Fix **C1** first: it is red now, it is roughly a ten-line change, and it removes a defect instead of
exempting one. Then amend Step 1 for the ordering contradiction (C2) and re-derive the three pinned
figures. Then add the two missing controls to Steps 5–6 (C3) **before** the bar is frozen.

Depth: deep · Rounds: 0 (second opinion unavailable) · Stop reason: single-model fallback.
