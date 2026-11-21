# Is The Markov-Chain Defeat A Property Of The Bacterial Substrate Or Of Generators

## Lens status
Repo-index: absent — `repo-index-status.sh` returned nothing, so research was direct file reading with no Layer 2 degradation.
Tools: `sentrux` present (plan is non-structural, so the gate will skip), `gh` present, `codex` ABSENT — the pre-exec review will be single-voice, a named weakness for a plan whose core risk is a cross-substrate comparison that can be made to say what one wants.
Baseline tests: 12 failed / 9054 passed / 191 skipped, all 12 in four documented absent-data families.
Predecessors read in full before this plan was designed: `wiki/eukaryotic_arm_result_HepG2_2026-11-01.md` and `..._K562_2026-11-01.md`.

## Problem Statement

The generator's one robust measured failure is that it loses to a **3-mer Markov chain**: against its
own seed-matched null it is more distinguishable from natural sequence at **10 of 10 seeds**, mean gap
**+0.1186** (`wiki/glm_generator_paired_sweep_2026-11-20.md`). That finding survived the paired-noise
sweep intact while the companion "fine-tuning made it worse" delta did not.

Everything measured about it is on ONE substrate: *E. coli* K-12 upstream-of-CDS windows, 150 bp,
4,340 of them, and they are PROXIES for promoters rather than mapped transcription start sites.

**Meanwhile the project has already measured that this substrate is the weaker one.** The eukaryotic
arm (frozen bar 2026-11-01, run locally in 38 minutes) asked whether the bacterial ceiling was the
method's or the substrate's, and returned **`SUBSTRATE_LIMITED`** on both cell types: a conv encoder
reaches median Spearman **0.4853** on human HepG2 lentiMPRA against **0.3234** on bacterial tiles,
with the human gain's paired-bootstrap CI **[+0.2399, +0.2536]** lying entirely above the bacterial
anchor **+0.0122**. K562 agrees (per-fold conv 0.5106-0.5357 against GC 0.2164-0.2444).

So the question this plan answers is narrow and cheap: **is the Markov-chain defeat a property of the
bacterial substrate, or of generators in general?**

If a 3-mer chain is a *weaker* opponent on human regulatory sequence — relative to that substrate's
own noise floor — then the 0.1186 gap is partly an artifact of having measured on the substrate the
project has already shown to be impoverished, and the generator arm should move before it is judged.
If the chain is an *equally strong or stronger* opponent there, the defeat is about generators and
moving substrate will not help.

**THE METHODOLOGICAL CONSTRAINT, taken from the predecessor rather than invented.** The K562 memo
states: *"The comparison that crosses the substrate boundary is the **within-substrate gain**, not raw
Spearman. A within-substrate contrast is the only quantity that crosses a substrate boundary. Raw
Spearman does not: the panels have different noise floors and different assays."* Our two substrates
also differ in LENGTH (230 bp vs 150 bp), which changes a k-mer discriminator's job directly. So this
plan must NOT compare a human distinguishability number to the bacterial 0.6679. The decision
quantity is each substrate's **null-minus-self-control** margin, which is within-substrate by
construction.

Non-goals: no neural generator is run here (deliberate — see below); no training of anything; no
change to the frozen bacterial surface; no claim about whether a generator would do better on human.

**Why no neural generator in this plan.** Adding one would require a new base model (the current
generator is `GENERator-v2-prokaryote-1.2b-base`, prokaryote-locked) and would confound the substrate
question with a model-family change. The Markov null and the uniform control are enough to
characterise the OPPONENT, they need no GPU, no Kaggle and no training, and they answer the substrate
question on its own. A generator arm is the follow-on, gated on this result.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

- Bacterial, from the committed paired-sweep grid: Markov null mode-both **0.6679** (seed 7),
  self-control ceiling mode-both **0.5865** -> **null-minus-self-control margin +0.0814**. Kmer mode:
  null 0.5990, ceiling 0.5447 -> margin **+0.0543**.
- Bacterial generator, for context only and NOT part of this plan's comparison: zero-shot 0.7787,
  worse than its own seed-matched null at 10 of 10 seeds, mean gap +0.1186.
- Bacterial substrate: 4,340 upstream-of-CDS PROXY windows, 150 bp.
- Human substrate: **159,795 measured pairs** (HepG2 61,463 + K562 98,332) at **230 bp**, CC-BY-4.0,
  authors' 10 folds — **PRESENT on this host** at `data/mpra_human/` (~72 MB, fetched 2026-11-16),
  both cells loading. No fetch is required.
- Human replicate ceiling: **NOT DERIVABLE** — no per-replicate columns in the deposit
  (`wiki/human_frame_shape_2026-11-16.json`). The bacterial side's 0.941 has no human counterpart
  and must not be borrowed as one.
- The label question is **already screened and CLEARED** for this candidate
  (`wiki/lentimpra_gate_screen_2026-10-27.json`, `intended_layer: L4_forward_continuous`), with the
  screen's own contract: a PASS bounds the LABEL only, not reachability, regime fit, or worth-doing.
- The published incumbents are **NOT available as a comparator**:
  `wiki/lentimpra_metric_resolution_2026-10-28.md` is `UNRESOLVED_TABLE_IDENTITY_UNVERIFIED`,
  attempted and refused.
- Human null-minus-self-control margin: **NEVER MEASURED.** That absence is this plan's subject.
- Human encoder performance, measured 2026-11-01 and NOT re-measured here: conv 0.4853 (HepG2 median)
  vs GC 0.2362; K562 per-fold conv 0.5106-0.5357.
- Test suite: 12 failed / 9054 passed / 191 skipped; the 12 are cyp4f2/GeT-RM 7,
  essentiality-ladder 3, README-quickstart 1, pneumo cps DB 1.

## Verification Signal

- The fetch completes and every file passes the publisher's MD5; a byte mismatch ABORTS rather than
  proceeding with a partial file.
- `usable_cells()` reports HepG2 and K562 usable, and **WTC11 still refusing** — a sudden WTC11 pass
  would mean the loader's gene-keyed guard had broken, not that the data improved.
- `alphabet_and_length_ok(seqs, 230)` passes on the loaded natural set, so the sequences really are
  230 bp over ACGT before any statistic is computed.
- Both the Markov null and the uniform control are scored, and the **uniform control lands high**
  (bacterially 0.9641); a uniform control near chance would mean the discriminator is broken and the
  run is void regardless of the null's value.
- `self_control_ceiling` is re-derived on the human set per feature mode, never inherited from the
  bacterial 0.5865.
- The reported decision quantity is a within-substrate margin for each substrate, and the memo
  contains **no direct comparison of raw human and bacterial distinguishability levels**.
- Frozen five byte-unchanged; `prospective_lock.verify_lock` returns `ok=True`.
- No regression against the 12-failure baseline, reconciled by CAUSE not count.

## Implementation Steps

### Step 1: Verify the already-present substrate and record its provenance
Files: wiki/human_falsifier_substrate_2026-11-21.json
Depends on: none

**What changes:**
- NO FETCH. The substrate is already at `data/mpra_human/` (6 files, ~72 MB, 2026-11-16). This step
  verifies what is there and records its provenance, so the measurement rests on a checked artifact
  rather than on the assumption that a five-day-old download is intact.
- Re-verify each file against the committed `_fetch_manifest.json` (and `_zenodo_record.json` for the
  deposit identity), recording a PER-FILE verdict rather than a single aggregate pass. A mismatch
  ABORTS: a truncated xlsx that still opens is the documented worst case, and the sibling lesson is
  that a truncated gzip decompresses its opening blocks perfectly and raises only at the end.
- Record `usable_cells('data/mpra_human')` per cell, and pass `data_dir` explicitly everywhere
  downstream — never relying on `DEFAULT_DIR`, which points at the detached drive and is exactly what
  made an earlier draft of this plan believe the data was absent.
- The D:-rule question is CLOSED, not deferred: 72 MB already sits on C:, so there is no destination
  decision to surface.

**Test strategy:**
- Assert `usable_cells('data/mpra_human')` reports HepG2 and K562 usable and WTC11 still refusing at
  its 0.95 join bar.
- Assert the artifact records a per-file manifest verdict, not one aggregate boolean.
- Assert the loader is called with an explicit `data_dir` and that no code path in this plan reads
  `DEFAULT_DIR`.

### Step 2: Pre-register the cross-substrate reading
Files: wiki/human_falsifier_bar_2026-11-21.json, tests/test_human_falsifier_bar.py
Depends on: none

**What changes:**
- A frozen bar, written before any human number exists, declaring the decision quantity as the
  **within-substrate margin** `margin = null_distinguishability - self_control_ceiling`, computed
  separately per substrate and per feature mode, with the bacterial values quoted as the committed
  anchors (mode-both +0.0814, kmer +0.0543).
- An explicit prohibition, in the bar itself: **the raw human distinguishability level may NOT be
  compared to the bacterial 0.6679.** The substrates differ in length (230 vs 150 bp), assay, kingdom
  and library construction, and the predecessor memo states that only a within-substrate contrast
  crosses that boundary. The bar records this so a later reader cannot quietly make the easier,
  invalid comparison.
- Outcomes, ordered, first match wins: `INSTRUMENT_INVALID` (the uniform control does not land high,
  or `alphabet_and_length_ok` fails, or the self-control cannot be derived); `CHAIN_IS_WEAKER_ON_HUMAN`
  (human margin exceeds the bacterial margin by more than the self-control's own seed spread — the
  3-mer opponent is relatively weaker there, so the bacterial defeat is partly substrate-bound);
  `CHAIN_IS_COMPARABLE` (the two margins differ by less than that spread — the defeat is about
  generators, and moving substrate will not relieve it); `CHAIN_IS_STRONGER_ON_HUMAN` (human margin
  lower by more than the spread — moving substrate would make the generator's task harder).
- The tolerance is the **self-control's own seed spread on the human set**, derived in the run rather
  than asserted here, because the magnitude is exactly what the run measures and freezing a constant
  now would be guessing at the scale.
- A statement of what this design does NOT bound: one human panel, two cell types, 230 bp only, no
  neural generator, and a k-mer-plus-positional discriminator that cannot see long-range organisation
  so any low score remains necessary-not-sufficient.
- A two-sided forecast with reasoning.

**Test strategy:**
- Assert the bar parses, carries all four outcomes in a declared order, defines the decision quantity
  as a difference of two measured quantities, and names the bacterial anchors.
- Assert the bar contains the explicit prohibition on the raw cross-substrate comparison.
- Assert no outcome condition is decided by a constant other than the two committed bacterial anchors.

### Step 3: A human natural-set path for the falsifier
Files: dna_decode/glm/human_natural.py
Depends on: Step 1

**What changes:**
- A new module providing `human_natural_set(cell, data_dir, *, fold=None, limit=None) -> (list[str],
  dict)`: loads lentiMPRA via `load_cell`, returns the 230 bp sequences plus a stats dict recording
  cell, n, length, fold selection and join rate.
- **The fit/eval split uses the authors' published 10-fold assignment**, not a split invented here, so
  the Markov null is fitted and scored on halves the publisher defined. This mirrors the bacterial
  harness's fit/eval discipline without importing its genome-shaped assumptions.
- The module REFUSES rather than degrading: an absent `data_dir`, an unusable cell (WTC11), or any
  sequence not exactly 230 bp over ACGT raises, carrying the loader's own reason string.
- It does NOT modify `scripts/glm_generator_falsifier.py` or `dna_decode/glm/corpus.py` — the
  bacterial path stays byte-identical, so nothing this plan does can alter a committed bacterial
  number.

**Test strategy:**
- With the substrate absent, assert the module raises and the message names `fetch_lentimpra.py`
  (so the failure is actionable, matching the loader's existing behaviour).
- With the substrate present, assert length 230 for every sequence, the fold selection round-trips,
  and WTC11 refuses.
- The bacterial-harness-untouched assertion lives in Verification item 2 rather than here, so this
  step's Files: line stays honest: it does not modify those modules, and naming them as step files
  would create a dependency edge that does not exist.

### Step 4: Tests for the human path
Files: tests/test_human_natural.py
Depends on: Step 3

**What changes:**
- The Step 3 strategies as real tests, skipping cleanly when the substrate is absent (the repo's
  established pattern for gitignored data) rather than failing and adding to the 12-failure set.
- A test that the bacterial anchors used by the bar match the committed grid artifact, so the
  comparison cannot drift from the numbers it cites.
- A NON-VACUITY check on the refusal path: assert the module raises for a 150 bp sequence list, since
  a length guard that never fires is not a guard.

**Test strategy:**
- `uv run pytest tests/test_human_natural.py tests/test_human_falsifier_bar.py -q` green.
- Each new load-bearing assertion confirmed to FAIL when its subject is reverted.

### Step 5: Measure the human margin locally and commit the raw numbers
Files: scripts/human_falsifier_margin.py, wiki/human_falsifier_margin_2026-11-21.json
Depends on: Step 2, Step 4

**What changes:**
- A script that, per cell and per feature mode: builds the natural set, fits `MarkovGenerator(order=3)`
  on the authors' fit folds, generates `n` sequences at 230 bp, scores `heldout_auroc` against the
  eval folds, derives `self_control_ceiling` on the eval folds, and also scores `UniformGenerator` as
  the instrument check.
- CPU-only and local. No GPU, no Kaggle, no training — the predecessor eukaryotic arm ran locally in
  38 minutes with a conv encoder, and this is strictly cheaper.
- Per-cell, per-mode records are written as the run proceeds, so a truncated run leaves a readable
  partial table rather than an average over an unknown denominator.
- The raw numbers are committed verbatim BEFORE any interpretation, as the predecessor plan
  established.

**Test strategy:**
- Assert the uniform control lands high on the human set; if it does not, the run is
  `INSTRUMENT_INVALID` by Step 2's frozen reading and no margin is reported.
- Assert every cell-by-mode cell is present with no gaps.

### Step 6: Apply the frozen reading and write the verdict
Files: wiki/human_falsifier_margin_2026-11-21.md
Depends on: Step 5

**What changes:**
- Compute each substrate's margin, apply Step 2's ordered rule with no post-hoc threshold, and report
  the human and bacterial margins side by side as MARGINS — never as raw levels.
- State plainly which of the two readings the result supports: that the Markov-chain defeat is partly
  an artifact of the bacterial substrate, or that it is a property of generators and substrate will
  not relieve it.
- Record the limits the bar named, plus one the predecessor insisted on: this is ONE substrate change,
  so it cannot say WHICH difference (length, assay, kingdom, library construction) carries any effect.
- If the result is `CHAIN_IS_WEAKER_ON_HUMAN`, say explicitly that this does NOT show a generator
  would succeed there — only that the opponent is relatively weaker — and that a generator arm is a
  separate, gated question.

**Test strategy:**
- A test asserting the memo's stated verdict is the one Step 2's rule yields from the committed
  numbers, so the prose cannot drift from the computation.
- A test asserting the memo contains no direct raw-level comparison between the two substrates.

## Execution Preview

5 waves, max parallelism 2, critical path 5 steps (Step 1 -> Step 6 via 3, 4, 5).

- Wave 0: Step 1 (fetch + feasibility), Step 2 (pre-registration) — disjoint files
- Wave 1: Step 3 (human natural-set module)
- Wave 2: Step 4 (tests)
- Wave 3: Step 5 (the local measurement)
- Wave 4: Step 6 (verdict)

Sequential execution is expected for the predecessors' reasons: gitignored data makes a worktree's
absent-data failure set diverge from the baseline this plan reconciles against, and the fetched
substrate is a host-global precondition rather than a per-worktree one.

## Risk Flags

- **The invalid comparison is the easy one, and it is the plan's central risk.** Comparing human
  distinguishability directly to the bacterial 0.6679 is one subtraction away at all times, is
  intuitive, and is wrong — the substrates differ in length, assay, kingdom and construction.
  Mitigated by freezing the prohibition into the bar itself and by a test asserting the memo contains
  no raw-level comparison. Stated as a risk rather than assumed solved, because a prose rule is weaker
  than a computation.
- **Length is a direct confound on a k-mer discriminator.** 230 bp gives a 4-mer discriminator more
  positions to work with than 150 bp, which can move distinguishability for reasons unrelated to
  substrate quality. The within-substrate margin absorbs part of this (the self-control is derived at
  the same length) but not necessarily all of it. Unresolved and named.
- **The lentiMPRA sequences carry 15 nt adaptors at each end** (*"15nt 5' adaptor - 200nt element -
  15nt 3' adaptor"*). Constant flanks are trivially learnable and could inflate or deflate
  distinguishability depending on whether generated sequences reproduce them. Whether to strip to the
  200 nt element is an Open Question, not silently decided.
- **The substrate is absent and must be fetched**, and the standing rule puts heavy files on a drive
  that is not attached. Step 1 surfaces the destination as a decision rather than writing to C: on
  its own judgement.
- **This plan deliberately contains no generator**, so it cannot say anything about whether generation
  would succeed on human. That is a scope choice to avoid confounding substrate with model family;
  the risk is that a reader takes a favourable margin as a green light for generation, which Step 6
  is required to forbid explicitly.
- **`codex` is absent**, so the pre-exec review is single-voice on a plan whose core risk is a
  cross-substrate comparison that can be made to say what one wants.
- **Pre-existing, non-blocking:** 12 suite failures in four absent-data families, none touched here.
- **Pre-existing and plan-relevant:** the human metric from the 2026-11-01 arm sits inside a criterion
  whose other ten metrics are bacterial, and `ceiling_verdict` takes a max over metrics. It is
  harmless today (0.4853 < bacterial best 0.5076) and a guard
  (`test_the_human_arms_metric_is_not_read_as_best`) fires if it stops being harmless. This plan adds
  no metric to that criterion, deliberately.

## Open Questions

- Strip the 15 nt adaptors and measure on the 200 nt element, or keep the full 230 bp as the publisher
  shipped it? Keeping them matches the measured activity labels and the encoder arm's geometry;
  stripping removes a constant flank that a discriminator can trivially exploit. Measuring BOTH is
  affordable here since the run is CPU-only, and may be the better answer — but it doubles the table
  and needs the bar to name which is primary.
- Which cell is primary, HepG2 or K562? They have different n (61,463 vs 98,332) and the encoder arm
  reported both. Reporting both is cheap; declaring one primary matters only if they disagree, and the
  bar should say what happens then rather than leaving it to the writer.
- ~~Where should the fetched substrate live?~~ **CLOSED** — it is already on C: at
  `data/mpra_human/`, 72 MB, fetched 2026-11-16. No decision needed.
- Given there is NO derivable human replicate ceiling, is a margin against the self-control ceiling
  sufficient on its own? The bacterial side can express a result as a fraction of 0.941; the human
  side cannot, so the two substrates' results will not be expressible on a common normalised scale
  even as margins. This plan's comparison does not require it, but a later reader may want it and
  cannot have it from this deposit.
- Should the bacterial margin be recomputed at 230 bp to isolate the length confound, by re-running
  the bacterial null on longer windows? It would separate length from substrate, which this plan
  otherwise cannot do — but it touches the bacterial harness, which this plan deliberately leaves
  byte-identical. Named rather than silently scoped out.

## Verification

1. `uv run pytest tests/test_human_natural.py tests/test_human_falsifier_bar.py -q` green, each new
   load-bearing test confirmed to fail when its subject is reverted.
2. `git diff` shows `dna_decode/glm/corpus.py` and `scripts/glm_generator_falsifier.py` unchanged.
3. Frozen five byte-unchanged; `prospective_lock.verify_lock` returns `ok=True`.
4. Full suite reconciled by CAUSE against the 12-failure baseline, not by count.
5. Every fetched file's publisher MD5 verified and recorded per file.
6. The uniform control lands high on the human set before any margin is read.
7. The memo reports margins, never raw cross-substrate levels, and a test enforces it.

## Save-time amendments

Captured at: 2026-11-21
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads
ONLY `## Implementation Steps` for executable work; nothing here is an instruction.

- no generator by design — characterise the Markov opponent only, so substrate is not confounded with model family; decision quantity is the within-substrate margin, never a raw cross-substrate level

- **Both clauses above are load-bearing and both are grounded in the plan body.** "No generator by
  design" is why Step 5 runs only `MarkovGenerator` and `UniformGenerator`: adding a neural generator
  would require a new base model, since the current one is `GENERator-v2-prokaryote-1.2b-base`, and
  that would confound "is this substrate better" with "is this model family better". "Within-substrate
  margin, never a raw cross-substrate level" is taken verbatim in spirit from the predecessor memo,
  which states that *"a within-substrate contrast is the only quantity that crosses a substrate
  boundary"* — and the two substrates here also differ in LENGTH (230 bp vs 150 bp), which changes a
  k-mer discriminator's job directly.

- **THE BODY WAS CORRECTED BEFORE SAVING, which deviates from save-verbatim deliberately.** The
  version emitted to chat had a Step 1 that fetched the lentiMPRA substrate. That fetch is
  unnecessary: the data has been on disk at `data/mpra_human/` since 2026-11-16 (6 files, ~72 MB,
  both cells loading). Persisting a plan whose first step instructs a pointless 72 MB download is
  worse than deviating from verbatim, so Step 1 was reframed from a fetch into a verification of
  what is already present, keeping its integrity discipline against the committed
  `_fetch_manifest.json`. Step numbering and the dependency graph are unchanged, so the wave
  structure is identical.

- **Five further corrections landed in the same pass**, all from a user challenge to double-check
  the human side: the Pre-Change Baseline now says PRESENT rather than absent; WTC11's refusal
  reason is the live join-rate bar (*"joins at 0.0000 (0 of 23064 measured rows), below the 0.95
  bar"*) rather than the gene-keying an earlier draft paraphrased from a docstring; the D:-rule Open
  Question is CLOSED rather than deferred; the already-CLEARED gate screen
  (`wiki/lentimpra_gate_screen_2026-10-27.json`) and the refused published-metric identity
  (`wiki/lentimpra_metric_resolution_2026-10-28.md`, `UNRESOLVED_TABLE_IDENTITY_UNVERIFIED`) are
  recorded in the baseline; and the measured absence of any human replicate ceiling
  (`wiki/human_frame_shape_2026-11-16.json`) is recorded along with the instruction that the
  bacterial 0.941 must not be borrowed as one.

- **How the "absent" error happened, recorded because it recurred three times in one session.**
  `usable_cells()` called with no argument defaults to `DEFAULT_DIR = D:/dna_decode_cache/mpra_human`
  on a detached drive, so it reports every cell unusable. The same shape of error had already
  produced "server1 cannot run this" (it was powered off) and "Kaggle is auth-gated" (the CLI was
  missing; credentials were valid). `CLAUDE.md` prescribes `scripts/data_inventory.py` before
  claiming any data is missing; it was skipped each time, and the committed
  `wiki/data_inventory.md` is itself stale (Oct 9, no `mpra_human` row).

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,reuse-missing -->
