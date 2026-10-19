# Promoter-Activity Substrate Ladder and the Mechanistic Baseline — Technical Plan

> ## ⚠ CORRECTION 2026-10-13 — READ BEFORE ACTING ON ANY STEP
>
> **Two load-bearing numbers in this plan are wrong, and the body below is left UNEDITED as the
> historical record of what was planned on a false premise.** Do not treat the step bodies as current.
>
> **1. The sigma70 rung is 10,898 usable, NOT 3,821.** Every occurrence of 3,821 / 35.1% / "the 65% gap
> is absent sequences" / "8 x 8 x 8 x 8" / "UP_element collapses to 1" in the body is VOID. The real
> shape, measured through `dna_decode.glm.expression.load_pairs`: **10,898 joined at rate 1.00**, zero
> dropped, **3 x 8 x 8 x 8 x 8 = 12,288 cells with 10,898 filled**, UP_element has **three** levels.
> 3,821 is exactly what a join without hyphen/underscore normalisation returns — the failure
> `expression.py` was written to prevent, enforces a `MIN_JOIN_RATE = 0.95` floor against, and
> test-pins as `rep.naive_join == 3821`.
>
> **2. The fragment-rung ceiling is 0.7928 (LB) / 0.8014 (M9), NOT 0.4583 / 0.4731.** Those are
> `single_replicate_agreement`; the ceiling on a sequence model is the Spearman-Brown step-up,
> square-rooted. **This reverses the plan's reasoning:** the frag rung costs ~15% of headroom for 7.3x
> the data, not ~50%, so "a 0.4583-ceiling substrate" is not the cost the plan priced in.
>
> **Step 5 is VOID and its artifact has been DELETED.** It shipped `dna_decode/glm/constructed.py` +
> `tests/test_glm_constructed.py`, a second sigma70 loader built on the naive join, duplicating a
> canonical loader that already existed with the fix. Both files are gone (neither was ever committed).
> The one original contribution — `irrelevant_feature_control` — now lives in
> `dna_decode/glm/expression.py`, generalised to all five element columns instead of Background only.
> **The canonical loader is `dna_decode.glm.expression.load_pairs`. There is no other.**
>
> Corrected substrate figures: `wiki/mpra_substrate_census_2026-10-12.{md,json}` (see its
> `corrections_2026_10_13` block). Open Question "is the 65% recoverable" is CLOSED — nothing was ever
> missing.

> Finish and verify two in-flight measurements, settle whether a 346-parameter biophysical model genuinely beats the learned encoder or merely re-encodes GC, split the ceiling criterion onto the substrate it was computed from, and register the two substrates found on disk with their leakage controls fixed in advance.

## Lens status
Inputs: conversation (2026-10-12 directives: scope to 150 bp tiles, run the 1/10/100/150 ladder, pursue same-organism-different-mutations, ladder to longer sequences, research industry practice)
Degradations: repo-index unconfigured (direct file reading only); sentrux absent (no architecture gate); `project-rules.md` absent (Project Rules Check omitted); no `graphify-out/` (reuse survey from filesystem + direct grep); **two measurements still RUNNING at plan time — every number from them is marked provisional and Steps 1-2 exist to finalize them**

## Problem Statement

The GLM ceiling criterion has been stuck at `blocked:untried-arm` because a declared 4-decade-span arm could not run on a 150 bp substrate. Four things discovered on 2026-10-12 change the problem rather than resolve it as stated:

1. **Span is the wrong lever, now twice-disconfirmed.** Four spans from 1.10 to 2.18 decades sit within noise of each other (provisional, 4 of 10 seeds), and the field moves the *opposite* way: LegNet, the DREAM-2022 winner at Spearman 0.979 on 150 bp promoters, deliberately caps its receptive field at 79 bp on a 150 bp input with a first-layer kernel of 7 and no dilated convolutions.
2. **A width-1 branch is provably information-free.** Measured: shuffling a sequence while preserving base counts changes its pooled output by 3e-8, and its max-pool half is exactly constant. The 2.18-decade ladder reaches its span only by adding a branch that re-encodes GC.
3. **A 346-parameter biophysical model appears to beat every learned arm.** The Salis-lab Promoter Calculator, never trained on this substrate, scores a provisional median 0.3856 against encoder-multi 0.3234 and GC 0.3112, winning 4/4 paired seeds against both. **This is measured on a non-random 15-34% subset in file order and is NOT yet a result.**
4. **Both substrates the user proposed were already cached.** A 321,123-fragment library at 48-499 bp, and a designed sigma70 mutational library. I had told the user the longer-sequence work needed data we did not have.

The size lever the literature ranks first is **unavailable here**: the peak-tile file holds 46,713 rows and 44,106 are already in use.

## Codebase Context

- `scripts/glm_encoder_gate.py` — the arm runner. `leave_peak_out(tiles, seed=s) -> (train, test, n_held_peaks)` is the evaluation protocol (per seed ~608 held peaks, ~11,200 test tiles). Arms are additive: without `--pyramid` the three published arms are byte-identical. Gained `--extra-widths W,W,W` (repeatable) on 2026-10-12, which names each arm by its MEASURED span so two width sets cannot collide in one artifact.
- `dna_decode/glm/encoder.py` — `DEFAULT_WIDTHS=(6,20,75)` (1.10 dec), `PYRAMID_WIDTHS=(3,12,48,150)` (1.70 dec), `validate_widths` refuses a kernel wider than the sequence. `encoder_feature_width(w, channels=32) = 2*32*len(w)`.
- `dna_decode/glm/genomewide.py` — `load_peak_tiles(cache_dir, allow_download)` + `noise_ceiling(tiles)`. **Cache root is `D:/dna_decode_cache/mpra`, NOT `D:/dna_decode_cache`** — the loader raises `GenomeWideDataError` on the parent dir, which cost a failed call during research.
- `dna_decode/glm/` holds 8 modules (`baselines`, `corpus`, `edit`, `encoder`, `expression`, `falsifier`, `generate`, `genomewide`); a new substrate loader belongs here, and 6 `tests/test_glm_*.py` files set the convention.
- **Substrate census, measured 2026-10-12** (all three in `D:/dna_decode_cache/mpra/`, space-delimited not tab):
  - peak tiles — 46,713 rows, 44,106 after the published filters, all 150 bp, **replicate ceiling 0.941**
  - `frag-rLP5` LB + M9 — **321,123 fragments, 48-499 bp (median 244, p90 290), 98.8% longer than 150 bp, replicate ceiling 0.4583**
  - `GSE108535` sigma70 — 10,898 designed promoters; **3,821 join to a sequence (35.1%), all exactly 150 bp**; normalization recovers **zero** more, so the 65% gap is absent sequences, not a name-format fix. Usable subset is crossed 8 x 8 x 8 x 8 (Minus35 x Spacer x Minus10 x Background) but **UP_element collapses to 1**; all 8 Minus10 variants still appear on >1 background, so background-cancellation survives. Zero label degeneracy (10,898 distinct values in 10,898 rows).
- **Promoter Calculator** — installed ISOLATED at `D:/tools/promcalc_env` because it pulls sklearn 1.9.1 against the project's pinned 1.8.0; verified the project venv still reports sklearn 1.8.0 / numpy 2.4.4 after install. API: `Promoter_Calculator()` class, `.run(seq, TSS_range=[0,len])` then `.output()` -> `Forward_Predictions_per_TSS` / `Reverse_Predictions_per_TSS`, each a dict of dicts carrying `Tx_rate`, `dG_total` and 8 per-element energies. ~0.26 s/tile single-core, 73 forward candidates per 150 bp tile.
- **Verified against the primary source, correcting the research agent:** the LaFleur/Hossain/Salis preprint states "346 interaction energy parameters" verbatim, R2 0.45-0.65 on external in-vivo sets — and contains **NO Urtecho/Kosuri comparison and NO CNN baseline**. The research's flagged-unverified claim of "R2=0.60 on the Urtecho data beating the CNNs" is not in the paper; its 0.60 is a genome-integrated modular-promoter set. So scoring it on endogenous tiles is a NEW measurement, not a reproduction.
- **Prior-work replay** (`wiki/decisions-log.md`): two HIGH-salience lessons bind here. *"Before proposing a direction as novel or a negative as closed, re-derive it from the artifacts — three separate strategic claims this session were confidently wrong with the grounding on disk"* (I repeated this today with the frag library). And *"when a signal targets a SUBSET, the validation denominator MUST be that subset"* — which is exactly the live risk in the provisional 0.3856.

### Reusable-Code Survey
- `scripts/glm_encoder_gate.py::leave_peak_out` — the published evaluation protocol; any new baseline must be scored on these exact per-seed splits or its number is not comparable to 0.3234.
- `dna_decode/glm/genomewide.py::load_peak_tiles` + `noise_ceiling` — reuse for every substrate rung rather than re-parsing; `noise_ceiling` is also the right tool to compute the frag library's own ceiling.
- `dna_decode/glm/encoder.py::fit_encoder` + `validate_widths` — the sigma70 rung is 150 bp, so the existing encoder runs unchanged; no new architecture.
- `~/.claude/skills/soraya/scripts/arms_ledger.py` — `measure_arm` / `declare_arm`; criteria are keyed by their text and there is **no rename or move API**, which constrains how Step 4 can split.
- `~/.claude/skills/soraya/scripts/ce_ledger.py::append_action` — the enum-validated AC9 row writer already used for this family (rows 4-6).
- Searched for shared helpers: `graphify-out/GRAPH_REPORT.md` (absent), `src/`, `lib/`, `utils/`, `common/`, `helpers/` (all absent), `dna_decode/glm/`, `scripts/`, `tests/`.

## Pre-Change Baseline

Captured 2026-10-12. **Two figures are PROVISIONAL because their runs are still executing — this is why Steps 1-2 exist.**

- Published arms on 44,106 tiles, leave-peak-out, 10 seeds: `gc` **0.3112** · `encoder-w20` 0.3163 · `encoder-multi` **0.3234** (best) · `encoder-pyramid` 0.3167 · shuffled nulls -0.0240 / -0.0198. Every encoder arm's paired CI vs GC **includes zero**.
- Replicate ceiling **0.941**; best measured is 34% of it.
- PROVISIONAL (4 of 10 seeds): `encoder-span2p18dec` 0.3206, `encoder-span1p88dec` 0.3220 — both within noise of the above.
- PROVISIONAL (14,885 of 44,106 tiles scored, **non-random, file order**): Promoter Calculator `fwd_sum_tx` median **0.3856**, `fwd_max_tx` 0.3804, `best_strand_tx` 0.3709; paired vs GC **+0.0676 (4/4 seeds)**, vs encoder-multi **+0.0390 (4/4 seeds)**. `fwd_min_dG` is the exact negative of `fwd_max_tx`, so dG and Tx_rate are one feature.
- Arms ledger: criterion `ceiling: best achievable prediction of bacterial promoter activity from sequence, reported against the replicate-derived 0.941 ceiling` -> `blocked:untried-arm`; 3 measured, 2 planned, best `conv-3width-13k-params-150bp`.
- Suite **6,409 passed / 14 skipped / 0 failed**. Frozen five byte-unchanged; prospective lock (2026-08-31) verifies `ok=True`, zero drift. Paths cross-checked as equal to the lock manifest's own list, so the frozen check is not vacuous.
- `scripts/glm_encoder_gate.py` is MODIFIED (uncommitted) with the `--extra-widths` flag; `dna_decode/glm/encoder.py` is hash-identical to HEAD (`fe71bd4`).

## Verification Signal

- Span ladder completes 10/10 seeds and the artifact reports, for each of the 4 spans, a paired mean vs GC with a CI. **Expected: every CI includes zero.** A span arm that separates from GC would falsify the twice-disconfirmed reading and must be reported loudly rather than smoothed.
- Promoter Calculator reaches 44,106/44,106 with `err=0`, and the full-set median is compared against the 14,885-row provisional. **If the median moves by more than 0.02, the non-random subset was carrying the result** and the provisional figure is retracted in the artifact, not footnoted.
- The hybrid test returns a paired `GC+dG` vs `GC-alone` delta over 10 seeds with a sign count. A delta whose CI includes zero means the mechanistic feature **re-encodes composition**; a positive delta means it supplies an orthogonal axis. Both outcomes are publishable and the artifact states which fired.
- `arms_ledger.all_ceiling_verdicts()` reports the peak-tile criterion scoped to its own substrate, and a second criterion carrying the frag library with **its own measured 0.4583 ceiling, never 0.941**.
- `dna_decode/glm/constructed.py` loads exactly **3,821** sigma70 promoters, every sequence 150 bp, and `held_out_background_split` produces folds where no background appears in both train and test — asserted, not assumed.
- The irrelevant-feature control scores **at or below chance** on the same splits. A control that scores well means the split leaks and the rung is not yet usable.
- Suite still 6,409+ with zero regressions; frozen five byte-unchanged by hash; prospective lock still `ok=True`.

## Implementation Steps

### Step 1: Finish the span ladder and record all four spans as measured arms
Files: wiki/glm_span_ladder_2026-10-12.md, wiki/glm_span_ladder_2026-10-12.json, soraya_runs/arms-ledger.json
Depends on: none

**What changes:**
- Let the running 10-seed job complete; it already emits its own JSON via `--out`.
- Write the memo reporting all four spans (1.10 / 1.70 / 1.88 / 2.18 decades) with paired means vs GC and CIs, plus each arm's shuffled null.
- Record the two new arms via `arms_ledger.measure_arm` on the EXISTING criterion with their evidence paths. This is additive and needs no authority: measuring an arm is not killing one.
- State in the memo that the 2.18-decade rung reaches its span only via an order-blind width-1 branch, with the 3e-8 shuffle measurement as evidence, so a future reader cannot read 2.18 as four real scales.
- Do NOT retire the declared 4-decade arm here — that is Step 4's gated question.

**Test strategy:**
- Assert the artifact carries 10 seeds per arm and that each shuffled null is near zero; a null that is not near zero invalidates its paired arm.
- Re-read the ledger and assert exactly two arms were added, that no existing arm's metric changed, and that `best` is recomputed by the code rather than asserted.
- Confirm the three published arms in the run reproduce their committed values, which is the reconcile check that the `--extra-widths` flag stayed additive.

### Step 2: Finish the Promoter Calculator sweep and test it against its three named failure modes
Files: wiki/promcalc_vs_learned_2026-10-12.md, wiki/promcalc_vs_learned_2026-10-12.json
Depends on: none

**What changes:**
- Complete the sweep to 44,106 with `err=0`; the runner is resumable by design, so a kill loses nothing.
- Then run the three checks that could each kill the provisional result, and report each explicitly:
  - **Subset bias.** Recompute the median on the full set and diff it against the 14,885-row figure. The scored subset was file order, not random, so this is the first-order risk.
  - **Orientation.** Forward-only vs best-of-both-strands. Forward already leads (0.3856 vs 0.3709), which is biologically expected for an orientation-defined assay; confirm it holds at full n.
  - **Reduction choice.** `max_tx` vs `sum_tx` vs `min_dG` across all seeds. `min_dG` is the exact negative of `max_tx`, so report them as ONE feature and do not present three correlated numbers as three results.
- Score on the SAME `leave_peak_out` per-seed splits as the encoder. An untrained model has every tile effectively held out, but the pairing must still be like-for-like or the comparison flatters whichever arm got the larger evaluation set.
- Record the primary-source correction: the paper makes no Urtecho and no CNN comparison, so this is a new measurement.

**Test strategy:**
- Assert `ok + err == 44,106` and that err rows are listed with their exception text, never silently dropped.
- Assert the joined prediction count equals the tile count; a short join must name the shortfall rather than report on partial coverage.
- Re-derive one tile's `Tx_rate` through the isolated env independently of the sweep and assert it matches the stored value, so a join-order bug cannot pass.
- Confirm the project venv still reports sklearn 1.8.0 / numpy 2.4.4 after all isolated-env work.

### Step 3: Settle the question that actually matters — does mechanism ADD to composition?
Files: wiki/promcalc_vs_learned_2026-10-12.md, wiki/promcalc_vs_learned_2026-10-12.json
Depends on: Step 2

**What changes:**
- Fit two models on the identical per-seed `leave_peak_out` splits: GC alone, and GC + the single biophysical feature. Report the paired delta with a sign count over 10 seeds.
- Keep it to two features. The literature on this exact substrate says "keep models simple, especially if training data is limited", a 2-feature model is what the Urtecho authors found competitive, and a larger model would confound the question with capacity.
- Also report the correlation between GC and the biophysical feature directly. A high correlation with no additive gain is the clean "re-encodes composition" verdict; low correlation with a gain is the orthogonal-axis verdict.
- Write whichever verdict fires. A null here is a real result: it would mean a 346-parameter mechanistic model and a one-number composition baseline are measuring the same thing on this substrate.

**Test strategy:**
- A shuffled-label control on the same splits must collapse to zero; a 2-feature fit is small enough to be checked this way cheaply.
- Assert the GC feature used is the gate's own, parsed from its log rather than reimplemented, so the baseline cannot drift from the one that produced 0.3112.
- Assert the fit never sees test rows, by construction from the split object rather than by inspection.

### Step 4: Split the ceiling criterion onto the substrate it was computed from — AUTHORITY-GATED
Files: soraya_runs/arms-ledger.json, project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 1, Step 2, Step 3

**What changes:** the user has chosen the split ("best prediction on 150-letter bacterial promoter tiles makes most sense to me"), with one amendment carried from discussion: scope by **substrate and ceiling**, not by length, because the frag library is also bacterial promoter data at a different ceiling and "150 letters" would quietly stand in for "the 0.941-ceiling library".

- Scope the current criterion to the peak-tile substrate and its 0.941 ceiling, carrying the now-measured arms.
- Declare a second criterion for the longer-context rung, stamped with the frag library's **own measured 0.4583 ceiling**. Quoting 0.941 against frag data would be the shared-key error one level up.
- `arms_ledger` has **no rename or move API** and criteria are keyed by text, so the mechanics must be chosen deliberately: either declare the new criteria and leave the old text as a superseded shell, or edit the sidecar directly. Surface the choice; do not silently hand-edit a Soraya-owned store.
- Append one AC9 row naming the option taken and the resulting `all_ceiling_verdicts()`.

**Why this is still gated despite the user's choice:** the split redefines what "finished" means for a criterion, and the 4-decade arm's disposition rides on it. R6 says an arm dies only by its own executed evidence. Steps 1-3 may make that retirement evidence-backed rather than a bookkeeping act, which is why this step sits last.

**Test strategy:**
- Before: snapshot `all_ceiling_verdicts()`. After: assert the peak-tile criterion's `best` is whatever the code computes from the measured set, and that no measured arm's metric or evidence changed.
- Assert the frag criterion carries 0.4583 and that the string 0.941 appears nowhere in it.
- Back up `arms-ledger.json` and diff after, so the change is reviewable and revertible.

### Step 5: Register the sigma70 constructed-variation substrate with its leakage controls fixed in advance
Files: dna_decode/glm/constructed.py, tests/test_glm_constructed.py
Depends on: none

**What changes:**
- New loader returning the **3,821** joined sigma70 promoters with sequence, expression, and the four design axes. Refuse rather than pad if a row has no sequence: 65% do not join, and that is absent data, not a formatting problem (normalization recovers zero).
- `held_out_background_split` — folds where no `Background` value appears in both train and test. This is the leakage-proof split the literature demands, and it is available here precisely because all 8 Minus10 variants appear on more than one background.
- An **irrelevant-feature control** built in from the start, because a random forest on deliberately irrelevant features scored 0.75 AUC on a saturation-mutagenesis benchmark and "reproduced the prediction quality of the best challenge submissions". The control must score at chance or the split leaks.
- Record in the module docstring that `UP_element` collapses to 1 level in the usable subset, so the design is 4-way crossed, not 5.
- Loader and splits only. **No model is trained in this step** — fitting before the controls are proven is how the cited benchmark was fooled.

**Test strategy:**
- Assert exactly 3,821 rows load and every sequence is 150 bp.
- Assert the background split is a true partition and that no background crosses the boundary — by construction from the axis values, not by sampling.
- Assert the irrelevant-feature control scores at chance on the real splits, and prove the control is non-vacuous by planting a leaked feature and confirming it scores well.
- Assert the loader REFUSES a variant row with no joined sequence rather than emitting a placeholder.

### Step 6: Write the substrate census so the frag library cannot be lost again
Files: wiki/mpra_substrate_census_2026-10-12.md, wiki/mpra_substrate_census_2026-10-12.json
Depends on: none

**What changes:**
- One artifact recording all three rungs with measured, not estimated, properties: row counts, length distributions, replicate ceilings, label-degeneracy checks, and the exact cache paths.
- State plainly that the project held a 321,123-fragment 48-499 bp library while a declared arm sat blocked on "no longer-context substrate", and that the blocker was a claim about what had been looked for. The literal 4-decade arm still needs ~10,000 bp and frag tops out at 499, so the blocker was narrower than stated rather than simply wrong.
- Record the size-lever wall: 46,713 rows is the whole peak-tile file and 44,106 are in use, so the literature's first-ranked lever cannot be pulled on the 0.941-ceiling substrate.
- Record that these files are space-delimited, since a tab-assuming parse returns zero rows and reads as an empty dataset.

**Test strategy:**
- A test asserting every cache path the census names resolves on disk, so the census cannot promise a file that is not there.
- Re-derive one figure per rung from the raw file and assert it matches the census, so the artifact is checkable rather than narrative.
- Assert the census's ceiling figures are each attached to their own rung, with no rung inheriting another's.

## Execution Preview

Wave 0 (4 parallel):  Step 1 — Finish the span ladder and record all four spans, Step 2 — Finish the Promoter Calculator sweep and test its three failure modes, Step 5 — Register the sigma70 substrate with leakage controls, Step 6 — Write the substrate census
Wave 1 (1 parallel):  Step 3 — Settle whether mechanism adds to composition
Wave 2 (1 parallel):  Step 4 — Split the ceiling criterion (AUTHORITY-GATED)

Critical path: Step 2 → Step 3 → Step 4 (3 waves)
Max parallelism: 4 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — **The headline result is provisional on a NON-RANDOM subset and may not survive Step 2.** The 0.3856 median rests on 14,885 of 44,106 tiles taken in file order, which for genomic tiles is plausibly a contiguous region. Per-seed test sets were ~1,600 rows against the protocol's ~11,200. Mitigation: Step 2 diffs the full-set median against the provisional one and the Verification Signal pre-commits to retracting rather than footnoting a move over 0.02. Do not cite 0.3856 outside this plan until Step 2 lands.
- Severity: high — **Step 4 is authority-gated and must not auto-execute.** It redefines a criterion's terminal condition and decides a declared arm's fate. `/execute-plan` reaching Wave 2 without a recorded decision must stop. Mitigation: Steps 1-3 are additive and may convert the retirement into an evidence-backed act.
- Severity: medium — **`arms_ledger` has no criterion rename or move API**, and criteria are keyed by their text. A split therefore cannot be done cleanly through the public surface; the options are a superseded shell or a direct sidecar edit, and the latter bypasses the enum-validated writer on a Soraya-owned store. Named, not resolved — Step 4 surfaces the choice rather than picking silently.
- Severity: medium — **The sigma70 rung is a third of what it first appeared.** 3,821 of 10,898 variants join to a sequence; I quoted 10,898 in discussion before checking. Normalization recovers zero, so the gap is absent data. `UP_element` also collapses to 1 level, making the design 4-way not 5-way crossed. The plan uses 3,821 throughout.
- Severity: medium — **The literature's first-ranked lever is unavailable and that is measured, not assumed.** Training-set size improves log-linearly and is unsaturated even at the replicate ceiling, but 46,713 rows is the entire peak-tile file. More data means leaving the 0.941 substrate, which is what makes the frag rung's 0.4583 ceiling a real cost rather than a footnote.
- Severity: medium — **External tool surfaces, all confirmed by direct execution or primary source:** `Promoter_Calculator()` + `.run(seq, TSS_range)` + `.output()` with `Forward_Predictions_per_TSS`/`Reverse_Predictions_per_TSS` dicts carrying `Tx_rate`/`dG_total` — verified by running it on a real tile. `--extra-widths` — verified via `--help` and a smoke run that named arms by span. `leave_peak_out` signature — read from source. `arms_ledger.measure_arm`/`declare_arm`/`kill_arm` and `ce_ledger.append_action` — verified via `inspect.signature`. The 346-parameter figure and the absence of a CNN/Urtecho comparison — read from the preprint, correcting a research claim rather than inheriting it.
- Severity: medium — **A tab-assuming parse of these files returns ZERO rows and reads as an empty dataset.** It cost one failed analysis today and only surfaced because the script crashed. Step 6 records the delimiter; every loader in Steps 2 and 5 must split on whitespace.
- Severity: low — **The isolated env is load-bearing.** The Promoter Calculator pulls sklearn 1.9.1 against the project's pinned 1.8.0; installing in-project would have moved sklearn under 6,409 tests. Verified post-install that the project still reports 1.8.0. Any future step touching that package must keep `D:/tools/promcalc_env` separate.
- Severity: low — **`python3` does not exist on this host**; every command uses `python` / `uv run`. The cache root for `load_peak_tiles` is `D:/dna_decode_cache/mpra`, not its parent.
- Severity: low — **Restructuring applied:** Step 3 was split out of Step 2 so the additive-versus-GC question gets its own verdict rather than being folded into a comparison table; Step 4 was moved to depend on Step 3 as well as 1 and 2, because what counts as the criterion's `best` arm may change; Steps 5 and 6 were kept independent of the measurement chain so a blocked Wave 2 cannot strand them.
- Severity: low — **Concurrent sessions hold uncommitted edits** to pgx, certification-capstone and feature-record files in this tree. No step may stage them; each step stages only its own named files.

## Open Questions

- **Step 4's mechanics:** given no rename API, should the old criterion text remain as a superseded shell (leaving three criteria, one permanently closed) or should the sidecar be edited directly to rewrite the key? The first is honest but noisy; the second is clean but bypasses the validated writer.
- **Does the frag rung deserve its own criterion at all,** or is a 0.4583-ceiling substrate better treated as a training-size resource for the 0.941 criterion via transfer? The plan assumes a separate criterion per the user's ladder framing, but training on frag and evaluating on tiles is a different and possibly stronger design that nothing here tests.
- **Is the 65% of sigma70 variants without sequences recoverable from another GSE108535 supplementary file?** The plan proceeds on 3,821 rather than blocking, but 10,898 would be a materially better rung and only the barcode-mapping file was checked.
- **Does the biophysical feature belong in the deployed decoder** if Step 3 shows it adds to GC? That is a surface-change question touching what the project ships, not an executor call, and it is deliberately out of scope here.

## Verification

1. `uv run pytest tests/ -q` — expect 6,409 plus any tests added by Step 5, zero failures, and compare against the recorded baseline rather than a remembered one.
2. Re-run `scripts/glm_encoder_gate.py --seeds 10 --no-download` WITHOUT `--extra-widths` and confirm the three published arms and `reconcile_ok` are unchanged, proving the flag stayed additive.
3. Assert the Promoter Calculator sweep reports `ok + err == 44,106` with every error listed, and that the full-set median is stated beside the provisional one with an explicit kept-or-retracted verdict.
4. Confirm the hybrid verdict in the artifact names which outcome fired, and that its shuffled-label control collapsed to zero.
5. `arms_ledger.all_ceiling_verdicts()` — each criterion carries its own ceiling, no rung inherits another's, and 0.941 appears nowhere in the frag criterion.
6. `dna_decode/glm/constructed.py` loads exactly 3,821 rows at 150 bp; the background split is a proven partition; the irrelevant-feature control scores at chance and is shown non-vacuous by a planted leak.
7. Frozen five byte-unchanged **by hash** (`git hash-object`, not byte count — `git checkout --` rewrites line endings while `git status` still reads clean); prospective lock re-verified `ok=True` with zero drift.
8. Every cache path named in the census resolves on disk, and one figure per rung is re-derived from the raw file.

## Save-time amendments

Captured at: 2026-10-12
Source: `/save-plan` arguments

AUDIT NOTES ONLY. `/execute-plan` reads only `## Implementation Steps` for executable work; the
bullets below are provenance for human readers, not instructions. If an amendment changes a step's
file list, dependencies or structure, re-run `/technical-plan` before `/execute-plan` rather than
relying on this section to carry it.

- scope the criterion by substrate+ceiling rather than by sequence length
- sigma70 rung is 3,821 usable not 10,898
- biophysical-baseline result held provisional pending full-set recompute

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,severity-high mode=conversation -->
