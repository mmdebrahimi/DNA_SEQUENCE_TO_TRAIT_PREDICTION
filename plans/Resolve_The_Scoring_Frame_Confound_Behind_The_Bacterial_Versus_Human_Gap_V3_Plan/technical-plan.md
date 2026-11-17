# Resolve The Scoring-Frame Confound Behind The Bacterial-Versus-Human Gap (v3)

## Lens status

Execution mode recommendation: **sequential** (see Risk Flags R3 — the substrate lives under gitignored `data/`, so a worktree has no inputs).
gh available: yes · sentrux available: yes (plan is non-structural — gate will skip) · codex: **absent on this host**
Repo-index: absent — no degradation, direct file reading used throughout
Wiki replay: 3 HIGH-salience entries consulted (tile-substrate-wall 2026-10-21, promcalc regeneration 2026-10-22, zero-shot-vs-supervised 2026-10-23) plus the 2026-11-16 vacuous-statistic retraction written this session
Supersedes: v1 (saved then halted before Step 1, no file touched) and v2 (emitted, never saved). v1 is **untracked** — see R12 before overwriting it.
Review: a single-agent adversarial pass over v2's five open questions resolved four of them; this body carries those resolutions. External review tooling was unreachable (no `codex` binary, local model endpoint down), so the review is one voice, not two.

## Problem Statement

The eukaryotic arm reported a conv-over-GC gain of **+0.2471** (HepG2) and **+0.3016** (K562) against **+0.0122** on bacterial promoter tiles, and `README.md:389` states that learning from sequence buys *"roughly 20–25x more on human regulatory DNA than on bacterial promoters"*. The arm's own first honest limit states it **cannot** say which of four differences — assay, element length, kingdom, library construction — carries the effect.

A scratch measurement on 2026-11-16 (`wiki/glm_bacterial_scoring_frame_2026-11-16.md`) proposed a fifth candidate: that the bacterial scoring frame is detection-dominated, so the two numbers measure different tasks. **Part of that measurement has since been retracted, and this plan is a rewrite around the retraction rather than a continuation of it.** What fell and what stands:

- **Retracted (committed, pushed, then corrected the same day).** The claim that a single binary split of the label recovers **0.8194** of the 0.941 ceiling — "87% of the achievable ranking signal" — and the inference that the pooled Spearman is *overwhelmingly a detection metric*. Spearman is Pearson on ranks and a threshold indicator is 1 on exactly the top *k* ranks, so `rho = sqrt(3k(n-k)/(n^2-1))`, which at `n=44,106 / k=14,913` equals 0.8194 **exactly**; uniform, Gaussian, Pareto and deliberately-bimodal samples all return 0.8194 to four decimals. The statistic measured the split proportion and said nothing about the label.
- **Weakened, and repairing it is what this plan is for.** The reproducibility contrast — pooled replicate agreement **0.7945**, within actives **0.8849**, within inactives **0.3927** — is computed on subsets *defined by a threshold on the label*, so range restriction attenuates the within-group figures by construction. The committed memo already carries that caveat. It is suggestive, not decisive.
- **Unaffected, and these are what the plan builds on.** Two direct measurements, neither a derived statistic. (1) `active` is a **pure threshold on `expn_med` at 0.5778 with zero overlap either side** — a measured fact about the file — and 29,193 of 44,106 tiles (66%) sit below it. (2) The **sign flip is unanimous**: GC scores +0.3112 pooled and **−0.0464** within actives, while a 4-mer ridge scores +0.1438 pooled and +0.1320 within actives, so the margin moves from −0.1793 / 0 of 10 seeds pooled to **+0.1820 / 10 of 10** within actives. A margin between two arms under the *same* restriction is not a design artifact: range restriction attenuates a between-group predictor and a within-group predictor differently, and that differential is the claim.

**What the retraction changes about the plan.** v1's Step 1 shipped a `label_bimodality(values)` helper whose headline output was "best single-threshold-versus-label Spearman and the share of a supplied ceiling it recovers" — the retracted quantity, generalised into a reusable function and pointed at the human label in Step 4 as the plan's decisive measurement. Building it would have produced a clean side-by-side table of two numbers that each encoded only a split proportion. It is replaced by a **quantile-binned reproducibility curve** with a selection-attenuation null, which measures what the retracted statistic was reaching for and can come out flat.

**The replacement is selection-independent by construction, and that is newly possible because the deposit carries a nested replicate design nobody here had read.** Measured 2026-11-16: `GSE144621` publishes **four** sub-replicate columns — `RNA_exp_sum_{1,2}_{1,2}` — all numeric on 44,106 of 44,106 tiles with zero NaN, and each group summary `RNA_exp_sum_g` is **exactly the mean of its own two sub-columns** (44,106 of 44,106, both groups). Within-group agreement is 0.7283 (group 1) and 0.7230 (group 2) against the cross-group 0.7945 the published ceiling is derived from. So the curve can bin on **group 1** and measure agreement **inside group 2** — the binning variable and the scored pair then share no measurement, which is what makes the curve a statement about the data rather than about the binning. Binning on one group and scoring the mirror group is a free replication.

**What the review changed, beyond v2.** Four of v2's five open questions are now closed in the body rather than carried, and two of them closed against v2's stated reasoning rather than for it:

- The **published verdict needs no reinterpreting**, so the authority call v2 reserved for it does not exist. `verdict_detail` reads *"the CI for the human gain lies entirely ABOVE the bacterial gain 0.0122"* — a comparison of two measured numbers, true whatever either measures. The frame finding lands inside two limits the artifact already states rather than contradicting any claim it makes. What remains is a narrower, genuinely open call, and it is the plan's single Open Question.
- **Keeping the strongest arm's active-frame score a diagnostic is still right, but not for the reason v2 gave.** The ledger already mixes frames (see Codebase Context), so "preserve commensurability" is not available as an argument. The reason is intrinsic: a range-restricted score is not a candidate for *best achievable prediction* of the population the criterion names.
- **Re-framing the substrate cannot be settled before the measurement, but the decision rule can be fixed now**, so it moves into Step 3's frozen bar as a pre-committed recommendation instead of staying an open question.
- **The retracted number stays in place with its correction attached**, with one new rule in Step 5 closing a provenance hazard the review located.

**Non-goals.** Unchanged. This plan does not re-run any human model, does not retune any bacterial arm, does not edit `scripts/glm_encoder_gate.py::PREREGISTERED`, does not edit any arm's placement in the ledger, and does not touch the eukaryotic arm's artifacts or `verdict` field.

## Codebase Context

**`dna_decode/glm/genomewide.py` (559 LOC)** — owns the substrate. `Tile` carries `seq`, `expression`, `peak`, `strand`, `category`, `active`, `rep1`, `rep2`, where `expression` is `expn_med` and `rep1`/`rep2` are `RNA_exp_sum_1`/`_2`. `load_peak_tiles(cache_dir, *, allow_download)` returns `(list[Tile], LoadStats)` and `ensure_cached` raises `GenomeWideDataError` with the `external_root` report when the configured volume is detached — it refuses to auto-substitute. `leave_peak_out(tiles, *, frac=0.25, seed)` returns `(train, test, n_held_peaks)` and holds out **whole peaks**. `noise_ceiling(frags)` is the Spearman-Brown step-up and returns `single_replicate_agreement` + `reliability_of_2rep_mean_spearman_brown` + `ceiling_for_a_sequence_model`; it takes a list of objects carrying `.rep1`/`.rep2`, so it **cannot** accept an arbitrary column pair — that fact shapes Step 1's interface. `spearman(a, b)` is the shared tie-handling wrapper and returns 0.0 on NaN. Feature builders `gc_features`, `kmer_features(seqs, k)` and `positional_onehot` are all present.

**The loader keeps the unique key, so the sub-replicate columns can be joined without editing a shared module.** `Tile.peak` is built as `f"{peak_start}_{peak_end}_{strand}"`, so `(tile.seq, tile.peak)` reconstructs the tuple `(variant, peak_start, peak_end, strand)` — **measured unique: 44,106 distinct keys, zero duplicates**, on the `category=="tile"` + `len==150` set. The sequence alone is **not** a key: **688 sequences recur, up to 3×**, so a join on `variant` would silently mis-pair 688 tiles' replicates. The loader also `.upper()`s the sequence, which any joining reader must match.

**Substrate is LOCAL and verified.** `data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz` reproduces every published anchor: `n=44,106` after `category=="tile"` + `len==150`, `14,913` active / `29,193` inactive, single-replicate `0.7945`, ceiling `0.9410`, GC `0.3110` against the published `0.3112`, and **GC seed 0 exactly `0.3000`** — the committed reconcile anchor, whose source is `wiki/glm_genomewide_oracle_2026-10-07.json`. The default cache (`D:/dna_decode_cache/mpra`) is detached; pass `data/mpra` explicitly. The label is **not** a function of the two group summaries — `expn_med` matches their mean on **0 of 44,106** rows while rank-correlating with it at 0.9413 (and with each group at 0.8921 / 0.8936).

**`scripts/glm_encoder_gate.py` (413 LOC)** — the published conv gate. The active flag is already plumbed: `_active_labels(test)` at L142, `auroc(scores, labels)` with tie handling at L120, and `main()` prints `{n_active} active / {n_inactive} inactive`; the flag feeds the secondary AUROC framing and **never restricts the Spearman**. `PREREGISTERED` (L65) is frozen and carries `reconcile_anchor = {"arm": "gc", "seed": 0, "expected": 0.3000, "tolerance": 0.01, "source": "wiki/glm_genomewide_oracle_2026-10-07.json"}`, `generalises_margin 0.05`, `null_max_abs 0.05`, `registered_seeds list(range(10))`, and six verdicts including `INDETERMINATE_RECONCILE_FAILED` and `INDETERMINATE_PROTOCOL_DRIFT`. `ridge_gc(train, test)` returns `(spearman, predictions)`, fitting `Ridge(alpha=1.0)` on `StandardScaler`-scaled `gc_features` against `np.log1p(expression)`; its docstring records that it was re-implemented once already and that the reconcile gate is what proves fidelity. The script's default artifact path is `wiki/glm_learned_representation_<date>.json`. **Artifact names corrected from v1:** the published result lives at `wiki/glm_encoder_gate_result_2026-10-08.md` with its JSON sidecar at `wiki/glm_learned_representation_2026-10-08.json` — v1 cited a `wiki/glm_encoder_gate_2026-10-08.json` that **does not exist**.

**`dna_decode/glm/encoder.py` (439 LOC)** — `fit_encoder(train, test, *, widths=DEFAULT_WIDTHS, protocol=REGISTERED_PROTOCOL, seed, device, shuffle_labels, group_key=tile_peak_key)` returns a dict with `spearman`, `predictions`, `n_parameters`, `epochs_ran`, `device`, and the protocol actually used. `DEFAULT_WIDTHS = (6, 20, 75)` = the published 13,217-param arm. `shuffle_labels` permutes **training** targets only — the validity null. `tile_peak_key` (L209) REFUSES a record with no `.peak`. Early stopping runs against a peak-blocked inner split; `test` is never seen until the final score.

**Committed per-tile predictions exist, so the strongest arm needs no retraining.** `data/promoteratlas/asgiven_arm_pertile_2026-11-02.json` (1.83 MB, `rung: bins10_mean_asgiven`) and `data/promoteratlas/per_tile_predictions_2026-11-01.json` (2.04 MB) each hold 10 seeds of `{tile_indices, pred}`. Lengths match the committed split exactly — seed 0 → 11,198, seed 1 → 11,018 — against `wiki/asgiven_arm_result_2026-11-02.json`'s `gates/split_equality/per_seed` (`held_out_peaks: 608` every seed). Committed per-seed Spearmans exist for both the as-given rung (`rung/per_seed`, median **0.5086**) and the atlas comparator (`comparator/atlas_per_seed` = `{0: 0.5461, 1: 0.5499, 2: 0.4983, 3: 0.4986, 4: 0.5075, 5: 0.5163, 6: 0.5078, 7: 0.4933, 8: 0.5016, 9: 0.5210}`, median **0.5076**).

**The arms ledger already mixes frames, which is why Step 5's annotation is a disclosure and not a repair.** `soraya_runs/arms-ledger.json` holds **exactly one** criterion — *"ceiling: best achievable prediction of bacterial promoter activity from sequence, reported against the replicate-derived 0.941 ceiling"* — with **16 arms** (14 `measured`, 1 `planned`, 1 `killed`). One of the measured arms is `substrate-change-constructed-eukaryotic-or-human-pan` at **0.4853**: a *human* lentiMPRA score filed under a criterion whose own name says *bacterial* and whose anchor is the *bacterial* 0.941 ceiling. It does not currently set `best` (incumbent `asgiven-position-resolved-promoteratlas-150bp` at 0.5086) but sits **0.0233** below doing so. See R11 — this is pre-existing and deliberately not fixed here.

**The reader-facing surfaces Step 6 must touch, with exact locations.** `README.md:389-391` carries the headline — *"roughly 20–25x more on human regulatory DNA than on bacterial promoters"*, closing *"the bacterial substrate was holding the method back, not the method"* — **with no caveat in its own paragraph**. The caveat opens 13 lines later at `README.md:402-403` (*"CAVEAT ADDED 2026-11-16 … so the 20-25x ratio is suspect"*) and closes at `README.md:413-415` with *"Read the ratio as unresolved, not as 20-25x"*. **v1 and v2 both said this caveat lives in "the README's roadmap step 4"; it does not.** Roadmap step 4 is `README.md:362` and carries no caveat at all — only *"First measured step taken 2026-11-01 on human regulatory DNA; see below."*

**One README sentence is already false, independent of this plan's outcome.** `README.md:413-414` attributes the unmeasured human check to *"detached volume + an uninstalled dependency"*. **`openpyxl` 3.1.5 is installed**, and `scripts/fetch_lentimpra.py` accepts `--dest`, so neither stated blocker holds — the check is unrun, not blocked. Step 6 corrects this regardless of which verdict fires.

**`dna_decode/glm/lentimpra.py` (444 LOC)** — `DEFAULT_DIR = D:/dna_decode_cache/mpra_human` (detached). `load_cell(cell, data_dir, *, min_join_rate)` joins design to activity and **refuses below the join rate**; `load_activity` returns `{sequence_id: (activity, fold)}` from Table S6; `usable_cells(data_dir)` reports per-cell usability with a reason. `LentiRecord` carries `name`, `seq`, `activity`, `fold`, `cell`. Activity is `log2(RNA/DNA)` and goes negative. **WTC11 refuses** at join rate 0.0000.

**`scripts/fetch_lentimpra.py` already exists and must be reused, not reimplemented.** CLI: `--dest` (default the detached `D:` path, overridable), `--select tables|all` (`tables` is the ~74 MB four-table subset, not the 3.5 GB + 4.33 GB dumps), `--list`, `--yes`. It verifies **every byte against the publisher's Zenodo MD5**, and the MD5 — not the size — is the integrity guard.

**The retracted number's provenance surface, measured.** `0.8194` occurs as a coincidental digit string in five unrelated repo files — `data/promoteratlas/{asgiven_arm_pertile_2026-11-02,per_tile_predictions_2026-11-01,positional_arm_pertile_2026-11-01}.json`, `wiki/gemme_holdout_tables.json`, `dna_decode/pigment/hirisplex_coefficients.json` — plus the memo, the plan and `wiki/decisions-log.md`. **No `cell_registry` contract cites the memo**, so the hazard in R13 is latent rather than live.

**Dependency reality, measured not assumed:** `openpyxl` **3.1.5 is PRESENT**, `torch` 2.7.1+cu118, `sklearn` 1.8.0, `scipy` 1.17.1. `C:` has **33.8 GB free of 511 GB (93.4% full)** — ample for a 74 MB fetch.

### Reusable-Code Survey

Candidates found and adopted rather than rebuilt:
- `dna_decode/glm/genomewide.noise_ceiling` — the ceiling computation, reused verbatim on `Tile` subsets (the published 0.9410 comes through it). It accepts only `.rep1`/`.rep2`-shaped objects, so Step 1 factors the Spearman-Brown step-up into one pure function **pinned by test to agree with it to the decimal**, rather than keeping a second divergent copy of the arithmetic.
- `dna_decode/glm/genomewide.spearman` — the shared tie-handling wrapper, used for every correlation in the new module so the curve cannot disagree with the published numbers on ties alone.
- `dna_decode/glm/genomewide.leave_peak_out` — the split, reused unchanged so split discipline is inherited, not reinvented.
- `dna_decode/glm/genomewide.{gc_features, kmer_features}` and `encoder.fit_encoder` — the arms, reused unchanged so the comparison stays architecture-matched to the published +0.0122.
- `scripts/glm_encoder_gate.{ridge_gc, auroc, _active_labels}` — imported, not copied.
- `dna_decode/fba/nulls.py` — not applicable numerically (its `curveball_shuffle` is for binary matrices) but adopted as the **pattern**: the null lives in a pure, separately-tested function beside the statistic it controls, never inline in the reporting script.
- `dna_decode/glm/positional_readout.{_check_bins, bin_means}` — read for the repo's existing bin-count refusal conventions; the binning there is over sequence POSITION rather than label rank, so the helpers do not transfer and Step 1 writes its own with the same refusal discipline.
- `scripts/fetch_lentimpra.py` — the MD5-verified fetcher, invoked as a subprocess with `--dest` overridden.
- `scripts/asgiven_arm_score.py` + `wiki/asgiven_arm_result_2026-11-02.json` — the committed per-seed values that become Step 2's reconcile anchor.

Searched: `graphify-out/GRAPH_REPORT.md` (absent), `dna_decode/glm/`, `dna_decode/eval/`, `dna_decode/fba/nulls.py`, `scripts/glm_*`, `scripts/*arm*`, `tests/test_glm_*`. No existing `reproducibility_curve` / `quantile_bin` / `rank_bin` helper anywhere in `dna_decode/`, `scripts/` or `tests/`.

## Pre-Change Baseline

Published, on the pooled frame (`wiki/glm_encoder_gate_result_2026-10-08.md` + `wiki/glm_learned_representation_2026-10-08.json`): GC **+0.3112** (2 params), `encoder-w20` **+0.3163** (2,657), `encoder-multi` **+0.3234** (13,217), `encoder-shuffled` **−0.0240**. Paired margin `encoder-multi` − GC = **+0.0125 mean, 6/10 seeds, 95% CI [−0.0099, +0.0349]** (includes zero); frozen verdict `WEAK`. Best arm overall **0.5086**, pooled margin over GC **+0.1964**. Ceiling **0.941**.

Scratch, 2026-11-16, 4-mer-ridge proxy, same splitter and seeds: pooled GC **+0.3112** / 4-mer **+0.1438** / margin **−0.1793** / 0 of 10 seeds; active-only GC **−0.0464** / 4-mer **+0.1320** / margin **+0.1820** / 10 of 10 seeds. Replicate agreement pooled **0.7945** (ceiling 0.9410), active **0.8849** (0.9690), inactive **0.3927** (0.7510) — **all three within-group figures carry a selection-attenuation caveat and are the quantity Step 1's curve replaces.** The threshold-versus-label Spearman **0.8194** is **RETRACTED as vacuous** and is not a baseline for anything.

Nested replicate design, measured 2026-11-16 and not previously read: `RNA_exp_sum_{1,2}_{1,2}` all numeric on 44,106 of 44,106 tiles; `RNA_exp_sum_g == mean(sub-columns)` on 44,106 of 44,106 for both groups; within-group Spearman **0.7283** (group 1) and **0.7230** (group 2); cross-group **0.7945**.

Human side (`wiki/eukaryotic_arm_result_{HepG2,K562}_2026-11-01.json`): HepG2 conv **0.4853** / GC **0.2362** / gain **+0.2471** `[+0.2399, +0.2536]`; K562 conv **0.5237** / GC **0.2200** / gain **+0.3016** `[+0.2916, +0.3100]`; `verdict` **`SUBSTRATE_LIMITED`**, `verdict_detail` *"the CI for the human gain lies entirely ABOVE the bacterial gain 0.0122"*, five `honest_limits`. **No human replicate ceiling exists** — both artifacts carry the *bacterial* anchor (0.941) as context.

Ledger: 1 criterion, 16 arms, `best` 0.5086, runner-up-by-value 0.5076, human arm 0.4853 filed under the bacterial criterion.

Test suite: **13 failed / 8,717 passed / 204 skipped in 479.13 s**, the 13 being the standing absent-data set (`test_claude_md_growth_ceiling` 1, `test_cyp4f2_getrm_concordance` 7, `test_essentiality_report_card_ladder` 3, `test_pneumo_selection_probe` 1, `test_readme_quickstart_output` 1).

## Verification Signal

This plan has worked when all of the following hold:

1. **Step 1's curve is proven to be a statement about the data, not the design.** Each of three unrelated synthetic distributions, generated at the *same* n, bin count, global pair correlation and binning-variable correlation, returns a curve whose trend statistic falls inside the null band — and a planted heteroscedastic sample returns one that exceeds the null's maximum over the registered number of draws. A curve that cannot come out flat is rejected before it is pointed at real data.
2. **Step 1's selection-independent join is exact or it fails.** The sub-replicate reader's key set is asserted identical to the loaded tiles' `(seq, peak)` key set, both sides asserted duplicate-free, and a deliberately shifted join raises rather than scoring a mis-paired vector.
3. **Step 2 reconciles.** Pooled per-seed Spearman recomputed from the committed per-tile predictions matches the committed per-seed values within 0.001 for all 10 seeds, for whichever arm each file holds — and the arm identity is established **by that reconcile, not by filename**. A mismatch exits non-zero and reports nothing.
4. **Step 3 reconciles on the frame that has a published anchor before changing frames.** GC seed 0 on the **pooled** frame reproduces `0.3000 ± 0.01`, and `encoder-multi` reports 13,217 parameters. Only then is the active frame scored.
5. **Step 3's active-frame validity null is clean** — the shuffled-label arm scores `|spearman| < 0.05` on the active frame. A new frame does not inherit the old frame's null.
6. **Step 3's framing recommendation is fixed before the run and never read backwards.** The five verdicts each map to exactly one framing recommendation inside the frozen bar, and a test asserts the mapping is total and that no recommendation can be produced without a verdict.
7. **Step 4 either measures or refuses by name.** Either the human label's shape is reported with its join rate, or the script exits non-zero naming the specific wall (fetch failure, MD5 mismatch, join rate below threshold). If the deposit carries no replicate columns, the human **reproducibility curve is reported UNAVAILABLE by name** — never approximated from a single measurement, and never substituted by the bacterial ceiling.
8. **Step 6's two verdict-independent README corrections land regardless of outcome** — the false dependency attribution at `README.md:413-414` is corrected, and the headline at `README.md:389-391` carries its caveat in its own paragraph rather than 13 lines later.
9. **No new test failures.** `uv run pytest tests/ -q` shows the same 13 absent-data failures and nothing else, plus the new tests passing.
10. **Nothing frozen moved.** The frozen five + `uv.lock` + `wiki/finetune_arm_acceptance_bar.json` byte-identical; `prospective_lock.verify_lock` returns `ok=True` on the 2026-08-31 manifest; `scripts/glm_encoder_gate.py::PREREGISTERED`, `dna_decode/glm/genomewide.py`, `wiki/glm_learned_representation_2026-10-08.json` and both `wiki/eukaryotic_arm_result_*.json` byte-identical.
11. **Every new guard is proven non-vacuous** — each new test file is shown to fail when its subject is reverted, and the mutation is confirmed to have applied.

## Implementation Steps

### Step 1: Selection-independent frame module, with the vacuous statistic removed
Files: dna_decode/glm/active_frame.py, tests/test_glm_active_frame.py
Depends on: none

**What changes:**
- New module holding every pure function the three measurement steps share, so the frame logic has one definition. `dna_decode/glm/genomewide.py` is **not** edited.
- `spearman_brown_ceiling(r) -> dict` — the step-up arithmetic as one pure function taking a bare correlation, because `noise_ceiling` accepts only `.rep1`/`.rep2`-shaped objects while the curve needs arbitrary column pairs. One definition, pinned against `noise_ceiling` rather than forked from it.
- `read_sub_replicates(path, *, category, length) -> dict[tuple[str, str], tuple[float, float, float, float]]` — a narrow reader for the four `RNA_exp_sum_{1,2}_{1,2}` columns only, keyed on the verified-unique `(upper(variant), f"{peak_start}_{peak_end}_{strand}")` tuple and applying the caller's filter predicate verbatim. Asserts the key is duplicate-free on its own side and reports its row accounting. It does not return a label and cannot be used to re-derive one.
- `join_sub_replicates(tiles, subs) -> list[SubReplicatePair]` — joins by `(tile.seq, tile.peak)` and **raises on any key mismatch in either direction**, because the sequence alone is not a key (688 recurrences, up to 3×) and a silent partial join would mis-pair replicates while still producing a plausible curve. Reports the matched count so a 100% join is asserted rather than assumed.
- `reproducibility_curve(bin_by, pair_a, pair_b, *, n_bins=10, binning_is_independent) -> dict` — equal-sized rank bins on `bin_by`; per bin it returns `n`, `agreement` (via the shared `spearman`), the implied ceiling, and the bin's value range. `binning_is_independent` is a **required** argument the caller must state, and the function additionally refuses when `bin_by` is the same object as either scored array or is element-wise identical to one. Returns a `trend` statistic (Spearman of bin index against bin agreement) alongside the per-bin table.
- `attenuation_null(bin_by, pair_a, pair_b, *, n_bins, n_draws, seed) -> dict` — the control the retraction demands. Draws synthetic samples with the **design held fixed** (same n, same bin count, same global `pair_a`–`pair_b` correlation, same `bin_by`-to-pair correlation, Gaussian copula with the observed marginals) and returns the null distribution of the `trend` statistic as mean / p95 / **max over `n_draws`**, in the repo's established style. A curve whose trend falls inside the null band measures range restriction rather than the data, and the returned dict states that in a `reads_as` field.
- `shape_summary(values) -> dict` — what survives of v1's `label_bimodality`: modal-bin share and distinct-value count, which are the repo's own shipped `assay_degeneracy` quantities and are not determined by the design. **The best-single-threshold-versus-label Spearman and its share-of-ceiling are deleted**, with a module-level comment recording the closed form so nobody reintroduces them.
- `threshold_audit(tiles) -> dict` — carried over from v1 and still sound: it **measures** whether `active` coincides with a clean cut on `expression`, returning the implied cut, the count on each side, and `overlap_either_side`. It never assumes a threshold, and on non-zero overlap it says so and callers must treat the flag as a label.
- `frame_decomposition(tiles) -> dict` — same shape as v1 with one addition: the active and inactive blocks carry **`selection_attenuated: true`** and a pointer to the curve, because both subsets are defined by a threshold on the label and their within-group agreement is attenuated by construction.
- `restrict_predictions(tile_indices, preds, tiles, keep) -> tuple[list, list]` — carried over from v1; raises on any index out of range or any length mismatch, so a mis-aligned prediction file fails loudly instead of scoring a shifted vector.
- Every returned dict carries `conditions_on_the_outcome: bool`, true for any block defined by a threshold on the label.

**Test strategy:**
- `reproducibility_curve` + `attenuation_null`: three unrelated distributions (uniform, Gaussian, heavy-tailed) at identical design parameters each land inside the null band; a planted heteroscedastic sample (agreement rising with the binning variable) exceeds the null's max. This is the design-held-fixed test the retraction requires, and it runs before any real data.
- `reproducibility_curve` refuses when the binning array is the scored array, and when `binning_is_independent` is not supplied.
- `spearman_brown_ceiling` agrees with `noise_ceiling` to the decimal on the same synthetic replicate pair.
- `join_sub_replicates` raises on a dropped key, an extra key, and a deliberately shifted key; a synthetic case with a duplicated sequence but distinct peaks joins correctly, pinning that the key is the pair and not the sequence.
- `shape_summary` has no threshold-Spearman field at all, asserted by key-set equality so a reintroduction fails.
- `threshold_audit` recovers a planted clean cut with zero overlap; with one mislabelled tile it reports non-zero overlap and does not claim a cut.
- `frame_decomposition` reproduces a hand-computed Spearman-Brown ceiling and stamps `selection_attenuated` on both conditioned blocks.
- Real-data test gated on `data/mpra/...`: reproduces `44,106 / 14,913 / 29,193`, pooled ceiling `0.9410`, the `0.5778` cut with zero overlap, the four sub-columns numeric on all 44,106 rows, `RNA_exp_sum_g == mean(subs)` on all rows for both groups, and within-group agreement `0.7283 / 0.7230`. Skips with a reason when the file is absent.

### Step 2: Re-score the committed arms on the active frame, with no retraining
Files: scripts/active_frame_rescore.py, tests/test_active_frame_rescore.py
Depends on: Step 1

**What changes:**
- Reads `data/promoteratlas/asgiven_arm_pertile_2026-11-02.json` and `data/promoteratlas/per_tile_predictions_2026-11-01.json`, loads tiles from `data/mpra` via the explicit-cache path, and re-scores each stored prediction vector on the active subset. **No model is trained and the `D:` volume is not needed.**
- **Reconcile gate runs first and is blocking.** For each file and seed, recompute the pooled Spearman from the stored predictions and compare against the committed per-seed values in `wiki/asgiven_arm_result_2026-11-02.json` (`rung/per_seed` for the as-given rung, `comparator/atlas_per_seed` for the atlas arm). Arm identity is **assigned by whichever committed series matches**, and if neither matches within tolerance the script exits non-zero with `INDETERMINATE_ARM_IDENTITY_UNRESOLVED`.
- Also cross-checks each seed's `n_test` against `gates/split_equality/per_seed`, so a prediction file from a different split cannot pass on Spearman agreement alone.
- Reports, per arm: pooled per-seed and median (the reconciled values), active-frame per-seed and median, the paired per-seed margin against `ridge_gc` re-scored on the same active test rows, win count, and a paired bootstrap CI on the margin.
- Also reports the margin against `ridge_gc` within each of Step 1's rank bins on the held-out rows, so the arm-versus-GC contrast is read across the label's range rather than at one chosen cut. This is the graded form of the binary active/inactive split and costs nothing extra: it re-uses the same stored predictions.
- **Stamps `not_a_ledger_candidate: true` with its reason on every active-frame figure it emits.** The reason is intrinsic, not bookkeeping: an active-frame score is range-restricted on the outcome, so it is not a candidate for *best achievable prediction* of the population the ledger's criterion names. This closes v2's Open Question 2 in the artifact rather than in prose.
- Writes `wiki/active_frame_rescore_<date>.json` under its own schema. Does not write into any existing arm artifact and does not touch the ledger.

**Test strategy:**
- Unit: a synthetic prediction file whose pooled recompute matches a synthetic committed series resolves to that arm; one matching neither exits with the indeterminate code; one whose `n_test` disagrees fails even when Spearman agrees.
- Unit: the reconcile tolerance is asserted to reject a 0.01 perturbation, so it is not so loose as to be vacuous.
- Unit: the per-bin margin table is asserted to have equal-sized bins and to sum its row counts back to the held-out total, so a binning bug cannot silently drop rows.
- Unit: every active-frame figure in the emitted artifact carries `not_a_ledger_candidate`, asserted structurally so a later field addition cannot ship without it.
- Real-data test gated on the prediction files: asserts the as-given series reconciles to median `0.5086` and the atlas series to `0.5076`.
- Non-vacuity: perturb one stored prediction vector and confirm the gate fires.

### Step 3: Run the published conv on the active frame, under its own frozen bar
Files: scripts/glm_active_frame_gate.py, tests/test_glm_active_frame_gate.py
Depends on: Step 1

**What changes:**
- A **separate script with its own frozen pre-registration**, importing `ridge_gc`, `auroc`, `_active_labels` from `scripts/glm_encoder_gate.py` and `fit_encoder` / `DEFAULT_WIDTHS` from `dna_decode.glm.encoder`. `scripts/glm_encoder_gate.py::PREREGISTERED` is **not edited** — the active frame breaks its reconcile anchor by construction (GC seed 0 is ≈ −0.04 there, not 0.3000), so sharing the dict would retroactively rewrite a pre-registration and force `INDETERMINATE_RECONCILE_FAILED`.
- **Two-phase reconcile, and the ordering is load-bearing.** Phase 1 scores GC and `encoder-multi` on the **pooled** frame and asserts the published anchor (GC seed 0 `0.3000 ± 0.01` per `wiki/glm_genomewide_oracle_2026-10-07.json`, and 13,217 parameters). Only after that passes does Phase 2 score the active frame. Validating on the frame that *has* a published anchor is the only way to prove the data and split are right before the frame changes.
- Arms on the active frame: `gc`, `encoder-multi` (widths `(6, 20, 75)`, the published 13,217-param arm) and `encoder-shuffled` as the active-frame validity null. Registered seeds `range(10)`, matching the published list.
- Phase 3 reports the conv-versus-GC margin within each of Step 1's rank bins, on the same fitted models and the same held-out rows, so the headline does not rest on the single 0.5778 cut. The frozen bar still grades on the active-frame margin; the binned curve is reported beside it as the shape of that margin and is explicitly **not** a second bar.
- Pre-committed outcomes: `FRAME_EXPLAINS_THE_GAP` (active-frame margin ≥ +0.15 and ≥ 8 of 10 seeds), `FRAME_CONTRIBUTES` (margin > +0.05, CI excludes zero), `FRAME_DOES_NOT_EXPLAIN` (margin ≤ +0.05 or CI includes zero), `INDETERMINATE_NULL_NOT_CLEAN`, `INDETERMINATE_POOLED_RECONCILE_FAILED`.
- **The bar additionally carries a total `framing_recommendation` map, one entry per verdict, fixed before the run.** This is where v2's Open Question 3 goes: whether future bacterial arms should be scored on actives only, or detection and strength declared as two criteria, cannot be settled by judgement before the measurement — but the decision *rule* can be, and leaving it open invites reading the framing off the result. The map is a **recommendation only**; adopting it edits a ledger criterion and stays a user decision. The load-bearing entry is the pessimistic one: under `FRAME_DOES_NOT_EXPLAIN` the recommendation is **keep the pooled frame and declare nothing**, because re-framing onto a subset with no measured headroom would score every future arm on a frame strictly worse than today's.
- The frozen bar is stamped `not_a_blind_pre_registration: true` with the reason — the 2026-11-16 scratch measurement motivated it, so the thresholds were set with a directional expectation in view, following the existing field of that name in `wiki/asgiven_arm_result_2026-11-02.json`. It is also stamped `frozen_after_the_mechanism_analysis: true`, which the repo's own standing rule requires: a bar frozen *before* the analysis that says what a correct fix can achieve is mis-specified, and the two bars mis-specified that way this quarter were both frozen too early.
- Pins `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` / `NUMEXPR_NUM_THREADS` to 1 and records the per-seed wall time, per the recorded BLAS-oversubscription lesson.
- Writes `wiki/glm_active_frame_gate_<date>.json`.

**Test strategy:**
- Unit on the verdict function: each of the five outcomes is reachable from a synthetic arm dict, and `INDETERMINATE_POOLED_RECONCILE_FAILED` takes precedence over any graded verdict.
- Unit: the `framing_recommendation` map is total over the five verdicts, asserted by key-set equality against the verdict enum so adding a verdict without a recommendation fails; and a recommendation cannot be produced without a verdict.
- Unit: the binned margin table cannot change the verdict — the verdict function is called with two different binned tables and the same arm medians, and returns the same verdict both times.
- Unit: the script's import of `PREREGISTERED` is read-only — a test asserts `scripts/glm_encoder_gate.py` is byte-identical to `HEAD` after the module is imported, and that the new script defines its own distinct pre-registration dict.
- Unit: a null arm scoring 0.2 forces `INDETERMINATE_NULL_NOT_CLEAN` regardless of the candidate margin.
- Smoke run on a small tile subset with `--seeds 1` to prove the end-to-end path without the full CPU cost.
- Thread-pinning is asserted by reading the environment the script sets, not by timing.

### Step 4: Fetch the human tables locally and measure the human label's shape
Files: scripts/human_frame_shape.py, tests/test_human_frame_shape.py
Depends on: Step 1

**What changes:**
- Invokes the existing `scripts/fetch_lentimpra.py --select tables --dest data/mpra_human --yes` as a subprocess. The fetcher's MD5 verification is the integrity guard and is **not** reimplemented; the `--dest` override exists precisely so the detached `D:` default can be bypassed.
- Loads HepG2 and K562 via `lentimpra.load_cell(cell, data_dir="data/mpra_human")` and reports each cell's join rate. WTC11 is expected to refuse at join rate 0.0000 and that refusal is recorded, not worked around.
- The human-versus-bacterial comparison is **two separate things, reported separately and never merged into one headline.** (a) `shape_summary` on each cell's `activity` and on the bacterial `expression` — modal-bin share and distinct-value count only. (b) `reproducibility_curve` on the human label **if and only if** the deposit carries per-replicate columns; the script enumerates the columns it found and records `human_reproducibility_curve: null` with the reason when it does not. There is no fallback, no single-measurement proxy, and no inherited bacterial curve.
- Reports explicitly whether a human ceiling is derivable. If replicate columns exist, computes it via `spearman_brown_ceiling`; if not, records `human_ceiling: null` with the reason, because the eukaryotic artifacts currently carry the *bacterial* 0.941 as context and that substitution must not be silently inherited.
- **Does not invent a human `active` flag.** lentiMPRA publishes no such call, so nothing here thresholds a label this plan chose.
- Writes `wiki/human_frame_shape_<date>.json`.

**Test strategy:**
- Unit: the subprocess invocation is asserted to carry `--select tables` and the overridden `--dest`, so a future edit cannot silently pull the 7.8 GB dumps or write to the detached volume.
- Unit: a simulated fetch failure, a simulated MD5 failure and a join rate below threshold each exit non-zero with distinct codes; none returns a partial shape measurement.
- Unit: with no replicate columns present, `human_reproducibility_curve` and `human_ceiling` are both `null` with reasons, and a test asserts the artifact never carries `0.941` in any human field.
- Unit: the artifact's shape keys are asserted to be exactly `shape_summary`'s key set, so a threshold-Spearman field cannot reappear on the human side.
- Real-data test gated on `data/mpra_human/` being present; skips with a reason otherwise.

### Step 5: Write the resolution artifact and make the frame caveat machine-visible
Files: wiki/glm_scoring_frame_resolution.md, wiki/glm_scoring_frame_resolution.json, soraya_runs/arms-ledger.json, tests/test_scoring_frame_artifact.py
Depends on: Step 2, Step 3, Step 4

**What changes:**
- One artifact consolidating Steps 2–4 into an answer to the original question, with the verdict taken from Step 3's frozen bar and the Step 2 upper bound and Step 4 human comparison beside it.
- States plainly which of the eukaryotic arm's four named differences this resolves and which it leaves open, and carries forward the five honest limits from the 2026-11-16 memo that still apply — in particular that the active frame conditions on the outcome and that its absolutes are not comparable to pooled ones.
- Carries a `supersedes` block naming the retracted statistic, its closed form, and the fact that this plan's curve is the replacement — so a reader who encounters the 0.8194 anywhere in git history lands on the correction rather than re-deriving it. Every curve trend is reported **with its null band in the same block**, never as a bare number.
- **Records that the published eukaryotic `verdict` is untouched and why, rather than leaving it as an open authority call.** `verdict_detail` is a comparison of two measured numbers and stays true whatever this plan finds; the frame result lands inside `honest_limits[0]` (which difference carries the effect — now five candidates) and strengthens `honest_limits[1]` (the 0.0122 anchor treated as a fixed point) from *uncertain* to *frame-dependent*. The artifact states both relationships explicitly so the next reader does not re-open a settled question.
- **Adds a `never_cite_for_provenance` field naming `wiki/glm_bacterial_scoring_frame_2026-11-16.json`.** The repo's contract-number audit verifies a cited number by whole-number presence in the cited artifact, so a contract citing that memo would find the retracted `0.8194` and verify it as provenance. No contract cites it today; this field plus a test keeps it that way. The memo itself is **not** edited and the retracted value is **not** deleted — deleting the evidence of the error is what lets it be re-derived.
- Adds **one** annotation key to the `ceiling:` criterion in `soraya_runs/arms-ledger.json` recording that bacterial metrics under that criterion are measured on a detection-dominated frame, pointing at the new artifact. The writing script asserts every pre-existing field unedited and **refuses to write unless the change is augment-only**. It does **not** add an arm, does not create a second criterion, and does not move the mis-filed human arm (R11).
- Does **not** edit any arm's `metric`, `status`, `evidence` or `verdict`, and does not touch the eukaryotic arm artifacts.

**Test strategy:**
- Artifact-shape tests: required keys present, every cited `wiki/` path exists, the JSON parses (the recorded trap — an invalid hand-edited JSON artifact was committed and cited once).
- A test asserting no curve trend statistic appears in the artifact without its null band in the same block, and that the `supersedes` block names the retracted quantity.
- A test asserting no `cell_registry` contract cites any path listed in `never_cite_for_provenance`, which both pins today's clean state and fires if a future contract adds the citation.
- Ledger augment-only test: arm count 16 before and after, key set identical, every `metric` / `status` value identical, exactly one new key added at the criterion level and none at arm level.
- A test asserting the artifact does not state a bare ratio without the frame caveat in the same section.
- Non-vacuity: delete the new annotation and confirm exactly the new tests fail while the pre-existing ledger tests pass.

### Step 6: Correct the reader-facing surfaces, two of them regardless of the verdict
Files: README.md, dna_decode/eval/regime.py, tests/test_readme_frame_caveat.py
Depends on: Step 5

**What changes:**
- **Verdict-independent correction 1 — a false statement.** `README.md:413-414` attributes the unmeasured human check to *"detached volume + an uninstalled dependency"*. `openpyxl` 3.1.5 is installed and `scripts/fetch_lentimpra.py` takes `--dest`, so neither blocker holds; the check is unrun, not blocked. Replace with what is actually true, naming the Step 4 path that reaches it. This lands whichever verdict fires.
- **Verdict-independent correction 2 — adjacency.** The headline at `README.md:389-391` states the 20–25× ratio and closes *"the bacterial substrate was holding the method back, not the method"* with **no caveat in its own paragraph**; the caveat opens at `README.md:402` and the retraction note at `README.md:413`. A reader who stops at 391 gets the uncaveated inference the frame work questions. Move a one-sentence pointer into the headline paragraph. Also corrects v1's and v2's own mislocation — the caveat is **not** in roadmap step 4 (`README.md:362`, which carries none).
- **Verdict-dependent part.** Replace the interim *"Read the ratio as unresolved, not as 20-25x"* with the measured outcome under whichever of Step 3's pre-committed verdicts fired, in the wording fixed in that bar before the run.
- Updates the `constructed → molecular` row in `dna_decode/eval/regime.py` only if Step 3 returns `FRAME_EXPLAINS_THE_GAP` or `FRAME_CONTRIBUTES`, since that is the row whose cited evidence the bacterial-versus-human contrast informs. If `FRAME_DOES_NOT_EXPLAIN`, the row is left unchanged and the artifact records that the frame hypothesis was tested and rejected.
- Does **not** act on Step 3's `framing_recommendation`; that is surfaced for the user, not applied.

**Test strategy:**
- A new guard asserting the README never again attributes the human check to an uninstalled dependency, and that the 20–25× headline paragraph contains a caveat marker. Proven non-vacuous by reverting the sentence and confirming it fails.
- A test asserting the README no longer contains the interim "unresolved" wording once a verdict has been recorded, and that it cites the new artifact.
- Existing README guards must stay green: `tests/test_advertised_commands.py` + `tests/test_advertised_env_vars.py` (148 tests) resolve every advertised command through its real parser.
- `tests/test_regime_boundary.py` must still make all three historical scope-compressions fail loudly; any row edit is checked not to weaken those.

## Execution Preview

| Wave | Steps | Parallelism |
|---|---|---|
| 0 | Step 1 | 1 |
| 1 | Step 2, Step 3, Step 4 | 3 |
| 2 | Step 5 | 1 |
| 3 | Step 6 | 1 |

Total waves: **4**. Max parallelism: **3**. Critical path (toolkit): **Step 1 → Step 2 → Step 5 → Step 6**. The toolkit's path is topological; the **wall-clock** long pole is Step 3 in the same wave — 10 seeds × 3 arms of CPU conv training, though on 14,913 active tiles rather than 44,106, so roughly a third of the published run's cost.

Wave 1 has no intra-wave file overlap: Steps 2, 3 and 4 each own one script and one test file, and all three only read Step 1's module.

Step 1 is the single wave-0 bottleneck and is deliberately **not** split: the null exists only to control the curve, and the reader exists only to feed that curve selection-independent inputs. Separating them would let the curve land in one step and its control in another, which is the shape that produced the retraction being repaired here.

## Risk Flags

**R0 — v3 exists because v1's decisive statistic was vacuous, and the replacement must not repeat the class.** The defence is structural rather than careful: `reproducibility_curve` ships with `attenuation_null` in the same module, the null holds the design fixed, it is exercised on three unrelated distributions before any real data is read, and a planted heteroscedastic case must exceed the null's max. If the curve cannot come out flat on unrelated data, it is measuring the design and Step 1 fails before Steps 2–4 can cite it.

**R1 — The active frame conditions on the outcome, and no step may hide it.** `active` is a threshold on `expression` itself, so every active-frame number is range-restricted. Step 1 stamps `conditions_on_the_outcome` on every such block and `selection_attenuated` on the within-group reproducibility figures; Step 2 stamps `not_a_ledger_candidate` on every active-frame figure it emits. The finding rests on the *differential* between GC and the learned arms under the *same* restriction, not on any absolute — which is why Steps 2 and 3 survive the retraction intact while the reproducibility contrast did not. Mitigated, not eliminated.

**R2 — `[unverified]` whether the lentiMPRA deposit carries replicate-level data.** If it does not, no human reproducibility curve and no human ceiling can be computed, and the human-versus-bacterial comparison is limited to distributional shape. Step 4 must record both as `null` with reasons rather than inherit the bacterial 0.941, which both eukaryotic artifacts currently carry as context. This is the one place the plan's decisive comparison can come back unavailable, and it must come back unavailable rather than approximated.

**R3 — Sequential execution recommended over parallel worktrees.** The substrate (`data/mpra/`) and the committed predictions (`data/promoteratlas/`) sit under gitignored `data/`, so a worktree would start with no inputs and an unreconcilable baseline — the same reasoning that made the previous two plans in this arc sequential.

**R4 — A frozen pre-registration must not be retrofitted.** Step 3 creates a new bar rather than extending `glm_encoder_gate.PREREGISTERED`, and a test pins that file byte-identical to `HEAD`. Extending it would both break its own reconcile anchor and rewrite a commitment after seeing results.

**R5 — Step 3's bar is not a blind pre-registration and says so twice.** The 2026-11-16 scratch result motivated the thresholds (`not_a_blind_pre_registration: true`), and the bar is frozen *after* the mechanism analysis (`frozen_after_the_mechanism_analysis: true`) because the two mis-specified bars this quarter were both frozen before it.

**R6 — BLAS oversubscription will mis-time Step 3 if threads are unpinned.** The recorded instance cost a 2.34× throughput error and a wrong projection. Step 3 pins all four thread variables to 1 and records per-seed wall time.

**R7 — Do not report the 4-mer proxy and the conv as the same thing.** The scratch margin (+0.1820) came from a 4-mer ridge, not the published conv, and the repo has a recorded instance of a zero-shot result being reported as replicating a supervised one. Steps 2, 3 and 5 must label each number with the arm that produced it.

**R8 — Disk.** `C:` is 93.4% full with 33.8 GB free. Step 4 adds ~74 MB and `--select tables` must not drift to `--select all` (7.8 GB); a test pins the flag.

**R9 — The sub-replicate join is the new silent-failure surface.** The sequence is not a unique key (688 recurrences, up to 3×), so a join written on `variant` alone would mis-pair 688 tiles' replicates and still produce a plausible curve. Step 1 joins on `(seq, peak)`, asserts both key sets duplicate-free and identical, and a deliberately shifted join is tested to raise.

**R10 — v1 and v2 both cited surfaces that are wrong.** `wiki/glm_encoder_gate_2026-10-08.json` does not exist (real paths: `wiki/glm_encoder_gate_result_2026-10-08.md` and `wiki/glm_learned_representation_2026-10-08.json`; the 0.3000 anchor's source is `wiki/glm_genomewide_oracle_2026-10-07.json`), and the README caveat Step 6 edits is **not** in roadmap step 4. Both corrected here. Every path and line reference in this body was verified against the working tree this session.

**R11 — PRE-EXISTING, not introduced here: the ledger already mixes frames, and one arm is mis-filed.** `soraya_runs/arms-ledger.json` has one criterion naming *bacterial* promoter activity against the *bacterial* 0.941 ceiling, and `substrate-change-constructed-eukaryotic-or-human-pan` (a human lentiMPRA score, 0.4853) sits inside it. It does not set `best` today (incumbent 0.5086) but is 0.0233 away. This means "preserve commensurability" is **not** available as the reason for keeping an active-frame score out — the reason is R1's intrinsic one. Moving the human arm is a ledger edit beyond this plan's augment-only discipline and is left to the user; Step 5's augment-only test pins arm count at 16 so this plan cannot silently change it either way.

**R12 — v1 is UNTRACKED, so saving v3 in place destroys it irrecoverably.** `git ls-files --error-unmatch` fails on `plans/Resolve_The_Scoring_Frame_Confound_Behind_The_Bacterial_Versus_Human_Gap_Plan/technical-plan.md`; the whole folder shows as `??`. The repo's save-once rule prefers an in-place rewrite plus an index `Status:` update over a new file, and that is the right call — but it assumes git holds the predecessor, and here it does not. **Commit v1 before overwriting it.** v1 is the only artifact showing what the vacuous statistic propagated into, and this session's auto-memory refers to it. Supporting facts: no `.claude/execute-plan-state/` file exists for this plan (12 files there, none matching), so in-place rewrite breaks no resume; the index entry keeps its single `[plan_file: ...]` key and gains `**Status:** candidate` (47 of 86 entries carry one; v1's does not). Not `superseded` — that value requires a `**Superseded by:** <basename>` and an in-place rewrite produces none.

**R13 — The retracted number is a latent provenance hazard, which is why it is disclosed rather than deleted.** The contract-number audit verifies a cited number by whole-number presence in the cited artifact, so a contract citing `wiki/glm_bacterial_scoring_frame_2026-11-16.json` would find `0.8194` and verify it. No contract cites it today — only the plan, `README.md:415` and `wiki/decisions-log.md`. Step 5's `never_cite_for_provenance` field plus its test hold that line. Note separately that `0.8194` already occurs as a coincidental digit string in five unrelated repo files, which is exactly why whole-number matching alone cannot be trusted as semantic provenance.

## Open Questions

1. **Should the two eukaryotic artifacts gain a sixth `honest_limits` entry naming the scoring-frame confound?** *For:* a reader who reaches `wiki/eukaryotic_arm_result_{HepG2,K562}_2026-11-01.json` first will not see the caveat, and the repo's own standing rule is that a correction must reach every surface restating the claim. *Against:* those files are dated records of what was known on 2026-11-01, the claim they state is not falsified by anything here, and editing a published artifact to graft in a later finding is the pattern that makes dated artifacts untrustworthy. The plan's current position is the **forward pointer only** (Step 5), leaving both artifacts byte-identical. This is a user authority call and the only one left.

## Verification

- `uv run pytest tests/ -q` → same 13 absent-data failures, no new ones, plus the new test files passing.
- Each new test file proven non-vacuous by reverting its subject and confirming only that file fails — with the mutation confirmed to have applied, per the recorded `MUTATION-DID-NOT-APPLY` trap.
- Step 1's null exercised on three unrelated distributions at identical design parameters, all landing inside the band, plus one planted heteroscedastic case exceeding the null's max over the registered draws.
- Step 3's `framing_recommendation` map asserted total over the verdict enum.
- Frozen five + `uv.lock` + `wiki/finetune_arm_acceptance_bar.json` byte-identical to `HEAD`; `prospective_lock.verify_lock` on the 2026-08-31 manifest returns `ok=True` with 5 pinned files.
- `scripts/glm_encoder_gate.py`, `dna_decode/glm/genomewide.py`, `wiki/glm_learned_representation_2026-10-08.json`, both `wiki/eukaryotic_arm_result_*.json` and `wiki/glm_bacterial_scoring_frame_2026-11-16.json` byte-identical to `HEAD`.
- `soraya_runs/arms-ledger.json` diffed before and after: 16 arms before and after, identical key set, no `metric` or `status` changed, exactly one new criterion-level key and none at arm level.
- Step 2's reconcile reproduces as-given median `0.5086` and atlas median `0.5076`; Step 3's pooled phase reproduces GC seed 0 `0.3000 ± 0.01` and 13,217 parameters.
- `uv run python scripts/project_status.py` still exits 0 and reports the registry unchanged at 128 cells.

## Save-time amendments

Captured at: 2026-11-16
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads
ONLY `## Implementation Steps` for executable work, so nothing here is an executable instruction.
If an amendment changes a Step contract (file lists, dependencies, structure), re-run
`/technical-plan` before `/execute-plan` rather than relying on this section.

- v3: folds in the review resolutions
- the published eukaryotic verdict needs no reinterpreting (verdict_detail is a CI comparison that survives the frame finding)
- the ledger already files a human arm under the bacterial criterion so commensurability was never the reason to keep an active-frame score out (R11)
- the framing decision becomes a pre-committed map inside Step 3's bar
- Step 5 gains a never_cite_for_provenance rule
- Step 6 gains two verdict-independent README corrections including one currently-false sentence
- and R12 records that v1 is untracked
<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,test-strategy-leak -->
