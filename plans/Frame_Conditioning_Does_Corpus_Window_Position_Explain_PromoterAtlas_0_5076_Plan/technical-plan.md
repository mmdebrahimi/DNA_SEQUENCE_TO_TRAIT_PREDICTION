# Frame Conditioning — Does Corpus-Window POSITION Explain PromoterAtlas's 0.5076? — Technical Plan

> The linear experiments showed that removing the detection direction does not cost the atlas channel; this asks the same question from the other side — hand the baseline the exact positional rule that built the pretraining corpus, and see whether the gain survives.

## Lens status
Inputs: conversation (candidate 7 of `project_state/glm-learned-representation-2026-10-07.md`, drafted earlier as Steps 4 and 7 of an emitted-but-unsaved plan).
Degradations: repo-index absent (`repo-index-status.sh` not on this host) — direct file reading only, no Layer 1/2 summaries; `project-rules.md` absent — Project Rules Check omitted; `DESIGN.md` absent — no UI surface; graphify absent — Reusable-Code Survey built from direct grep; `codex` absent; `sentrux` and `gh` present. `wiki/plans-index.md` carries no entry for this domain and `wiki/decisions-log.md` has no entry key mentioning glm/promoter/tile/coordinate, so the replay contributed lessons (below) rather than prior decisions.

## Problem Statement

The exposure investigation closed its linear half on 2026-11-01: removing a 32-dimensional discriminative subspace drives control-vs-tile separability from 0.9156 to 0.6823 while costing the atlas channel only 0.0080 of attributable performance, and a nonlinear probe buys just +0.0099 over a linear one, so the "nonlinear survives" caveat is tight. That is a **subtractive** argument — take the detection signal away and the number holds.

This plan runs the **additive** argument, which is the one a sceptic asks for: give a cheap baseline the exact positional rule that generated the pretraining corpus and see whether the atlas gain is absorbed. It is now sharply specifiable, because the rule is known rather than guessed — `extract_promoter_regions(gbs, min_intergenic_length=5, regulatory_region_length=200)` at `src/promoter_atlas/utils/genomics.py:6` takes the 200 bp immediately upstream of each annotated gene start, strand-aware, gated on a 5 bp minimum intergenic distance.

Measured at full scale for this plan, and the numbers enable the design rather than merely permitting it:

- A faithful replica of that rule over the cached MG1655 annotation yields **3,496 windows covering 625,255 bases = 13.47% of the genome**.
- **28,622 of 43,061 uniquely-located tiles (66.47%) overlap such a window**; 8,239 (19.13%) sit entirely inside one. So corpus membership is a feature with real variance, not a near-constant — which is what makes the absorption test possible at all.

Three cells, in increasing strength:

- **C1 — does the encoder KNOW corpus membership?** Predict the in-window flag from the 128-dim embedding, held out. This is the positional restatement of the already-measured 0.9156 separability.
- **C2 — does handing the frame to the baseline absorb the gain?** Add the frame features to the GC baseline and re-measure the atlas margin over that strengthened baseline.
- **C3 — the stratified cell, and the strongest.** Within in-window tiles only, and within out-of-window tiles only, does atlas still beat GC? This controls for frame membership **by construction** rather than by regression, so it does not depend on the frame being linearly modellable.

**Non-goal:** this is not a new arm. A positional or fine-tuned PromoterAtlas variant is a NEW ARM under R6 (that is candidate 8); this plan measures a limit on the existing measured arm and must not move `arms.native150` or the ceiling verdict.

## Codebase Context

**The coordinate trap is confirmed at scale and is the plan's central hazard.** The substrate's reference is `U00096.2`; the cached annotation is `GCF_000005845.2` / `NC_000913.3`. Measured on 2,000 tiles: the stated `start` column agrees with the genome at **63/2000 = 3.15%**, a coincidence rate. Located by direct sequence search instead, the top stated-vs-located offsets are **+1978 (387), +1977 (374), +1976 (294), +777 (286), +776 (85), 0 (83)** — offset *families*, i.e. cumulative indel drift between assembly versions, not a constant shift that could be corrected. **Any frame feature must locate by sequence and ignore `start`/`end` entirely.**

**Locating is cheap and nearly complete.** Sweeping both strands of the 4,641,652 bp genome against a dict of tile sequences takes **5 s** and resolves **43,386 of 43,417 distinct sequences (0.9993)**, of which **43,061 are UNIQUE (0.9918)**, 325 are ambiguous (>1 occurrence) and 31 are not found. The sweep direction matters: indexing 44k tile keys and sliding the genome is fast, whereas 44k independent `str.find` calls over 4.6 Mb is not.

**A pre-existing duplication fact the incumbent measurement also carries.** The 44,106 scored rows map to **43,417 distinct sequences**, so 689 rows share a sequence with another row. `leave_peak_out` splits on PEAK, so a sequence duplicated across two peaks can land on both sides of a split. This is a property of the substrate and of every arm on this rung, not of this plan, and it is reported rather than fixed here.

**The extraction rule has two subtleties that must be replicated, not corrected.** In `extract_promoter_regions`, the `if start < regulatory_region_length: continue` guard sits *before* `last_end_forward = end`, so a skipped gene does **not** update the running `last_end`; and the reverse-strand branch tests `last_end_reverse - end >= min_intergenic_length` where `last_end_reverse` holds the previous gene's **start**. Both look like defects. The question this plan asks is *what the corpus contained*, not *what it should have contained*, so the replica reproduces them and says so.

**Annotation parsing is already available and verified.** `dna_decode.data.annotations.parse_gff3` returns 9,523 rows from the cached GFF3, of which **4,340 are `CDS` on the single contig `NC_000913.3`** (strand 2,222 `-` / 2,118 `+`), all **4,340 carrying a `locus_tag`** — which matters because the extractor filters on `"locus_tag" in feature.qualifiers`. GFF3 is 1-based inclusive while Biopython's `location.start` is 0-based and `.end` exclusive, so the replica converts explicitly.

**`genomewide.py` and `expression.py` must not be edited.** `scripts/kaggle/build_promoteratlas_tile_arm_kernel.py` embeds `genomewide.py` as base64+gzip and `scripts/promoteratlas_arm_dryrun.py` asserts byte-identity with the repo copy. Both files were verified byte-unchanged at the end of the preceding run and must stay so; all new code goes in new modules.

**The harness and its cached products already exist from the preceding run.** `dna_decode/glm/atlas_embed.py` (`load_atlas`, `masked_mean_pool`, `ridge_spearman` mirroring `fit_predict_ridge` exactly), `dna_decode/glm/tile_controls.py`, `dna_decode/glm/exposure_verdict.py`, and the gitignored embedding cache `data/promoteratlas/tile_embeddings_native150_n44106.npy` — so no step here re-embeds, and the 15.7-minute pass is not repeated.

**`DEFAULT_CACHE` is on a detached volume.** `genomewide.py:53` sets `DEFAULT_CACHE = Path("D:/dna_decode_cache/mpra")` and `D:` is absent; every call passes `cache_dir` explicitly with `allow_download=False`.

**Test infrastructure.** pytest, configured in `pyproject.toml`. Real-data tests skip when gitignored inputs are absent — the pattern `tests/test_tile_controls.py` uses. `tests/test_contract_number_audit.py` uses `wiki/` as its decoy pool, so float dumps stay in gitignored `data/promoteratlas/`.

### Reusable-Code Survey
- `dna_decode/data/annotations.py` — `parse_gff3` (verified: 4,340 CDS rows, all with `locus_tag`), `_load_genome_dict`, `_revcomp`; read-only reuse for the window replica
- `dna_decode/glm/atlas_embed.py` — `ridge_spearman` (the incumbents' exact protocol), `masked_mean_pool`, `load_atlas`; built in the preceding run, reused unchanged
- `dna_decode/glm/exposure_verdict.py` — the validity-first grader shape and its `CHECK_ORDER`/reachability test pattern, to be mirrored not imported
- `dna_decode/glm/genomewide.py` — `load_peak_tiles`, `leave_peak_out`, `spearman`, `gc_features`; read-only, must not be edited
- `scripts/promoteratlas_subspace_removal.py` — the `auroc` mid-rank implementation (ties averaged, because `sorted()`-order tie-breaking has silently shifted medians in this repo) and the seeded control pattern
- `data/promoteratlas/tile_embeddings_native150_n44106.npy` — the cached embedding pass; reused so nothing re-embeds
- None additionally — searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/`, `lib/`, `utils/`, `helpers/` (none exist at repo root), `dna_decode/glm/`, `dna_decode/data/`, `scripts/`

## Pre-Change Baseline

Every figure below is re-derived from a committed artifact or measured by the probe written for this plan; none is carried from conversation prose.

- From `wiki/promoteratlas_tile_arm_2026-11-01.json`: `arms.native150.median` **0.5076**, `arms.padded200.median` **0.4983**, `gc_anchor.median` **0.3112**, `gc_anchor.seed0` **0.3000** (`tol` 0.01), `n_seeds` **10**, `n_tiles_scored` **44106**, `determinism_max_abs_delta` **0.00e+00**.
- From `wiki/promoteratlas_control_separation_2026-11-01.json`: separability embedding **0.9156**, GC-only **0.7641**; atlas minus 1-D separation direction **0.5095**; atlas minus a random direction **0.5075**; verdict `EXPOSURE_EXCLUDED_LINEAR`.
- From `wiki/promoteratlas_subspace_removal_2026-11-01.json`: at k=32, atlas **0.4956**, residual separability **0.6823**, random-subspace atlas **0.5036**, attributable drop **0.0080**.
- From `wiki/promoteratlas_nonlinear_bound_2026-11-01.json`: linear probe **0.9176**, nonlinear **0.9275**, gain **+0.0099**, verdict `LINEAR_CAVEAT_TIGHT`.
- From `data/promoteratlas/local_reconcile_2026-11-01.json`: verdict `RECONCILED`; atlas per-seed **10/10 exact**, GC per-seed **10/10 exact**, max \|Δ\| **0.00e+00**; per-tile correlation with the arm's committed predictions median **r = 1.000000**.
- Measured for this plan (to be re-emitted by code in Steps 1–2): stated-coordinate agreement **63/2000**; distinct tile sequences **43,417** of 44,106 rows; located **43,386 (0.9993)**; uniquely located **43,061 (0.9918)**; ambiguous **325**; not found **31**; replicated corpus windows **3,496** covering **625,255 bp (13.47%)**; uniquely-located tiles overlapping a window **28,622 (66.47%)**, fully inside **8,239 (19.13%)**.
- Inherited margin: **0.02**, from `wiki/promoteratlas_arm_acceptance_bar.json` via `wiki/promoteratlas_exposure_acceptance_bar.json`. Not re-derived, not invented.

## Verification Signal

The run succeeds when the Step-5 frozen bar returns one of its named outcomes with every validity check passing first, and the arm's exposure limit gains a positional clause. Named outcomes:

- **FRAME_EXPLAINS** — handing the frame to the baseline drops the atlas margin below the inherited 0.02, **or** the stratified margin collapses within either stratum. Corpus-window position was carrying the gain.
- **FRAME_PARTIAL** — the margin shrinks materially but stays above 0.02. The limit becomes a quantified positional discount.
- **FRAME_EXCLUDED_CONDITIONAL** — the margin survives both the strengthened baseline and both strata. Bounded: see the pre-commitments.
- **INDETERMINATE_HARNESS_DRIFT** — the cached embedding pass no longer reproduces `arms.native150` 0.5076.
- **INDETERMINATE_LOCATE_RATE** — unique-locate rate below the floor stated in Step 2.
- **INDETERMINATE_FRAME_UNINFORMATIVE** — the frame features carry no expression signal on their own. Then "adding them does not absorb the gain" is vacuously true and says nothing.
- **INDETERMINATE_BASELINE_NOT_STRENGTHENED** — adding the frame did not improve the baseline's own training fit, so the baseline never used what it was given.

Two gates sit in front of all of it: the cached embeddings must still reproduce 0.5076, and the frame must be shown informative before any absorption claim is read.

## Implementation Steps

### Step 1: Faithful replica of the corpus extraction window
Files: dna_decode/glm/corpus_window.py, tests/test_corpus_window.py
Depends on: none

**What changes:**
- `corpus_windows(cds_rows, genome_length, *, regulatory_region_length=200, min_intergenic_length=5)` returns the ordered list of `(start0, end0_exclusive, strand, locus_tag)` windows, reproducing `extract_promoter_regions` exactly — including that a gene skipped by the `start < regulatory_region_length` guard does NOT update the running `last_end`, and that the reverse-strand branch compares against the previous gene's START.
- Both quirks are documented in the module as DELIBERATE REPLICATION of the authors' code, with the reasoning: the question is what the corpus contained, not what it should have contained. A "corrected" window would measure a rule nobody trained on.
- `window_coverage(windows, genome_length) -> bytes` gives a per-base mask; `overlap_fraction(mask, pos0, length)` gives a tile's overlap.
- CDS rows are filtered to those carrying a non-empty `locus_tag`, because the extractor filters on `"locus_tag" in feature.qualifiers`.
- Explicit 1-based-inclusive GFF3 to 0-based-half-open conversion at the single boundary where it happens, asserted by test rather than assumed.
- Pure: standard library plus the already-parsed annotation table. No torch, no network.

**Test strategy:**
- Unit on a synthetic genome with planted CDS features: forward window placement, reverse-strand window placement, the `min_intergenic_length` rejection, and the `start < 200` edge rejection.
- A test that pins the `last_end`-not-updated-on-skip quirk by constructing a case where correcting it would change the window set — so a future "tidy-up" fails loudly.
- A test that pins the reverse-strand previous-START comparison the same way.
- Real-data test, skipped without the cached GFF3: asserts **3,496** windows covering **625,255** bases, the figures measured for this plan.

### Step 2: Locate every tile by sequence search, and refuse the ambiguous
Files: dna_decode/glm/tile_locate.py, tests/test_tile_locate.py
Depends on: none

**What changes:**
- `locate_sequences(seqs, genome) -> dict[str, LocateResult]` with `pos0`, `strand`, `n_occurrences`, built by indexing the query sequences and sliding both strands of the genome once — measured at 5 s for 43,417 queries against 4.6 Mb, where per-query search would not be viable.
- A sequence occurring more than once is **AMBIGUOUS, excluded, and COUNTED** — never resolved by taking the first hit. Same for not-found.
- `LOCATE_FLOOR = 0.95`: `locate_report` REFUSES a verdict below it. Measured 0.9918 unique, so the floor has real headroom and is not fitted to the observation.
- The module NEVER reads the substrate's `start`/`end` columns, and a test asserts that by source inspection — measured agreement is 3.15%, so using them would produce a fully self-consistent and entirely wrong frame.
- Emits a cached locate artifact into gitignored `data/promoteratlas/` so Step 6 does not re-sweep.

**Test strategy:**
- Unit on a synthetic 2 kb genome: forward location, reverse-complement location with the position normalised back to forward coordinates, a planted duplicate excluded and counted, and a absent sequence counted as not-found.
- Non-vacuity: a test that the reverse-strand path is genuinely exercised, by locating a sequence present ONLY as a reverse complement.
- Real-data test, skipped without the cached genome: asserts uniquely-located **43,061** of **43,417** distinct, ambiguous **325**, not-found **31**, and that the top stated-vs-located offsets form families rather than a single constant.
Test-strategy-ignore: data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz

### Step 3: Frame features on located coordinates
Files: dna_decode/glm/frame_features.py, tests/test_frame_features.py
Depends on: Step 1, Step 2

**What changes:**
- `frame_features(located, windows, mask, cds_rows) -> (names, matrix)` emitting, per uniquely-located tile: `in_window_any` (overlaps a corpus window), `in_window_full` (entirely inside), `overlap_fraction`, `signed_distance_to_next_cds_start` on the matching strand, and `inside_cds`.
- Both `in_window_any` and `in_window_full` ship because their variance differs sharply (66.47% vs 19.13%) and neither is obviously the right granularity; the bar grades on the named primary and reports the other.
- Strand handling is explicit: distance is measured toward the gene the window belongs to, so a minus-strand window's "downstream" runs the other way along the genome.
- REFUSES to emit for a tile that is not uniquely located, rather than imputing — excluded tiles are returned as a counted set.
- Pure and torch-free, so it imports without the ML extra.

**Test strategy:**
- Unit on the synthetic genome from Steps 1–2: a tile fully inside a window, one straddling an edge (fraction strictly between 0 and 1), one far from any gene, one inside a CDS, and the sign convention on both strands.
- Arithmetic reconciliation: emitted rows plus excluded rows equals the input count exactly. An earlier abstention arm in this repo turned true positives into an apparent bias through a quietly reduced denominator, and that defect was caught by this arithmetic rather than by reading code.
- Real-data test, skipped without the cache: `in_window_any` fires on **28,622** tiles and `in_window_full` on **8,239**, and neither is constant.

### Step 4: The grader — validity-first, every outcome reachable
Files: dna_decode/glm/frame_verdict.py, tests/test_frame_verdict.py
Depends on: none

**What changes:**
- A pure `grade(FrameInput) -> FrameVerdict` over the seven named outcomes, in a fixed `CHECK_ORDER` with the four indeterminates evaluated BEFORE the three substantive outcomes.
- Mirrors `dna_decode/glm/exposure_verdict.py` in shape rather than importing it, because the inputs differ; the shared discipline is the ordering plus reachability.
- Each outcome's detail string carries its own bound, so a downstream quoter cannot drop it — the pattern that keeps `EXPOSURE_EXCLUDED_LINEAR` from being restated as "exposure excluded".
- Depends on nothing: it reads a dataclass of medians, so it is Wave 0 and can be written and proven before any measurement exists.

**Test strategy:**
- One constructed input per outcome, held as DATA so the reachability test cannot silently cover fewer outcomes than the order declares.
- The reachable set must EQUAL `CHECK_ORDER`, not merely be contained in it. Run 1 of the parent arm emitted its strongest verdict over a failing gate whose predicate could never fire; containment would permit exactly that.
- Ordering tests: harness drift beats every substantive outcome; an uninformative frame beats `FRAME_EXPLAINS`; a non-strengthened baseline beats both.

### Step 5: Freeze the frame-conditioning acceptance bar
Files: scripts/freeze_frame_bar.py, wiki/promoteratlas_frame_acceptance_bar.json
Depends on: Step 1, Step 2, Step 3, Step 4

**What changes:**
- Emitted BY CODE, after feasibility is established (window count, locate rate and feature variance all known) and BEFORE any conditioning number exists. Specifying a threshold before asking what the substrate permits produced four amended thresholds in a single session in this repo.
- Inherits the **0.02** margin; asserts it is still present in the parent bar rather than restating it.
- Fixes the primary frame granularity (`in_window_any`) in advance, so the 66%/19% choice cannot be made after seeing which one gives a better answer.
- Records a forecast with explicit probabilities over the named outcomes. The standing record is **2 of 7** and the errors have been pessimistic; the forecast exists to be scored, not deferred to.
- Pre-commitments: (a) the frame must be shown INFORMATIVE before any absorption claim is read, and the response to an uninformative frame is a better positional feature, never a weaker test; (b) `FRAME_EXCLUDED_CONDITIONAL` is bounded — it excludes the corpus-window rule as replicated, at 200 bp with a 5 bp intergenic gate, on one organism's annotation, and may never be written as "position excluded"; (c) the stratified cell's baseline must be re-measured ON THE SAME STRATUM, never compared to the full-cohort 0.5076, because stratifying changes both n and the peak composition; (d) `arms.native150` keeps its measured value whatever the outcome.

**Test strategy:**
- Schema test: every outcome carries a numeric predicate, `verdict_check_order` lists each exactly once and matches the Step-4 grader's `CHECK_ORDER`, and the forecast probabilities are present and sum sanely.
- A test that the bar names its primary granularity, so the post-hoc choice is structurally unavailable.
Test-strategy-ignore: wiki/promoteratlas_arm_acceptance_bar.json, wiki/promoteratlas_exposure_acceptance_bar.json

### Step 6: The experiment — three cells, graded
Files: scripts/promoteratlas_frame_conditioning.py, wiki/promoteratlas_frame_conditioning_2026-11-01.json
Depends on: Step 3, Step 5

**What changes:**
- Reuses the cached embedding pass; no re-embedding, no GPU, no network.
- **Anchor first:** reproduces `arms.native150` 0.5076 on the committed splits before any new number is read, exactly as the reconcile gate did.
- **C1** — predict `in_window_any` from the 128-dim embedding, held out over the committed splits, reporting AUROC beside the GC-only equivalent. The positional restatement of the measured 0.9156 separability.
- **C2** — fit the identical ridge on `gc`, on `gc + frame`, and on `atlas`, and report the atlas margin over each. The strengthened baseline must show a better TRAINING fit than plain GC, else the run returns `INDETERMINATE_BASELINE_NOT_STRENGTHENED`.
- **C3** — restrict to in-window tiles, then to out-of-window tiles, and re-measure BOTH atlas and GC within each stratum, reporting the within-stratum margin. This controls for frame membership by construction, so it does not require the frame to be linearly modellable.
- Excluded-tile accounting is explicit in the artifact: ambiguous, not-found and non-unique counts reconcile with the scored count.
- Grades through the Step-4 grader against the Step-5 bar, validity checks first.
- Per-tile vectors and the locate cache stay in gitignored `data/promoteratlas/`; the wiki artifact carries medians only. A 102,755-decimal dump broke two contract-number guards on 2026-11-01.

**Test strategy:**
- Unit: the strengthened baseline is genuinely strengthened — adding a column cannot lower the training fit — so a "no absorption" result cannot come from a baseline that ignored the feature.
- Unit: stratum sizes plus excluded counts reconcile with the full cohort exactly.
- Integration: the unprojected anchor reproduces 0.5076, catching a harness that drifted since the preceding run.
Test-strategy-ignore: wiki/promoteratlas_frame_acceptance_bar.json, data/promoteratlas/tile_embeddings_native150_n44106.npy

### Step 7: Record the verdict on every surface that restates the claim
Files: wiki/promoteratlas_frame_result_2026-11-01.md, wiki/promoteratlas_tile_arm_2026-11-01.json, project_state/glm-learned-representation-2026-10-07.md, tests/test_promoteratlas_frame_claims.py
Depends on: Step 6

**What changes:**
- Memo ordered DECISIONS NEEDED, then CORRECTIONS, then findings, with the verdict and its bound in the first screen.
- The machine-readable surface is updated BEFORE the prose: the arm artifact's `corpus_selection_exposure_resolved` block gains a `frame_conditioning` sub-block, additively, with `arms{}` asserted byte-identical after the write.
- The exposure memo's bounded limit gains a positional clause; the original text is preserved and amended, never deleted.
- Retires candidate 7 in the G-A candidate table with a `[done:<id>]` marker from `candidate_replay.done_marker`, and appends one Action-Log row through the sanctioned cheap path with a class from the 7-value enum. The five pre-existing done rows are LEFT ALONE: a 2026-10-19 session measured that the ranker's "journaled done but still listed as pending" advisory is unclearable by any safe action, and redoing that is foreclosed work.
- Scores the Step-5 forecast and updates the standing record from 2 of 7.
- Writes no float dump into `wiki/`.

**Test strategy:**
- Assert every number in the memo appears in the artifact it cites — the contract-number-audit discipline.
- Assert the `FRAME_EXCLUDED_CONDITIONAL` prose never appears without its bound, the guard shape that keeps the linear bound attached to the exposure claim.
- Run `scripts/contract_number_audit.py` and assert the verdict is still `NOTHING_TO_ADJUDICATE` and the decoy-pool decimal count has not materially grown.
- Run the affected-file subset only, never the full suite while a torch job holds a model — that dies with an access violation that reads like a test failure.
- Re-pin `tests/test_arms_ledger_stack_arm.py` ONLY if its subject actually moved. It should NOT: this is a limit on an existing arm, so `best` and the ceiling verdict must be unchanged, and a passing guard is the evidence.
Test-strategy-ignore: wiki/promoteratlas_frame_conditioning_2026-11-01.json, wiki/promoteratlas_exposure_result_2026-11-01.md, scripts/contract_number_audit.py, tests/test_arms_ledger_stack_arm.py

## Execution Preview

Wave 0 (3 parallel):  Step 1 — corpus-window replica, Step 2 — tile location, Step 4 — the grader
Wave 1 (1):           Step 3 — frame features
Wave 2 (1):           Step 5 — freeze the bar
Wave 3 (1):           Step 6 — the experiment
Wave 4 (1):           Step 7 — record

Critical path: Step 1 → Step 3 → Step 5 → Step 6 → Step 7 (5 waves)
Max parallelism: 3 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — **The substrate's stated coordinates are wrong and using them yields a self-consistent, entirely wrong frame.** Measured agreement 63/2000 = 3.15%, with offsets in families (+1976/+1977/+1978, +776/+777) rather than a constant. Resolved by locating via sequence search only, plus a Step-2 source-inspection test asserting the stated columns are never read.
- Severity: high — **A "corrected" window replica would measure a rule nobody trained on.** The authors' extractor carries two apparent defects (a skipped gene not updating `last_end`; the reverse branch comparing against the previous gene's start). Step 1 replicates both and pins each with a test that fails if a future edit tidies them.
- Severity: medium — **`genomewide.py` is embedded base64+gzip in the arm kernel with a byte-identity assert in `scripts/promoteratlas_arm_dryrun.py`.** Any edit silently invalidates the dry run. All new code lands in new modules; verified by an empty diff at the end of every wave.
- Severity: medium — **The stratified cell changes n and peak composition, so its baseline must be re-measured within the stratum.** Comparing a within-stratum atlas number to the full-cohort 0.5076 would manufacture a difference. Pinned as a Step-5 pre-commitment, not left to the writer.
- Severity: medium — **689 of 44,106 rows share a sequence with another row** and `leave_peak_out` splits on peak, so a duplicated sequence can straddle a split. This is a property of the substrate and of every arm on this rung, including the incumbent 0.5076. Reported in the artifact; deliberately NOT fixed here, since changing the split would make every number incomparable with the committed tables.
- Severity: medium — **`FRAME_EXCLUDED_CONDITIONAL` is bounded by construction.** It excludes the corpus-window rule *as replicated* — 200 bp, 5 bp intergenic gate, one organism's annotation — and does not exclude positional explanations generally. Fixed as a pre-commitment before the run.
- Severity: medium — **325 ambiguous and 31 not-found tiles are excluded, and exclusions have produced a false finding in this repo before** via a quietly reduced denominator. Step 3 reconciles emitted plus excluded against the input count exactly.
- Severity: low — **`genomewide.DEFAULT_CACHE` points at the detached `D:` volume**; every call passes `cache_dir` with `allow_download=False`.
- Severity: low — **`torch` reports `cuda False` on this host**, so all compute is CPU. No step re-embeds (the 15.7-minute pass is cached), so the plan needs no GPU and no Kaggle push.
- Severity: low — **External tool surfaces.** `plan-toolkit.py check|waves|validate` verified by direct invocation during this run; `python3` is the Windows Store alias stub on this host and `python` must be used. `sentrux` resolves; `codex` is absent and unused. No step in this plan touches a network or service surface, so no external schema is left unconfirmed.
- Severity: low — **No restructuring was applied during wave grouping.** Wave 0's three steps touch six disjoint paths and no intra-wave overlap exists. The one ordering choice worth naming: Step 4 (the grader) was placed in Wave 0 rather than beside the bar, because it is a pure function over a dataclass and proving its reachability before any measurement exists is strictly better than proving it after.

## Open Questions

- **Which frame granularity is primary — `in_window_any` (66.47%) or `in_window_full` (19.13%)?** Step 5 fixes `in_window_any` in advance so the choice cannot be made after seeing the answer, but the decision is genuinely arguable: `any` has far more variance while `full` matches the corpus construction more literally. Both ship in the artifact; only the primary grades.
- **If the verdict is FRAME_EXPLAINS, does `arms.native150` stay at 0.5076?** This plan says yes — a measured metric is not retro-edited and the discount lives in the honest limit. That is the conservative reading, and it leaves `best` on the ceiling criterion pointing at a number the same artifact discounts. That is an acceptance-surface question rather than a technical one, and it is the same open question the exposure plan raised; it has not been settled.
- **Should the 689 sequence-duplicate rows be excluded from the rung's splits?** Doing so would change every committed number on this rung, so it is not in scope here. Whether it is worth a separate re-measurement of the incumbent tables is a user call.

## Architecture Impact

**Expected structural changes:**
- Modules added: `dna_decode/glm/corpus_window.py`, `dna_decode/glm/tile_locate.py`, `dna_decode/glm/frame_features.py`, `dna_decode/glm/frame_verdict.py` — four new leaves in an existing package. Modules modified: none under `dna_decode/`.
- Dependency direction: `corpus_window`, `tile_locate` and `frame_verdict` depend on the standard library (plus the already-parsed annotation table) and nothing else; `frame_features` depends on the first two. All are consumers; scripts depend on `dna_decode/glm/`, which is the existing direction, and nothing under `dna_decode/` gains a dependency on a script.
- Boundary changes: none. All four new modules are torch-free, so the only torch dependency on this path stays isolated in the pre-existing `atlas_embed`.

**Risk areas:**
- `dna_decode/glm/genomewide.py` is the shared core abstraction in play and is read-only in every step, for the kernel-identity reason above.
- `dna_decode/glm/expression.py` gains no new callers here (Step 6 goes through `atlas_embed.ridge_spearman`, which mirrors it) and is not modified; the ridge protocol must stay byte-identical or every comparison loses its anchor.
- Potential cycles: none available. The four new modules form a two-level DAG with no back edge.

**Sentrux gates:**
- Baseline before Wave 0: yes — `sentrux` resolves on this host.
- Compare after each merged wave: no. Compare once after Wave 1, when all four new modules exist; per-wave comparison on a purely additive change is noise.
- Acceptable delta: a minor score movement reflecting four added leaf modules. A new cycle, a forbidden-layer violation, or any new edge into `genomewide.py` is not.
- Regression-response stance: investigate. On an additive, cycle-free change a flagged regression most likely means an unintended import, so the first action is to locate the new edge rather than accept or revert it.

**Anti-drift checks:**
- `dna_decode/glm/genomewide.py` and `dna_decode/glm/expression.py` diff empty at plan end.
- `scripts/promoteratlas_arm_dryrun.py` still passes its embedded-payload byte-identity assert.
- All four new modules import in an environment without torch, proving the isolation held.
- No new import of a script module appears anywhere under `dna_decode/`.
- The frozen five plus `uv.lock` are byte-unchanged and the active 2026-08-31 prospective lock re-verifies.

## Verification

1. **Frozen-surface and lock check, before and after:** the frozen five plus `uv.lock` diff empty; `genomewide.py` and `expression.py` diff empty; the active prospective lock re-verifies.
2. **Window replica:** 3,496 windows covering 625,255 bases on the cached annotation, and both replicated quirks pinned by tests that fail if corrected.
3. **Locate rate:** 43,061 of 43,417 distinct sequences uniquely located, 325 ambiguous and 31 not-found counted, and the unique rate at or above the 0.95 floor.
4. **Feature variance:** `in_window_any` fires on 28,622 tiles and `in_window_full` on 8,239; neither is constant.
5. **Bar reachability:** every outcome reachable by a constructed input, and the reachable set EQUALS the grader's `CHECK_ORDER`.
6. **Anchor:** the experiment reproduces `arms.native150` 0.5076 before reading any new number.
7. **Non-vacuity:** the frame is shown informative, and the strengthened baseline shows a better training fit than plain GC, before any absorption claim is graded.
8. **Accounting:** emitted plus excluded equals the input count; stratum sizes plus exclusions equal the full cohort.
9. **Claim propagation:** the arm artifact, the exposure memo, the project-state ledger and the new memo all carry the same verdict and the same bound, machine-readable surfaces first, with `arms{}` asserted byte-identical.
10. **Decoy-pool hygiene:** `scripts/contract_number_audit.py` still returns `NOTHING_TO_ADJUDICATE` and `wiki/` has gained no float dump.
11. **Test accounting:** the affected-file subset passes with a stated before/after count, run while no torch job holds a model, and `tests/test_arms_ledger_stack_arm.py` passes UNCHANGED.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,severity-high -->
