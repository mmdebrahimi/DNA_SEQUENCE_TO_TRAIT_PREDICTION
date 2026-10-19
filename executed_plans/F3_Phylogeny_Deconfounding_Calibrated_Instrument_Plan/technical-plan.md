# F3 — Phylogeny De-confounding: a Calibrated Instrument, Not Another Arm

## Lens status

Conversation-driven dispatch. Lenses applied: reuse-survey (graphify absent; filesystem + grep over
`dna_decode/eval/` and `scripts/*_arm.py`), external-surface verification (`cv_ridge_gp` /
`greedy_representative_clusters_from_matrix` read in source, not inherited from a reference),
PER replay over `wiki/decisions-log.md` (HIGH slots + module overlap on `dna_decode/eval/`).
No `project-rules.md`, no `AUDIT_REPORT.md`, no `DESIGN.md` at root — `## Project Rules Check` omitted.
`sentrux` not installed.

## Problem Statement

Every organism-level g→p number in this repo is produced under **random k-fold CV**. On a natural
population that measures *taxonomic* signal and *functional* signal together, and the project has already
paid for this once: the Arabidopsis FT10 arm learned population structure rather than flowering-time
biology, and the de-confounded within-group r² went **negative** (`wiki/phase2_arabidopsis_result_2026-06-12.md`).
The corrected regime map (`wiki/organism_gp_regime_correction_2026-08-29.md`) names the discriminating
variable as **population design**, not organism complexity — which means the repo's defence against the
failure is currently *substrate choice*, not an instrument.

`dna_decode/eval/genomic_prediction.py` says so in its own module docstring:

> *"On natural-population genotypes the same code would need within-group de-confounding (see the
> Arabidopsis closed negative)."*

The gap is named and unbuilt. `cv_ridge_gp` takes no `groups` argument; `_kfold_indices` is a plain
`rng.permutation`. So the de-confounding discipline exists as prose in one docstring and as a bespoke
script for one unrelated axis (`scripts/crossaxis_lineage_deconfound.py`, Mash clades for AMR cohorts),
and the four constructed-panel arms (`yeast_bloom_gp_arm`, `bxd_gp_arm`, `arabmagic_gp_arm`,
`dgrp_learned_decoder`) cannot ask the question at all.

**What this plan builds is an instrument, and it does not ship until it is calibrated against a
known-answer control pair.** The deliverable is *not* another panel result. It is: a grouped-CV seam
threaded through the existing GP helper, a genotype-derived grouping module, and a **double control**
that can fail in both directions — over-aggression (destroys signal where there is no structure to
remove) and insensitivity (fails to respond where structure is measured to exist).

**Non-goals.** No new substrate is fetched. No published number is re-scored. No learned model is
trained. No frozen file is touched (`amr_rules.py` and the four other lock-pinned files are untouched;
`dna_decode/eval/genomic_prediction.py` is not in the prospective-lock surface).

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

Captured before any edit, and recorded in the artifact so the delta is auditable:

- `uv run pytest tests/ -q` → current pass/skip/fail counts (CLAUDE.md records the suite as fully green:
  3875 passed / 11 skipped / 0 failed as of 2026-08-26; **re-measured, not quoted**).
- `cv_ridge_gp` / `cv_model_gp` current public signatures — no `groups` parameter. Pinned by a test added
  in Step 2 asserting the *new* parameter is keyword-only and defaults to `None`, so every existing
  caller's behaviour is byte-identical.
- The four panel arms' committed artifacts (`wiki/yeast_bloom_gp_arm_*.json`,
  `wiki/arabmagic_gp_arm_*.json`, …) are the **unchanged-output anchor**: Step 2 must leave
  `groups=None` numerically identical, verified by re-running one arm and diffing its JSON.
- `dna_decode/eval/regime.py` row count + the live `organism_transfer` census (CLAUDE.md: 8 of 9 regimes
  carry no held-out-organism evidence) — re-derived from the table, not restated.

## Verification Signal

The plan succeeds when **all** of the following hold:

1. `groups=None` is a provable no-op: one existing panel arm re-run produces a byte-identical artifact.
2. A grouped run's permutation null is demonstrably computed under the **grouped** folds (asserted
   directly, by capturing the `folds` object the null consumes — not inferred from the number).
3. The inner alpha split is group-aware, pinned by a test that **fails** on the positional-slice version.
4. The grouping module refuses, rather than returns, on: fewer groups than splits, a group larger than
   `n/n_splits`, >50% singletons (the over-split guard), and a non-finite distance.
5. The calibration harness emits one of its **six pre-registered verdicts**, every one of them named and
   committed to a frozen dict *before* the run, and at least two of them adverse to the instrument.
6. Whatever the verdict, `dna_decode/eval/regime.py` gains the de-confounding protocol's status and the
   `wiki/` artifact exists with its own limits section.

**A `blocked:` verdict is a success of this plan, not a failure of it** — the same shape as
`INDETERMINATE_INSUFFICIENT_RUNGS` in the essentiality transfer ladder, which was that plan's
pre-committed outcome.

## Implementation Steps

### Step 1: Substrate integrity and structure census (W0 gate)
Files: scripts/f3_substrate_census.py, tests/test_f3_substrate_census.py
Depends on: none

**What changes:**
- New read-only probe over the four cached panels (`bloom`, `arabmagic`, `bxd`, `dgrp`). For each, resolve
  its genotype + phenotype paths from the owning arm script (single source of truth — do not re-hardcode).
- **Integrity, not presence.** Parse the genotype matrix **in full** and record: row count, per-row field-width
  agreement, marker count, individual count, allele-coding alphabet, NaN/missing fraction, and — for any
  gzipped member — a full decompress. Emit `integrity: intact | unusable:<reason>`. A panel that is present
  but truncated/short-rowed/empty is `unusable`, **named**, and excluded.
  *Rationale, from this repo's own record:* `file exists` is not `file intact` — 23 cached GPS assemblies
  passed an existence check while 11 were 199-byte HTTP 403 pages and 10 were truncated gzip. A
  head-only read is the wrong guard: a truncated gzip decompresses its opening blocks perfectly and
  raises only at the end.
- **Structure census.** For each intact panel compute, from the genotype matrix alone:
  `pc1_variance_fraction` (leading eigenvalue share of the centred genotype covariance — **no free
  parameter**, the primary statistic) and, as corroboration, `n_clusters` /
  `largest_cluster_fraction` / `singleton_fraction` from Step 3's grouping module.
- Emit `wiki/f3_substrate_census_<date>.json`. Exit 0 always (a census, not a gate on the suite), but
  the harness in Step 5 **refuses** on a census that reports no intact panel.

**Test strategy:**
- Synthetic fixtures for every `unusable` branch: truncated gzip (must be caught by full decompress and
  **not** by a 2048-byte head read — pin both, so the weak guard is a failing test), ragged row widths,
  zero data rows, all-missing column.
- `pc1_variance_fraction` on two constructed matrices with known structure: an i.i.d. matrix (low) and a
  two-block matrix (high). Assert the ordering, not a value.
- A panel whose paths are absent reports `unusable:absent` and does not raise.

### Step 2: Grouped-fold seam in the GP helper
Files: dna_decode/eval/genomic_prediction.py, tests/test_genomic_prediction_grouped.py
Depends on: none

**What changes:**
- Add `_grouped_fold_indices(groups, k, seed) -> list[np.ndarray]`: assign whole groups to folds
  (descending group size, greedily to the lightest fold — deterministic), returning the same
  list-of-test-index-arrays shape `_kfold_indices` returns, so **every downstream consumer is unchanged**.
- Thread a keyword-only `groups: Sequence[Hashable] | None = None` through `cv_ridge_gp` and
  `cv_model_gp`. `None` → existing `_kfold_indices` path, bit-for-bit.
- **Refusals** (raise, never degrade): `len(groups) != len(y)`; fewer distinct groups than `n_splits`;
  any single group holding > `1/n_splits` of the rows after NaN-filtering (a dominant group makes
  grouped CV effectively leave-a-fraction-out — the recorded Mash finding where a 57.8% clade made
  LOMO-clade-out into leave-58%-out).
- **Apply the `groups` mask alongside the `y`-NaN mask.** `cv_*_gp` already drops `~isnan(y)` rows;
  `groups` must be filtered identically or the labels silently misalign with the rows.
- **Fix the inner alpha leak** in `_cv_predictive_r`: when `groups` is supplied, build the inner
  validation split by whole groups drawn from the outer training groups, not by positional slice.
- Add `fold_scheme: str` (`"random"` | `"grouped"`) and `n_groups: int | None` to `GPResult` +
  `as_dict()`, so an artifact can never be ambiguous about which scheme produced it.

**Test strategy:**
- **No-op proof:** `cv_ridge_gp(X, y, groups=None)` and `cv_model_gp(...)` equal the pre-change values to
  full float precision on a fixed seed (golden values captured in the test, not recomputed from the
  implementation under test).
- **Group purity:** for every fold, no group appears in both train and test. Checked on an uneven group
  size distribution.
- **The null inherits the folds:** monkeypatch `_cv_r_with` to capture the `folds` argument on every call
  and assert the null's folds are *the same object* as the measured run's. This is the property the
  whole design leans on; it is asserted directly, not inferred.
- **Inner-leak regression:** a fixture where the positional inner split necessarily mixes a group across
  inner train/validation; the test fails against the old slice and passes against the group-aware split.
- Each refusal raises with the offending quantity in the message.
- Cross-check fold assignment against `sklearn.model_selection.GroupKFold` on a benign case (equal group
  sizes), documenting in the test *why* GroupKFold is not the implementation.

### Step 3: Genotype-derived grouping module
Files: dna_decode/eval/panel_structure.py, tests/test_panel_structure.py
Depends on: none

**What changes:**
- `ibs_distance_matrix(X, ids) -> DistanceMatrix` — pairwise allele-sharing distance (fraction of markers
  at which two individuals differ, missing-aware), returning `clonality.DistanceMatrix` so the existing
  clusterer consumes it directly.
- `derive_threshold(distance_matrix) -> float` — the threshold is **derived from the distance
  distribution** (a low percentile of the off-diagonal distances), never a swept grid.
  *Rationale, from this repo's own record:* the Mash threshold sweep minimised
  `variance_ratio = intra/inter`, which → 0 as clusters → singletons, so MINIMIZE was monotonically
  biased toward over-splitting; extending the grid downward picked 0.002 → 101 clades / 56.5% singletons.
  The qualifying criteria guarded only against a dominant cluster, never against over-splitting.
- `panel_groups(X, ids) -> PanelGroups` — derive threshold, call
  `greedy_representative_clusters_from_matrix`, then enforce **both** guards and **refuse**:
  `largest_cluster_fraction > 0.60` (dominant-cluster) **and** `singleton_fraction > 0.50`
  (over-split — the guard the Mash work lacked).
- `pc1_variance_fraction(X) -> float` — the parameter-free structure statistic Step 1 leads with.
- Pure, offline, no Docker, no network. Mash is **not** used: these are dense genotype matrices, and IBS
  on the markers is the direct measurement where Mash would be a sketch of a proxy.

**Test strategy:**
- IBS on hand-computable 3×4 matrices, including a missing-data case.
- Two-block synthetic panel: clusters recover the blocks; an i.i.d. panel refuses with
  `singleton_fraction > 0.50`; a panel with one 70% block refuses with the dominant-cluster guard.
- **Chaining regression:** a deliberate A~B~C chain where single-linkage merges all three and greedy-rep
  does not — asserting the reused primitive's property is actually relied on.
- Determinism: identical output across two calls and under input row permutation.

### Step 4: Calibration harness with pre-registered verdicts
Files: scripts/f3_deconfound_calibration.py, tests/test_f3_deconfound_calibration.py
Depends on: Step 1, Step 2, Step 3

**What changes:**
- **Control-pair selection rule, frozen in code before any run:**
  - **Negative control** = a panel whose structure is randomised **by construction** — a bi-parental
    segregant cross (`bloom`). This is a structural property of the design, not a measurement.
  - **Positive control** = among intact panels, the one with the **highest measured**
    `pc1_variance_fraction`, required to clear a pre-registered structure floor. If none clears it →
    `blocked:no-structured-panel`, with the unblock named (fetch the 1001G natural-population genotype
    matrix) and **no metrics emitted**.
- Per (panel, trait): run `cv_ridge_gp` twice on identical rows — `groups=None` and
  `groups=panel_groups(...)` — each with its own `n_perm=100` null, which is grouped-or-random
  automatically via the `folds` seam.
- **Paired per trait, never a difference of aggregate medians.** Per-trait
  `drop = r_random - r_grouped`; a sign test across that panel's traits. Report
  `n_traits`, `n_dropped`, `median_drop`, and the sign-test p for **each** control separately.
- **The verdict criterion is differential, and asserts no magnitude.** Six states, ordered, first match
  wins:
  | verdict | condition |
  |---|---|
  | `blocked:no-intact-panel` | census reports no intact panel |
  | `blocked:no-structured-panel` | no intact panel clears the structure floor |
  | `INDETERMINATE_UNDERPOWERED` | either control has too few traits for a sign test |
  | `INVALID_OVER_AGGRESSIVE` | negative control **stops** beating its own grouped null |
  | `INVALID_INSENSITIVE` | positive control's drop is **not** greater than the negative control's |
  | `VALID_DIFFERENTIAL_RESPONSE` | negative control keeps beating its own grouped null **and** the positive control's drop exceeds it |
- **Why the criterion is differential and not "the positive control falls to its null":** a natural panel
  may carry structure **and** real functional signal, and a *correct* de-confounder removes the former
  and leaves the latter. A "falls to its null" bar would therefore condemn a working instrument. What a
  structure-de-confounder must do is respond **more** on a structured panel than on one where structure
  is randomised by construction — so the criterion is the contrast between the two controls, which is
  what makes them a pair rather than two runs.
- **This is instrument calibration, not an experiment with a foregone outcome.** The distinction matters
  because this repo has just recorded the opposite failure (a revival trigger whose experiment could not
  surprise anyone). Here both adverse verdicts are genuinely reachable: the de-confounder could destroy
  the bi-parental panel's signal (the `groups` partition on a cross is near-arbitrary, so an
  over-aggressive partition is a live possibility, not a strawman), and it could fail to respond on the
  structured panel (if IBS clustering does not recover the structure that `pc1_variance_fraction`
  detects). Using a known-answer pair to *calibrate an instrument* is a positive+negative control pair;
  it is not the same act as citing a foregone result as new information.

**Test strategy:**
- A synthetic end-to-end fixture pair — one simulated bi-parental cross, one simulated structured panel
  with a *planted* structure-only phenotype — driving the harness to `VALID_DIFFERENTIAL_RESPONSE`
  offline, with no cached panel present.
- A fixture driving **each** adverse verdict: signal-destroying partition → `INVALID_OVER_AGGRESSIVE`;
  structure-blind partition → `INVALID_INSENSITIVE`. A verdict state with no test that reaches it is
  decoration.
- The verdict function is pure and tested on hand-built inputs; ordering is pinned (a case satisfying two
  conditions returns the earlier one).
- `blocked:` states emit **no** metrics (pinned by asserting the keys are absent, not empty).
- The structure floor and the sign-test minimum are module constants read by the test, never restated in
  it (the shipped-default-is-a-promise failure: a CLI default silently overrode a constant lowered five
  days earlier).

### Step 5: Run the calibration, emit the artifact
Files: wiki/f3_substrate_census_2026-10-09.json, wiki/f3_deconfound_calibration_2026-10-09.md, wiki/f3_deconfound_calibration_2026-10-09.json
Depends on: Step 4

**What changes:**
- Run Step 1's census, then Step 4's harness, on the real cached panels.
- The memo leads with the **verdict and its scope**, and carries a limits section stating at minimum:
  the two controls have **different traits**, so the between-control comparison is two-sample, not paired
  (the within-control tests are paired); the positive control's structure is **measured**, while the
  negative control's randomisation is **by construction** — a deliberate asymmetry in evidence strength;
  IBS-greedy-rep grouping is one grouping of many; and `n_perm=100` bounds the resolution of every
  `beats_null`.
- If the verdict is `blocked:`, the memo says what is blocked and what unblocks it, and **no instrument
  claim is made**. CLAUDE.md is not edited to claim a capability the run did not establish.

**Test strategy:**
- Manual: the emitted JSON parses (the hand-edited-JSON trap — `wiki/fba_regulatory_conditional_test_2026-08-12.json`
  was invalid JSON for a day, cited as evidence, and nothing read it).
- The existing `test_every_fba_wiki_artifact_is_parseable_json`-class guard is extended to the new
  artifacts, so a hand-edit cannot ship unparseable.
- Re-run determinism: two runs at the same seed produce identical verdicts and metrics.

### Step 6: Register the protocol on the trust surface
Files: dna_decode/eval/regime.py, tests/test_regime_boundary.py, wiki/decisions-log.md, CLAUDE.md
Depends on: Step 5

**What changes:**
- Record the de-confounding protocol's calibration status in `regime.py` alongside the regime rows, so a
  future proposal that reaches for grouped CV finds out whether the instrument is calibrated **and on
  which panels** — not merely that a function exists.
- Extend the regime screen so a natural-population proposal returns the calibration status rather than
  only `REQUIRES_DECONFOUNDING`.
- Append a `wiki/decisions-log.md` entry.
- CLAUDE.md: a single bullet, written **from the artifact**, carrying the verdict, the control pair, the
  differential criterion, and the named limits. If the verdict is `blocked:`, the bullet says so.
  Under-claiming is as much a trust-surface falsehood as over-claiming — but so is claiming a calibrated
  instrument on a blocked run.

**Test strategy:**
- `test_regime_boundary.py` gains a case asserting the three historical compressions still fail loudly
  (unchanged behaviour) **plus** one asserting the new status is derived from the artifact, not restated
  as a literal — the `contract_number_audit` discipline.
- The existing CLAUDE.md word-ceiling guard must stay green **unmodified** — this step does not edit it.
  The bullet is budgeted to fit; if it does not, an existing stale bullet is pruned rather than the
  ceiling raised. (Raising the ceiling to fit new prose is how the guard stops being a guard.)

## Execution Preview

| Wave | Steps | Parallelism |
|---|---|---|
| 0 | Step 1, Step 2, Step 3 | 3 |
| 1 | Step 4 | 1 |
| 2 | Step 5 | 1 |
| 3 | Step 6 | 1 |

Total waves: 4. Max parallelism: 3. Critical path: Step 1 → Step 4 → Step 5 → Step 6 (4 steps, per the
toolkit's deterministic `waves` output; Steps 2 and 3 are equally long but the toolkit resolves to Step 1).

Wave 0 has no intra-wave file overlap: Step 1 writes `scripts/f3_substrate_census.py`, Step 2 writes
`dna_decode/eval/genomic_prediction.py`, Step 3 writes `dna_decode/eval/panel_structure.py`, each with its
own test file. Step 1 *calls* Step 3's module at runtime for its corroborating cluster statistics, so if
executed in parallel its census must degrade to `pc1_variance_fraction` only when `panel_structure` is
absent — or Step 1 moves to Wave 1. **Prefer moving Step 1 to Wave 1** in sequential execution; the
dependency is declared as `none` because the integrity half (the gate) has no dependency, and that half
is what blocks Step 4.

## Risk Flags

- **The positive control may not exist on disk.** Measured: there is no cached natural-population
  Arabidopsis genotype matrix; `arabmagic` is a constructed multi-parent panel. `bxd` (recombinant inbred
  lines) and `dgrp` (inbred lines from one natural *Drosophila* population, with documented cryptic
  relatedness) are the live candidates, and **their integrity is unverified** — which is precisely why
  Step 1 is a gate and the control pair is selected by measurement. If none clears the structure floor,
  the honest outcome is `blocked:no-structured-panel`.
- **`arabmagic` is constructed, so a structured-panel result from it is weaker than a natural-population
  one.** If it becomes the positive control, the memo must say the structure is *residual family
  structure in a constructed cross*, not natural population structure, and the instrument's calibration
  is correspondingly scoped.
- **IBS-greedy-rep is one grouping.** A failure to respond could be the grouping's blindness rather than
  the instrument's insensitivity. Step 4's `INVALID_INSENSITIVE` verdict text must name this alternative
  explicitly rather than concluding against the de-confounder.
- **The inner-alpha fix changes the ridge path for grouped runs only**, but `_cv_predictive_r` is shared.
  The no-op proof (Step 2) is what bounds the blast radius; if it fails, stop and re-derive rather than
  adjusting the golden values.
- **`n_perm=100` on four panels × several traits × two schemes** is the cost centre. Each permutation
  re-runs a full k-fold ridge. If wall-clock becomes a problem, reduce **panels**, never `n_perm` — the
  null's resolution is the instrument.
- **No frozen file is touched**, but Step 6 edits `CLAUDE.md` and `regime.py`. Verify the prospective lock
  still verifies after Step 6 (`verify_lock`), since the five pinned files must be byte-unchanged — they
  are not in this plan's file list, so this is a confirmation, not a risk.
- `sentrux` not installed → no architecture gate available for this change.

## Open Questions

1. **Which panel becomes the positive control is deliberately left to Step 1's measurement**, not decided
   here. The alternative — naming it now — is what produced the falsified premise this plan corrects.
2. **Is `blocked:no-structured-panel` an acceptable terminal for this plan, or should it trigger fetching
   the 1001G genotype matrix?** Fetching is a substrate acquisition (large download, new cache, new
   provenance) and is scoped **out**; whether to authorise it is a user call, not an executor one.
3. **Should the calibrated grouped-CV path become the *default* for natural-population panels?** That
   would change what existing arms report and is an authority decision; this plan only builds and
   calibrates the instrument.
4. **The structure floor's value.** It is a module constant with a stated rationale, but it is
   **asserted, not derived** — the same honest status as the `>60%` dominant-source bar elsewhere in this
   repo. Deriving it would need more than four panels.

## Architecture Impact

- **Boundaries:** `dna_decode/eval/panel_structure.py` is new and depends only on `numpy` +
  `eval.clonality`. `eval.genomic_prediction` gains no new import (the `groups` plumbing is pure numpy);
  it does **not** import `panel_structure`, so the grouping policy and the CV mechanics stay separable —
  a caller supplies `groups` from any source.
- **Direction of dependency:** `scripts/` → `dna_decode/eval/` only. No `eval` module imports from
  `scripts/`.
- **Frozen surface:** untouched. The five lock-pinned files (`amr_rules.py`, `calibrated_amr_rules.json`,
  `mic_tiers.py`, `shipped_decoder_surface.py`, `cohort_manifest.py`) are not in any step's file list.
- **Gates:** `n/a — sentrux not installed`.

## Verification

1. `uv run pytest tests/ -q` — pass count ≥ baseline, 0 new failures.
2. `uv run pytest tests/test_genomic_prediction_grouped.py tests/test_panel_structure.py tests/test_f3_substrate_census.py tests/test_f3_deconfound_calibration.py -q` — all green.
3. **No-op proof:** re-run one existing panel arm; `git diff --stat` on its `wiki/*.json` artifact shows
   no change. (Not `git diff --quiet` on a path that may not exist — a guard on a missing path passes
   vacuously.)
4. `uv run python scripts/f3_substrate_census.py` → census JSON parses; every panel carries an
   `integrity` verdict.
5. `uv run python scripts/f3_deconfound_calibration.py` → one of the six verdicts, with metrics present
   iff the verdict is not `blocked:`.
6. `uv run python -m scripts.prospective_lock_validate` → lock still verifies (confirmation that Step 6
   touched nothing pinned).
7. `uv run python scripts/contract_number_audit.py` → still `NOTHING_TO_ADJUDICATE`; any number the
   CLAUDE.md bullet quotes is present in the cited artifact.
8. CLAUDE.md word-ceiling guard — still green, unmodified.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
