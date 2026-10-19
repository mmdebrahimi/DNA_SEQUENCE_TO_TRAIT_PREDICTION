# Mission Reconciliation + A Native Clonality Route (v2, decisions answered) — Technical Plan

> v1 with D1-D7 answered: the crosscheck bar is now DERIVED from the committed lineage table rather than asserted, the arm-kill option is withdrawn because code refuses it, and the versioning recommendation is reversed to `git init` in place on measured evidence.

## Lens status
Inputs: conversation (v1 plan at `plans/Mission_Reconciliation_A_Native_Clonality_Route_Plan/technical-plan.md` + the answered open questions D1-D7)
Degradations: sentrux unavailable (Architecture Impact gates record `n/a`) · repo-index script absent at `$HOME/.claude/scripts/repo-index-status.sh` (Step 1 research was direct file reading) · `python3` on PATH is the WindowsApps stub; the toolkit ran through `C:\ProgramData\miniforge3\envs\dna_decode\python.exe` · `sourmash` CLI surface `[unverified]` (wheel verified, package not installed) · D: and E: both offline, so Step 5 cannot produce a number today · `git check-ignore` is UNRUNNABLE against the live rules because the working copy is not a repository, which is itself what Step 9 fixes

## Problem Statement

v1 ended with eight unresolved questions. All eight are now answered, and **three answers contradict v1's own text** — this revision exists to make the plan say what the evidence says.

**What changed, and why each is not a preference:**

- **The crosscheck bar was wrong and the committed data shows it.** v1 proposed ARI >= 0.80 + cluster count within 25% + identical G8 verdict. Reading all 20 cell-threshold pairs out of `wiki/provdisjoint_lineage_metrics.json`: minority effective-lineage counts are `2, 4, 5, 6, 7, 9, 11, 13, 14, 14, 14, 15, 16, 16, 17, 19, 21, 23, 23, 24` against a bar of 3 — **19 PASS, exactly 1 TRIP, median 14**. So "identical verdict per cell" passes almost automatically and is nearly vacuous as a gate; **all the discriminating power sits in one pair** (`Klebsiella x ciprofloxacin @ 0.005`, minority 2). `cluster count within 25%` binds only at m=2 where verdict-identity already binds harder, and permits 18-30 at m=24 where the gate does not care. `ARI >= 0.80` is unanchored: nothing in `g8_dedup_collapses_balance` consumes partition similarity, so ARI 0.75 with every verdict identical is fine while ARI 0.85 with the TRIP flipped is a disaster. Step 4 now freezes a bar derived from the gate's own arithmetic.
- **v1 offered an option the code refuses.** Its D3 and its index entry both offered "kill one of the three unrun arms with scope". `arms_ledger.kill_arm` requires the arm's OWN executed evidence and refuses redundancy reasons; its refusal text says a no-evidence kill is "an explicit user authority decision". That is R6 working. Withdrawn from the plan.
- **The versioning recommendation is reversed, on a measured fact.** v1 leaned toward committing the live copy into `C:\Users\Laptop Portable\sg8\`. **The live `.gitignore` carries 4 `!soraya_runs/` negations; sg8 carries 3 — it is missing `!soraya_runs/constraint-expiry.json`.** Committing into sg8 under sg8's rules would silently drop the R9 constraint registry (10 registered scope cuts, the rail whose whole purpose is outliving the session that wrote it) — reproducing verbatim the defect the live `.gitignore` comment says was fixed on 2026-10-10. `git init` in place also fixes the anchoring bug this session hit live: `git rev-parse --show-toplevel` resolves to the home directory, where a *different* project's `wiki/plans-index.md` sits, and the v1 save nearly appended a dna_decode entry to it.

**Also settled, without changing v1's shape:** D1 ratified but split into three edits of different standing (two corrections, one genuine choice); D2 recommended-accept and left as the user's only affirmative `user:acceptance-bar` call; D4 adopt, recording explicitly that adoption is NOT evidence; D5 hold, and reword the mission so step 4 reads as *carried by a declared arm* rather than absent; D7 ratified.

Scope boundary unchanged: this plan does **not** execute the G8 measurement (that is the saved CABBAGE_G8 plan, whose Steps 2/3/5/6 are unbuilt — `scripts/cabbage_lineage_metrics.py` does not exist).

## Codebase Context

**The clonality seam is a two-field dataclass, which is what makes a second producer additive.** `dna_decode/eval/phylogeny.py:39` defines `DistanceMatrix(strain_ids: list[str], matrix: np.ndarray)` with an n x n shape check in `__post_init__`. `compute_mash_distances(...) -> DistanceMatrix` (line 94) is one producer. `dna_decode/eval/clonality.py:59` consumes only that: `greedy_representative_clusters_from_matrix(distance_matrix, threshold) -> dict[str, int]`, greedy-representative and chaining-resistant. No consumer changes.

**The G8 verdict is code-owned and needs no edit.** `rejection_gates.py:51` holds `MIN_EFFECTIVE_LINEAGES = 3`; `g8_dedup_collapses_balance` (line 246) returns `INSUFFICIENT_DATA` with no count, `TRIP` below the bar, `PASS` at or above.

**The committed lineage table is the bar's substrate, and it is fully populated.** `wiki/provdisjoint_lineage_metrics.json` (18,537 B) holds 10 cells x 2 thresholds, each carrying `effective_lineage_N_R` / `effective_lineage_N_S` plus a `cluster_weighted` block with Wilson CIs, `n_discordant` and `n_cluster_abstain`. Every cell reads `partial: false`, `n_genomes_missing: 0`, `lineage_tier_emitted: true`. **The three near-bar pairs are**: `Klebsiella x ciprofloxacin @ 0.005` (R 2 / S 18, minority **2**, the only **TRIP**), `Escherichia_coli_Shigella x ciprofloxacin @ 0.005` (R 4 / S 21, minority 4), `Escherichia_coli_Shigella x gentamicin @ 0.005` (R 5 / S 9, minority 5).

**Half the G8 plan landed already.** `dna_decode/eval/cabbage_lineage.py` (269 L) ships `SCHEMA = "cabbage-lineage-metrics-v1"`, `SIDECAR_NAME`, a `FORBIDDEN_SIDECAR_NAME` guard against writing `provdisjoint_lineage_metrics.json`, plus `cluster_structure`, `min_effective_lineages`, `weighted_point` and `assert_weighted_renderable`; `tests/test_cabbage_lineage.py` is 187 L.

**Care disk pre-check passes, measured:** `C:` is 476 G with **66 G free (87% used)** — adequate for a `sourmash` install into the conda env on a host CLAUDE.md warns runs 85-99% full.

**The `.gitignore` is already correct here and its correctness is untested.** It uses `soraya_runs/*` (trailing `/*`, because git cannot re-include a file whose parent directory is excluded) with four negations: `standing-positions.json`, `negative-results.json`, `arms-ledger.json`, `constraint-expiry.json`. Those negations have **never been verified by git** — there is no repository to run `check-ignore` against.

**Test infrastructure.** Flat offline-first pytest; real-data tests skip when a gitignored dataset is absent. Optional-dependency unavailability has a shipped precedent in `dna_decode/forward/prosst_scorer.py::StructureMethodUnavailable`.

### Reusable-Code Survey
- `dna_decode/eval/phylogeny.py` — owns the `DistanceMatrix` the new producer must return; `compute_mash_distances` is the shape to mirror.
- `dna_decode/eval/clonality.py` — the unchanged consumer; reused verbatim, which is what makes the producer additive.
- `dna_decode/eval/cabbage_lineage.py` — shipped G8 sidecar owner; Step 5 reuses `cluster_structure` rather than recomputing structure metrics.
- `dna_decode/eval/rejection_gates.py` — owns `MIN_EFFECTIVE_LINEAGES` and the G8 verdict; imported by both Step 4 and Step 5, never restated.
- `dna_decode/forward/prosst_scorer.py` — the shipped lazy-import + named-unavailability pattern copied for `SourmashUnavailable`.

(Searched: `graphify-out/GRAPH_REPORT.md` — absent; `src/`, `lib/`, `utils/`, `common/`, `helpers/` — none exist in this layout; `dna_decode/eval/`, `scripts/`, `tests/`.)

## Pre-Change Baseline

- **Mission predicates met: 1 / 10**, via `big_idea.parse_big_idea` + `os.path.exists` + exit status. Per gap: `amr-label-substrate-validated` 0/3 · `deconfounder-calibrated` 0/2 · `generator-beats-its-baseline-panel` 1/3 (the `test-exit-0` passes: 45 passed, 1 skipped) · `trait-to-edit-end-to-end` 0/2. Sigma budget 12, none released.
- **The 20-pair lineage table** above: 19 PASS / 1 TRIP, minority counts `2..24`, median 14, three pairs within +/-2 of the bar. This is the comparison set for Step 5.
- **Arms on the gap-3 ceiling: 7 measured, 3 `planned`.** Incumbent `promcalc-plus-conv-rank-average-hybrid-150bp` rho **0.4254** (10/10 seeds, CI95 [0.0238, 0.0409] excludes zero) against a 0.941 replicate ceiling. Verdict `blocked:untried-arm`.
- **`learned-model-is-the-product`:** `times_raised: 3`, `status: "open"`, `resolved_by: ""`, three recorded quotes.
- **`.gitignore` negation counts: live 4, sg8 3** (sg8 missing `!soraya_runs/constraint-expiry.json`).
- **`C:` free: 66 G of 476 G.**
- **Versioning: 0 files** under `pythonprojects/dna_decode` tracked by any git repository; `git rev-parse --show-toplevel` returns `C:/Users/Laptop Portable`.
- **sha256 of the five prospective-lock files** (`dna_decode/eval/amr_rules.py`, `calibrated_amr_rules.json`, `mic_tiers.py`, `shipped_decoder_surface.py`, `cohort_manifest.py`) and of `wiki/decoder_validation_report_card.json` and `wiki/cabbage_g8_acceptance_bar.json`.
- **Test suite:** `6409 passed / 14 skipped / 0 failed` (CLAUDE.md, 2026-10-12); re-measure through the conda interpreter before Wave 0.
- **`sourmash` not installed**; only wheel availability verified (`Would install sourmash-4.9.4`, `py3-none-win_amd64`).

## Verification Signal

- **Predicates met rises 1 -> 4 of 10**, `deconfounder-calibrated` reaching 2/2 so its budget-2 allocation releases. Gap 1's repointed first predicate must still read **open** — it now names `wiki/cabbage_lineage_metrics.json`, which does not exist yet; a flip to met would mean Step 1 pointed it at the wrong file.
- **The frozen crosscheck bar names three specific cell-threshold pairs and one TRIP cell by name**, and is computable with no drive, no `sourmash` and no network — so the bar half lands even while Step 5 is blocked.
- **`git check-ignore` confirms all four `soraya_runs/` sidecars are re-included** and that `data/`, `*.h5` and dated run directories are excluded. This assertion is currently **unrunnable**; making it runnable is the point of Step 9.
- **`compute_sourmash_distances` output is accepted unmodified by `greedy_representative_clusters_from_matrix`** — asserted on a synthetic fixture, no Mash, no network.
- **The five prospective-lock files are byte-identical** at every wave end (`verify_lock` -> `ok=True`), and `wiki/decoder_validation_report_card.json` is byte-identical.
- **Full suite still `0 failed`** on the conda interpreter.

## Implementation Steps

### Step 1: Repoint the four unsatisfiable predicates, as three edits of different standing
Files: project_state/trait-to-edit/big-idea.md, tests/test_mission_gap_predicates.py
Depends on: none

**What changes:**
- **Correction (gap 2, two predicates).** Add the `_2026-10-10` suffix the real artifacts carry. Both exist; `f3_deconfound_calibration_1001g_2026-10-10.json` reads `verdict: VALID_DIFFERENTIAL_RESPONSE`, `smoke: false`. Nothing about the bar's intent changes — the gap is already satisfied and merely unregisterable.
- **Correction (gap 1, first predicate).** `wiki/cabbage_g8_effective_lineages.json` becomes `wiki/cabbage_lineage_metrics.json`, the artifact the saved CABBAGE_G8 plan actually writes (its lines 55/63/77, schema `cabbage-lineage-metrics-v1`). Leaving it means that gap can never close on any successful measurement.
- **Genuine choice (gap 1, third predicate).** `wiki/cabbage_scored_cell_2026.json` becomes `wiki/cabbage_scored_cell_escherichia_gentamicin_2026-10-11.json` — the one cell of four whose own `source_concentration_of_scored_set` reads `source_diverse` (10 database labels, largest share 0.3534, 21 countries). **Record the flattering-selection risk in the gap text**: naming the best of four flatters, and it is only defensible because Step 2 puts the other three in the prose. A stricter alternative (name the weakest cell, or require all four) is noted inline as a one-line amendment the user can take.
- **Unchanged, deliberately:** gap 1's `cabbage_bioproject_diversity.json`. It is not a naming problem — `database`-level diversity already ships per-cell and the gap's own text correctly says `database` is an overlapping aggregator label, not a BioProject. No producer exists.
- Re-run `big_idea.parse_big_idea` + `big_idea.validate` so a malformed predicate list fails here.

**Test strategy:**
- Assert exactly 4/10 met, gap 2 at 2/2, and gap 1's repointed first predicate still open.
- `big_idea.validate(parse_big_idea(text))` raises nothing.
- Assert every `file-exists` argument either resolves on disk or is named by a saved plan under `plans/` — the check that would have caught all five misnames, including the one v1 found only by reading the G8 plan.

### Step 2: Correct four stale mission-prose claims
Files: project_state/trait-to-edit/big-idea.md, tests/test_mission_prose_claims.py
Depends on: Step 1
Test-strategy-ignore: scripts/contract_number_audit.py

**What changes:**
- **F3 (D-nothing, pure fact).** Replace `**NOT CALIBRATED** (blocked:positive_control_not_fold_compatible)` with the measured state: CALIBRATED, `VALID_DIFFERENTIAL_RESPONSE`, +0.0441 natural against 0.0000 bi-parental, carrying the artifact's own limit verbatim — **bounds the instrument only; one positive and one negative panel; validates no downstream organism-level prediction.**
- **Gap-1 diversity caveat (D1 support).** One sentence recording that `database`-level diversity already ships per-cell with 3 of 4 cells reading `source_concentrated`, and that **all four are 96-98% USA** (largest-country shares 0.9637 / 0.5844 / 0.9775 / 0.9597) — a limit no BioProject count removes. This is what makes Step 1's cell choice honest rather than flattering.
- **Modality framing (D7).** Keep the generator-bottleneck claim; correct the implied lever. On this substrate a 346-parameter biophysical model never trained on it scores rho 0.3968, above every purely-learned arm; the only tile-substrate arm whose paired CI excludes zero is the orthogonal hybrid at 0.4254, whose shuffled-conv control **hurts** (-0.1232, 0/10) so the gain is channel content. Second independent replication of rank-average-orthogonal-modalities (first: ESM2+GEMME+ProSST, 84-90% of proteins). State explicitly that it does **not** overturn the G-A encoder-gate negative — conv versus one GC number still has a CI including zero, a different question.
- **Step 4 of the north star (D5, NEW in v2).** Change "Step 4 (climb to animals) is deliberately absent" to record that step 4 is **CARRIED by the already-declared `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` arm**, naming it. "Deliberately absent" reads as closed; "carried by a declared, unrun arm" is accurate and findable, and it answers `substrate-is-not-the-subject` without spending Sigma budget on a duplicate gap.
- Leave "What must NOT be re-proposed" untouched; no edit here reopens a closed negative.

**Test strategy:**
- Assert `NOT CALIBRATED` and `positive_control_not_fold_compatible` are gone and `VALID_DIFFERENTIAL_RESPONSE` is present — a prose-relationship guard on the same pattern as the existing essentiality over-claim guard, which pins a relationship between adjacent clauses rather than a substring.
- Assert the step-4 paragraph names the arm id verbatim, so a future rename of the arm breaks the test rather than silently orphaning the cross-reference.
- Assert every number added appears in a cited `wiki/` artifact, via the existing `scripts/contract_number_audit.py` whole-number discipline.
- Re-run `big_idea.validate`; Step 1's predicate list must be unchanged.

### Step 3: A sourmash-backed DistanceMatrix producer behind an optional extra
Files: dna_decode/eval/sourmash_distance.py, tests/test_sourmash_distance.py, pyproject.toml
Depends on: none

**What changes:**
- `compute_sourmash_distances(strain_genomes: dict[str, Path], *, ksize: int = 21, scaled: int = 1000) -> DistanceMatrix`, importing `DistanceMatrix` from `dna_decode.eval.phylogeny` so the dataclass keeps exactly one definition.
- Lazy import raising a named `SourmashUnavailable`, mirroring `prosst_scorer.StructureMethodUnavailable`. No module-level import, so the test file stays collectable without the package.
- **The similarity-to-distance conversion is one named, tested helper, not an inline expression.** `ani_to_distance(ani: float) -> float` returning `1 - ani` so units match what the fixed `(0.001, 0.005)` thresholds mean. The docstring states plainly that FracMinHash at `scaled=1000` is not Mash's bottom-sketch and that threshold semantics are assumed until Step 5 measures them.
- **Verify the real sourmash surface before writing against it** (`sourmash sketch dna --help`, `sourmash compare --help`, or the Python API) and record the verified surface in the module docstring. `[unverified]` at planning time.
- `pyproject.toml`: add `lineage = ["sourmash>=4.9"]` to `[project.optional-dependencies]`, following `dev` / `quantize`. Not a default dependency. Install into an isolated check and re-run the full suite before committing the extra — `sourmash` pulls its own numeric stack and the Kaptive lesson is that an unpinned scientific package can move `numpy` under the whole suite. Disk headroom measured adequate (66 G free).

**Test strategy:**
- `ani_to_distance`: 1.0 -> 0.0, 0.999 -> 0.001, monotone decreasing, rejects out-of-range.
- Synthetic three-genome fixture in a tmp dir: returns a `DistanceMatrix`, passes the shape check, symmetric, zero diagonal, and two near-identical sequences come out closer than an unrelated one.
- **The seam test:** feed the result straight into `greedy_representative_clusters_from_matrix` and assert a `dict[str, int]` covering every input id — the assertion that proves drop-in compatibility.
- `SourmashUnavailable` (not `ImportError`, not bare `Exception`) when absent; real-sourmash cases skip via `importlib.util.find_spec`.

### Step 4: Freeze the DERIVED crosscheck bar, offline, before any number exists
Files: wiki/sourmash_vs_mash_crosscheck_bar.json, scripts/freeze_crosscheck_bar.py, tests/test_crosscheck_bar.py
Depends on: none

**What changes:**
- A generator that reads `wiki/provdisjoint_lineage_metrics.json`, computes each pair's minority effective-lineage count and its `g8_dedup_collapses_balance` verdict through the **imported** gate, and writes the frozen bar. **No drive, no `sourmash`, no network — this step is runnable today.** Splitting it out of v1's single Step 4 is what makes the bar landable while the measurement is blocked.
- The bar's three gate clauses, **derived from the gate's own arithmetic rather than asserted**:
  1. **Verdict identity** across all 20 pairs. Recorded as NECESSARY BUT NEARLY VACUOUS, with the reason in the artifact: 19 of 20 PASS with median minority 14 against a bar of 3, so a badly wrong clustering still reproduces 19/20.
  2. **Exact equality of the minority count** on every pair whose Mash count is within +/-2 of the bar. The generator computes that set rather than hardcoding it; at freeze time it resolves to exactly three pairs — `Klebsiella x ciprofloxacin @ 0.005` (2), `Escherichia_coli_Shigella x ciprofloxacin @ 0.005` (4), `Escherichia_coli_Shigella x gentamicin @ 0.005` (5).
  3. **The TRIP must reproduce.** `Klebsiella x ciprofloxacin @ 0.005` must come back TRIP. A method that converts the project's only real negative into a PASS manufactures lineage diversity — the most dangerous direction, because it would license `sourmash` on CABBAGE and report false lineage independence. **A bar that cannot fail on its one real negative is theatre.**
- **Reported, never gated:** ARI, cluster count, largest-cluster fraction, discordant fraction, and the full 20-pair table. v1 had ARI and cluster-count as gate clauses; they are demoted with the reason recorded — ARI is unanchored to anything the gate consumes, and cluster-count-within-25% binds only where verdict-identity already binds harder.
- Stamps `frozen_before_any_agreement_number_existed: true`, `derived_from` (the lineage artifact path + its sha256), and a `non_vacuity` block stating the 19/20 base rate so no reader mistakes a 19/20 headline for a result.
- **Naming the three pairs is not leaking future information** — they are computed from already-committed Mash numbers, so pinning them is pinning against known data, and the artifact says so.

**Test strategy:**
- Assert the generator derives the near-bar set rather than hardcoding it: feed a synthetic lineage artifact whose counts are all large and assert the near-bar set comes back empty; feed one with a count of 3 and assert it is included.
- Assert the verdict is computed through the imported `g8_dedup_collapses_balance`, not a local reimplementation (patch `MIN_EFFECTIVE_LINEAGES` to 10 and assert the frozen verdicts move).
- Assert the bar refuses to write if `provdisjoint_lineage_metrics.json` has any cell with `partial: true` or `n_genomes_missing > 0` — a bar derived from an incomplete table is a bar derived from nothing.
- Assert `non_vacuity` records a base rate and that the generator raises if the computed base rate is 20/20 (a table with no TRIP cannot calibrate a TRIP-reproduction clause).

### Step 5: Run the crosscheck measurement against the frozen bar
Files: scripts/sourmash_vs_mash_crosscheck.py, tests/test_sourmash_vs_mash_crosscheck.py
Depends on: Step 3, Step 4

**What changes:**
- For each cell in `wiki/provdisjoint_lineage_metrics.json`, re-cluster the same assemblies with `compute_sourmash_distances` at both fixed thresholds and evaluate the frozen bar's three clauses.
- **Compares partitions, not matrices, which resolves a circularity.** Comparing raw distances would require Mash to run — the blocked thing. The committed partition is a Mash product already on disk.
- Reuses `cabbage_lineage.cluster_structure` for the reported structure metrics rather than recomputing them.
- Emits `wiki/sourmash_vs_mash_crosscheck_2026-10-18.{md,json}`. **A cohort whose assemblies cannot all be resolved emits `INSUFFICIENT_DATA` with `n_genomes_missing` and no agreement number** — the same incompleteness rule the G8 bar commits to. With D: and E: offline this is today's expected outcome and is a valid result.
- Exit codes: 0 = all three clauses met, 1 = a clause failed, 2 = `INSUFFICIENT_DATA`, 3 = `sourmash` unavailable.
- Reads the bar from disk and **raises if `frozen_before_any_agreement_number_existed` is false** or if the bar's `derived_from` sha256 no longer matches the lineage artifact — a bar derived from a since-changed table is not a commitment.

**Test strategy:**
- Offline unit tests on synthetic partitions: identical partitions satisfy clause 1; a partition that splits the TRIP cell's minority class into 4 fails clause 3 while still satisfying clause 1 — **the test that proves clause 3 is not redundant with clause 1.**
- Assert the runner refuses an agreement number when any assembly is missing (manifest with a deliberate hole).
- Assert a lineage-artifact sha256 mismatch raises rather than scoring.
- Assert exit 2 on an empty cohort and exit 3 with `sourmash` absent, and that neither writes a bar-met verdict.

### Step 6: Amend the G8 frozen bar's distance source, carrying the measured agreement
Files: wiki/cabbage_g8_acceptance_bar.json
Depends on: Step 5

**What changes:**
- `method_commitments.distance_source` gains the sourmash route with its justification: the three clause results, the artifact path, and the explicit statement that FracMinHash is not bottom-sketch MinHash.
- **The forbidden-proxy clause is preserved verbatim.** The portal's pre-computed AMRFinderPlus tree stays forbidden — the decoder's prediction is a function of those same determinants, so a determinant-derived partition is circular. This widens the distance METHOD only.
- Add `amended: 2026-10-18`, `amended_because`, `original_distance_source`, `agreement_evidence`. Relabel the existing `frozen_before_any_lineage_number_existed` / `frozen_before_any_assembly_was_fetched` as properties of the ORIGINAL freeze rather than letting the amendment silently inherit them.
- **Runs only on Step 5 exit 0.** Exit 1 means the correct outcome is to leave the bar unamended and record that sourmash does not substitute — a finding. Exit 2 or 3 amends nothing.
- Record that CABBAGE's four minority classes are 1018, 61, 627 and 67, and that a 60-isolate Klebsiella R class collapsed to 2 at 0.005 — so the 61 and 67 classes have a real chance of TRIPPING. The bar already pre-commits a TRIP as a valid substrate finding; accepting sourmash is not assuming G8 passes.

**Test strategy:**
- Assert the amended file parses, still imports `MIN_EFFECTIVE_LINEAGES` by reference, and still contains the forbidden-proxy clause (string-level guard).
- Assert `agreement_evidence` names a `wiki/` file that exists.
- Assert the five prospective-lock files are byte-identical before and after. The G8 bar is **not** lock-pinned, so this amendment cannot retire the active lock — verified, not assumed.

### Step 7: Declare the stack arm; record that a no-evidence kill is user authority
Files: soraya_runs/arms-ledger.json
Depends on: none

**What changes:**
- `--declare` `orthogonal-modality-stack-promcalc-plus-N-channels` at `--cost-tier cheap` against the existing ceiling criterion. Rationale in the entry: the incumbent best is a 2-channel rank-average winning 10/10 with a shuffled-conv control that hurts, so the gain is channel content; whether it **stacks** to 3+ orthogonal channels is untested, and it is the exact shape of the published ESM2+GEMME+ProSST result the 2-channel version echoed.
- **The kill option v1 offered is WITHDRAWN.** `kill_arm` requires the arm's own executed evidence and refuses redundancy reasons; its refusal text routes a no-evidence kill to explicit user authority. Record in the entry that the ceiling stays `blocked:untried-arm` while any arm is `planned`, and that reaching `plateau` by declining to declare would buy a verdict through scope reduction — the 2026-07-09 defect in new clothes.
- Leave the three existing unrun arms untouched, `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` included: Step 2 now cross-references it as the carrier of north-star step 4, so deleting it would orphan that reference.

**Test strategy:**
- `arms_ledger.py --show` before and after; assert the seven measured arms' metrics and evidence strings are byte-unchanged and only the new entry differs.
- Assert the ceiling verdict after the change is `blocked:untried-arm` naming four arms.
- Assert `kill_arm` on an unrun arm with a redundancy reason still raises — a regression guard on the refusal this step depends on.
- Assert `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` is present, so Step 2's cross-reference cannot silently dangle.

### Step 8: Adopt the three-times-raised standing position, recording that adoption is not evidence
Files: soraya_runs/standing-positions.json
Depends on: Step 2

**What changes:**
- Set `learned-model-is-the-product` to `status: "adopted-user-was-right"` (the `ADOPTED` constant at `standing_positions.py:73`), with `resolved_by` pointing at the Step 2 mission text, `resolved_at`, and `tested_exactly`.
- **`tested_exactly` must state that no experiment tested the claim.** Adoption is because the user holds the position and the mission now encodes it. The tempting evidence — `glm_headroom_allrt_2026-09-30`'s supervised 0.8141 leave-one-study-out against 0.5121 zero-shot, 7 of 8 RT drugs passing versus 0 of 8, on a catalog-negative subset — supports the learned DIRECTION but is **DISCRIMINATIVE**, while the position is about **GENERATION**. Those are different quantities and conflating them is this project's recorded orientation-error 0. Citing it as the resolving evidence is exactly what `proxy_warning()` exists to catch.
- Record that adopting this does not contradict Step 2's modality correction: the winning hybrid **contains** a learned channel, so the measurement says a learned channel pays when combined with orthogonal non-learned ones, not that learning is the wrong bet.
- Direct JSON edit — the CLI exposes `--record` only, with no `--adopt` or `--resolve` verb, and the script lives outside this project at `C:\my_skills\skills\soraya\scripts\standing_positions.py`.
- Leave `substrate-is-not-the-subject` at `times_raised: 2`, `status: "open"`. Step 2's step-4 rewording addresses it; it is not resolved by this plan.

**Test strategy:**
- `standing_positions.py --show`: the 3x-open flag no longer fires; `substrate-is-not-the-subject` still listed open.
- Assert `resolved_by` names a file that exists AND that the referenced text is present in it.
- Assert the three recorded `quotes` are byte-unchanged — adoption must not rewrite the user's own words.
- Assert `tested_exactly` does not cite a wiki artifact as the resolving measurement (a string-level guard against the proxy trap this step is warning about).

### Step 9: git init in place, point at the existing origin, and verify the ignore rules
Files: .gitignore, scripts/live_vs_sg8_delta.py, tests/test_live_vs_sg8_delta.py
Depends on: none

**What changes:**
- A read-only diagnostic enumerating every path where the live copy and `C:\Users\Laptop Portable\sg8\` differ, classified live-only / sg8-only / content-differs, writing a manifest so the reconcile is done against a list. Expected live-only: 26 `wiki/*.json`, 14 `scripts/*.py`, plus `dna_decode/sklearn_compat.py`, `tests/test_sklearn_compat.py`, `tests/test_advertised_env_vars.py`.
- **`git init` in place, not a commit into sg8 — reversing v1 on measured evidence.** The live `.gitignore` carries 4 `!soraya_runs/` negations and sg8 carries 3, missing `!soraya_runs/constraint-expiry.json`; committing under sg8's rules would silently drop the R9 constraint registry. In-place init also makes `git rev-parse --show-toplevel` resolve to the project rather than the home directory, which is the bug that nearly wrote a dna_decode plan entry into another project's `wiki/plans-index.md` this session.
- Sequence, each step non-destructive: `git init` -> `git remote add origin https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION.git` -> `git fetch origin main` -> **`git reset origin/main`** (MIXED: sets the index to origin/main and leaves every working file untouched, so `git status` then renders the full drift as unstaged changes).
- **`git reset --hard` is FORBIDDEN in this step and must not appear in any command.** It would destroy the 682 live artifacts and 507 scripts, which exist in no other copy at their current state. Same for `git checkout .` and `git clean -fd`.
- **First post-init act is verifying the ignore rules, which have never been checked by git**: `git check-ignore -v` must report all four `soraya_runs/` sidecars as RE-INCLUDED and `data/`, a `*.h5`, and a dated `soraya_runs/<run-id>/` path as excluded. The trailing-`/*` trick the file's own comment describes is asserted by comment and untested until now.
- Nothing is committed in this step. Staging and the first commit are a separate decision once the delta manifest has been read.

**Test strategy:**
- Unit-test the delta classifier on a synthetic pair of trees (live-only / sg8-only / content-differs / identical), with no reference to the real paths.
- Assert the diagnostic writes only its own artifact: run it twice and confirm `git status` in the home repo and in `sg8` is unchanged.
- After init, assert `git check-ignore -q soraya_runs/constraint-expiry.json` exits NON-zero (i.e. not ignored) and that `git check-ignore -q data/processed/x.parquet` exits zero — the non-vacuity pair for the negation block.
- Assert `git rev-parse --show-toplevel` now returns the project directory, pinning the anchoring fix.
- Assert the forbidden commands appear nowhere in the step's scripts (`grep` guard for `reset --hard`, `clean -fd`, `checkout .`).

## Execution Preview

```
Wave 0 (5 parallel):  Step 1 — Repoint the four unsatisfiable predicates
                      Step 3 — A sourmash-backed DistanceMatrix producer
                      Step 4 — Freeze the DERIVED crosscheck bar (offline, runnable today)
                      Step 7 — Declare the stack arm
                      Step 9 — git init in place + verify the ignore rules
Wave 1 (2 parallel):  Step 2 — Correct four stale mission-prose claims
                      Step 5 — Run the crosscheck measurement
Wave 2 (2 parallel):  Step 6 — Amend the G8 frozen bar's distance source
                      Step 8 — Adopt the standing position
```

Critical path: Step 3 -> Step 5 -> Step 6 (3 waves). Co-critical: Step 1 -> Step 2 -> Step 8.
Max parallelism: 5 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode. **Today that fallback applies** — the working copy is not a repository. Step 9 is what changes it, and Step 9 is itself in Wave 0, so a first run of this plan is sequential by construction.

## Risk Flags

- Severity: critical — **`git reset --hard`, `git clean -fd` and `git checkout .` in Step 9 would destroy 682 artifacts and 507 scripts that exist in no other copy at their current state.** Only `git reset origin/main` (mixed) is permitted; the step carries a grep guard and the test suite asserts the forbidden forms appear nowhere.
- Severity: high — **Clause 1 of the frozen bar is nearly vacuous and the artifact must say so.** 19 of 20 pairs PASS with median minority 14 against a bar of 3. A reader seeing "19/20 verdicts matched" would read a broken instrument as validated. This is the project's own recorded trap (the PEAR WT-row at mode-share 0.2678, the vacuous HIV censoring filter, the vacuous ResFinder `POINT` filter): a screen whose pass rate is ~100% by construction is not a control. Mitigated by the `non_vacuity` block, clause 3, and the test that proves clause 3 is not redundant with clause 1.
- Severity: high — **FracMinHash is not bottom-sketch MinHash and the threshold semantics are assumed until Step 5 measures them.** The `(0.001, 0.005)` thresholds were inherited from Mash output. Step 6 is pre-committed not to run on a non-zero Step 5 exit.
- Severity: high — **External tool integrations.** `sourmash`: wheel availability **verified** (`pip install --dry-run --no-deps --only-binary=:all: sourmash` printed `Downloading sourmash-4.9.4-py3-none-win_amd64.whl.metadata` and `Would install sourmash-4.9.4` on py3.13.16/AMD64); its **CLI and Python API surface is `[unverified]`** — no flag, subcommand, output format or ANI accessor confirmed, because the package is not installed. Step 3 verifies the surface before writing against it. `mash`: confirmed unavailable natively — v2.3 ships only Linux64 and OSX64 tarballs, bioconda publishes no `win-64`, and both remaining routes need the same virtualization layer `wsl --status` reports as BIOS-gated. `git`: surfaces used are `init` / `remote add` / `fetch` / `reset` (mixed) / `check-ignore -v` / `rev-parse --show-toplevel`, all standard; `gh` present but unused.
- Severity: high — **Steps 1 and 2 both write `project_state/trait-to-edit/big-idea.md`.** Resolved by making Step 2 depend on Step 1, placing them in different waves. No intra-wave file overlap remains.
- Severity: medium — **Step 5 cannot produce a number while D: and E: are offline**, and `data/` currently resolves to an empty directory. Expected outcome is exit 2 `INSUFFICIENT_DATA`, which blocks Step 6 by design. The drive has proved intermittent (present ~28 minutes then dropped with bus-level USB faults), so Step 5 should not be launched as a long run until it is stable. **Step 4 was split out of v1's Step 4 precisely so the bar lands regardless.**
- Severity: medium — **Amending a frozen bar after discovering a block is the shape of the mis-specified-bar failures this project recorded twice.** Mitigated by: the amendment changes the METHOD, its justification comes from a separately-frozen bar, the circularity rail is preserved verbatim, and `original_distance_source` keeps the prior commitment legible. It remains a bar change and belongs to the user.
- Severity: medium — **New optional dependency.** `lineage = ["sourmash>=4.9"]` is additive and not default-installed, but `sourmash` pulls its own numeric stack and the Kaptive lesson is that an unpinned scientific package can move `numpy` under the whole suite. Install isolated, re-run the suite, then commit the extra. Care disk pre-check **passes as measured**: 66 G free on `C:`.
- Severity: medium — **Transitive imports not visible in the `Files:` lists.** Step 3 imports `DistanceMatrix` from `phylogeny.py`; Step 4 imports `g8_dedup_collapses_balance` + `MIN_EFFECTIVE_LINEAGES` from `rejection_gates.py`; Step 5 additionally imports `cluster_structure` from `cabbage_lineage.py`. None of those files is modified by any step, so the edges are read-only — but a change to `DistanceMatrix`'s shape contract breaks Step 3 silently.
- Severity: low — **`standing_positions.py` has no adoption verb.** Step 8 is a direct JSON edit; adding the verb would be a skill-code change outside this plan's scope.
- Severity: low — **Restructuring applied, v1 -> v2.** v1's Step 4 split into Step 4 (freeze the bar, offline) + Step 5 (measure, drive-gated), so the bar lands while the measurement is blocked; its ARI and cluster-count clauses demoted from gates to reported diagnostics with the reason recorded; the arm-kill option withdrawn (code refuses it); the versioning recommendation reversed with the `.gitignore` negation count as evidence; a fourth prose edit added to Step 2 for the north-star step-4 rewording. Nine steps, three waves, max parallelism 5 (v1: eight steps, three waves, max parallelism 4).

## Open Questions

- **D2 is the one decision still genuinely outstanding** — accept `sourmash` as an admissible G8 distance source, conditional on Step 5 clearing the Step 4 bar? It is the only item the ask gate classified affirmatively (`user:acceptance-bar`: "redefines what counts as done"). Recommended: accept. Declining leaves gap 1 blocked on a BIOS setting indefinitely.
- **Does the derived three-clause bar replace the asserted one?** Step 4 freezes it as written unless amended now; amending after seeing the number is the failure this project has recorded twice. The derivation is in Problem Statement and Codebase Context; clause 3 (the TRIP must reproduce) is the one carrying the discriminating power.
- **Is the gap-1 cell choice right?** Step 1 names the one `source_diverse` cell of four. A stricter alternative is to name the weakest cell or require all four. Flagged because choosing the best of four flatters, and only Step 2's prose edit keeps it honest.
- **A no-evidence arm kill, if you want the ceiling to reach `plateau` this horizon.** `kill_arm` routes it to explicit user authority; it is not available to the executor. Not recommended — `blocked:untried-arm` is the honest state.
- **Staging and the first commit after Step 9's `git init`.** Step 9 deliberately stops at a mixed reset and a read delta manifest. What gets committed first, and on which branch, is a separate call.

## Architecture Impact

**Expected structural changes:**
- Modules added: `dna_decode/eval/sourmash_distance.py` (new leaf). Modified: none. Deleted: none. New scripts: `scripts/freeze_crosscheck_bar.py`, `scripts/sourmash_vs_mash_crosscheck.py`, `scripts/live_vs_sg8_delta.py`.
- Dependency direction: `sourmash_distance` -> `eval/phylogeny` (the dataclass only) + the optional third-party `sourmash`. `freeze_crosscheck_bar` -> `eval/rejection_gates`. `sourmash_vs_mash_crosscheck` -> `sourmash_distance`, `eval/clonality`, `eval/cabbage_lineage`, `eval/rejection_gates`. **Every edge points from new code into existing code; nothing existing gains an import of anything new.**
- Boundary changes: one new external-tool boundary (`sourmash`) behind a lazy import and a named `SourmashUnavailable`; one new repository boundary (Step 9's in-place `git init` with the existing origin). No layer relocation.

**Risk areas:**
- `DistanceMatrix` becomes a two-producer abstraction. Intended, but its shape contract is now load-bearing for two callers.
- `rejection_gates.g8_dedup_collapses_balance` gains a second consumer (Step 4's bar generator alongside the CABBAGE screen). Read-only, and Step 4's test patches `MIN_EFFECTIVE_LINEAGES` to prove the import is live rather than reimplemented.
- No shared utility is modified: `clonality.py`, `rejection_gates.py`, `cabbage_lineage.py`, `phylogeny.py` are all read-only here.
- Cycles: none possible. `sourmash_distance` imports `phylogeny`, which gains no import of it, so no A-to-B-to-A can form.

**Sentrux gates:**
- Baseline before Wave 0: no — `n/a, sentrux not installed` (tool check found only `gh`).
- Compare after each merged wave: no — `n/a, sentrux not installed`.
- Acceptable delta: `n/a, sentrux not installed`. The substituted structural checks are the no-reverse-import and no-cycle assertions below.
- Regression-response stance: `n/a` for sentrux. For the substituted checks, a failure is investigate-and-fix before the wave counts as merged, never accept-and-document.

**Anti-drift checks:**
- `grep -n "sourmash" dna_decode/eval/{phylogeny,clonality,rejection_gates,cabbage_lineage}.py` returns nothing.
- `grep -rn "class DistanceMatrix" dna_decode/` returns exactly one line, `phylogeny.py:39`.
- `import sourmash` appears at no module top level, so `python -c "import dna_decode.eval.sourmash_distance"` succeeds with the package absent.
- `grep -rn "MIN_EFFECTIVE_LINEAGES *=" dna_decode/ scripts/` returns exactly one definition, `rejection_gates.py:51` — no step restates the bar.
- `grep -rnE "reset --hard|clean -fd|checkout \." scripts/live_vs_sg8_delta.py` returns nothing.
- The five prospective-lock files match their pre-plan sha256 and `verify_lock` returns `ok=True`.
- `wiki/decoder_validation_report_card.json` is byte-identical — this plan scores no cell.
- Full suite on the conda interpreter still reports `0 failed`; skip count rises by at most the number of sourmash-gated tests.

## Verification

1. Re-run the predicate evaluator (`big_idea.parse_big_idea` + `os.path.exists` + exit status). Expect **4 of 10 met**, `deconfounder-calibrated` 2/2, gap 1's repointed first predicate still **open**.
2. `uv run python scripts/freeze_crosscheck_bar.py` with no drive and no `sourmash` — writes `wiki/sourmash_vs_mash_crosscheck_bar.json` naming 3 near-bar pairs, 1 TRIP cell, and a `non_vacuity` base rate of 19/20.
3. `python -c "import dna_decode.eval.sourmash_distance"` with `sourmash` absent succeeds; `compute_sourmash_distances` raises `SourmashUnavailable`.
4. `pytest tests/test_sourmash_distance.py tests/test_crosscheck_bar.py tests/test_sourmash_vs_mash_crosscheck.py tests/test_mission_gap_predicates.py tests/test_mission_prose_claims.py tests/test_live_vs_sg8_delta.py -q` passes, real-sourmash cases skipped when absent.
5. `pytest tests/ -q` through `C:\ProgramData\miniforge3\envs\dna_decode\python.exe` reports `0 failed`.
6. `sha256sum` the five prospective-lock files against baseline — identical; `verify_lock` returns `ok=True`. Same for `wiki/decoder_validation_report_card.json`.
7. `arms_ledger.py --show` — seven measured arms byte-unchanged; ceiling verdict `blocked:untried-arm` naming four arms.
8. `standing_positions.py --show` — the 3x-open flag no longer fires; `substrate-is-not-the-subject` still open; quotes byte-unchanged.
9. `scripts/sourmash_vs_mash_crosscheck.py` exits 0 with all three clauses met, or exits 2 `INSUFFICIENT_DATA` naming `n_genomes_missing`. **Exit 2 is expected while D: and E: are offline and is not a failure.**
10. After Step 9: `git rev-parse --show-toplevel` returns the project directory; `git check-ignore -q soraya_runs` style checks confirm all four `soraya_runs/` sidecars re-included and `data/` + `*.h5` + dated run dirs excluded; `git status` in the home repo and in `sg8` unchanged by the diagnostic.

## Save-time amendments

Captured at: 2026-10-18
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads ONLY
`## Implementation Steps` for executable work — amendments here are NOT executable instructions. If an
amendment changes a Step contract (file lists, dependencies, structure), re-run `/technical-plan` before
`/execute-plan`.

**This save REPLACED v1 in place** rather than creating a `_2` folder, per the save-once rule in
`rules/planning-pipeline.md`: a v1 -> v2 rewrite prefers an in-place edit plus an index `Status:` update,
because a second folder is exactly the superseded artifact that rule exists to keep out of
`plans-index.md`. v1 is preserved beside this file as `technical-plan.md.v1bak`.

- v2: crosscheck bar DERIVED from the 20-pair lineage table (19 PASS/1 TRIP) replacing the asserted ARI+count clauses
- arm-kill option withdrawn (kill_arm refuses without the arm's own evidence)
- D6 reversed to git init in place because sg8's gitignore is missing the constraint-expiry negation
- Step 4 split so the bar lands offline

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high -->
