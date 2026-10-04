# Transfer Verdict Shape On The Existing De-Confound Package (F1, Phase 0, v3)

## Lens status

- **Dispatch:** conversation-driven — v3 lands the four Open-Question resolutions into the v2 body, which is already on disk at `plans/Transfer_Benchmark_And_Regime_Transfer_Axis_F1_Phase_0_Plan/technical-plan.md`. Per the save-once rule (a v1→v2→v3 rewrite after review is an **in-place edit**, not a new plan folder) this body replaces that file again; `## Repo grounding` and `## Save-time amendments` are preserved on save.
- **New since v2, found during this plan's own Step 1 and it changes a step:** `within_group_r2`'s `min_n=30` is a **tighter floor than `MIN_SCORED=10`**, and it governs the cell the PASS condition reads. Measured: at 20 members per group `within_group_r2` returns `(nan, 0 groups used)` while `cv_r2` on the identical data returns a confident-looking `−0.4230`. A k-grid derived from `MIN_SCORED` alone would therefore produce runs whose pooled cell reports a number and whose verdict cell is silently undecidable.
- **Repo-index:** absent (`update-index` CLI not installed) — direct file reading, no degradation.
- **Tools:** `gh` present · `codex` present · **`sentrux` absent** → no architecture gate.
- **Optional context:** `wiki/` present → `plans-index.md` + HIGH-salience `decisions-log.md` replayed. `AUDIT_REPORT.md` / `DESIGN.md` / `project-rules.md` absent → Project Rules Check omitted.
- **Scope gate:** $0, laptop-only, offline. No money, no new labels, no GPU, no network. Touches **no frozen file**.

## Problem Statement

Every learned number in this repo was measured with the held-out unit being a **position**, a **study**, or a **protein** — never an **organism**, and `dna_decode/eval/regime.py` cannot express that. The one-virus and one-species limits live in free-text `note` fields where they go stale, which is the exact failure the module was built to stop.

The recorded failure mode is specific. On Arabidopsis FT10 the embedding scored spearman 0.213/0.221/0.226 against a structure-only baseline of **0.484**, and within-group r² of **−0.173/−0.110/−0.128** against structure-only **+0.038** — it transferred phylogenetic nearest-neighbour lookup. **It did win one metric** (global r² −0.046/−0.036/−0.031 vs −0.449), the metric most polluted by population structure. Any evaluation that read the winning metric would have scored it a success.

**v1 scope was measured unachievable, and this plan narrows it deliberately rather than silently.** v1 proposed a falsifier with a real transfer-POSITIVE arm. Verification found that **every genotype-phenotype panel on disk is a single constructed population with no internal organism axis** — Bloom BYxRM, BXD mouse, Arabidopsis 19-parent MAGIC — so a held-out-organism fold has nothing to hold out and would silently degrade to K-fold. The rescue (hold out a whole panel/species) is also invalid: the three loaders emit panel-specific marker matrices with no shared feature ontology or orthology map, so a zero would prove *feature-space incompatibility*, not a transfer floor, and labelling it `held_out_organism` would encode a false biological reading into the benchmark.

**So the deliverable is the verdict SHAPE plus an explicit, first-class `TRANSFER_POSITIVE_UNMEASURED` stamp** — not a green falsifier. That stamp is the point: it is what stops a later session from assuming transfer works. F1 still owes everyone the verdict shape; no gLM number is interpretable before it exists.

**This is settled by measurement, not pending ratification.** v2 filed the no-positive-arm choice as an open question needing sign-off; that was wrong. The alternatives are not *ship with a positive arm* vs *ship without one* — they are *ship without one* vs *fabricate one* (invent a grouping inside a single cross, or relabel feature-space incompatibility as biology). There is no decision left in that. **One genuine fork does survive and it is the only question this plan asks:** whether an instrument that can currently measure only an absence is worth building at all. That is a scope call, recorded in Open Questions.

**Non-goals:** no backbone, no pretraining, no corpus, no embedding inference, no new de-confounding mathematics anywhere.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

- **Test suite:** `uv run pytest tests/ -q --collect-only` → **5428 collected** (measured live; unchanged since v1 — only markdown has been edited). The last full run reported 5424 passed / 12 skipped / 0 failed, which does not reconcile with 5428; `/execute-plan`'s baseline pass establishes the authoritative figure and the regression check uses **that**.
- **`eval/regime.py`:** 8 `Regime` rows, 3 axes, **zero rows able to express which split produced their number**. 13 tests in `tests/test_regime_boundary.py`. Imports in 0.005 s.
- **`eval/transfer.py` / `eval/transfer_gauntlet.py`:** do not exist.
- **G2 replay inputs:** not machine-readable (markdown table only).
- **Arabidopsis committed numbers** (seeds 42/7/99) — global r² −0.046/−0.036/−0.031 vs structure-only −0.449; spearman 0.213/0.221/0.226 vs 0.484; within-group r² −0.173/−0.110/−0.128 vs +0.038.
- **`dna_decode.deconfound` consumers today:** 8 scripts + its own `__init__` + 1 test file. Adding `eval/` consumers must not change its surface.
- **Frozen surface:** 5 files under `wiki/prospective_lock_manifest_2026-08-31.json`; `verify_lock` currently `ok=True`.
- **scipy:** present (1.17.1) but **not declared in `pyproject.toml`** — it arrives transitively via `scikit-learn>=1.4`. Pre-existing; this plan increases reliance on it.

## Verification Signal

1. **The falsifier discriminates on the two arms it can, and refuses to claim the third.** A1 POSITIVE on Bloom; B NEGATIVE on the synthetic structure-only input; the run emits `transfer_positive: "TRANSFER_POSITIVE_UNMEASURED"` as a **first-class verdict value**, never an omission or a pass.
2. **Arm C cannot contribute to PASS.** A test asserts PASS is reachable with Arm C absent, and that Arm C alone never produces PASS.
3. **The replay regression pins the hard case.** Arm C must return NEGATIVE on the Arabidopsis numbers *even though the candidate WINS global r²* — proving the verdict is driven by the de-confounded metric, not by whichever metric the candidate happens to win.
4. **Skipped groups are reported, never dropped silently.** For any fold set where `cv_r2` would skip a group, the harness names it and the count appears in the artifact.
5. **Four-cell metric matrix, not one score.** `(zero_shot_group_r2, kshot_pooled_r2, kshot_within_group_r2, kshot_group_offset_baseline)`; transfer-positive requires beating the offset baseline on the **within-group** cell; a pooled-only win reports `GROUP_OFFSET_LEARNED`.
   - **The verdict refuses when its own cell is undecidable, and never substitutes the pooled one.** `WITHIN_GROUP_UNSCORABLE` on `isnan(within_group)` or `n_groups_used == 0`, asserted on the measured 6×20 case where the pooled cell still reads `−0.4230`.
   - **The k grid is derived and refuses rather than degrades.** `cap = min_group_size − max(MIN_SCORED, WITHIN_GROUP_MIN_N)`; `cap < 1` raises; `WITHIN_GROUP_MIN_N` is passed explicitly to `within_group_r2` at every call site rather than inherited from its default.
6. **No new de-confounding math.** An AST scan asserts `eval/transfer.py` and `eval/transfer_gauntlet.py` import **no** `sklearn.*` or `scipy.*` symbol; every statistic arrives via `dna_decode.deconfound`.
7. **Provenance recorded, not promised.** The gauntlet report carries the qualname of every imported primitive plus the sha256 of `deconfound/deconfound.py`, so drift between the replay and the live arms is detectable from the artifact.
8. **The regime axis is augment-only, verified by diff** — all 8 regimes × both targets, with and without the new argument; only `conditions` may differ.
9. **Import posture preserved:** `python -c "import dna_decode.eval.transfer"` stays under ~0.05 s (no eager scipy/sklearn).
10. **Nothing frozen moved;** `verify_lock` re-verifies `ok=True`; `wiki/decoder_validation_report_card.json` unchanged (27 cells, identical key set). **No regressions** vs the baseline pass.

## Implementation Steps

### Step 1: Split the regime axis into an unordered fact and a narrow ordered field
Files: dna_decode/eval/regime.py, scripts/regime_map.py, tests/test_regime_boundary.py
Depends on: none

**What changes:**
- `SPLIT_UNITS = frozenset({"position", "study", "protein", "organism", "clade", "condition", "none"})` — an **unordered** vocabulary. No comparability implied between members.
- `ORGANISM_TRANSFER_LEVELS = ("unmeasured", "held_out_organism", "held_out_clade")` — the **only** ordered field, answering one question: does this row carry organism-transfer evidence.
- `Regime` gains `split_unit: frozenset[str]` (**no default** — an empty frozenset means "no split was performed", a different claim from "unknown") and `organism_transfer: str = "unmeasured"` (**defaulted**, because `unmeasured` is simultaneously the fail-closed value and the true value for all 8 rows today).
- Values read off the committed artifacts, never invented: `constructed_molecular` → `{"protein"}`; `constructed_molecular_zeroshot` → `{"position"}`; `natural_molecular_supervised_blindspot` → `{"study"}`; `constructed_organism_condition_switch` → `{"condition"}`; the rest → `{"none"}`. **Every row keeps `organism_transfer="unmeasured"`.**
- `organism_transfer_is_unmeasured() -> tuple[str, ...]` (all 8 today) and `split_units_for(key) -> frozenset[str]`.
- `screen_proposal(..., claims_organism_transfer: str | None = None)` is **augment-only**: when a claim exceeds the matched regime's `organism_transfer`, append a CONDITION. It never touches `verdict`, `reason`, `evidence` or `artifact`. A value outside `ORGANISM_TRANSFER_LEVELS` → `ScreenResult(None, "UNKNOWN", ...)`, mirroring the existing bad-axis handling. `classify_regime`'s 3-tuple key is unchanged.
- `scripts/regime_map.py` renders **both** new fields as printed columns plus the headline `N of 8 regimes carry no held-out-ORGANISM transfer evidence`. `as_dict()` carries them into `wiki/learned_regime_map.json`; a field reaching only the JSON is not a disclosure.
- **Exit code stays 0 when all 8 rows read `unmeasured`, deliberately.** A non-zero exit means *the tool failed*; all-unmeasured is a correct description of the world. The existing non-zero exit is reserved for a genuine integrity break — a cited artifact that does not resolve. Conflating "the evidence is absent" with "the map is broken" would make a truthful map look broken on every run, and this repo has the cost on record: *"a permanently-red suite trains readers to discount red, which is how a genuine new failure gets missed."* The headline line in stdout carries the loudness instead.

**Test strategy:**
- `split_unit` is required (constructing a `Regime` without it raises `TypeError`); `organism_transfer` defaults to `unmeasured`; every value is in its vocabulary.
- **No ordering exists over `SPLIT_UNITS`** — a test asserts the module exposes no comparison helper over split units, so the v1 ladder cannot reappear.
- `organism_transfer_is_unmeasured()` returns all 8 today, with the assertion written as a tripwire on the count changing rather than a frozen 8.
- Augment-only diff over 8 regimes × `("replace", "blind_spot_complement")` with and without the new argument; non-vacuous (at least one case gains a condition).
- Printed stdout contains both columns and the headline count equals `organism_transfer_is_unmeasured()`.

### Step 2: Make the G2 replay inputs machine-readable, and verify rather than trust them
Files: wiki/phase2_arabidopsis_result_2026-06-12.json, tests/test_g2_replay_fixture.py
Depends on: none

**What changes:**
- New committed sidecar holding the three metrics × three seeds × two arms transcribed from the markdown table, plus `source: "wiki/phase2_arabidopsis_result_2026-06-12.md"` and the verbatim table rows it was transcribed from.
- The alternative — hard-coding twelve numbers inside the falsifier — would make the replay an assertion that can silently drift from its source. This follows the repo's established **verified-not-trusted** binding pattern.

**Test strategy:**
- A test **re-parses the markdown table** and asserts every JSON value matches it; a transcription error fails loudly.
- Asserts the structural fact the replay turns on: the embedding arm **beats** structure-only on `global_r2` while **losing** on both `spearman` and `within_group_r2`. If that shape is ever absent, Arm C is not testing what it claims.
- Asserts `within_group_r2` is negative for all three seeds (the signature of the recorded failure).

### Step 3: Transfer harness — group folds that report skips, k-shot bookkeeping, power check
Files: dna_decode/eval/transfer.py, tests/test_eval_transfer.py
Depends on: none

**What changes:**
- New module; numpy at module scope, **`dna_decode.deconfound` imported lazily inside each function** (preserving the measured 0.002–0.005 s `eval/` import posture).
- `held_out_group_folds(group_of, *, allow_unassigned=False) -> list[TransferFold]`; raises `UnassignedMemberError` on an id with no group. Fold membership is a pure function of the group key, so disjointness is **by construction**.
- `fold_report(folds)` surfaces **structure, not just count**: n_folds, largest-group fraction, singleton fraction.
- `skipped_groups(folds, min_train=5, min_test=1) -> tuple[str, ...]` — **precomputes exactly which groups `cv_r2` will silently `continue` past** (`deconfound.py:30-31`), so the harness can name them instead of quietly scoring a subset of organisms. Reports, never repairs.
- `k_shot_curve(X, y, group_of, *, ks=None, fit_predict, seed=0) -> KShotCurve` — k shots drawn from the held-out group, excluded from the scored set by construction; the same draw reused across every method compared (paired).
- **The k grid is DERIVED, not asserted.** `ks=None` computes `cap = min_group_size - max(MIN_SCORED, WITHIN_GROUP_MIN_N)` then `[0] + [2**i for i in … if 2**i <= cap]`. `MIN_SCORED = 10`; `WITHIN_GROUP_MIN_N = 30` **mirrors `within_group_r2`'s own `min_n` and is passed explicitly at the call site**, never left to that function's default — a default in another module silently governing this plan's PASS condition is the drift class this plan exists to avoid. `cap < 1` **raises**: a group too small to give one shot and still satisfy the within-group floor cannot carry a k-shot curve, and returning an empty curve would be the fail-open shape. An explicit `ks=` is still accepted for reproducing a published grid.
- This also matters substantively: with group sizes of 198–1008 a grid stopping at 16 would very likely never reach saturation, so the v2 fixed tuple risked reporting the uninteresting left end of the curve as if it were the curve.
- `power_check(curve_a, curve_b) -> str` returns `UNDERPOWERED_METHODS_NEVER_DIFFERED` when no k cell yields a different prediction set; callers must consult it **before** any effect size, since a zero-difference run necessarily has a zero gain.

**Test strategy:**
- k shots never intersect the scored set at any k; unassigned id raises by default and buckets only under `allow_unassigned=True`.
- `skipped_groups` flags a 3-member group against `min_train=5` and the flagged set matches what a direct `cv_r2` call actually drops (asserted against observed NaNs, so the prediction is checked against behaviour rather than against the docstring).
- Two identical methods → `UNDERPOWERED_METHODS_NEVER_DIFFERED`, returned before any gain is computed.
- Paired determinism: same seed → identical shot indices across methods.
- **Derived grid:** groups of 80 → `cap = 50` → grid tops out at 32; groups of 35 → `cap = 5` → grid `[0,1,2,4]`; groups of 30 → `cap = 0` → **raises**. The 30-member boundary is asserted against the measured `within_group_r2` behaviour (`nan`, `n_groups_used=0`), not against the literal constant.
- `import dna_decode.eval.transfer` does not pull `sklearn` or `scipy` into `sys.modules`.

### Step 4: Fail-closed leakage audit
Files: dna_decode/eval/transfer.py, tests/test_eval_transfer.py
Depends on: Step 3

**What changes:**
- `leakage_audit(train_ids, test_ids, *, resolver=None) -> LeakageAudit` with `overlap`, `unresolved`, `incomplete`, `reason`. Exact-id overlap fails immediately; an optional `resolver` may collapse aliases and **any** resolver failure sets `incomplete=True`.
- `incomplete=True` is treated as a possible false-independence claim, not a harmless offline fallback — the posture `provenance_disjoint_validate.py` already takes. Consumers refuse to emit a transfer number unless `--allow-incomplete-leakage-audit`, which stamps `leakage_audit_degraded=True` into the artifact.
- Imports nothing from the FROZEN `eval/cohort_manifest.py` at module scope and never writes to it; a caller may pass one in.

**Test strategy:**
- Overlap → fail; disjoint → pass; a raising resolver → `incomplete=True` with a reason.
- A consumer handed an incomplete audit refuses by default, succeeds with the override, and the override is visible in the record.
- A test asserts `cohort_manifest` is absent from the module's imports.

### Step 5: Gauntlet of five continuous controls plus the four-cell verdict matrix
Files: dna_decode/eval/transfer_gauntlet.py, tests/test_transfer_gauntlet.py
Depends on: Step 3

**What changes:**
- Five controls, every statistic delegated to `dna_decode.deconfound` (lazily imported): `shuffled_label_null` → `permutation_null` (within-group, the correct null for this estimand); `structure_only_baseline` → `cv_r2` on one-hot group membership; `nearest_neighbour_sequence` → 1-NN by `numpy.argmin` on a supplied distance matrix, scored with `r2`; `pca_only_baseline` → **numpy SVD** to reduce `X`, then `cv_r2` (numpy, so the module needs no sklearn of its own); `held_out_clade` → `cluster_from_distance` groups through `cv_r2`.
- `run_gauntlet(...) -> GauntletReport` carrying each control's score, `n_controls_run`, the **four-cell metric matrix**, `skipped_groups`, and `primitive_provenance` (qualname per imported primitive + sha256 of `deconfound/deconfound.py`).
- `assert_all_controls_ran()` raises below 5 — and, because 4 of 5 controls are imported, the **load-bearing** guard is `assert_same_evaluation_basis()`: identical group vector, identical scored-id set, identical residualization convention and identical k-shot exclusions across every control. Five labels existing proves far less than five controls having been scored on the same items.
- `verdict(matrix, report) -> str` returns `BEATS_ALL_CONTROLS` / `LOSES_TO_<control>` / `GROUP_OFFSET_LEARNED` / `WITHIN_NOISE` / **`WITHIN_GROUP_UNSCORABLE`**. **Transfer-positive requires beating the group-offset baseline on the within-group cell**; a pooled-only win is `GROUP_OFFSET_LEARNED`.
- **`WITHIN_GROUP_UNSCORABLE` is a required verdict, not an error path.** `within_group_r2` returns `(nan, 0)` whenever every group falls under `min_n`, and the pooled cell remains a confident-looking number on the identical data. The verdict function must therefore refuse on `isnan(within_group)` or `n_groups_used == 0` and **never fall back to the pooled cell** — that fallback is how a `GROUP_OFFSET_LEARNED` run would be reported as a pass.
- `min_n` is passed **explicitly** (`WITHIN_GROUP_MIN_N`) at every `within_group_r2` call site, and `n_groups_used` is carried into `GauntletReport` so the reader can see how many groups actually contributed.
- A guard refuses routing a continuous score into `clade_baseline.validation_gate` (AUROC-shaped).

**Test strategy:**
- Signal-free input loses to the shuffled-label null; pure-structure input loses to `structure_only_baseline`; within-group signal beats both. Fixtures follow `test_deconfound_package.py:23`'s construction, with the negative being that fixture minus its signal term.
- `assert_all_controls_ran()` raises on a 4-control report; `assert_same_evaluation_basis()` raises when one control is scored on a different id set — **the non-vacuity pin that matters**, and tested by mutating the basis rather than the count.
- **AST scan:** the module imports no `sklearn.*` / `scipy.*` symbol. Non-vacuous: PCA goes through numpy SVD precisely so this guard can hold.
- `GROUP_OFFSET_LEARNED` fires on an input where pooled wins and within-group does not.
- **`WITHIN_GROUP_UNSCORABLE` fires on 6 groups × 20** (the measured case: within-group `nan`/0 groups, pooled `−0.4230`), and a test asserts the verdict is NOT `LOSES_TO_*` or any pooled-derived value there — the fallback is closed by test, not by comment.
- `n_groups_used` appears in the report for both the scorable and the unscorable case.

### Step 6: Falsifier with two live arms, one regression arm, and an explicit unmeasured verdict
Files: scripts/transfer_benchmark_falsifier.py, tests/test_transfer_benchmark_falsifier.py
Depends on: Step 1, Step 2, Step 4, Step 5

**What changes:**
- `PREREGISTERED` bar as a module constant **above any data loading**: PASS requires **A1 positive AND B negative**, and emits `transfer_positive="TRANSFER_POSITIVE_UNMEASURED"` whenever no live held-out-organism substrate was supplied.
- **A1 — scoring-core control** (not transfer evidence): Bloom BYxRM through the gauntlet, stamped `is_transfer_evidence=False`.
- **B — synthetic negative, full machinery**: g groups with a per-group phenotype offset and **zero** within-group signal, run through folds → leakage audit → gauntlet. Deterministic, offline. This is the only arm that exercises the live machinery on a known negative.
- **C — verdict-replay regression** (not an arm that can pass): the Step 2 sidecar through `verdict`, stamped `mode="verdict_replay"`, `machinery_exercised=False`, `counts_toward_pass=False`. Must return NEGATIVE **despite the candidate winning global r²**.
- **DROPPED in v3 — the cross-panel control is not built.** v2 carried an optional `--incompatible-feature-space-control` running yeast/mouse/MAGIC together, barred from PASS. It is removed: its result is **knowable without running it**. The three panels use three mutually unmappable identifier *schemes* (`27915_chr01_27915_T_C` coordinate+allele / `rs31443144` dbSNP / `MN1_29291`) with no orthology or liftover resource anywhere in this plan, so "the model scores ~0" restates the input rather than measuring anything. A number with no information content that can be quoted out of context is a liability, and this repo has the precedent on record — a filter shipped as an active control while removing nothing, later retracted. If a shared cross-organism representation ever arrives, that is a real arm on a real substrate; a vacuous stand-in in the meantime only gives a future reader something to misquote.
- Refusals: missing substrate, absent sidecar, or a failed leakage audit → exit 3 `INDETERMINATE` naming the arm.
- Emits `wiki/transfer_benchmark_falsifier_<date>.{md,json}` under the `transfer_benchmark_` prefix (namespace verified free).
- Exits: 0 = A1+B discriminated, 1 = bar not met, 3 = refused/indeterminate.

**Test strategy:**
- B and C run offline with no `D:` dependency and produce the pre-registered verdicts.
- **PASS is reachable with Arm C absent**, and Arm C alone never yields PASS — the loophole closed by test.
- Arm A1 skipped when `D:/dna_decode_cache/bloom/` is absent, and its skip forces `INDETERMINATE`; a test asserts no pass on one arm.
- `TRANSFER_POSITIVE_UNMEASURED` is present in the emitted record on every run of this plan's configuration — it is a value, not an omission.
- **No cross-panel flag exists** — a test asserts the CLI rejects `--incompatible-feature-space-control`, so the dropped control cannot be quietly reinstated without the plan changing.
- Frozen-surface assertion **after asserting each of the 5 paths exists** — `git diff --quiet` on a missing path exits 0 and tests nothing.
- Emitted JSON is re-parsed in-test.

## Execution Preview

| Wave | Steps | Parallel | Files |
|---|---|---|---|
| 0 | Step 1, Step 2, Step 3 | 3 | `eval/regime.py` + `scripts/regime_map.py` · `wiki/...json` · `eval/transfer.py` (all disjoint) |
| 1 | Step 4, Step 5 | 2 | `eval/transfer.py` · `eval/transfer_gauntlet.py` (disjoint) |
| 2 | Step 6 | 1 | `scripts/transfer_benchmark_falsifier.py` |

- **Total waves:** 3 · **max parallelism:** 3 · **critical path:** Step 3 → 4 → 6 (3 steps).
- Steps 3 and 4 share `eval/transfer.py`, hence the explicit dependency rather than co-scheduling.
- Fewer waves and a shorter critical path than v1 (4 waves / 4 steps) because the de-confounding primitives are imported rather than built.
- No network, Docker, GPU or money. Every step reversible.

## Risk Flags

1. **F1 ships with no transfer-positive arm.** Mitigated only by making the absence a first-class verdict value rather than silence. **Not eliminated:** an instrument that has never returned a positive on real data is half-validated, and the synthetic positive could pass for reasons that will not hold on real genotypes.
2. **The transcription in Step 2 is a hand-copy of a markdown table.** Mitigated by the re-parse test, which turns it from trusted to verified. Residual: if the markdown itself is wrong, both sides agree and the test passes.
3. **`skipped_groups` predicts another function's internal behaviour.** It re-implements `cv_r2`'s `tr.sum() < 5 or te.sum() < 1` condition, so a change inside `deconfound` would desynchronise it. Mitigated by testing the prediction against observed NaNs rather than against the constant, so a drift fails loudly — but it is duplicated knowledge, which is a cost this plan otherwise refuses. **The same shape applies to `WITHIN_GROUP_MIN_N = 30` mirroring `within_group_r2`'s `min_n`:** mitigated by passing it explicitly (so the two cannot diverge silently at the call site) and by asserting the 30-member boundary against measured behaviour rather than the literal, but it is a second mirrored constant.
4. **`permutation_null` has no seed parameter** (per-draw seed = loop index), so the null cannot be varied across runs. Deterministic and reproducible, but a seed-sensitivity check is impossible without touching `deconfound`, which this plan does not do.
5. **scipy is undeclared in `pyproject.toml`** and arrives transitively through `scikit-learn>=1.4`. Pre-existing — `deconfound` already relies on it and ships — but this plan increases that reliance. Plan-relevant, not plan-caused; declaring it is a one-line change deliberately left out of scope.
6. **Eight split-unit assignments are my reading of committed artifacts.** A mis-declaration would make the map confidently wrong in a new way. **A human should spot-check the eight** before they are treated as authoritative.
7. **`sentrux` absent** → no architecture gate on a change adding two `eval/` modules. Recorded, not mitigated.
8. **Baseline count does not reconcile** (5428 collected vs 5424+12 reported). `/execute-plan`'s baseline pass is authoritative.
9. **Scope creep into F3.** Finding a substrate with a shared cross-organism representation is the real unblock and is explicitly out of scope here.
10. **The synthetic arm's group sizing is now load-bearing, where v2's was free.** With the within-group floor at 30, Arm B's groups must satisfy `size ≥ 30 + max(ks)`; the `test_deconfound_package.py:23` template (4 × 80) clears it with `cap = 50`. A future edit shrinking the synthetic groups would push the verdict cell to `WITHIN_GROUP_UNSCORABLE` rather than failing — visible, but only if the reader notices the verdict value changed.
11. **`min_n = 30` is `deconfound`'s choice, not a derived statistical threshold**, and this plan adopts it rather than justifying it. Raising or lowering it would move which groups contribute to the only cell PASS reads. Left untouched deliberately — amending a tested shared primitive to suit a new consumer is the drift this plan refuses — but it means the PASS condition inherits an unexamined constant.

## Open Questions

**One, and it is a scope call rather than a technical one.**

1. **Is F1 worth building when it can currently measure only an absence?** Everything downstream needs the verdict shape — no gLM number is interpretable without it — and the `TRANSFER_POSITIVE_UNMEASURED` stamp is itself the artifact that prevents the recorded failure of a later session assuming transfer works. Against that: the instrument will not return a positive on real data until a substrate with a shared cross-organism representation exists, which is F3 work with no committed date. Declining to build a measuring device that currently measures an absence is a legitimate position; proceeding is my recommendation.

**v2's other four open questions are resolved in this body, not deferred** — recorded here so they are not re-opened:
- *Ratify the three decisions* — all three ratified. (a) is **not** an authority call: the measurement settled it, and the real fork it concealed is question 1 above. (b) now rests on the measured 8-consumer eager-import fact (Codebase Context), not on preference. (c) stands as drafted.
- *Cross-panel control* — **dropped** (Step 6), with a test asserting the flag does not exist.
- *k grid* — **derived** from the two measured floors (Step 3), replacing the asserted tuple.
- *`regime_map.py` exit code* — **stays 0** (Step 1), on the recorded red-suite precedent.

## Verification

```bash
# Per-step unit tests (offline, no D:, no network)
uv run pytest tests/test_regime_boundary.py tests/test_g2_replay_fixture.py \
              tests/test_eval_transfer.py tests/test_transfer_gauntlet.py \
              tests/test_transfer_benchmark_falsifier.py -q

# The headline fact, derived not asserted
uv run python -c "from dna_decode.eval import regime as R; print(R.organism_transfer_is_unmeasured())"

# Both new fields must RENDER, not just sit in JSON
uv run python scripts/regime_map.py

# The derived k grid and the two floors (expect a raise at cap < 1)
uv run python -c "
from dna_decode.eval.transfer import derive_k_grid, MIN_SCORED, WITHIN_GROUP_MIN_N
print('floors:', MIN_SCORED, WITHIN_GROUP_MIN_N)
print('size 80 ->', derive_k_grid(80)); print('size 35 ->', derive_k_grid(35))
try: derive_k_grid(30); print('BUG: should have raised')
except Exception as e: print('size 30 -> raises:', type(e).__name__)"

# Import posture: no eager scipy/sklearn under eval/
uv run python -c "
import sys, dna_decode.eval.transfer
assert not [m for m in sys.modules if m.startswith(('sklearn','scipy'))], 'eager heavy import'
print('import posture OK')"

# Falsifier: offline arms (expect A1 skipped -> INDETERMINATE without D:)
uv run python scripts/transfer_benchmark_falsifier.py --offline-only

# Full run (needs D:/dna_decode_cache/bloom/)
uv run python scripts/transfer_benchmark_falsifier.py
# 0 = A1+B discriminated · 1 = bar not met · 3 = refused/indeterminate

# Frozen-surface + lock (paths asserted to EXIST before diffing)
uv run python -c "
from dna_decode.eval import prospective_lock as P
p,m = P.resolve_active_lock(); print(m['lock_date'], P.verify_lock(m).ok)"
git diff --quiet HEAD -- dna_decode/eval/amr_rules.py dna_decode/data/calibrated_amr_rules.json \
  dna_decode/data/mic_tiers.py dna_decode/data/shipped_decoder_surface.py dna_decode/eval/cohort_manifest.py

# deconfound surface unchanged (8 script consumers + 1 test must not break)
uv run pytest tests/test_deconfound_package.py -q
git diff --quiet HEAD -- dna_decode/deconfound/

# Report card unchanged (27 cells)
uv run python scripts/build_validation_report_card.py && git diff --stat wiki/decoder_validation_report_card.json

# Full suite against /execute-plan's baseline
uv run pytest tests/ -q
```

**Manual verification required:** the eight `split_unit` assignments (Risk Flag 6) and the Step 2 transcription against the markdown table (Risk Flag 2).

## Repo grounding

### Captured by: brainstorm @ 2026-10-03
- Files read: `dna_decode/eval/regime.py`, `scripts/regime_map.py`, `tests/test_regime_boundary.py`,
  `dna_decode/eval/genomic_prediction.py`, `dna_decode/eval/clonality.py`, `dna_decode/eval/cv.py`,
  `dna_decode/eval/clade_baseline.py`, `dna_decode/eval/cohort_manifest.py`, `dna_decode/eval/phylogeny.py`,
  `dna_decode/fba/nulls.py`, `dna_decode/deconfound/{__init__,deconfound,scorecard}.py`,
  `tests/test_deconfound_package.py`, `scripts/{yeast_bloom_gp_arm,bxd_gp_arm,arabmagic_gp_arm,yeast_growth_decoder}.py`,
  `pyproject.toml`, `wiki/phase2_arabidopsis_result_2026-06-12.md`, `wiki/glm_phase5_decomposition_2026-10-01.md`.
- Key claims:
  - `dna_decode/deconfound/` is a canonical, exported, tested package (`tests/test_deconfound_package.py:17
    test_public_api_imports` pins the public surface) already providing `cv_r2(..., groups=)` (leave-one-group-out,
    `deconfound.py:24`), `within_group_r2` (group-centered residual r2, `:43`), `cluster_from_distance`
    (`:70`), `group_centered_spearman` (`:79`) and group-aware `permutation_null` (`:138`). The plan's
    Steps 2/3/5 would RE-IMPLEMENT the fold generator and 4 of the 5 gauntlet controls; the Reusable-Code
    Survey missed the package entirely.
  - `within_group_r2` is the metric family that produced the Arabidopsis per-seed numbers Arm C replays, so a
    second implementation in a new `transfer_gauntlet.py` would create two truth surfaces for the metric the
    replay is compared against.
  - EVERY genotype-phenotype panel on disk is a SINGLE constructed population with no internal organism axis:
    Bloom BYxRM (`D:/dna_decode_cache/bloom/`), BXD mouse (`D:/dna_decode_cache/bxd/`), Arabidopsis 19-parent
    MAGIC (`D:/dna_decode_cache/arabmagic/`, 677 lines x 1260 markers x 8 traits). Arm A therefore has nothing
    to hold out and would silently degrade to K-fold.
  - Cross-PANEL hold-out is NOT a held-out-organism arm: the three loaders produce panel-specific marker
    matrices with no shared feature ontology, orthology projection or coordinate map, so a zero result proves
    feature-space incompatibility rather than a transfer floor.
  - `phylogeny.cluster_by_ani:232` is union-find (single-linkage, chaining-prone), but
    `deconfound.cluster_from_distance:75` is AVERAGE-linkage -- so the plan's chaining argument for
    greedy-representative clustering argues against the wrong incumbent.

## Save-time amendments

Captured at: 2026-10-03
Source: `/save-plan` arguments

**Audit-notes-only contract.** Provenance for human readers. `/execute-plan` reads ONLY
`## Implementation Steps` for executable work; nothing here is an instruction. An amendment that
changes a Step's file list, dependencies or structure means the executable body is out of sync with
the amendment intent — re-run `/technical-plan` before `/execute-plan`.

- v2 supersedes v1: no transfer-positive arm exists on disk, so TRANSFER_POSITIVE_UNMEASURED is a first-class verdict
- de-confounding maths stays solely in dna_decode/deconfound and an AST guard enforces it
- ladder replaced by unordered split_unit plus a narrow organism_transfer
- G2 replay inputs become a verified machine-readable sidecar
- skipped_groups surfaces cv_r2's silent group drops

### Captured at: 2026-10-03 (v3)

- v3 lands the four Q resolutions: cross-panel control DROPPED (result knowable without running it), k grid DERIVED from two measured floors, regime_map exit stays 0 on the red-suite precedent, gauntlet home settled by the measured 8-consumer eager import
- new in v3 — within_group_r2's min_n=30 is the binding floor, not MIN_SCORED=10, so WITHIN_GROUP_UNSCORABLE becomes a required verdict and the pooled-cell fallback is closed by test

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
