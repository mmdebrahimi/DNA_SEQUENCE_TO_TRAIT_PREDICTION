## Lens status
Dispatch: conversation-driven (no argument). Subject resolved from this session: of the three decisions surfaced by `/soraya --advance`, the drive is physical and the portfolio fork is strategic — regenerating promcalc is the only technical build, and the executed tile-chain plan explicitly scoped it out ("does not run the pre-registered comparison, regenerate promcalc, install its toolchain, or touch the frozen acceptance bar"), calling it "a reachable-but-user-gated regeneration path".
Repo-index: absent — direct file reading only, no degradation.
sentrux: n/a — not installed. codex: n/a — not installed. gh: present.
graphify: absent — `### Reusable-Code Survey` uses filesystem + import-graph reading.
Wiki replay: `wiki/plans-index.md` + `wiki/decisions-log.md` read. Two prior entries bind this plan — `Make_The_Tile_Chain_S_Substrate_Legible_And_Classify_Its_Wall_Plan/` (executed 2026-10-21; left promcalc `data-absent`, user-gated) and `Tile_Substrate_Selection_Value_And_The_Untried_Hybrid_Arm_Plan/` (executed; produced the current incumbent arm and the shared reconcile anchors).
Scope: restore the promcalc input and run the EXISTING frozen pre-registered comparison. No new arms, no new acceptance bar.

## Problem Statement

The tile chain's two substrate dependencies are now cleanly separated and exactly one is left. As of 2026-10-22 the chain gets past the tile load and fails on only promcalc:

```
tiles 44106 | ceiling 0.9410
--in-dir D:/dna_decode_cache/promcalc is not reachable [absent].
```

The MPRA half was closed by re-fetching ~100 MB of GSE144621 to `data/mpra/`. The promcalc half is the per-site mechanistic decomposition — Stage B's output — which lives only on a USB volume that detached permanently and whose weights-bearing sibling killed four sweep runs.

Three things make regeneration the right move now rather than later, and all three are measured rather than argued:

1. **It needs neither the drive nor a GPU.** Stage B is CPU-only and reads the TSV Stage A writes from the now-local MPRA cache.
2. **The cost is known, not estimated.** 0.26 s/tile single-core; the recorded resume ran at 785 tiles/min on 7 workers (5,132 tiles in 6.5 min), consistent with the published 74 min for the full 44,106. Output 140 MB gzipped.
3. **The install is measured safe.** `uv pip install --dry-run promotercalculator` resolves 7 packages and would install exactly one — `promotercalculator==1.2.4` — leaving sklearn 1.8.0 untouched. This REFUTES the sklearn-1.9.1 conflict that both Stage A's and Stage B's docstrings still assert as their reason for existing.

What this unblocks: Stage C's own pre-registered comparison (does the full per-site decomposition beat the single `max_tx` scalar reduced from it?) has never run on this host, and the arms-ledger ceiling criterion is `blocked:untried-arm` with `orthogonal-modality-stack-promcalc-plus-N-channels` among the four never-run arms.

**Non-goals, stated so scope cannot drift:** this plan does not add an arm, does not touch `PREREGISTERED` in `scripts/promcalc_decomposition_arm.py`, does not restage the 4.33 GiB GENERator weights, does not touch the sweep, and does not edit any of the frozen five or `uv.lock`.

## Codebase Context

**The three-stage architecture, read from source.** The process boundary is a FILE, deliberately:

| stage | script | env | reads | writes |
|---|---|---|---|---|
| A | `scripts/promcalc_export_tiles.py` | project | `genomewide.load_peak_tiles` (owner) | `tiles.tsv` + `substrate.json` |
| B | `scripts/promcalc_decompose.py` | separate interpreter | `tiles.tsv` | `decomposition.csv.gz` (+ `part-NN` shards) |
| C | `scripts/promcalc_decomposition_arm.py` | project | the decomposition | `wiki/promcalc_decomposition_arm_*.{md,json}` |

**Stage B's contract, verified by reading the AST rather than the prose.** Module-level imports are stdlib only (`argparse csv gzip json os sys time multiprocessing pathlib`); `promoter_calculator` is imported LAZILY at line 58 inside the worker (`from promoter_calculator import Promoter_Calculator`); the only two `dna_decode` occurrences are a docstring sentence and the `D:/dna_decode_cache/promcalc` path string — neither is an import. So the stated boundary holds today and **nothing enforces it**: `grep -rln promcalc_decompose tests/` returns nothing. Stage B has zero test coverage.

**Install name differs from import name** — `pip install promotercalculator`, `import promoter_calculator`. Already recorded as a trap in the toolchain probe; a naive install of the import name 404s.

**Flags, confirmed by reading argparse (not inherited from prose).**
- Stage A: `--cache-dir` (default `genomewide.DEFAULT_CACHE`), `--out-dir` (default `D:/dna_decode_cache/promcalc`), `--no-download`, `--limit`.
- Stage B: `--in-dir` (default `D:/dna_decode_cache/promcalc`), `--workers` (default `cpu_count()-1`), `--smoke`, `--limit`, `--resume`.
- Stage C: `--in-dir` (default `D:/dna_decode_cache/promcalc`), `--mpra-cache` (added 2026-10-22), `--smoke`, `--seeds` (10), `--out`.

**Resume already exists and is shard-based, because a full run was killed once.** `shard_paths()` enumerates shards EXPLICITLY rather than globbing `decomposition*`; `--resume` skips tile_ids already present across the main file and every shard and writes a new `decomposition.part-NN.csv.gz`, because a gzip stream cannot be safely appended. The recorded cause of that kill was session teardown at 88.4% (38,974 of 44,106).

**The acceptance bar is frozen in code and must not be touched.** `PREREGISTERED` in Stage C carries `margin: 0.02`, `min_seeds_winning: 8`, `registered_seeds: 0..9`, `control: feature_permutation`, `control_max_abs: 0.05`, and `margin_is` recording that +0.02 is DERIVED from this substrate's measured effect sizes. The verdict ladder returns `INDETERMINATE_RECONCILE_FAILED` / `INDETERMINATE_CONTROL_NOT_CLEAN` / `NO_GAIN` / `DECOMPOSITION_ADDS` / `WEAK`.

**The reconcile block IS this plan's validity anchor, and it already exists.** `PREREGISTERED["reconcile"]` holds `gc_seed0` expected `0.3000` ±`0.01` (sourced to `glm_encoder_gate.PREREGISTERED.reconcile_anchor`) and `max_tx_median` expected `0.3968` ±`0.02` (sourced to `wiki/promcalc_vs_learned_2026-10-12.json`). A regenerated table that reproduces `max_tx` median within tolerance IS the published substrate; one that does not causes Stage C to refuse rather than score a different substrate. No new bar is needed and none is proposed.

**A recorded measurement trap that directly governs Step 5.** From `wiki/promcalc_decomposition_arm_2026-10-13.md`: a 90-second incremental count gave 138 tiles/min and therefore a 5.0 h projection, while the actual resume ran at 785 tiles/min — because measuring progress meant repeatedly decompressing a growing 140 MB gzip while 7 workers competed for 8 cores. *"Observing the sweep slowed the sweep, and the slowed rate was then quoted as the sweep's rate."*

**Stale rationale in two shipped docstrings.** Stage A: *"`promoter_calculator` pulls sklearn 1.9.1 while this project is pinned at 1.8.0, so it lives in an ISOLATED interpreter at `D:/tools/promcalc_env`"*. Stage B: *"installing it in the project venv would move sklearn under all 6,409 tests."* Both reasons are refuted by the dry-run. The separate environment is still correct, but for a DIFFERENT reason — see Step 3.

**`uv.lock` is a frozen reproducibility unit**, named in the 2026-06-13 freeze alongside the frozen five. That is what rules out `uv add promotercalculator`, independently of the sklearn question.

### Reusable-Code Survey

Searched: `scripts/promcalc_*.py`, `scripts/glm_tile_predictions.py`, `dna_decode/glm/genomewide.py`, `dna_decode/data/external_root.py`, `scripts/glm_run_lock.py`, `run_glm_anchor_detached.bat`, `tests/test_glm_tile_substrate.py`, `tests/test_substrate_loader_ownership.py`, `tests/test_promcalc_toolchain_probe.py`. (`graphify-out/GRAPH_REPORT.md` absent; no `src/`, `lib/`, `utils/` candidate dirs — this repo is a flat package plus `scripts/`.)

Reuse candidates, all preferred over new code:

1. **`scripts/glm_tile_predictions.py::require_external_dir`** — the shipped refusal shape for an unreachable configured directory: `resolve_external_root(...)`, return on `configured`, else `SystemExit` carrying `describe()`. Step 1 reuses this rather than writing a second resolver.
2. **The `--mpra-cache` pattern landed 2026-10-22** (flag defaulting to `None`, fallback to the module default resolved AT POINT OF USE, `except GenomeWideDataError → SystemExit(...) from None`). Step 1 applies the same shape to the promcalc directory.
3. **`dna_decode/glm/genomewide.ensure_cached`** — already carries both the volume-reachability refusal and (as of 2026-10-22) the file-named-as-directory refusal. Stage A inherits both through `load_peak_tiles`; no new validation needed on the MPRA side.
4. **`scripts/glm_run_lock.py run`** — PID-liveness one-job-at-a-time wrapper that holds the lock for the CHILD's lifetime. Step 5 reuses it; the documented failure it avoids is a launcher's `acquire` process exiting immediately so its recorded pid is already dead.
5. **`run_glm_anchor_detached.bat` launch discipline** — `Start-Process ... -WindowStyle Hidden` from PowerShell, plus the sentinel convention (`EXIT=-1` = ran and crashed, absent sentinel = killed before the exit line). Step 5 reuses the convention; it does NOT reuse the `.bat` itself, which is sweep-specific.
6. **`tests/test_glm_tile_substrate.py::_calls` / `_callers_of`** — the derived-caller-set AST helpers added 2026-10-22. Step 2 reuses `_calls`-shaped AST walking rather than regexing imports.
7. **`scripts/promcalc_decompose.py::shard_paths` / `done_tile_ids`** — existing resume machinery. Step 5 uses `--resume` as-is; no new resume logic.

Not reused, deliberately: `scripts/promcalc_toolchain_probe.py` is a read-only diagnostic and stays one — it must not become a provisioner.

## Pre-Change Baseline

Measured 2026-10-22, this host:

- **promcalc input: ABSENT.** `scripts/glm_tile_predictions.py --mpra-cache data/mpra` reaches the promcalc check and refuses: `--in-dir D:/dna_decode_cache/promcalc is not reachable [absent]`, `discovered candidates: none on any resolving drive`. The 88.4%-complete shard set from the killed run is on the detached volume and is NOT recoverable.
- **MPRA substrate: PRESENT and validated.** `data/mpra/` holds 100,602,188 bytes across three files; `load_peak_tiles(cache_dir='data/mpra')` → 46,713 rows, 44,106 150-bp tiles, noise ceiling 0.9410 — three figures matching values recorded before the re-download.
- **Stage C has never run on this host.** Last run is `wiki/promcalc_decomposition_arm_2026-10-13.{md,json}`, produced on the machine that had D:.
- **Stage B has zero tests.** `grep -rln promcalc_decompose tests/` → no matches.
- **Project env:** sklearn 1.8.0, Python 3.11.17 under `uv`.
- **Suite:** 14 failed / 6920 passed / 186 skipped — the 14 being the standing data-absent set (cyp4f2 7, essentiality-ladder 5, pneumo cps 1, README-quickstart 1).
- **Arms ledger:** the `ceiling:` criterion is `blocked:untried-arm`; 7 measured arms; best `promcalc-plus-conv-rank-average-hybrid-150bp`; 4 planned-and-never-run, one of which (`orthogonal-modality-stack-promcalc-plus-N-channels`) is gated on exactly this plan's output.
- **Frozen surface:** the frozen five byte-unchanged vs HEAD; active prospective lock `wiki/prospective_lock_manifest_2026-08-31.json` verifies.
- **Disk:** C: 58.5 GB free of 475.9 GB. D: does not resolve.

## Verification Signal

What makes this plan succeeded-or-failed, decided before running it:

1. **The regenerated table reproduces the published substrate.** Stage C's frozen reconcile gate passes on BOTH anchors: `max_tx` median within 0.02 of 0.3968, `gc` seed 0 within 0.01 of 0.3000. This is the load-bearing signal — it is what distinguishes "regenerated the same decomposition" from "computed a different one that happens to run". A failure here is a REAL finding (likely `promoter_calculator` version drift) and must be reported as such, never tuned away.
2. **Stage C emits a verdict from the frozen ladder** — `DECOMPOSITION_ADDS` / `WEAK` / `NO_GAIN` — rather than `INDETERMINATE_RECONCILE_FAILED` or `INDETERMINATE_CONTROL_NOT_CLEAN`. **Any of the three substantive verdicts is a successful plan**; `NO_GAIN` answers the pre-registered question negatively and is not a plan failure.
3. **Coverage is complete, not partial.** Stage B's reported coverage over the whole substrate is 44,106/44,106 tiles, and the site count is ~73 forward candidates per tile (the published figure).
4. **The project environment is untouched.** `uv run python -c "import sklearn; print(sklearn.__version__)"` still prints `1.8.0`, and `git diff --stat uv.lock pyproject.toml` is empty.
5. **No regression.** Suite still 14 failed / 6920+N passed / 186 skipped, with the same 14 names.
6. **Stage B's stated boundary is now enforced** — the new guard fails if Stage B gains a `dna_decode` import or a module-level third-party import.

## Implementation Steps

### Step 1: Resolve the promcalc directory at point of use instead of a hardcoded D: default
Files: scripts/promcalc_export_tiles.py, scripts/promcalc_decompose.py, scripts/promcalc_decomposition_arm.py
Depends on: none

**What changes:**
- All three stages hardcode `D:/dna_decode_cache/promcalc` (`DEFAULT_OUT` in A, `DEFAULT_DIR` in B and C). Replace the bare constant with resolution AT POINT OF USE, mirroring the `--mpra-cache` shape shipped on 2026-10-22: an explicit flag stays authoritative; when the flag is absent, fall back to the module default and, if that is unreachable, refuse with the `external_root.describe()` report rather than an `OSError`.
- Stage A and C are project-env and may import `dna_decode.data.external_root`. **Stage B may NOT** — it is the isolated stage and importing `dna_decode` would break the boundary Step 2 enforces. Stage B therefore gets a self-contained reachability check (`Path.exists()` on the directory or an ancestor) and a refusal message that names the flag to pass, with NO `dna_decode` import.
- Keep the D: path as the recorded default in all three, so a host that HAS the volume behaves bit-identically. This is parameterisation, not relocation.

**Test strategy:**
- Unit: each stage's `--in-dir`/`--out-dir` default is unchanged (the D: string) — pinning that this is additive.
- Unit: an unreachable directory yields a refusal naming the flag, not an `OSError`/`FileNotFoundError` traceback — asserted per stage.
- Unit (Stage B specifically): the refusal path is reachable without importing `dna_decode`, exercised by parsing the module and asserting the check uses only stdlib.
- Manual: `uv run python scripts/promcalc_decomposition_arm.py --in-dir Z:/nope --mpra-cache data/mpra` prints a refusal and exits non-zero.

### Step 2: Guard Stage B's stated import boundary, which has zero tests today
Files: tests/test_promcalc_decompose.py
Depends on: none

**What changes:**
- New test file. Stage B's docstring asserts *"This module therefore imports NOTHING from `dna_decode`... Keep it that way"* and nothing enforces it — the same shape as the `--help` over-claims this repo has repeatedly found: a stated invariant with no mechanism.
- AST-based guards: (a) module-level imports are stdlib-only; (b) no `dna_decode` import appears anywhere in the module, at any nesting depth; (c) the `promoter_calculator` import is LAZY (inside a function, not module level) so the project env can import the module to test it without the package installed; (d) the import name is `promoter_calculator` while the distribution is `promotercalculator` — pinned so the 404-on-naive-install trap stays recorded next to the code.
- Non-vacuity: synthetic source containing a `dna_decode` import, and synthetic source with a module-level `promoter_calculator` import, must each trip the respective detector.
- Reuse `_calls`-style AST walking rather than regexing import lines; a regex cannot tell a docstring mention from an import, and this module contains exactly that ambiguity (two `dna_decode` occurrences, neither an import).

**Test strategy:**
- The file IS the test. Must pass against Stage B as it stands today (the boundary currently holds), proven non-vacuous by the synthetic cases above.
- Guard-the-guard: assert the discovered module-level import set is non-empty, so a parse failure cannot read as "no third-party imports".

### Step 3: Provision the Stage-B interpreter on a reachable path, and correct the refuted isolation rationale
Files: scripts/promcalc_decompose.py, scripts/promcalc_export_tiles.py, wiki/promcalc_regeneration_runbook_2026-10-22.md, wiki/promcalc_toolchain_probe_2026-10-22.json
Depends on: Step 1, Step 2

**What changes:**
- Create a dedicated venv on a reachable path (`uv venv` + `uv pip install promotercalculator` into it). The documented interpreter `D:/tools/promcalc_env/Scripts/python.exe` is on the detached volume and is gone.
- **Rewrite the isolation rationale in both docstrings, because the stated reason is measured FALSE and the conclusion is nonetheless unchanged.** The sklearn-1.9.1 conflict does not exist (`--dry-run`: one package, sklearn 1.8.0 untouched). A separate environment is still correct for a different and now-primary reason: **`uv.lock` is a frozen reproducibility unit**, so `uv add` is not available, and an imperative `uv pip install` into the project env is silently reverted by the next `uv sync`. A separate venv is near-free and keeps the boundary Step 2 enforces. Record the refutation beside the new reason rather than deleting the old claim — the same discipline the wall memo's amendments use.
- Emit a fresh `wiki/promcalc_toolchain_probe_2026-10-22.json` by re-running the read-only probe, because the committed 2026-10-21 record still carries the refuted `isolation_reason` string verbatim. Leave the 10-21 artifact in place, dated.
- Write the runbook: the exact commands, the interpreter path, the measured cost (785 tiles/min on 7 workers), and the polling trap from Step 5 stated up front.

**Test strategy:**
- Real-surface: the provisioned interpreter resolves `from promoter_calculator import Promoter_Calculator` and reports its `__version__` / distribution version; recorded in the runbook. This is the R3 "first integration test is the real surface" check — a green mock proves nothing about a lazily-imported third-party package.
- Real-surface: `uv run python -c "import sklearn; print(sklearn.__version__)"` in the PROJECT env still prints `1.8.0` AFTER the install, and `git diff --stat uv.lock pyproject.toml` is empty. Belt-and-braces: the dry-run predicted this, and a prediction is not a measurement.
- Unit: a docstring guard asserting neither Stage A nor Stage B still claims the sklearn-1.9.1 conflict as its reason (the over-claim pattern), and that both name `uv.lock`'s frozen status instead.

### Step 4: Stage A export plus a bounded Stage B smoke before committing an hour
Files: wiki/promcalc_regeneration_runbook_2026-10-22.md
Depends on: Step 3

**What changes:**
- Run Stage A against the restored local MPRA cache, writing `tiles.tsv` + `substrate.json` to the resolved promcalc directory. Stage A's own `EXPECTED_TILES = 44106` / `EXPECTED_CEILING = 0.941` asserts fire here, so a drifting loader fails before Stage B starts.
- Run Stage B with `--smoke` (and/or `--limit`) on a bounded slice, which writes `*_smoke.*` and cannot collide with the full run's shard namespace.
- Verify the smoke CSV header equals `SITE_FIELDS` exactly (`tile_id, tss, dG_total, dG_10, dG_35, dG_ext10, dG_spacer, dG_disc, dG_ITR, dG_UP, dG_bind, Tx_rate, hex35_start, hex10_start, spacer_len`) and that the per-tile site count is ~73 forward candidates. A header mismatch means the installed `promoter_calculator` returns different fields — i.e. version drift — and that is far cheaper to learn here than after an hour.
- Record the observed smoke throughput, measured WITHOUT decompressing the output (see Step 5).

**Test strategy:**
- Manual, with recorded output: the smoke run's header, site-per-tile count, and elapsed time land in the runbook.
- Unit: a pure assertion that the committed `SITE_FIELDS` tuple matches the header the runbook records, so the two cannot drift apart silently.
- Gate: if the header or the site count disagrees with the published shape, STOP and report version drift. Do not proceed to Step 5.

### Step 5: Full Stage B regeneration, detached, resume-safe, and not polled
Files: wiki/promcalc_regeneration_runbook_2026-10-22.md
Depends on: Step 4

**What changes:**
- Run all 44,106 tiles at `--workers 7`. Expected ~56-74 min from the recorded 785 tiles/min; output ~140 MB gzipped.
- **Launch detached via `Start-Process -WindowStyle Hidden`, not from a session-parented shell.** The recorded cause of the previous full run's death at 88.4% was session teardown, and three sweep runs died the same way this week. Wrap in `scripts/glm_run_lock.py run` so one job runs at a time and the lock is held for the CHILD's lifetime.
- **Do NOT poll progress by decompressing the output.** This is the recorded trap: doing so turned a ~1 h job into a 5.0 h projection and the slowed rate was then quoted as the sweep's rate. Progress comes from the worker's own stderr/stdout line count or a wall-clock checkpoint, never from re-reading the growing gzip.
- On any death, re-launch with `--resume`, which skips tile_ids present across the main file and every shard and writes a new `part-NN`. Resume is existing machinery; this plan adds none.
- Write a sentinel on exit, reusing the established convention so the failure mode is legible: `EXIT=-1` means ran-and-crashed, an ABSENT sentinel means killed before the exit line.

**Test strategy:**
- Manual, gated: Stage B's returned record reports coverage 44,106/44,106 over the whole substrate (not over a resume slice — the script already distinguishes these).
- Manual: the output directory holds the main gzip plus any shards, and `shard_paths()` enumerates exactly those (never a bare `decomposition*` glob).
- Negative: deliberately re-run with `--resume` after completion and confirm it reports "nothing to resume", proving the done-set logic sees the full table.

### Step 6: Validate against the frozen reconcile anchors, then run the pre-registered comparison
Files: wiki/promcalc_decomposition_arm_2026-10-22.md, wiki/promcalc_decomposition_arm_2026-10-22.json
Depends on: Step 5

**What changes:**
- Run Stage C over the regenerated decomposition with `--seeds 10` (the registered seed set) and `--mpra-cache data/mpra`.
- The reconcile gate fires FIRST and is the plan's validity anchor: `gc` seed 0 within 0.01 of 0.3000, `max_tx` median within 0.02 of 0.3968. Reproducing `max_tx` is what establishes that the regenerated table IS the published substrate rather than a different computation that merely runs.
- Then the frozen verdict ladder reports `DECOMPOSITION_ADDS` / `WEAK` / `NO_GAIN` against `margin: 0.02` and `min_seeds_winning: 8`, with the feature-permutation control checked against `control_max_abs: 0.05`.
- **`PREREGISTERED` is not edited.** If the result is `NO_GAIN`, that is the answer to the pre-registered question and it ships as the answer. Re-sizing a frozen bar after seeing the result is the recorded failure this project has twice had to undo.
- The artifact records the regeneration provenance explicitly: the installed `promoter_calculator` version, that the input was regenerated locally rather than inherited, and that the original shard set was lost with the volume.

**Test strategy:**
- Unit: the emitted verdict equals the frozen `_verdict` function applied to the artifact's own recorded inputs — so the written verdict cannot disagree with its own numbers.
- Unit: the artifact cites both reconcile anchors with their recorded expected values and sources, and does not invent a band.
- Unit: the artifact records the `promoter_calculator` version, since its absence from the 2026-10-12 record is the reason version drift is an open risk at all.

### Step 7: Record the wall-class change and close the ledger loop
Files: wiki/tile_chain_substrate_wall_2026-10-21.md, wiki/tile_chain_substrate_wall_2026-10-21.json, project_state/glm-umbrella-2026-10-07.md
Depends on: Step 6

**What changes:**
- Append amendment 5 to the wall memo and `amendment_5` to its sidecar: the promcalc dependency moves off `data-absent`, with the regeneration provenance and the reconcile outcome. The memo exists to record which dependency has which wall class, so leaving it stale after the second one resolves would make it wrong — the failure mode it was written to prevent.
- Do NOT edit the frozen `dependencies` array; amend, as amendments 3 and 4 did.
- Append one Action Log row to the umbrella via `/project-state --append-action --class run-tests`, now that the Step-4 path gate is unblocked (all 15 ledgers were corrected 2026-10-22).
- Update the arms-ledger note for the ceiling criterion: this plan does not retire it — `orthogonal-modality-stack-promcalc-plus-N-channels` becomes RUNNABLE, which is not the same as run. State that distinction explicitly so the next `--advance` does not read the criterion as settled.

**Test strategy:**
- Unit: amendment 5 pins the wall class from the memo's own closed vocabulary (`code-closable` / `data-absent` / `external-wall`) — the first draft of amendment 3 invented a fourth spelling and this guard exists because of it.
- Unit: the memo prose and the sidecar both carry amendment 5 (both or neither, as amendment 4 is pinned).
- Manual: `git diff` shows no change to `PREREGISTERED`, to the frozen five, or to `uv.lock`.

## Execution Preview

Waves (topological, from the `Depends on` edges):

| wave | steps | parallelism |
|---|---|---|
| 0 | Step 1, Step 2 | 2 |
| 1 | Step 3 | 1 |
| 2 | Step 4 | 1 |
| 3 | Step 5 | 1 |
| 4 | Step 6 | 1 |
| 5 | Step 7 | 1 |

Total waves: 6. Max parallelism: 2. Critical path: 1 → 3 → 4 → 5 → 6 → 7 (6 of 7 steps).

Intra-wave file overlap: none. Wave 0's two steps touch disjoint files (three `scripts/promcalc_*.py` vs one new `tests/` file).

**Recommended execution mode: SEQUENTIAL, overriding the parallelism.** Wave 0 offers only 2-way parallelism, and `/data/*` is gitignored with force-added files, which makes a worktree's absent-data failure set unreconcilable against the parent — the same reason the 2026-10-21 tile-chain plan ran sequentially. The one real cost in this plan is Step 5's wall-clock hour, which no amount of worktree parallelism shortens.

## Risk Flags

1. **`promoter_calculator` version drift — the single biggest risk, and it is `[unverified]`.** The version used for the 2026-10-12 run is recorded NOWHERE; only today's PyPI `promotercalculator==1.2.4` is known. If the model's output changed, `max_tx` median drifts and Step 6's reconcile gate FAILS. Mitigation is already built rather than added: the gate refuses instead of scoring a different substrate, and Step 4's header/site-count check catches a shape change in minutes rather than after an hour. Step 6 records the version so this risk cannot recur.
2. **Cumulative local disk against the standing D: rule.** This plan adds ~140 MB on top of the ~100 MB MPRA cache, so ~240 MB total on C: (58.5 GB free). Two orders of magnitude below the GB-scale assets the rule targets, and `data/` already holds 4.0 GB here — but it is now the SECOND such judgement in one day and should be visible as a running total, not re-argued each time.
3. **Step 5 is a ~1 h detached job on a host that killed four this week.** Three of those four were the USB paging fault, which is structurally removed here (nothing is memory-mapped off the drive; the drive is not involved). The fourth cause — session teardown — is NOT removed, which is why Step 5 mandates `Start-Process` plus `--resume`.
4. **`NO_GAIN` is a likely and acceptable outcome.** The incumbent `max_tx` scalar beat every learned arm, and the margin is deliberately set above the learned arms' CI-including-zero band. A negative answer to the pre-registered question is the plan succeeding, and Verification Signal 2 says so explicitly, so it cannot be retro-read as failure.
5. **Step 3 rewrites shipped docstrings.** Low blast radius (prose, no behavior) but it edits the stated reason two modules exist. The refutation is recorded beside the new reason rather than replacing it silently, so a future reader can see the correction rather than wondering which account is right.
6. **No `project-rules.md`** in this repo, so the `## Project Rules Check` section is omitted per the skill contract. `sentrux` is absent, so no architecture gate runs.
7. **Stage B gets tests for the first time (Step 2) and they are AST-only.** They pin the import boundary and the laziness of the third-party import; they do NOT exercise Stage B's decomposition logic, which needs the real package and is covered by Steps 4-5 at the real surface instead. Stated so the coverage is not over-read.

## Open Questions

1. **Where should the regenerated 140 MB live?** The plan assumes a gitignored repo-local directory (consistent with `data/mpra/`), but a different reachable location is equally valid and is the user's call if the running C: total matters.
2. **Should the lost 88.4% shard set be treated as recoverable?** The plan assumes not — it was on the detached volume. If that volume returns, resuming from those shards would cut Step 5 to roughly seven minutes, which changes the cost calculus but not the plan's shape.
3. **Does the `orthogonal-modality-stack-promcalc-plus-N-channels` arm run in this plan or the next?** Scoped OUT here deliberately: this plan restores the input and answers the EXISTING pre-registered question. Running a new arm is a separate pre-registration and should not ride along on a substrate-restoration plan.
4. **Should `promotercalculator` be pinned?** Unpinning is how drift reaches Step 6. Pinning it in the runbook is free; pinning it in `pyproject.toml` is not available, because `uv.lock` is frozen. Recording the installed version (Step 6) is the minimum; whether to also pin the runbook's install command to `==1.2.4` is a judgement the user may want to make.

## Verification

- Steps 1-3 and 7 are verified by the unit tests named in each step, run with `uv run pytest tests/ -q`.
- Steps 4-6 are verified at the REAL surface — the provisioned interpreter, the real export, the real decomposition, the real Stage C run — because every one of them depends on a lazily-imported third-party package whose behavior a mock cannot establish.
- The plan is complete when: the reconcile gate passes on both anchors; Stage C emits one of the three substantive verdicts; coverage is 44,106/44,106; sklearn in the project env is still 1.8.0 with `uv.lock` and `pyproject.toml` byte-unchanged; the frozen five and the active prospective lock still verify; and the suite shows the same 14 pre-existing failures and no new ones.

## Save-time amendments

Captured at: 2026-10-22
Source: `/save-plan` arguments

**Audit notes only.** `/execute-plan` reads ONLY `## Implementation Steps` for executable work; the
bullets below are provenance for human readers and are not executable instructions.

- subject resolved to promcalc regeneration (the only technical build of the three open decisions)
- separate venv RETAINED but its rationale corrected — the sklearn-1.9.1 conflict is refuted, the real reason is that uv.lock is a frozen reproducibility unit
- the reconcile gate adopted as the validity anchor rather than writing a new bar
- sequential execution recommended over the computed 2-way parallelism
- running the orthogonal-modality-stack arm and pinning promotercalculator left as open questions rather than folded into scope

### Superseded by measurement after the save (2026-10-22, same day)

The bullets above are preserved verbatim as the author wrote them. Two were overtaken by measurement
hours later, during the execution of Steps 3-5. Recorded here rather than by editing either the
bullets or the plan body, so the correction is visible instead of silently applied.

- **Amendment 2's "the sklearn-1.9.1 conflict is refuted" is TOO STRONG, and the plan body's own
  Problem Statement item 3 inherits the same over-claim.** Measured in a FRESH environment,
  `promotercalculator==1.2.4` genuinely DOES resolve `scikit-learn==1.9.1` (with numpy 2.4.6, scipy
  1.17.1, biopython 1.88). The project-env `--dry-run` that the refutation rested on stayed at
  sklearn 1.8.0 only because the project's existing pin already satisfied the resolver. Both facts
  are true and the original docstring rationale was closer to right than the correction was. The
  amendment's CONCLUSION stands untouched -- the separate venv is retained -- but it is now correct
  for BOTH reasons: the fresh-resolve conflict is real AND `uv.lock` is a frozen reproducibility unit.
- **The plan body's cost figures are wrong by more than 10x and must not be used for budgeting.**
  Problem Statement item 2 ("the cost is known, not estimated"), Step 5's "~56-74 min", and Risk Flag
  3's "~1 h detached job" all derive from the recorded 0.26 s/tile and 785 tiles/min in
  `wiki/promcalc_decomposition_arm_2026-10-13.md`. Measured on this host: **8.7 core-seconds per
  tile**, and the real run projects **~11 h at 65 tiles/min** on 7 thread-pinned workers (~26 h
  unpinned). 785 tiles/min would require ~0.5 core-seconds/tile, a ~17x gap no thread setting
  explains; left UNRESOLVED rather than reconciled by argument.
- **Risk Flag 1 (`promoter_calculator` version drift, marked `[unverified]`) MATERIALIZED.** Two
  independent breaks on 1.2.4, each fatal to every tile: the CLASS no longer re-exports at the
  package root, and each per-TSS value is now a `PromoCalcResults` object with no `.get`. Both were
  caught by Step 4's bounded smoke gate in 24 s rather than after the full run -- which is the
  specific thing that gate was put in the plan to do. The class itself is unchanged, so stage C's
  frozen reconcile anchors remain the right acceptance test.
- **Amendment 4 (sequential execution) was followed.** Steps 1-5 ran sequentially in the parent
  checkout; no worktree was created.
- Measurements, commands and honest limits: `wiki/promcalc_regeneration_runbook_2026-10-22.md` and
  its `.json` sidecar.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified -->
