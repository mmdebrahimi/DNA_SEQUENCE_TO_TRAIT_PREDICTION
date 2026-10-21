# Make the Tile Chain's Substrate Legible and Classify Its Wall — Technical Plan

> A pre-registered run has sat unexecuted since 2026-10-17 behind a substrate nobody had located; half the blocker turns out to be a hardcoded drive letter and half is genuinely absent data.

## Lens status
Inputs: conversation (the VOI question, and a correction to my own answer to it)
Degradations: repo-index absent (`scripts/repo-index-status.sh` not installed) — direct file reading only; sentrux unavailable; `python3` is a non-functional WindowsApps stub, toolkit run under `python`

## Problem Statement

`wiki/glm_tile_hybrid_primary_acceptance_bar.json` is a frozen pre-registration — authored 2026-10-17, `frozen_before_any_number_was_seen: true`, with an anchor-reconciliation refusal gate — and its own status field reads `PRE-REGISTRATION ONLY. No arm has been run.` It asks the strategic question the GLM arm is currently sitting on: substituted as primary scorer into the unchanged bar, does `hybrid_conv_maxtx` clear bar A at k=10? Its verdicts are `ORACLE_READY_AS_SCORING_HEAD` (de-prioritise scoring-head work relative to the generator) versus `SELECTION_IS_THE_BOTTLENECK` (prediction work IS justified). Either answer redirects effort.

It has not run because its substrate was never located. **I asserted in conversation that the substrate was present, and that was wrong** — I had found `E:/PythonProjects/dna_decode_cache/mpra` and taken it for the input the chain consumes. Measurement says there are **two** substrate dependencies with **different wall classes**:

- **The MPRA peak-tile data — present, misaddressed.** `scripts/glm_tile_predictions.py:282` calls `load_peak_tiles()` with **no arguments**, so it takes `genomewide.DEFAULT_CACHE = D:/dna_decode_cache/mpra` unconditionally. `D:` is a card-reader slot with no media, so `ensure_cached` reaches `d.mkdir(parents=True, exist_ok=True)` at `genomewide.py:120` and dies with an unhandled `PermissionError: [WinError 21] The device is not ready`. The file it wants, `GSE144621_peak_tile_expression_formatted_std.txt.gz`, **is on this host** at `E:/PythonProjects/dna_decode_cache/mpra`. `ensure_cached` takes a `cache_dir` parameter; nothing passes one.
- **The promcalc per-tile predictions — genuinely absent.** `--in-dir` defaults to `D:/dna_decode_cache/promcalc`; a probe of every attached drive finds no `promcalc` anywhere, and `scripts/promcalc_decompose.py` additionally wants a separate interpreter at `D:/tools/promcalc_env/Scripts/python.exe`.

The bar's own `execution_notes` say it is `not_runnable_on: the SHNEITZER laptop: no D: data tree, no uv, no Docker`. This host **is** `Shneitzer`, the `no uv` clause is stale (uv has run the whole suite today), and "no D: data tree" is right in substance but conflates a reachable-but-misaddressed dependency with an absent one.

**Scope is diagnosis and legibility, not execution.** Non-goals, each deliberate: this plan does **not** run the pre-registered comparison, does **not** regenerate promcalc, does **not** install the promcalc toolchain, and does **not** touch the frozen acceptance bar or the five sha256-pinned AMR files.

## Codebase Context

- `dna_decode/glm/genomewide.py` — `DEFAULT_CACHE = Path("D:/dna_decode_cache/mpra")` (`:53`); `ensure_cached(cache_dir=DEFAULT_CACHE, *, which, allow_download)` (`:113`), `load_fragments` (`:132`) and `load_peak_tiles(cache_dir=DEFAULT_CACHE, ...)` (`:318`) all take the directory as a parameter. `ensure_cached` calls `mkdir` **before** any reachability check, which is why an absent volume surfaces as a raw `WinError 21` rather than a diagnosis.
- `scripts/glm_tile_predictions.py` — `--in-dir` default `D:/dna_decode_cache/promcalc` (`:251`), `--out-dir` default `DEFAULT_OUT = D:/dna_decode_cache/promcalc/tile_predictions` (`:48`), writes `seed-NN.json.gz` + `manifest.json`. Line `:282` is `tiles_all, _ = load_peak_tiles()` — the un-parameterised call. Line `:298` embeds a `D:/tools/promcalc_env/...` hint string.
- `scripts/glm_tile_selection.py` — reads `seed-*.json.gz` from `--in-dir`, **refuses without `manifest.json`** ("without the manifest there is no record of WHICH table this is"), and calls `require_real(tables)` so a smoke table cannot score. Holds the frozen `PREREGISTERED` block (`:61`) and `verdict()` (`:108`) that the acceptance bar restates verbatim.
- `dna_decode/data/external_root.py` — exactly the machinery this needs, and already written: `probe()` lists every signature dir **with its entries** (a name match is not evidence of the right contents), `resolve_external_root(configured, *, allow_discovered=False, required=())` returns `(path, status)` across four states and **refuses to substitute by default**, and `describe()` renders a human report that never substitutes. Live here: one candidate, `E:/PythonProjects/dna_decode_cache`, holding `('amr_portal', 'mpra')`.
- Consumers of `external_root` today: `scripts/data_inventory.py` only. The tile chain does not use it.
- `~/.claude/skills/soraya/scripts/constraint_expiry.py` — reachable; `register(..., expiry_kind, expiry_arg)` over kinds `gpu-vram-gte` / `file-exists` / `test-exit-0` / `manual`, store at `soraya_runs/constraint-expiry.json` (schema 1, **11** constraints, none covering this).
- Test infrastructure: six `tests/test_glm_tile_*.py` files, **64 of which pass with the substrate absent** — they are pure/offline, so Step 1's changes are testable here. Established conventions: an explicit non-vacuity partner per guard, and a pre-committed verdict vocabulary for any probe.

### Reusable-Code Survey
- `dna_decode/data/external_root.py` — the resolver, the four-state status, and the never-substitute default. Step 1 reuses it rather than writing path handling; building a second resolver beside it is the duplication this repo has a recorded rail against.
- `dna_decode/glm/genomewide.py` — `ensure_cached` / `load_peak_tiles` already take `cache_dir`, so Step 1 threads an existing parameter rather than adding a mechanism.
- `~/.claude/skills/soraya/scripts/constraint_expiry.py` — the expiry register and its four kinds; Step 4 appends a 12th constraint rather than inventing a blocker record.
- `scripts/data_inventory.py` — the one existing `external_root` consumer; its call shape is the pattern Step 1 follows.
- None — searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/`, `lib/`, `utils/`, `common/` (absent; package root is `dna_decode/`) for a generic cache-path resolver other than `external_root`.

## Pre-Change Baseline

Measured 2026-10-21 at `bbb3612e`:

- `uv run python scripts/glm_tile_predictions.py` → unhandled `PermissionError: [WinError 21] The device is not ready: 'D:\dna_decode_cache\mpra'`, raised from `genomewide.py:120` in `ensure_cached`. A traceback, not a diagnosis.
- `uv run python scripts/glm_tile_selection.py` → `no seed tables in D:\dna_decode_cache\promcalc\tile_predictions; run scripts/glm_tile_predictions.py first`. Actionable in form and **misleading in fact**: running that script also fails, for an unrelated reason, and neither message says the volume is absent.
- `external_root.probe()` → exactly one candidate: `E:/PythonProjects/dna_decode_cache`, entries `('amr_portal', 'mpra')`. No `promcalc` on any attached drive.
- `E:/PythonProjects/dna_decode_cache/mpra` holds all three GSE144621 files, including the peak-tile table the chain needs.
- Tile-chain tests: 64 passed, substrate absent.
- Full suite: 14 failed / 6834 passed / 186 skipped; all 14 pre-existing absent-data or host-dependent.
- `soraya_runs/constraint-expiry.json`: schema 1, 11 constraints, **none** naming the tile chain's substrate.

## Verification Signal

1. With the substrate absent, `glm_tile_predictions.py` prints the `external_root` report — configured path marked ABSENT, the discovered candidate named with its contents — and exits non-zero. No raw `WinError` traceback.
2. **The load-bearing one.** With `--mpra-cache E:/PythonProjects/dna_decode_cache/mpra`, the script gets **past** the tile load and fails only on the absent promcalc input. That single run is what decides the wall classification, so the plan's central empirical claim stands or falls on it rather than on argument.
3. Behaviour is bit-identical when a configured path exists: `resolve_external_root` returns `configured` and no report is printed.
4. The toolchain probe emits one of its pre-committed verdicts, with the verdict vocabulary frozen before anything was read.
5. The wall memo names both dependencies with **distinct** classes and states what is needed and from whom; it quotes no number from the frozen bar as if measured here.
6. A 12th registered constraint exists with an auto-checkable `file-exists` expiry, and its stated limitation is the drive-letter one.
7. Full suite: 14 failed / ≥6834 passed plus the new cases, no committed-artifact churn, frozen five byte-unchanged.

## Implementation Steps

### Step 1: Route the tile chain's two substrate paths through external_root
Files: dna_decode/glm/genomewide.py, scripts/glm_tile_predictions.py, tests/test_glm_tile_substrate.py
Depends on: none

**What changes:**
- `ensure_cached` checks the directory is reachable **before** `mkdir`. On an absent/not-ready volume it raises `GenomeWideDataError` carrying `external_root.describe(cache_dir, required=('mpra',))`, so the failure names the configured path, its ABSENT status, and the discovered candidate with its entries. Today the same condition is an unhandled `WinError 21`.
- `scripts/glm_tile_predictions.py` gains `--mpra-cache` (default `genomewide.DEFAULT_CACHE`) and passes it to `load_peak_tiles(cache_dir=...)` at `:282`, which is the un-parameterised call that produces the current crash.
- Both `--mpra-cache` and `--in-dir` resolve through `resolve_external_root(..., required=(...))` and **refuse** on a non-`configured` status with the `describe()` report.
- **`allow_discovered` stays False.** The resolver must not silently adopt `E:` — a same-named root with different contents is the exact failure `external_root` was written about, and auto-substitution would convert a loud missing-ROOT into a quiet missing-DATASET. The operator passes the path; the tool only says what is there.
- Keep the behaviour bit-identical when the configured path exists, so nothing changes on the host that has `D:`.

**Test strategy:**
- Unit: an absent cache dir raises `GenomeWideDataError` whose message contains the configured path and the word ABSENT — not a `PermissionError`.
- Unit (non-vacuity): a cache dir that EXISTS takes the unchanged path and raises nothing, proving the guard is not refusing everything.
- Unit: `resolve_external_root` is called with `allow_discovered` absent-or-False — pinned, because a later "convenience" flip is the silent-substitution failure.
- Unit: injected candidates (the seam added to `describe`/`probe` earlier today) keep these tests hermetic, so the verdict does not track whichever volumes are attached.

### Step 2: Probe whether the promcalc toolchain can exist on this host
Files: scripts/promcalc_toolchain_probe.py, tests/test_promcalc_toolchain_probe.py
Depends on: none

**What changes:**
- A **read-only** probe that reports, in one artifact: whether `D:/tools/promcalc_env/Scripts/python.exe` resolves; what `scripts/promcalc_decompose.py` and `scripts/promcalc_export_tiles.py` import at module level; whether those imports are satisfiable by the current venv; free space on every attached volume; and whether the raw inputs `promcalc_export_tiles.py` needs are present.
- **Verdict vocabulary frozen before anything is read**, in the module: `ENV_PRESENT` / `ENV_ABSENT_DEPS_SATISFIABLE_HERE` / `ENV_ABSENT_DEPS_UNSATISFIABLE` / `INPUTS_ABSENT` / `INDETERMINATE_PROBE_FAILED`. A probe that picks its vocabulary after seeing the answer measures the analyst.
- **It installs nothing.** Reporting that an install is possible is a different claim from having done one, and this plan makes only the first.
- An `INDETERMINATE_PROBE_FAILED` path exists so a broken probe cannot be read as a clean negative — the repo has a recorded instance of a scorer erroring and printing a confident verdict from zero resolved calls.

**Test strategy:**
- Unit: each verdict is reachable from an injected filesystem/import state; the ladder's order is pinned.
- Unit (non-vacuity): a probe that resolves nothing returns `INDETERMINATE_PROBE_FAILED`, not `INPUTS_ABSENT` — "I could not look" and "I looked and it is not there" are different findings.
- Unit: the probe performs no writes and no installs — assert via a monkeypatched `subprocess` that refuses to be called.

### Step 3: Classify the wall and record both halves
Files: wiki/tile_chain_substrate_wall_2026-10-21.md, wiki/tile_chain_substrate_wall_2026-10-21.json
Depends on: Step 1, Step 2

**What changes:**
- One memo recording **two dependencies with two classes**, each with the measurement behind it: the MPRA peak-tile data as **code-closable** (present on `E:`, addressed at `D:`, and Verification Signal 2 is the proof), and the promcalc predictions as **data-absent** with Step 2's verdict deciding whether regeneration is a code path or an external one.
- States what is needed and from whom — the distinction that separates "write more code" from "a drive, a host, or an install must appear".
- Records the correction plainly: the substrate was asserted present in conversation on the strength of the wrong directory, and the raw MPRA data is not the input the chain consumes.
- Records that the frozen bar's `no uv` clause is **stale** while its `no D: data tree` clause is substantively right, so a future reader does not discard the whole note or trust all of it.
- **Quotes no number from the frozen bar as if measured here.** Its figures came from the two committed 2026-10-14 artifacts; this memo measures substrate reachability and nothing else.

**Test strategy:**
- Unit: the memo names both dependencies and assigns each a class from a closed vocabulary.
- Unit: the JSON sidecar is parseable — guarded because this repo committed an invalid hand-edited JSON artifact for a day and cited it as evidence.
- Unit: the memo contains none of the bar's headline figures, so it cannot be mistaken for a result.

### Step 4: Register the blocker with an expiry test
Files: soraya_runs/constraint-expiry.json
Depends on: Step 3

**What changes:**
- A 12th constraint via `constraint_expiry.register`, cid `glm-tile-chain-promcalc-substrate-absent`: what is given up (the frozen pre-registration cannot run, so the oracle-versus-selection fork stays unresolved), the rationale as written today, and an auto-checkable `file-exists` expiry.
- **The limitation is stated, not hidden.** A `file-exists` test on the configured `D:` path fires only if that drive letter returns; data arriving on another volume would not trip it, degrading this to a scheduled human check. That is the same honest limit the project already recorded for R9, where the compute that eventually arrived was cloud and a local-VRAM test would never have fired. Going from no expiry record at all to a dated one naming exactly what to re-open is the gain.
- A drive-independent form of the same check is written into the constraint's text as the complete test, since `test-exit-0` is not auto-run by `check_all`.

**Test strategy:**
- Unit: the store still parses and the pre-existing 11 constraints are field-identical afterwards — this is an append, and a diff is what proves it.
- Unit: the new constraint's `expiry_kind` is one of the module's four, and its `expiry_arg` is a path (not prose).

### Step 5: Correct the in-flight record
Files: NEXT.md, wiki/plans-index.md
Depends on: Step 3

**What changes:**
- `NEXT.md` gains the tile chain's real blocker, replacing nothing that is still true, under its own prune-don't-grow rule.
- Records that the `orthogonal-modality-stack-promcalc-plus-N-channels` arm stays `planned` in `soraya_runs/arms-ledger.json` and is **correctly** so — per R6 a declared arm cannot evaporate, and its blocker is substrate, not merit. **The arms ledger itself is not edited**: changing an arm's status to reflect a substrate problem would destroy the information that it is still owed.
- Names the one runnable high-VOI item that this diagnosis leaves standing — resuming Steps 4–6 of `plans/Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan/`, whose substrate (GLM weights on `E:/hf_cache`, MG1655 reference in `data/cache/refseq`) is verified present.

**Test strategy:**
- Unit: `NEXT.md` names the promcalc blocker and the arms ledger is byte-unchanged.
- Manual: read `NEXT.md` end to end for staleness, since its own rule is to prune rather than grow and a mechanical append defeats that.

## Execution Preview

Wave 0 (2 parallel):  Step 1 — route substrate through external_root, Step 2 — promcalc toolchain probe
Wave 1 (1):           Step 3 — classify and record the wall
Wave 2 (2 parallel):  Step 4 — register the expiry, Step 5 — correct the in-flight record

Critical path: Step 1 → Step 3 → Step 4 (3 waves)
Max parallelism: 2 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — Step 1 edits `dna_decode/glm/genomewide.py`, which `load_fragments` and `load_peak_tiles` share with the fragment-rung arms. A reachability check that is too eager would break a working path on the host that HAS `D:`. Mitigated by keeping the `configured`-status path bit-identical and pinning that with a non-vacuity test, but the regression would not show up on this host at all — it would only appear on the machine with the drive.
- Severity: medium — The whole plan's value rests on Verification Signal 2, a single run proving the MPRA half is code-closable. If that run fails for a third, unanticipated reason, Step 3's classification is unsupported and the memo must say so rather than classify on argument.
- Severity: medium — `external_root.probe()` scans every drive root. On a host with network or removable volumes that can be slow or raise; the module catches `OSError` per candidate, but the new call sites sit in a script's startup path where a hang would look like the script being broken. `[unverified]` on this host beyond the attached volumes.
- Severity: low — Step 4 writes a store owned by a skill outside this repo (`~/.claude/skills/soraya/scripts/constraint_expiry.py`). Its schema is not pinned by this repo's tests, so a skill-side schema change could invalidate the append. The parse-and-diff test bounds the damage to detection, not prevention.
- Severity: low — External tool integrations: `uv run` (verified in use all session), `pytest` (verified), `plan-toolkit.py` (verified by `--help`), `constraint_expiry.register` (signature read from source at `~/.claude/skills/soraya/scripts/constraint_expiry.py`, kinds `gpu-vram-gte`/`file-exists`/`test-exit-0`/`manual`). `python3` is a **non-functional** WindowsApps stub here — it prints a Store-install message and exits 0, so any script invoking `python3` fails open; the toolkit was run under `python`.
- Severity: low — Restructuring applied: the toolchain probe was split from the substrate routing so a probe failure cannot block the code-closable half, and the two record-correcting steps were separated from the memo that supplies their content so neither can cite a classification that does not exist yet.

## Open Questions

- **Is regenerating promcalc on this host worth doing at all, if Step 2 says it is possible?** It means installing a Promoter Calculator environment and recomputing per-tile predictions for 44,106 tiles across 10 seeds on a C:-tight host, and the frozen bar's anchor gate then demands the regenerated pipeline reproduce committed per-seed numbers within 0.01 on 9 of 10 seeds. That is a real risk of spending hours and landing on `REFUSED_ANCHOR_MISMATCH`. The alternative is waiting for the host that already has the data.
- **Should `allow_discovered=True` ever be used for the MPRA cache specifically?** The data on `E:` is the right file by name and size, and `CACHE_REQUIRED_SUBDIRS` exists to express exactly this kind of fitness check. Against: the module's founding argument is that a name match is not a contents match, and the one place this project adopted a discovered root by name is recorded as a hazard. Left as the operator passing `--mpra-cache` explicitly.
- Does the frozen bar's `execution_notes` deserve a correction in place, given its `no uv` clause is now false? It is a frozen pre-registration, and editing any part of one is a decision about what "frozen" means here.

## Verification

1. `uv run pytest tests/test_glm_tile_substrate.py tests/test_promcalc_toolchain_probe.py tests/test_glm_tile_predictions.py tests/test_glm_tile_selection.py -q` — green.
2. `uv run python scripts/glm_tile_predictions.py` with no arguments prints the `external_root` report and exits non-zero, with no `PermissionError` in the output.
3. `uv run python scripts/glm_tile_predictions.py --mpra-cache E:/PythonProjects/dna_decode_cache/mpra` gets past the tile load and fails only on the promcalc input — captured verbatim into the Step 3 memo.
4. `uv run python scripts/promcalc_toolchain_probe.py` prints one pre-committed verdict.
5. `soraya_runs/constraint-expiry.json` parses, holds 12 constraints, and the original 11 are field-identical.
6. Full suite reconciles against 14 failed / 6834 passed / 186 skipped with every moved count named; `git status` clean afterwards; frozen five byte-unchanged; prospective lock verifies at 2026-08-31.

## Save-time amendments

Captured at: 2026-10-21
Source: `/save-plan` arguments

Audit notes only. `/execute-plan` reads work from `## Implementation Steps`; nothing here is an
executable instruction. An amendment that changes a Step contract must go back through
`/technical-plan`.

- corrected my own VOI answer from the prior turn — I called the substrate present on the strength of the wrong directory, and the arm's real input (promcalc) is absent from every attached drive
- the blocker split into two dependencies with different wall classes once measured, so the plan became diagnosis-and-legibility rather than execution
- routing goes through the existing external_root rather than a second resolver, with allow_discovered deliberately left False
- regenerating promcalc, adopting a discovered root, and editing a frozen pre-registration all left as open questions rather than folded into scope

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high -->
