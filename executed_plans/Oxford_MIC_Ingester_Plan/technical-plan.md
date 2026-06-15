## Lens status
- **probe:** not applied (extends an already-probed arm; new surface is the ingestion layer).
- **brainstorm:** applied (pre-save, 2 rounds, grounded-verified) — 3 critical + 3 medium fixes folded below (exact-cohort manifest; crosswalk conflict taxonomy; run-scoped roll-up; operator-aware censoring; richer W0 probe; additive resolver method).
- **review:** not applied.

## Problem Statement
Build the INGESTION layer that turns a real cohort's MIC table into the `{BioSample: R/S}` labels the shipped+hardened external-revalidation arm consumes, for the Oxford pilot (`PRJNA604975`), and wire a one-command live run. SUCCESS = a single driver that re-validates the FROZEN decoder on Oxford, gated so the leakage/availability/powering verdicts cover EXACTLY the scored cohort and no stale/failed cell can publish. The decoder + the 5 FROZEN files (`amr_rules.py`, `mic_tiers.py`, `dna_decode/eval/cohort_manifest.py`, `build_validation_report_card.py`, `compute_lineage_metrics.py`) stay byte-unchanged. NOTE: the new per-run `cohort_manifest_external_<run_id>.json` ARTIFACT is unrelated to the frozen `cohort_manifest.py` MODULE.

## Codebase Context
- `dna_decode/data/external_mic_labels.py` — `build_drug_labels`/`write_labels`/`canonical_drug`/`parse_mic_token`/`tier_for_isolate`. parse_mic_token currently STRIPS the operator to a numeric bound (loses censor direction). Modified in Step 2 (operator-aware `MicValue`).
- `dna_decode/eval/biosample_resolver.py` — `runs_for_project` returns `(run, sample)` tuples (fields run_accession,sample_accession only). Step 1 ADDS `read_run_records_for_project(fields=...)` without touching the existing tuple contract/tests.
- `dna_decode/data/external_cohort_genomes.py` — `resolve_cohort_genomes` (crosswalk param + `ExternalKeyError` + `is_biosample_key`); the downstream dedup floor. Reused unchanged; the crosswalk's conflict detection (Step 4) sits UPSTREAM.
- `scripts/external_cohort_preflight.py` — `preflight(project, ...)` builds cohort BioSamples from `runs_for_project(project)` = the WHOLE project. Step 7 adds an EXACT-SET mode (reads the cohort manifest) + renames the project-wide path's artifact to `external_project_probe_*` (a coarse "worth ingesting?" diagnostic, NOT a scored gate).
- `scripts/external_cohort_revalidate.py` — scorer with `powering_gate` (exit 3 hard-fail / 1 degraded / 0). Writes the artifact before returning. Step 8 makes it REQUIRE the cohort manifest (drift guard) + stamp a `run_id` into the artifact filename.
- `scripts/build_external_validation_report.py` — `load_external_artifacts` blind-globs `external_validation_*.json`. Step 9 makes it refuse glob-all by default (require `--run-id`/`--artifacts`), skip `powering.hard_fail` cells, and publish degraded only with `--allow-degraded`.
- Pilot facts: Oxford `PRJNA604975` is a BROAD Gram-neg bacteraemia sequencing project (~3468 isolates) — the scored set is the MIC-tested E. coli SUBSET, hence the exact-cohort-manifest requirement. MIC-table schema + location UNVERIFIED (the W0 probe pins it).

### Reusable-Code Survey
- `external_mic_labels` (build_drug_labels/classify_tier path) — reused; Step 2 extends censoring, Step 5 feeds it.
- `biosample_resolver` parse/fetch pattern — reused; Step 1 adds a records accessor.
- `external_cohort_{preflight,revalidate}` + `build_external_validation_report` + `external_cohort_genomes` — reused; Steps 7/8/9 add manifest-scoping + run-scoping.
- `mic_tiers.classify_tier`/`breakpoints_for` (FROZEN) — read-only consumer.
- None new — searched: dna_decode/data/, dna_decode/eval/, scripts/, graphify-out/ (absent).

## Pre-Change Baseline
- Arm shipped+hardened: full suite 1161 passed (0 fail); 5 frozen files byte-unchanged; commits `0b1fba0..d0b4e12` local on main (unpushed).
- No ingestion layer; no live run; live preflight on `PRJNA604975` not yet run.

## Verification Signal
- W0 probe emits a schema/crosswalk-feasibility audit (cardinality, unique-key-by-candidate-field, dup MIC rows per isolate/drug, operator/censoring distribution per drug, MIC-key→BioSample resolution rate).
- `cohort_manifest_external_<run_id>.json` is the single handoff; preflight (exact-set), scorer, and roll-up all consume it; scorer drift-guards selected.tsv ⊆ manifest.
- Crosswalk hard-fails on 1-MIC-key→>1-BioSample + cross-field collisions, with field provenance; never silently collapses.
- Roll-up refuses glob-all by default, skips hard_fail cells, publishes degraded only with `--allow-degraded`.
- New tests green; existing 1161 unchanged; FROZEN 5 files byte-unchanged.

## Implementation Steps

### Step 1: Additive ENA read_run records accessor
Files: dna_decode/eval/biosample_resolver.py, tests/test_biosample_resolver.py
Depends on: none

**What changes:**
- Add `read_run_records_for_project(project, fields=("run_accession","sample_accession","sample_alias","secondary_sample_accession")) -> list[dict]` + a pure `parse_ena_read_run_records(tsv, fields)` parser. The existing `runs_for_project` `(run, sample)` tuple + `parse_ena_read_run` stay UNCHANGED (and their tests).

**Test strategy:**
- Unit (offline): parse multi-field TSV → list[dict]; missing optional column tolerated; existing tuple tests still pass.

### Step 2: Operator-aware MIC censoring
Files: dna_decode/data/external_mic_labels.py, tests/test_external_mic_labels.py
Depends on: none

**What changes:**
- Add `MicValue(value, operator, raw)` + `parse_mic_value(token) -> MicValue` (keeps `parse_mic_token` for back-compat). `tier_for_isolate` becomes operator-aware: a `>X`/`>=X` bound that itself meets the R tier → R; a `<X`/`<=X` bound that meets the S tier → S; otherwise the token is interval-censored → excluded. Add `CENSORED_EXCLUDED`/`CENSORED_HIGH_R`/`CENSORED_HIGH_S` to the bucket vocabulary so the censored disposition is visible.
- Prevents excluding a clearly-R `>2` cipro (was stripped to 2 → BORDERLINE). Existing non-censored behavior unchanged.

**Test strategy:**
- Unit: `>2` cipro → R (CENSORED_HIGH_R); `<=0.06` → S; mid-range `>1` → excluded (CENSORED_EXCLUDED); plain values unchanged; existing parse_mic_token tests still pass.

### Step 3: W0 empirical probe
Files: scripts/oxford_w0_probe.py, tests/test_oxford_w0_probe.py
Depends on: Step 1

**What changes:**
- New script. Given the project + a MIC-table path: fetch read_run records (Step 1) + read the MIC table; emit `wiki/oxford_w0_probe_<date>.json` with row cardinality, unique-key counts per candidate field (run/sample_alias/secondary/native), duplicate MIC rows per isolate/drug, operator/censoring distribution per drug, and the MIC-key→BioSample resolution rate. PURE summarize helpers split from the fetch/read for offline tests. This output PINS Step 4's key candidates + Step 5's column map.

**Test strategy:**
- Unit (offline fixtures): the summarize helpers over a fixture MIC table + fixture read_run records → expected cardinality/uniqueness/censoring/resolution counts.

### Step 4: Alias→BioSample crosswalk with conflict taxonomy
Files: dna_decode/data/external_crosswalk.py, tests/test_external_crosswalk.py
Depends on: Step 1

**What changes:**
- New module. `build_crosswalk(records) -> {native_key: BioSample}` over run/sample_alias/secondary/BioSample candidate fields. `resolve_keys(native_keys, records) -> {resolved, unresolved, conflicts}` with an explicit taxonomy: many runs→1 BioSample = OK; **1 MIC-native-key→>1 BioSample = HARD CONFLICT**; a candidate key colliding across fields to different BioSamples = **HARD CONFLICT**. Each conflict row carries `mic_key, candidate_field, candidate_value, resolved_biosample, source_row_id` (field provenance). Persist the crosswalk JSON.

**Test strategy:**
- Unit (fixtures): runs→one-BioSample OK; one-key→two-BioSamples = conflict; cross-field collision = conflict; unresolved reported; provenance fields present; persisted crosswalk round-trips.

### Step 5: Config-driven MIC-table ingester
Files: dna_decode/data/external_mic_ingest.py, tests/test_external_mic_ingest.py
Depends on: Step 2

**What changes:**
- New module. `ingest_mic_table(path, *, key_col, drug_cols, call_cols=None) -> {native_key: {canonical_drug: {"mics":[MicValue], "calls":set()}}}`. Config-driven column map (no hardcoded names; pinned from the W0 probe). Normalize drugs via `canonical_drug` (skip non-pilot); parse MIC cells via `parse_mic_value` (Step 2); collect categorical calls; emit an `unmapped_columns` report (loud, not silent). CSV/TSV tolerant.

**Test strategy:**
- Unit (fixtures): CSV+TSV; drug-column mapping + non-pilot skip; MicValue parsing incl censored; calls collected; unmapped-column report; missing key-col raises.

### Step 6: Label-emission driver + cohort manifest
Files: scripts/build_oxford_labels.py, tests/test_build_oxford_labels.py
Depends on: Step 4, Step 5

**What changes:**
- New script. Ingest (S5) → resolve native keys → BioSample via crosswalk (S4), ABORT on any HARD CONFLICT (no silent collapse) → drop+report unresolved → per pilot drug `build_drug_labels` → write `data/raw/oxford_extval_<drug>/selected_{strict,relaxed}.tsv` + `buckets_<drug>.json` (BioSample-keyed) AND the single `cohort_manifest_external_<run_id>.json` (rows: `mic_key, biosample, drug, tier, label, censor_meta, conflict_status`) — the exact scored-cohort definition consumed downstream. Emits a per-drug ingest summary.

**Test strategy:**
- Unit (fixtures + fake crosswalk): end-to-end → BioSample-keyed selected.tsv + manifest; hard-conflict aborts; unresolved reported not dropped-silently; keys pass `is_biosample_key`; manifest rows match labels.

### Step 7: Preflight exact-set mode + project-probe rename
Files: scripts/external_cohort_preflight.py, tests/test_external_cohort_preflight.py
Depends on: Step 6

**What changes:**
- Add `--cohort-manifest <json>` (EXACT-set mode): cohort BioSamples come from the manifest, so leakage + assembly-availability cover EXACTLY the scored set. The legacy `--project` path stays but its artifact is renamed `external_project_probe_<...>.json` + stamped `scored_gate: false` (a coarse diagnostic, never a scored-cohort gate). Verdict logic unchanged.

**Test strategy:**
- Unit: exact-set mode reads manifest BioSamples (not the whole project); project mode emits `external_project_probe_*` with `scored_gate=false`; existing preflight tests still pass.

### Step 8: Scorer requires the manifest + run-id scoping
Files: scripts/external_cohort_revalidate.py, tests/test_external_cohort_revalidate.py
Depends on: Step 6

**What changes:**
- Add `--cohort-manifest` (required for a live run) + `--run-id`. Drift guard: assert `selected.tsv` BioSamples ⊆ manifest BioSamples (raise on mismatch). Stamp `run_id` into the artifact + filename (`external_validation_<cohort>_<drug>_<run_id>_<date>.json`). powering_gate / exit codes unchanged.

**Test strategy:**
- Unit: drift guard raises on selected-vs-manifest mismatch; run_id in artifact + filename; existing powering/gate tests still pass.

### Step 9: Run-scoped roll-up
Files: scripts/build_external_validation_report.py, tests/test_build_external_validation_report.py
Depends on: Step 8

**What changes:**
- `load_external_artifacts` gains `--run-id`/`--artifacts <manifest>` scoping and REFUSES glob-all unless `--allow-unscoped-glob`. Skip any cell with `powering.hard_fail=true`; include `run_degraded` cells only with `--allow-degraded` (stamped, never headline). Render unchanged otherwise.

**Test strategy:**
- Unit: stale artifact beside a current failed artifact → only the scoped, non-hard-fail cells roll up; glob-all refused without the flag; degraded gated by --allow-degraded.

### Step 10: One-command live-run driver
Files: scripts/run_oxford_revalidation.py, tests/test_run_oxford_revalidation.py
Depends on: Step 6, Step 7, Step 8, Step 9

**What changes:**
- New orchestrator: W0 probe (advisory) → `build_oxford_labels` (S6) → `external_cohort_preflight --cohort-manifest` (abort != PASS unless `--allow-degraded`) → `external_cohort_revalidate` per drug (propagate exit 3/1) → roll-up ONLY IF every required drug exited 0 (or 1 under `--allow-degraded`) — driver gating is the PRIMARY invariant. Pure `plan_steps`/`worst_exit`/`roll_up_allowed` helpers split for offline tests; subprocess execution + run_id generation in `main`.

**Test strategy:**
- Unit: `worst_exit` precedence (3>1>2>0); `roll_up_allowed` blocks on any hard-fail/non-accepted-degraded drug; mocked dry-run aborts on non-PASS preflight; clean path returns 0. No network/Docker.

### Step 11: Docs + runbook
Files: CLAUDE.md, README.md, wiki/oxford_revalidation_runbook.md
Depends on: Step 10

**What changes:**
- Runbook (fetch table → W0 probe → driver → read report card) + the pinned MIC-table schema + gate/exit-code/run-id semantics + the cohort-manifest-vs-frozen-module naming caveat. CLAUDE.md gotcha + README command rows for the new scripts.

**Test strategy:**
- Docs only. Final full-suite run captured in Verification.

## Execution Preview
- Wave 0: Step 1, Step 2.
- Wave 1: Step 3 (S1), Step 4 (S1), Step 5 (S2).
- Wave 2: Step 6 (S4,S5).
- Wave 3: Step 7 (S6), Step 8 (S6).
- Wave 4: Step 9 (S8).
- Wave 5: Step 10 (S6,S7,S8,S9).
- Wave 6: Step 11 (S10).
- Total waves: 7. Max parallelism: 3 (Wave 1). Critical path: 1→4→6→8→9→10→11 (length 7).

## Risk Flags
- **MIC-table schema still unverified until the W0 probe runs against the real table** [unverified]: Steps 1-2-4 are schema-light; Step 5's column map + Step 6 MUST be pinned from the W0 probe output. Run the probe (+ live preflight) BEFORE executing Steps 5-6 for real — building the ingester against a guessed schema repeats the "build against real data" lesson.
- **Modifies 4 SHIPPED arm modules** [grounded]: external_mic_labels (S2), external_cohort_preflight (S7), external_cohort_revalidate (S8), build_external_validation_report (S9). All additive + their existing tests must stay green; the 5 FROZEN files are NOT among them.
- **`cohort_manifest_external_*` vs frozen `cohort_manifest.py`** [grounded]: same word, different things (per-run artifact vs the leakage-registry module). Documented in Step 11 to prevent confusion; do not route the artifact through the frozen module.
- **Preflight not yet run** [grounded]: if Oxford's scored BioSamples mostly lack free GCAs, the powering floor (≥10 R/≥10 S) won't clear — the exact-set preflight (S7) is the empirical gate.
- **sentrux absent** → architecture gate n/a.

## Open Questions
- `--allow-degraded` publishing: confirmed degraded cells roll up only with the flag (stamped, never headline); hard_fail never. Acceptable for the first Oxford clinical claim?
- Min strict per-class stays 10 (mirrors the project gate); W0 may justify a documented lower pilot threshold — decide at W0.
- Whether to pilot cipro-only first vs all of cipro/cef/gent (driver `--drugs`).

## Verification
- `uv run pytest tests/ -q --ignore=tests/test_models_foundation.py` → 0 regressions + new tests green.
- Manual (after fetching the real table): `scripts/run_oxford_revalidation.py --project PRJNA604975 --mic-table <path> --drugs ciprofloxacin` → W0 probe + exact-set preflight verdict + per-drug run-id-scoped artifact + run-scoped report card, with honest exit codes.
- Frozen invariant: git diff shows the 5 frozen files byte-unchanged.

## Save-time amendments

Audit-notes-only: this block is provenance for human readers. `/execute-plan` reads ONLY `## Implementation Steps`. If an amendment changes a Step's contract, re-run `/technical-plan` before `/execute-plan`.

Captured at: 2026-06-15
Source: `/save-plan` arguments

- W0 probe first
- operator-aware censoring
- crosswalk conflict taxonomy
- exact-cohort manifest
- run-scoped roll-up
- one-command driver

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified -->
