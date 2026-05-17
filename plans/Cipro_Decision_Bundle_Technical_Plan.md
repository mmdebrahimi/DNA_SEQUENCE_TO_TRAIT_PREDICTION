# Cipro Decision Bundle Plan — Technical Plan

> Implementation blueprint for the post-SUSPEND_CONDITION_4 decision bundle, scoped down per /review reductions: collapse Tier 0 Steps 1-2 into one script, drop Bakta + circular variants, defer mean+max preflight v3.

---

## Problem Statement

The post-SUSPEND_CONDITION_4 pipeline needs a tight Tier 0 implementation that produces a load-bearing verdict on whether Databricks burst spend is justified. The high-level plan at `plans/Cipro_Decision_Bundle_Plan.md` enumerates 8 steps across 4 tiers; /review surfaced 8 accepted reductions that collapse the workload onto code changes in 3 production files + 2 new script files + 1 contract fix.

The technical work is to:
1. Refactor existing code to expose reusable helpers (`_confidence_tier`, `gene_level_mutagenesis` aggregation kwarg, `load_bvbrc_ast` `method_filter=None`).
2. Ship the Tier 0 census + label-manifest combined script that consumes those helpers.
3. Wire a `--labels-from-manifest` flag into the Stage 1b runner so label-sensitivity reruns are code-enforced single-variant.
4. Defer Tier 1 (mean+max preflight v3) to a conditional follow-up plan; **drop** Bakta smoke + circular `label_mechanism_derived` variant entirely.

## Codebase Context

Verified findings from research:

- `scripts/cipro_mic_audit.py:54-84` defines `_confidence_tier(mics, calls)` — script-local helper that buckets MICs under CLSI + EUCAST simultaneously. Currently not importable from other modules.
- `dna_decode/data/ast_data.py:81-84,142` exposes `load_bvbrc_ast(path, organism=..., method_filter=DEFAULT_METHOD_FILTER)`. The check at line 142 (`if method_filter and method not in method_filter`) means an empty tuple `()` already bypasses filtering, but the contract is non-obvious. Eng review flagged this as a footgun for the census script which needs disk + etest rows.
- `dna_decode/interp/mutagenesis.py:95-149` (`gene_level_mutagenesis`) and `:170-216` (`saturation_mutagenesis`) both hardcode `aggregate_strain_features(matrix, "mean")` at lines 124, 142, 186, 204. `aggregate_strain_features` at `dna_decode/models/classifiers.py:56-86` already supports `"mean+max"` and `"max"` natively.
- `scripts/cipro_curated_baseline.py:80-90` defines `ABLATION_FEATURE_SETS` + emits `original_condition_4_verdict` + `amended_condition_4_verdict` at lines 402-407. No `given_suspended_gate` field yet.
- `scripts/stage1_n40_cipro.py:158-167` builds `labels_by_strain` directly from `s.ast_labels[drug_lower]`. CLI already supports `--aggregation` flag (line 521-526).
- `dna_decode/data/bvbrc_genome.py:72` defines `load_bvbrc_genome_metadata` — paired with `load_bvbrc_ast` it provides assembly_accession + MLST + contig_count + n50. `scripts/diagnose_bvbrc_mlst_gaps.py:26-103` already does AST × metadata × MLST layered counts at the full-BV-BRC universe level — its structure is the closest pattern for the census.
- Test infrastructure: `tests/test_data_ast.py`, `tests/test_interp_mutagenesis.py`, `tests/test_stage1_n40_cipro.py` exist; no test coverage for the cipro audit scripts (mic_audit, mechanism_audit, curated_baseline, mechanism_phenotype_merge). New test files needed.
- Prior-art lesson (HIGH salience, 2026-05-14): when /brainstorm or /review catches issues in a saved plan, edit the plan file in-place BEFORE /execute-plan rather than carrying corrections as a delta. This technical plan encodes the /review reductions into steps so they don't slip; **the source plan at `plans/Cipro_Decision_Bundle_Plan.md` still contains the un-reduced spec and should be updated separately** (flagged in Risk Flags).
- Prior-art lesson (HIGH salience, 2026-05-14): result label and decision label are separate fields. D4's relabeling estimand must remain "per-strain error concentration + rank-order stability" (the result), not "max AUROC" (a downstream interpretation).

## Pre-conditions (user-locked before execution)

Two user decisions per the /review synthesis MUST be locked before Wave 3 runs. They are NOT code-implementation steps but they gate the runtime interpretation:

- **PC1: D9 framing lock.** "Publish defensible negative/indeterminate" vs "ship working classifier." Affects whether Tier 3 cohort build is on the table.
- **PC2: D4 estimand — stored statistical context, NOT a bare fraction** (updated 2026-05-17 per /brainstorm Round 2). Use a Fisher exact label-stratified enrichment test with α=0.10 + enrichment ratio ≥ 1.25; report BOTH `strict_NOISY_only` AND `include_SUSPECT_S` noise-class definitions separately. The originally-proposed bare 60% threshold was BELOW the uniform-error null baseline (K/N = 26/38 = 68.4%) and would have biased the decision toward "label noise" regardless of evidence. See `wiki/cipro_decision_bundle_pre_conditions_2026-05-17.md` for the locked PC2 schema.

Steps 1-10 are not blocked on PC1/PC2; Step 11 + Step 12 are.

## Implementation Steps

### Step 0: Reconcile source plan in-place (blocking pre-Wave-0)
Files: plans/Cipro_Decision_Bundle_Plan.md
Depends on: none

**What changes:**
- Edit `plans/Cipro_Decision_Bundle_Plan.md` in place to absorb the /review reductions that this technical plan encodes. Per the 2026-05-14 HIGH-salience lesson, both plan files must be consistent before `/execute-plan` runs; otherwise worktree agents reading the source plan recreate the rejected scope.
- Specific edits:
  - Delete the source plan's Step 5 (Bakta 4-strain smoke) + Step 6 (mechanism completeness differential-test) from Tier 2 entirely.
  - Strip `label_mechanism_derived` + `label_suspect_s_as_r` + `label_exclude_no_mic` (standalone) from the D3 manifest variant list. Keep only `label_original` + `label_exclude_no_mic_borderline` (2 variants).
  - Strip the `--ignore-gate` CLI flag from D6 description. Replace with "JSON-field annotation only."
  - Re-class Tier 1 (Step 4 mean+max preflight v3) as "DEFERRED to a separate plan, conditional on the decision-cell `True_low_threshold` outcome from Step 11."
  - Delete D7 (Bakta selection) + D8 (mechanism completeness AMRFinder-first discipline) — both Tier 2 design decisions for dropped steps.
  - At the top of the source plan, immediately after the `>` summary line, add a `> **Executable scope:** see `plans/Cipro_Decision_Bundle_Technical_Plan.md`. This document is now a reconciled rationale.`
- Update `wiki/plans-index.md` entry for `Cipro_Decision_Bundle_Plan.md` if the key-decisions bullets reference dropped scope.

**Key details:**
- This step is a single-file edit + a small `wiki/plans-index.md` touch-up. Argued by the 2026-05-14 HIGH-salience lesson as the right discipline: edit the plan file in-place rather than carrying a delta document.
- After Step 0, the technical plan is the executable contract; the source plan is the rationale (D-decision narrative + open tradeoffs context that this technical plan doesn't repeat).

**Verification:**
- `grep -E "Bakta|label_mechanism_derived|label_suspect_s_as_r|--ignore-gate"` against `plans/Cipro_Decision_Bundle_Plan.md` returns 0 matches.
- `wiki/plans-index.md` Cipro_Decision_Bundle_Plan.md entry references only the surviving decisions.

### Step 1: Promote `_confidence_tier` helper + allow `method_filter=None` in ast_data
Files: dna_decode/data/ast_data.py, scripts/cipro_mic_audit.py
Depends on: Step 0

**What changes:**
- `dna_decode/data/ast_data.py` — add module-level `confidence_tier(mics: list[float], calls: list[str], breakpoints: dict | None = None) -> tuple[str, dict]`. Body is the function currently at `scripts/cipro_mic_audit.py:54-84`. Default `breakpoints` = `{"clsi_r": 2.0, "clsi_s": 0.5, "eucast_r": 1.0, "eucast_s": 0.25}`. Add module-level breakpoint constants (`CIPRO_CLSI_R`, etc.) for direct re-export.
- `dna_decode/data/ast_data.py:84,142` — change `method_filter: tuple[str, ...] = DEFAULT_METHOD_FILTER` to `method_filter: tuple[str, ...] | None = DEFAULT_METHOD_FILTER`; update line 142 check to `if method_filter is not None and method not in method_filter`. Empty tuple `()` still bypasses (backwards compatible); explicit `None` is now the documented "no filter" path.
- `scripts/cipro_mic_audit.py:54-84` — delete local `_confidence_tier`, import from `dna_decode.data.ast_data`. Callsite at line 116 updated.

**Key details:**
- Helper signature is policy-agnostic via the `breakpoints` dict. Callers can pass alternate CLSI / EUCAST / custom thresholds (used by the Step 7 census).
- `method_filter=None` semantics: documented in docstring; `()` retained as legal alias.
- No new public re-exports beyond `confidence_tier` + breakpoint constants.

**Test strategy:**
- Covered by Step 4.

### Step 2: `gene_level_mutagenesis` aggregation kwarg + `saturation_mutagenesis` shield
Files: dna_decode/interp/mutagenesis.py
Depends on: Step 0

**What changes:**
- `dna_decode/interp/mutagenesis.py:95` — add `aggregation: str = "mean"` kwarg to `gene_level_mutagenesis` signature.
- Lines 124, 142 — replace `aggregate_strain_features(full_matrix, "mean")` and `aggregate_strain_features(full_matrix[keep_mask], "mean")` with `aggregate_strain_features(..., aggregation)`.
- `dna_decode/interp/mutagenesis.py:170` — add `aggregation: str = "mean"` kwarg to `saturation_mutagenesis` signature.
- Lines 186, 204 — raise `NotImplementedError(f"saturation_mutagenesis does not yet support aggregation={aggregation!r}; only 'mean' is supported (Phase 2 work).")` when `aggregation != "mean"`. Keep current behavior unchanged for `aggregation == "mean"` (default).

**Key details:**
- Default = `"mean"` for both → all existing callers (preflight v1 / v2, viz, mutagenesis tests) preserved unchanged. Backwards-compatible refactor.
- `saturation_mutagenesis` shield is a latent-bug fence (eng-review finding): prevents silent mean-pool re-use if someone later passes `aggregation="mean+max"` while the per-position embedding loop still hardcodes mean.

**Test strategy:**
- Covered by Step 5.

### Step 3: Append `given_suspended_gate` field to curated baseline JSON
Files: scripts/cipro_curated_baseline.py
Depends on: Step 0

**What changes:**
- `scripts/cipro_curated_baseline.py` — append a new JSON-payload field `given_suspended_gate: "INFORMATIONAL_ONLY"` to the payload dict around line 405 (next to `amended_condition_4_verdict`). Single-line annotation only; **no CLI flag added** per /review (the curated baseline doesn't read the merge gate, so `--ignore-gate` would be a no-op).
- Update the markdown packet's "How to interpret" section to mention the field's meaning.

**Key details:**
- Field is a constant when SUSPEND_CONDITION_4 has fired upstream; the script itself doesn't dynamically read the merge gate state. Field's purpose is signaling to downstream readers, not branching control flow.
- Pure additive change to the JSON schema; existing consumers unaffected.

**Test strategy:**
- Covered by Step 6.

### Step 4: Tests for Step 1 (`confidence_tier` helper + `method_filter=None`)
Files: tests/test_data_ast.py
Depends on: Step 1

**What changes:**
- `tests/test_data_ast.py` — add tests:
  - `test_confidence_tier_high_r_extreme` — MIC list `[8.0, 16.0]` → tier `HIGH_R`.
  - `test_confidence_tier_high_s_extreme` — MIC list `[0.06, 0.125]` → tier `HIGH_S`.
  - `test_confidence_tier_borderline` — MIC list `[1.0]` → tier `BORDERLINE`.
  - `test_confidence_tier_ambiguous_clsi_eucast_disagree` — MIC `[0.75]` → CLSI=S, EUCAST=I/R; tier `AMBIGUOUS`.
  - `test_confidence_tier_no_mic` — empty mic list → `NO_MIC`.
  - `test_confidence_tier_conflict` — calls `["R", "S"]` → `CONFLICT`.
  - `test_confidence_tier_custom_breakpoints` — pass custom dict; verify thresholds applied.
  - `test_load_bvbrc_ast_method_filter_none_keeps_all_methods` — fixture with broth + disk + etest rows; `method_filter=None` keeps all; default keeps only broth.
  - `test_load_bvbrc_ast_method_filter_empty_tuple_equivalent_to_none` — backwards-compat regression guard.

### Step 5: Tests for Step 2 (`aggregation` kwarg + `saturation_mutagenesis` shield)
Files: tests/test_interp_mutagenesis.py
Depends on: Step 2

**What changes:**
- `tests/test_interp_mutagenesis.py` — add tests:
  - `test_gene_level_mutagenesis_default_aggregation_is_mean` — pin default; assert baseline_features.shape[1] == embedding_dim (NOT 2× for mean+max).
  - `test_gene_level_mutagenesis_mean_plus_max_aggregation` — pass `aggregation="mean+max"`; assert baseline_features.shape[1] == 2 * embedding_dim.
  - `test_gene_level_mutagenesis_unknown_aggregation_raises` — passes through to `aggregate_strain_features` which raises ValueError; assert propagation.
  - `test_saturation_mutagenesis_mean_pass_through_unchanged` — default mean call returns expected per-position deltas.
  - `test_saturation_mutagenesis_non_mean_raises_not_implemented` — `aggregation="mean+max"` raises `NotImplementedError` with the expected message substring `"Phase 2 work"`.

### Step 6: Tests for Step 3 (`given_suspended_gate` JSON field)
Files: tests/test_cipro_curated_baseline.py
Depends on: Step 3

**What changes:**
- NEW file `tests/test_cipro_curated_baseline.py`.
- 5 tests total — both the new `given_suspended_gate` field AND the existing 2-layer verdict logic at `scripts/cipro_curated_baseline.py:402-407` (pinning the gate-bearing behavior per /review C6):
  - `test_payload_emits_given_suspended_gate_field` — synthesize the smallest possible bundle (or mock the run_loso call) and assert the emitted JSON payload includes the new field with value `INFORMATIONAL_ONLY`.
  - `test_no_ignore_gate_cli_flag_present` — argparse introspection; assert `--ignore-gate` is NOT a registered argument (regression guard for the /review reduction).
  - `test_two_layer_verdict_original_pass_amended_fail` — synthesize a bundle where best all-feature AUROC ≥ 0.80 (ABSOLUTE_PASS) but no_POINT AUROC < `AMENDED_NO_POINT_GATE_AUROC` (0.773); assert `original_condition_4_verdict == "ABSOLUTE_PASS"` AND `amended_condition_4_verdict == "FAIL"` AND `pivot_trigger_condition_4_load_bearing == False`. POINT-dominated tautology case.
  - `test_two_layer_verdict_amended_no_point_pass` — synthesize where no_POINT AUROC ≥ 0.773; assert `amended_condition_4_verdict == "NO_POINT_PASS"` AND `pivot_trigger_condition_4_load_bearing == True`.
  - `test_two_layer_verdict_mechanism_only_pass` — synthesize where mechanism_only AUROC ≥ 0.80 but no_POINT < 0.773; assert `amended_condition_4_verdict == "MECHANISM_ONLY_PASS"`.

**Key details:**
- This is the first test file for `scripts/cipro_curated_baseline.py`; sets the test pattern for future cipro-script tests.
- Mock heavy dependencies (NT cache, refseq paths) — test only the JSON-emission + verdict-derivation path. The synthetic bundle can be a `FeatureBundle` dataclass with hand-crafted feature matrices that force known LOSO outcomes; alternatively, monkey-patch `run_loso` to return pre-baked `CVResult` objects.

### Step 7: New `scripts/cipro_feasibility_and_label_audit.py` (Tier 0 census + manifest)
Files: scripts/cipro_feasibility_and_label_audit.py
Depends on: Step 1

**What changes:**
- NEW file `scripts/cipro_feasibility_and_label_audit.py`.
- **CLI shape: argparse subparsers — two subcommands** (`census` + `manifest`), no shared mode flag. Pins the data-flow contract: census reads raw BV-BRC universe; manifest mode consumes the already-emitted mechanism-phenotype audit JSON and does NOT recompute MIC tiers.
  - `census`: required args `--ast-csv <path>`, `--metadata-csv <path>`. Optional: `--output-prefix <str>`, `--drug <str>` (default cipro for forward-compat with cef/tet). Output: `wiki/cipro_bvbrc_feasibility_census_<date>.{md,json}`. Census tiers each cipro row under 3 phenotype policies (HIGH extreme; CLSI-strict excluding I/borderline; EUCAST-strict excluding CLSI/EUCAST disagreement). Cross with assembly_accession non-null, MLST non-null, contig_count + N50 thresholds. Lineage-confounding sub-census: per policy, MLST count + largest-clade fraction + R/S separability by MLST alone. Method/source/testing-standard strata. Reports counts as a multi-tab JSON sidecar + markdown packet. **Output gate field:** `pass_path_a_gate: bool` — True iff any policy reaches ≥75R/≥75S decisive strains with MLST + downloadable assembly + assembly_quality pass.
  - `manifest`: required args `--cohort-parquet <path>`, `--mechanism-phenotype-json <path>`. Optional: `--output-prefix <str>`, `--drug <str>`. Output: `data/processed/cipro_label_manifest_<date>.parquet` with columns `strain_id`, `accession`, `mlst`, `original_label`, `mic_tier`, `mic_median`, `has_primary_cipro_mechanism`, `noise_class`, `label_original`, `label_exclude_no_mic_borderline`, `inclusion_flag_original`, `inclusion_flag_exclude_no_mic_borderline`. **Only 2 variants** per /review reduction; `label_mechanism_derived` + `label_suspect_s_as_r` deliberately omitted (circular / post-hoc-risky).
- Imports `confidence_tier` from `dna_decode.data.ast_data` (Step 1 helper) — census mode only.
- Imports `load_bvbrc_ast(method_filter=None)` from `dna_decode.data.ast_data` to keep disk + etest strata — census mode only.
- Imports `load_bvbrc_genome_metadata` from `dna_decode.data.bvbrc_genome` — census mode only.
- Manifest mode reads `--mechanism-phenotype-json` (today's `wiki/cipro_mechanism_phenotype_audit_2026-05-17.json`) for per-strain `noise_class` + `has_primary_cipro_mechanism` flags; reads `--cohort-parquet` for the strain set. **Does not touch the raw AST CSV.**

**Key details:**
- Reuse `_filter_by_assembly_quality` from `scripts/build_stage2_n150_cohort.py:75-85` rather than reimplementing.
- Reuse the layered counting pattern from `scripts/diagnose_bvbrc_mlst_gaps.py:26-103`.
- ≤350 LOC target across both subcommands.
- Manifest's `strain_id` set MUST equal the effective cohort exactly (1:1). The `inclusion_flag_<variant>` boolean is the only mechanism for filtering — missing rows are a data-integrity failure, not implicit exclusion (per /review M1).
- Error handling: raise on missing input files; warn-and-skip on rows with unparseable MIC; raise on missing required argparse args (subparsers enforce).

**Test strategy:**
- Covered by Step 9.

### Step 8: `stage1_n40_cipro.py` — `--labels-from-manifest --variant <col>` wiring + JSON sidecar
Files: scripts/stage1_n40_cipro.py
Depends on: Step 0

**What changes:**
- `scripts/stage1_n40_cipro.py:158-167` (`load_features`) — accept new args `labels_from_manifest: Path | None = None`, `variant_column: str | None = None`. When both provided:
  1. Load manifest parquet; raise informative error if `strain_id` set ≠ effective cohort strain_id set (strict 1:1; per /review M1).
  2. Override `labels_by_strain[s.strain_id]` from the manifest's `<variant_column>` value.
  3. Filter strains by `inclusion_flag_<variant_column>` column (rows where flag is False are excluded from LOSO).
  4. Raise informative error if after filtering, either class has fewer than 2 strains.
- CLI: add `--labels-from-manifest <path>` + `--variant <column-name>`. Both required-together; argparse should enforce. Allow only ONE `--variant` per invocation (D4 discipline — single variant per run, structurally enforced).
- **Refactor:** extract `build_parser() -> ArgumentParser` helper so tests can introspect arguments without invoking main. Drive-by cleanup of the existing parser-duplication pattern at `tests/test_stage1_n40_cipro.py:465-476`.
- Output filename: append `_<variant>` suffix to packet filename when variant is provided.
- **NEW JSON sidecar (per /review C4):** alongside the existing markdown packet, emit `wiki/stage1_n40_cipro_<aggregation>_<date>[_<variant>].predictions.json` with per-strain CVResult.folds data. Schema:
  ```
  {
    "drug": str,
    "aggregation": str,
    "variant": str | None,
    "strain_ids": [str],
    "per_strain": [
      {
        "strain_id": str,
        "mlst": str,
        "y_true": int,
        "nt_lr_score": float,
        "nt_xgb_score": float,
        "kmer_xgb_score": float,
        "fusion_score": float | null,
        "noise_class": str | null  // null when no manifest passed
      }
    ]
  }
  ```
- The sidecar consumes existing `CVResult.folds[i].y_true[0]` + `y_score[0]` already in memory at packet-write time. No new computation; pure additive serialization.

**Key details:**
- Pure additive change; default behavior (no flags) unchanged.
- Manifest must have `strain_id` + `<variant_column>` + `inclusion_flag_<variant_column>` columns; raise informative error on missing.
- Aggregation default still `mean+max` (from Stage 1b lock).
- The JSON sidecar is the upstream producer for the runtime decision artifact (Step 11). Worktree agents implementing this MUST emit the JSON in the same `write_packet` block so the sidecar lands atomically with the markdown.

**Test strategy:**
- Covered by Step 10.

### Step 9: Tests for Step 7 (census + manifest script)
Files: tests/test_cipro_feasibility_audit.py
Depends on: Step 7

**What changes:**
- NEW file `tests/test_cipro_feasibility_audit.py`.
- Tests:
  - `test_census_emits_three_policy_tables` — small synthetic AST CSV (10 strains across MIC range); assert JSON sidecar has 3 policy keys (`high_extreme` / `clsi_strict` / `eucast_strict`).
  - `test_census_high_r_count_matches_existing_mic_audit_for_n38_subset` — sanity regression: filter synthetic to the 7 known HIGH_R strain_ids, assert count matches today's audit.
  - `test_census_method_strata_includes_disk_and_etest` — fixture with disk + etest rows; assert both stratum counts non-zero (regression guard for the `method_filter=None` plumbing).
  - `test_census_pass_path_a_gate_false_on_undersized_cohort` — synthetic with 10 HIGH_R / 5 HIGH_S → `pass_path_a_gate=False`.
  - `test_census_pass_path_a_gate_true_on_sized_cohort` — synthetic with 80 HIGH_R / 80 HIGH_S → True.
  - `test_manifest_emits_exactly_two_label_variants` — assert columns include `label_original` + `label_exclude_no_mic_borderline` + corresponding inclusion flags; assert `label_mechanism_derived` + `label_suspect_s_as_r` NOT present (regression guard for the /review reduction).
  - `test_manifest_inclusion_flag_filters_no_mic_and_borderline` — strain with mic_tier in {NO_MIC, BORDERLINE, AMBIGUOUS} → `inclusion_flag_exclude_no_mic_borderline = False`.
  - `test_manifest_row_count_matches_cohort_size` — N=38 → 38 rows.

**Key details:**
- Use a small synthetic AST CSV fixture under `tests/fixtures/` if not already present.
- Don't load the real 2.5M-row BV-BRC CSV in tests.

### Step 10: Tests for Step 8 (manifest-driven stage1_n40_cipro + JSON sidecar)
Files: tests/test_stage1_n40_cipro.py
Depends on: Step 8

**What changes:**
- `tests/test_stage1_n40_cipro.py` — refactor + add tests.
- **Refactor:** existing local-parser tests at `tests/test_stage1_n40_cipro.py:465-476` should call `scripts.stage1_n40_cipro.build_parser()` (introduced in Step 8) instead of building a local `argparse.ArgumentParser()`. Drive-by cleanup of the duplicate-setup pattern per /review M2.
- New tests (use `stage1_n40_cipro.main(argv=[...])` or `build_parser().parse_args([...])`, NOT local parser):
  - `test_load_features_default_uses_cohort_labels` — regression guard: with no manifest flag, `labels_by_strain` derived from cohort.
  - `test_load_features_with_manifest_overrides_labels` — pass synthetic manifest with one variant flipped; assert `labels_by_strain` reflects the override.
  - `test_load_features_with_manifest_filters_inclusion_flag_false_strains` — assert excluded strains do not appear in returned strain_ids.
  - `test_cli_requires_variant_when_manifest_passed` — argparse rejects `--labels-from-manifest` without `--variant`.
  - `test_cli_rejects_unknown_variant_column` — synthetic manifest with only `label_original`; `--variant label_exclude_no_mic_borderline` raises informative error.
  - `test_output_filename_includes_variant_suffix` — packet path ends with `_label_exclude_no_mic_borderline.md`.
  - **Manifest schema enforcement (per /review M1):**
    - `test_load_features_rejects_manifest_with_missing_strain` — manifest's `strain_id` set is a strict subset of the effective cohort; assert raises informative error referencing the missing strain_id.
    - `test_load_features_rejects_manifest_with_extra_strain` — manifest contains a strain_id not in the effective cohort; assert raises informative error.
  - **Minimum-class-count enforcement (per /review C4/C6):**
    - `test_load_features_raises_when_variant_filter_leaves_fewer_than_two_strains_per_class` — synthetic manifest whose `inclusion_flag_<variant>` leaves 1 R and 5 S strains; assert raises with "fewer than 2 strains" substring before LOSO starts.
  - **Alignment preservation:**
    - `test_load_features_alignment_preserved_after_variant_filter` — strain_ids returned in same relative order regardless of which strains were filtered out; assert NT-vs-k-mer-vs-fusion alignment property holds (per the existing `compute_gate_outcome` strain_ids validation).
  - **JSON sidecar (Step 8's new artifact):**
    - `test_predictions_json_sidecar_emitted` — synthetic small run; assert `wiki/stage1_n40_cipro_*.predictions.json` exists alongside the markdown packet.
    - `test_predictions_json_schema_pins_per_strain_columns` — load the sidecar; assert each per_strain entry has `strain_id`, `mlst`, `y_true`, `nt_lr_score`, `nt_xgb_score`, `kmer_xgb_score`, `fusion_score`, `noise_class` keys.
    - `test_predictions_json_noise_class_null_when_no_manifest` — no manifest passed; assert all per_strain entries have `noise_class: None`.

### Step 10.5: NEW `scripts/cipro_error_audit.py` (gate-bearing producer)
Files: scripts/cipro_error_audit.py
Depends on: Step 8, Step 10

**What changes:**
- NEW file `scripts/cipro_error_audit.py`.
- Inputs (CLI): `--predictions-json <path>` (Step 8's per-strain JSON sidecar), `--manifest <path>` (Step 7's manifest), `--pre-conditions-md <path>` (the user-locked PC1/PC2 markdown).
- Logic: load the predictions JSON; for each strain, compute `error = abs(y_true - nt_lr_score >= 0.5)` (binary error). For each PC2 noise-class definition (`strict_NOISY_only`, `include_SUSPECT_S`), run a label-stratified Fisher exact (or Freeman-Halton) test conditional on the observed true-label distribution; compute observed_fraction, null_fraction, p_value, enrichment_ratio.
- Output: `wiki/cipro_per_strain_error_audit_<date>.{md,json}`. JSON schema:
  ```json
  {
    "predictions_path": str,
    "manifest_path": str,
    "pre_conditions_artifact_path": str,
    "pc2_test": "fisher_exact_label_stratified",
    "pc2_alpha": float,
    "pc2_enrichment_ratio_min": float,
    "total_errors": int,
    "results_by_definition": {
      "strict_NOISY_only": {
        "observed_fraction": float,
        "null_fraction": float,
        "p_value": float,
        "enrichment_ratio": float,
        "label_noise_concentration_established": bool
      },
      "include_SUSPECT_S": {
        "observed_fraction": float,
        "null_fraction": float,
        "p_value": float,
        "enrichment_ratio": float,
        "label_noise_concentration_established": bool
      }
    },
    "error_cluster_bucket": "established_under_at_least_one_definition" | "not_established",
    "per_strain": [
      {"strain_id": str, "error": bool, "noise_class": str}
    ]
  }
  ```
- Decision rule: `label_noise_concentration_established == True` iff (p_value < pc2_alpha) AND (enrichment_ratio >= pc2_enrichment_ratio_min). `error_cluster_bucket = "established_under_at_least_one_definition"` iff any definition's flag is True.
- Refuses to run if `--pre-conditions-md` file is missing (structural enforcement of PC1/PC2 lock).

**Test strategy:**
- Covered by Step 10.6 (new test file).

### Step 10.6: Tests for Step 10.5 (`cipro_error_audit.py`)
Files: tests/test_cipro_error_audit.py
Depends on: Step 10.5

**What changes:**
- NEW file `tests/test_cipro_error_audit.py`.
- Tests:
  - `test_fisher_exact_rejects_when_enrichment_clear` — synthesize predictions with high concentration on noisy strains (e.g., 14/15 strict-NOISY errors out of 25 noisy strains in N=38); assert p_value < 0.10 AND enrichment_ratio >= 1.25 AND `label_noise_concentration_established == True`.
  - `test_fisher_exact_fails_to_reject_at_null_baseline` — synthesize predictions with uniform errors (matching K/N null); assert `label_noise_concentration_established == False`.
  - `test_both_noise_class_definitions_reported_separately` — assert output JSON has both `strict_NOISY_only` AND `include_SUSPECT_S` keys with independent results.
  - `test_zero_total_errors_handled` — predictions all correct; assert no division-by-zero crash; emit `total_errors: 0` + appropriate null-result handling.
  - `test_refuses_without_pre_conditions_artifact` — missing pre-conditions file → RuntimeError with substring "PC1/PC2 pre-conditions".

### Step 10.7: NEW `scripts/cipro_decision_cell.py` (gate-bearing wrapper)
Files: scripts/cipro_decision_cell.py
Depends on: Step 10.5

**What changes:**
- NEW file `scripts/cipro_decision_cell.py`.
- Inputs (CLI): `--census-json <path>`, `--error-audit-json <path>`, `--pc1-framing <internal_closeout|publish|ship>`, `--pre-conditions-md <path>`. (Note: `--pc2-threshold` arg DROPPED — PC2 is now a structured statistical context, not a single float; the error-audit JSON already encodes the per-definition decision flags per the locked PC2 schema.)
- Logic:
  1. Verify `--pre-conditions-md` file exists. If absent, raise `RuntimeError("PC1/PC2 pre-conditions artifact not found at <path>; declare PC1 + PC2 in wiki/cipro_decision_bundle_pre_conditions_<date>.md before running this script.")`. Structural enforcement of PC1/PC2.
  2. Load `census-json` → extract `pass_path_a_gate: bool`.
  3. Load `error-audit-json` → extract `error_cluster_bucket: str` (= `"established_under_at_least_one_definition"` or `"not_established"`).
  4. Derive `decision_cell` from the 2x2 (census × error-audit-bucket):
     - `(True, established)` → `"True_high_threshold"`
     - `(True, not_established)` → `"True_low_threshold"`
     - `(False, established)` → `"False_high_threshold"`
     - `(False, not_established)` → `"False_low_threshold"`
  5. Compute `recommended_next_step` from the decision matrix conditioned on PC1 framing. Under `pc1_framing=internal_closeout`, no cell unlocks Databricks burst; the matrix only influences narrative tone for the closeout packet.
- Output: `wiki/cipro_decision_bundle_runtime_<date>.json` with schema:
  ```json
  {
    "pc1_framing": "internal_closeout" | "publish" | "ship",
    "pc2_test": "fisher_exact_label_stratified",
    "pc2_alpha": float,
    "pc2_enrichment_ratio_min": float,
    "pre_conditions_artifact_path": str,
    "census_pass_path_a_gate": bool,
    "error_audit_results_by_definition": {...},
    "error_cluster_bucket": str,
    "decision_cell": str,
    "recommended_next_step": str
  }
  ```
- **Step 12 reads BOTH `decision_cell` AND `pc1_framing`** (per /review C5). `decision_cell` does NOT encode PC1; Step 12 explicitly gates on both. Under `pc1_framing=internal_closeout`, Step 12 is a no-op (no relabel-LOSO authorized).

**Test strategy:**
- Covered by Step 10.8 (new test file).

### Step 10.8: Tests for Step 10.7 (`cipro_decision_cell.py`)
Files: tests/test_cipro_decision_cell.py
Depends on: Step 10.7

**What changes:**
- NEW file `tests/test_cipro_decision_cell.py`.
- Tests for each of the 4+ cells:
  - `test_decision_cell_true_high_threshold` — census pass + error bucket above → decision_cell = "True_high_threshold".
  - `test_decision_cell_true_low_threshold` — census pass + error bucket below → "True_low_threshold".
  - `test_decision_cell_false_high_threshold` — census fail + error bucket above → "False_high_threshold".
  - `test_decision_cell_false_low_threshold` — census fail + error bucket below → "False_low_threshold".
  - `test_decision_cell_ambiguous_band_suspends` — ambiguous_in_band → "AMBIGUOUS_*" prefix.
- Pre-conditions enforcement:
  - `test_raises_if_pre_conditions_artifact_missing` — `--pre-conditions-md` points to nonexistent file → RuntimeError with expected substring "PC1/PC2 pre-conditions".

### Step 11: Runtime — execute census + manifest + error audit + decision cell (manual; sequential)
Files: (runtime only — no source-tree changes)
Depends on: Step 7, Step 9, Step 10.5, Step 10.6, Step 10.7, Step 10.8, PC1, PC2

**What runs (in order):**
1. **Pre-conditions artifact:** user writes `wiki/cipro_decision_bundle_pre_conditions_<date>.md` declaring PC1 (publish vs ship) + PC2 (numeric threshold + ambiguous band if any). MUST exist before subsequent commands.
2. **Census:** `uv run python scripts/cipro_feasibility_and_label_audit.py census --ast-csv "C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv" --metadata-csv "C:/Users/Farshad/Downloads/BVBRC_genome (1).csv"` — emits `wiki/cipro_bvbrc_feasibility_census_<date>.{md,json}`.
3. **Manifest:** `uv run python scripts/cipro_feasibility_and_label_audit.py manifest --cohort-parquet data/processed/gate_b_n40_cipro_cohort.parquet --mechanism-phenotype-json wiki/cipro_mechanism_phenotype_audit_2026-05-17.json` — emits `data/processed/cipro_label_manifest_<date>.parquet`.
4. **Per-strain error audit:** `uv run python scripts/cipro_error_audit.py --predictions-json wiki/stage1_n40_cipro_mean-plus-max_2026-05-16.predictions.json --manifest data/processed/cipro_label_manifest_<date>.parquet --pc2-threshold <float>` — emits `wiki/cipro_per_strain_error_audit_<date>.{md,json}`.
   - **Pre-req:** Step 8's JSON sidecar must exist. If Stage 1b ran before Step 8 shipped, re-run Stage 1b under the patched script to produce the sidecar.
5. **Decision cell:** `uv run python scripts/cipro_decision_cell.py --census-json wiki/cipro_bvbrc_feasibility_census_<date>.json --error-audit-json wiki/cipro_per_strain_error_audit_<date>.json --pc1-framing <publish|ship> --pc2-threshold <float> --pre-conditions-md wiki/cipro_decision_bundle_pre_conditions_<date>.md` — emits `wiki/cipro_decision_bundle_runtime_<date>.json`.

**Decision matrix (computed by `cipro_decision_cell.py`; do NOT re-interpret by hand):**

| `decision_cell` | `pc1_framing` | recommended_next_step |
|---|---|---|
| True_high_threshold | publish | Step 12 (relabel-LOSO) to confirm noise hypothesis; close out as "publish defensible negative." |
| True_high_threshold | ship | Step 12 first, then Path A burst (label-noise rescue + cohort expansion). |
| True_low_threshold | publish | Tier 1 mean+max preflight v3 closeout; publish defensible architecture-bottleneck narrative. |
| True_low_threshold | ship | Tier 1 closeout, then Path A burst. |
| False_high_threshold | * | Cipro irreparable at N=38; **pivot to cef/tet smoke** (`plans/EP2_Cef_Tet_Smoke_Design_Plan.md`). |
| False_low_threshold | * | Cipro both noisy + NT-untestable; abandon cipro for Phase 1 EP1; close out as "publish defensible negative." |
| AMBIGUOUS_* | * | Re-lock PC2 numeric threshold; re-run Step 11. |

### Step 12: Runtime — conditional pre-declared label-sensitivity LOSO (manual)
Files: (runtime only)
Depends on: Step 8, Step 10, Step 11

**What runs:**
- Read `wiki/cipro_decision_bundle_runtime_<date>.json`. Step 12 fires IFF:
  - `decision_cell == "True_high_threshold"` AND
  - `pc1_framing in {"publish", "ship"}` (both framings authorize Step 12 in this cell; the difference is the downstream Path A burst).
- Command: `HF_HOME=D:/hf_cache uv run python scripts/stage1_n40_cipro.py --labels-from-manifest data/processed/cipro_label_manifest_<date>.parquet --variant label_exclude_no_mic_borderline`
- Primary estimand reported (NOT max AUROC): per-strain prediction direction stability vs original. Compare against PC2 threshold (re-using the same numeric bar).

## Execution Preview

```
Wave 0 (1 sequential, blocking):
  Step 0 — reconcile source plan in-place

Wave 1 (3 parallel):
  Step 1 — ast_data helper + filter (+ cipro_mic_audit consumer)
  Step 2 — mutagenesis aggregation + shield
  Step 3 — curated_baseline given_suspended_gate field

Wave 2 (5 parallel):
  Step 4 — tests for Step 1
  Step 5 — tests for Step 2
  Step 6 — tests for Step 3 (5 tests pinning 2-layer verdict)
  Step 7 — feasibility_and_label_audit.py (census + manifest subcommands)
  Step 8 — stage1_n40_cipro --labels-from-manifest + JSON sidecar

Wave 3 (2 parallel):
  Step 9  — tests for Step 7
  Step 10 — tests for Step 8 (uses main(argv); manifest-mismatch rejection)

Wave 4 (2 parallel):
  Step 10.5 — scripts/cipro_error_audit.py (consumes Step 8's JSON sidecar)
  Step 10.7 — scripts/cipro_decision_cell.py (gate-bearing wrapper)

Wave 5 (2 parallel):
  Step 10.6 — tests for Step 10.5
  Step 10.8 — tests for Step 10.7

Wave 6 (sequential, manual): Step 11 — census + manifest + error audit + decision cell
Wave 7 (sequential, manual, conditional): Step 12 — label-sensitivity LOSO

Critical path: Step 0 → Step 1 → Step 7 → Step 9 → Step 10.5 → Step 10.6 → Step 11 → Step 12 (8 waves; 6 code + 2 runtime)
Max parallelism: 5 agents (Wave 2)
```

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- **Source plan drift resolved structurally via Step 0.** Step 0 is a blocking pre-Wave-1 step that edits `plans/Cipro_Decision_Bundle_Plan.md` in place to match this technical plan's reduced scope. Per the 2026-05-14 HIGH-salience lesson, this is the right discipline: corrections land in the plan file, not as a delta document. All later steps `Depends on: Step 0`.
- **Tier 1 (mean+max preflight v3) deferred but not deleted.** The `gene_level_mutagenesis` aggregation kwarg ships in Step 2 so the refactor is ready. The runtime use is gated on Step 11's decision matrix (only fires in `True_low_threshold` cell). Document in CLAUDE.md gotchas so it doesn't get lost.
- **Pre-conditions PC1 + PC2 are NOW structurally enforced.** `scripts/cipro_decision_cell.py` (Step 10.7) refuses to write the runtime artifact unless `wiki/cipro_decision_bundle_pre_conditions_<date>.md` exists. Step 12 reads BOTH `decision_cell` AND `pc1_framing` from the runtime JSON; the matrix is machine-checked, not prose-interpreted.
- **Step 7's census script has no end-to-end smoke test in this plan.** Step 9 covers synthetic-fixture tests; real-data census runs in Step 11. If Step 11 produces unexpected counts (e.g., HIGH_R count diverges from today's audit), the cause may be silent parser drift on the raw CSV. Recommend running `scripts/cipro_mic_audit.py` on the N=38 cohort AS A REGRESSION CHECK after Step 7 lands, before firing Step 11.
- **Step 1 modifies TWO files** (`dna_decode/data/ast_data.py` + `scripts/cipro_mic_audit.py`); both listed in the `Files:` metadata. Step 1 verification must include a manual run of `cipro_mic_audit.py` to confirm the import + callsite refactor holds.
- **Step 8 ships the JSON sidecar that gates Step 10.5 + downstream.** If Stage 1b's existing markdown packet at `wiki/stage1_n40_cipro_mean-plus-max_2026-05-16.md` was emitted BEFORE Step 8 lands (it was — landed 2026-05-16), Step 10.5 has no sidecar to read. Mitigation: Step 11 runtime block 4 includes "re-run Stage 1b under the patched script to produce the sidecar" as an explicit pre-req for the error-audit step.
- **Restructuring applied.** Source plan listed Step 1 (census) + Step 2 (manifest) as separate Tier 0 entries; this technical plan merges into Step 7 (single script with `census` + `manifest` subcommands) per /review's accepted reduction. Source plan's Step 3 = this plan's Step 12 (renamed to clarify it's the conditional runtime, not a fresh coding step). Step 11's "small notebook / inline script" prose replaced with three concrete scripts (Step 10.5 + Step 10.7) + per-strain JSON sidecar (Step 8 extension).

## Verification

- **Wave 0 (Step 0):** `grep -E "Bakta|label_mechanism_derived|label_suspect_s_as_r|--ignore-gate"` against `plans/Cipro_Decision_Bundle_Plan.md` returns 0 matches. Source plan has the supersedes banner at the top.
- **Wave 1 + Wave 2 + Wave 3:** `uv run pytest tests/ -v` passes with new tests in Steps 4, 5, 6, 9, 10. Existing test count + new test count sum cleanly (no regression).
- **Step 1 verification:** `uv run python scripts/cipro_mic_audit.py` re-runs against the existing N=38 cohort and produces a packet identical (modulo timestamp) to `wiki/cipro_mic_audit_2026-05-17.md`. Tests at `tests/test_data_ast.py` pass with the new `confidence_tier` + `method_filter=None` cases.
- **Step 2 verification:** `uv run python -u scripts/cipro_attribution_preflight.py` (v1/v2 mean-pool path) re-runs and produces identical INCONCLUSIVE_MISS verdict. Backwards-compat preserved. Tests at `tests/test_interp_mutagenesis.py` pass.
- **Step 3 verification:** Re-running `cipro_curated_baseline.py` (with synthetic input) emits a JSON payload containing the new `given_suspended_gate: "INFORMATIONAL_ONLY"` field. Markdown packet's "How to interpret" mentions the field.
- **Step 7 verification:** `census` subcommand runs against the real BV-BRC CSV in under 5 min wallclock. Output sidecar has 3 policy tables. `manifest` subcommand produces a parquet with 38 rows × required columns. Manual sanity check: HIGH_R count under `high_extreme` policy matches today's MIC audit (7).
- **Step 8 verification:** Running `stage1_n40_cipro.py` with no manifest flags produces identical output to Stage 1b's existing verdict packet (modulo timestamp). The new `.predictions.json` sidecar exists alongside the markdown. Running with `--labels-from-manifest <synthetic.parquet> --variant label_original` produces same verdict. Running with `--variant label_exclude_no_mic_borderline` produces a smaller cohort + a `_label_exclude_no_mic_borderline.md` packet + matching `.predictions.json`.
- **Step 10.5 / 10.7 verification:** `cipro_error_audit.py` emits the `wiki/cipro_per_strain_error_audit_<date>.json` with `error_cluster_bucket` field; `cipro_decision_cell.py` refuses to run without the PC artifact and emits the runtime JSON with all required fields populated.
- **End-to-end gate (after Step 11):** `wiki/cipro_decision_bundle_runtime_<date>.json` contains `decision_cell` + `pc1_framing` + `recommended_next_step`. Step 12's runtime check is mechanical (`if decision_cell == "True_high_threshold"`). The user can answer "should Tier 3 (Path A burst) fire?" with grounded yes/no.
