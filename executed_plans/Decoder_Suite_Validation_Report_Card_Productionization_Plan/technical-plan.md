## Lens status
- **probe:** applied — Anchor-4 probe (this session) surfaced the 4 issues this plan addresses (over-claim non-cells, unknown powering, fragile leakage exclusion, aggregation-tier erosion).
- **brainstorm:** applied — pre-exec `/brainstorm` ran (class-d cross-cutting: shared `wiki/` namespace + multi-producer artifact). Findings folded into the Save-time amendments below.
- **review:** not applied.

## Problem Statement
The Anchor-4 v0 report card (`scripts/build_validation_report_card.py`, shipped this session) proved the standing-capability shape but has four productionization gaps the probe named: (1) leakage exclusion in `scripts/provenance_disjoint_validate.py` hardcodes only 3 of the 8 parquet cohorts in `data/processed/` — a maintenance trap that silently under-excludes as the suite grows; (2) the census (`scripts/ncbi_pd_provenance_census.py`) prints powering verdicts to stdout only, so `wiki/provdisjoint_census_results.json` had to be hand-written; (3) the report card's row set is the union of observed scored/censused/registry keys, so a shipped decoder that has never been censused (e.g. E. coli cef/tet/gent/meropenem/oxacillin) silently does not appear as a `NOT_CENSUSED` cell — the grid looks more complete than it is; (4) the shipped-decoder enumeration is partly hand-maintained. Goal: make leakage exclusion data-driven (an accession-manifest registry), make census powering self-persisting, and make the report-card grid auto-derive the full shipped-decoder surface — preserving the 6-state cell machine, the honest per-cell tier, and the no-aggregate-headline rule.

## Codebase Context
- `scripts/provenance_disjoint_validate.py` — Stage-2 scorer. `_prior_cohort_accessions(slug, exclude_self)` (lines 46-68) globs `data/raw/<slug>_*/selected.tsv` + a HARDCODED `_FLAGSHIP_PARQUET_COHORTS` list of 3 paths (lines 39-43). `select_disjoint(...)` consumes the exclusion set. Offline-safe (parquet load wrapped in try/except).
- `scripts/ncbi_pd_provenance_census.py` — Stage-1 census. `census_one(group, drug)` returns a dict `{other_R, other_S, powered, top_other_centers, ...}` keyed on `group` (NOT `organism`). `main()` prints only; writes nothing (confirmed).
- `scripts/build_validation_report_card.py` — v0 roll-up. `load_scored()` / `load_census()` / `load_registry()` + `classify()` (6-state machine) + `main()` unions scored+census+registry keys then appends OTHER_KINGDOM. No DRUG_RULE-driven row completeness.
- `dna_decode/eval/amr_rules.py` — `DRUG_RULE` dict (line 112, 6 drugs) + `load_calibrated_registry()` (line 144) + `rule_for()`. These are the authoritative shipped-bacterial-decoder sources.
- `dna_decode/data/{mic_tiers,antiviral_amr,antimalarial_amr,fungal_amr}.py` — CLI drug catalogs (`supported_drugs()`, antiviral peramivir/zanamivir, antimalarial artesunate/dihydroartemisinin/chloroquine, fungal anidulafungin) — the coverage source for the shipped-surface registry.
- `data/processed/*.parquet` — 8 cohorts (gate_a, gate_b, gate_b_mini_cef/tet/cohort, gate_b_n40_cipro, stage2_n150_cipro); only 3 are in the hardcoded exclusion list. `dna_decode/data/cohort.py::load_cohort` reads them (`.strains[].assembly_accession`).
- `data/raw/<slug>_*/selected.tsv` — per-cohort accession+label lists (tuning/calibration/validation roles encoded only by directory-name convention).

### Reusable-Code Survey
- `dna_decode/eval/amr_rules.py::DRUG_RULE` + `load_calibrated_registry` — reuse as the bacterial-decoder enumeration source (Step 4), no duplication.
- `dna_decode/data/cohort.py::load_cohort` — reuse to read parquet cohort accessions in the manifest builder (Step 1).
- `scripts/provenance_disjoint_validate.py::_prior_cohort_accessions` — the function the registry REPLACES; its glob + parquet logic is the seed for Step 1's builder.
- None — searched: scripts/, dna_decode/eval/, dna_decode/data/, graphify-out/ (absent) for any existing accession_manifest/cohort_role/prior_cohort concept; only provenance_disjoint_validate.py has it.

## Pre-Change Baseline
- v0 report card: 16 cells (3 SCORED, 3 POWERED_UNSCORED, 1 UNDERPOWERED, 2 ABSTAINS_BY_DESIGN, 6 NO_FREE_PHENOTYPE_SOURCE) at `wiki/decoder_validation_report_card.md`. (Post-Klebsiella-scoring: 6 SCORED.)
- Leakage exclusion: hardcoded to 3 parquet cohorts; 5 parquet cohorts in `data/processed/` are NOT excluded (gate_a, gate_b_mini_cef, gate_b_mini_tet, gate_b_mini_cohort, gate_b_n40_cipro).
- Census: stdout-only; `wiki/provdisjoint_census_results.json` was hand-written this session.
- Report-card rows: union of observed keys only; E. coli cef/tet/gent/meropenem/oxacillin do NOT appear (no NOT_CENSUSED completeness).
- Existing tests: `tests/test_build_validation_report_card.py` 7 green; `tests/test_amr_rules*.py` green.

## Verification Signal
- Manifest registry enumerates ALL parquet cohorts in `data/processed/` (assert >=8, vs current hardcoded 3) plus all `data/raw/*/selected.tsv`.
- `provenance_disjoint_validate.py` exclusion count for an E. coli run rises to include every parquet-cohort accession (no hardcoded list remains in the file).
- A census run writes/updates `wiki/provdisjoint_census_results.json` idempotently (re-running the same organism x drug updates in place, not duplicates).
- Report card renders the full shipped-decoder grid: E. coli x {cipro,cef,tet,gent,meropenem,oxacillin} present, with un-censused cells as `NOT_CENSUSED`.
- All new unit tests green; the existing 7 report-card tests + amr_rules tests unchanged (0 regressions).

## Implementation Steps

### Step 1: Accession-manifest registry module
Files: dna_decode/eval/cohort_manifest.py
Depends on: none

**What changes:**
- New module. `build_manifest(data_raw="data/raw", data_processed="data/processed") -> dict` scans every `data/raw/*/selected.tsv` AND every `data/processed/*.parquet`, recording each accession with its cohort name, `role` (heuristic from name: `*_provdisjoint_*`->validation, `stage2_n150*`->tuning, `gate_b*`->held-out, `*_cipro`/`*_<drug>`->calibration, else unknown), `source` (`selected_tsv`|`parquet`), `organism`, `drug` (parsed from cohort name where present).
- `prior_accessions(manifest, exclude_self_substr) -> set[str]` returns the union of accessions across all cohorts whose name does NOT contain `exclude_self_substr` — the data-driven replacement for the hardcoded `_FLAGSHIP_PARQUET_COHORTS`.
- Offline-safe: a parquet that fails to load is skipped with a recorded warning, not a crash (preserve current try/except posture).

**Test strategy:**
- Unit: synthetic data/raw + data/processed fixtures -> assert role/source/organism classification, union correctness, exclude_self filtering, and that a malformed parquet degrades to skip-with-warning.

### Step 2: Wire the validator to the registry
Files: scripts/provenance_disjoint_validate.py
Depends on: Step 1

**What changes:**
- Replace `_prior_cohort_accessions` body + delete the hardcoded `_FLAGSHIP_PARQUET_COHORTS` list; call `cohort_manifest.build_manifest()` + `prior_accessions(manifest, exclude_self_substr=f"{slug}_provdisjoint_")`.
- Preserve the existing print line (counts of excluded accessions) and the offline-safe fallback (data/raw glob still works if parquet load fails).
- Behavior-preserving for the 3 already-hardcoded cohorts; ADDITIVE for the other 5 (cross-organism accessions don't collide, so harmless where organisms differ; correct where they match).

**Test strategy:**
- Unit: monkeypatch build_manifest to a fixture; assert the exclusion set now contains all parquet-cohort accessions; assert a self-cohort substring is excluded from the exclusion set (not self-blocked).

### Step 3: Census self-persists powering JSON
Files: scripts/ncbi_pd_provenance_census.py
Depends on: none

**What changes:**
- After `census_one`, upsert the result into `wiki/provdisjoint_census_results.json` (schema `provdisjoint-census-results-v1`): match on `(organism, drug)`, update in place or append, keep `min_per_class` + `ecosystem_excluded` header. Preserve all existing stdout.
- Idempotent: re-running an organism x drug replaces its row, never duplicates.

**Test strategy:**
- Unit: call the upsert helper on a temp JSON; assert insert-then-update leaves exactly one row per (organism,drug) with the latest counts; assert header fields preserved.

### Step 4: Auto-derive the full shipped-decoder grid
Files: scripts/build_validation_report_card.py
Depends on: none

**What changes:**
- Add a `shipped_decoder_rows()` helper that yields the canonical row set: E. coli x every `DRUG_RULE` drug (via `amr_rules.DRUG_RULE`), plus a small declared canonical-organism map for drugs whose validated organism is not E. coli (meropenem->Klebsiella, oxacillin->Staphylococcus_aureus), plus `load_calibrated_registry()` keys, plus the existing OTHER_KINGDOM catalog.
- `main()` unions `shipped_decoder_rows()` with the observed scored+census+registry keys before classifying, so un-censused shipped decoders render as `NOT_CENSUSED` cells (grid completeness). Classification + 6-state machine + honest tier + no-aggregate-headline unchanged.

**Test strategy:**
- Unit: assert `shipped_decoder_rows()` includes E. coli x all 6 DRUG_RULE drugs + oxacillin->S. aureus + meropenem->Klebsiella; assert a drug with no census/score classifies `NOT_CENSUSED`; re-run existing 7 state-machine tests unchanged.

### Step 5: Tests + docs for the standing capability
Files: tests/test_cohort_manifest.py, tests/test_census_persist.py, CLAUDE.md
Depends on: Step 1, Step 2, Step 3, Step 4

**What changes:**
- New test files for the manifest registry (Step 1/2) and census persistence (Step 3); extend `tests/test_build_validation_report_card.py` for the auto-derived grid (Step 4).
- Document the standing workflow (census->persist->score->roll-up) + the append-incremental re-run model (new decoder -> one new cell; full-suite rebuild is ~1-2hr Docker per scored cell) + the accession-manifest registry as the leakage source-of-truth, in the CLAUDE.md decoder-suite gotchas.

**Test strategy:**
- Run the full suite (`uv run pytest tests/ -q`); assert 0 regressions + all new tests green.

## Execution Preview
- Wave 0 (parallel): Step 1, Step 3, Step 4 — independent (manifest module, census script, report-card script touch disjoint files).
- Wave 1: Step 2 (depends on Step 1).
- Wave 2: Step 5 (depends on 1-4).
- Total waves: 3. Max parallelism: 3 (Wave 0). Critical path: Step 1 -> Step 2 -> Step 5 (length 3).

## Risk Flags
- **Role-classification heuristic is name-based** — a cohort whose directory/parquet name doesn't match the role patterns lands as `unknown` role. Mitigation: `unknown` still contributes accessions to the exclusion union (fail-safe over-exclude), and the manifest is inspectable JSON; role precision is advisory, exclusion is the load-bearing output.
- **Census upsert writes a shared `wiki/` artifact** — multi-producer namespace (build_validation_report_card.py reads it). Pre-exec `/brainstorm` flagged silent-overwrite / schema-drift — addressed by the M1 amendment (group->organism normalizer + refuse capped/error rows).
- **Cross-organism over-exclusion** — including all parquet accessions for any organism's run over-excludes when organisms differ, but accessions don't collide across species so this only shrinks the available pool harmlessly; documented, not a correctness risk.
- **sentrux absent** — architecture gate n/a (not installed).

## Open Questions
- Canonical organism label for E. coli across census/registry/scored keys (`Escherichia_coli_Shigella` vs `Escherichia` vs `E_coli`) — pick one + a normalizer at implementation time.
- Whether `LABEL_CONFOUNDED` (M2) generalizes beyond oxacillin to future surrogate mismatches — leave open; add as needed.

## Verification
- `uv run pytest tests/ -q` -> 0 regressions, all new tests green.
- Manual: run census for one organism x drug -> confirm `wiki/provdisjoint_census_results.json` upserts one row.
- Manual: `--select-only` run of `provenance_disjoint_validate.py` for E. coli -> confirm excluded-accession count exceeds the old hardcoded-3-cohort count.
- Rebuild `scripts/build_validation_report_card.py` -> confirm E. coli x 6 DRUG_RULE drugs appear, un-censored ones as NOT_CENSUSED.

## Save-time amendments

Captured at: 2026-06-10
Source: `/save-plan` arguments

**Audit-notes-only:** these document the pre-exec `/brainstorm` findings. They are STRUCTURAL (they change Step 2/3/4 contracts + add a file + a 7th cell-state), so the Implementation Steps above are now out-of-sync with this intent. Re-run `/technical-plan` to fold these into the Steps before `/execute-plan`.

- C1 (Step 2, critical): leakage exclusion must use EXACT-self cohort identity (resolved `selected.tsv` path / cohort-dir), NOT the `exclude_self_substr` substring — the substring excludes ALL same-slug provdisjoint cohorts, so a re-run could reuse accessions from an earlier provdisjoint validation (false-independence). AND an incomplete accession-source load must produce `INCOMPLETE_MANIFEST` and FAIL CLOSED by default (no scoring, no artifact); opt-in `--allow-incomplete-manifest` stamps the degraded independence into the JSON+MD.
- C2 (Step 4, critical): replace ad-hoc `shipped_decoder_rows()` (E. coli x every DRUG_RULE drug over-enumerates E. coli oxacillin; OTHER_KINGDOM under-enumerates peramivir/zanamivir/artesunate/dihydroartemisinin/chloroquine) with a checked-in `dna_decode/data/shipped_decoder_surface.{py,json}` registry of `{engine, organism_scope, drug, phenotype_source_status, census_group?}`. Rows = registry ∪ observed scored/census keys. Add a coverage test asserting every CLI-routable drug (`supported_drugs() | supported_fungal_drugs() | antiviral | antimalarial`) appears in >=1 registry row. Surface = deployed-claim set, NOT a cross-product.
- M1 (Step 3, medium): add a `census_result_to_sidecar_row(result, date, min_per_class)` normalizer mapping `group`->`organism` (census_one returns `group`, sidecar+load_census use `organism`) AND refusing to upsert error/row-capped rows (never overwrite a prior good powering row with degraded data).
- M2 (Step 4, medium): add a 7th cell-state `LABEL_CONFOUNDED` for S. aureus x oxacillin (oxacillin AST is an unreliable mecA surrogate; cefoxitin is the CLSI surrogate) — distinct from `NOT_CENSUSED` so the honesty surface isn't misread.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
