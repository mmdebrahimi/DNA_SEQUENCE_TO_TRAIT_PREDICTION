# Constraint Hierarchy — Type A Laws, Clade Specialization, and Kill-Count Reporting — Technical Plan

> A non-frozen `constraints/` package: universal biological laws, two-axis (taxonomic × regime) specialization with explicit precedence, kill-count reporting so an inert filter can never pose as a control, and consolidation of five duplicated codon tables onto one clade-aware source.

## Lens status
Inputs: conversation
Degradations: repo-index unconfigured (no Layer-1/Layer-2 context; direct file reading used instead) · sentrux unavailable (Architecture Impact gates record `n/a`) · NCBI translation tables 2/4/6 are `[unverified]` — assignments must be fetched with a verbatim quote before the clade tables ship (see Risk Flags) · no `project-rules.md`, so the Project Rules Check section is omitted · no `DESIGN.md`, but this plan touches no UI

## Problem Statement

Phase 5's premise is that we will never have wet-lab characterisation for all life, so the model must carry constraints that hold for organisms it has never seen. Two gaps block that today, and both were measured during Step 1 rather than assumed.

**Gap 1 — the genetic code is hardcoded as universal, in five places.** `dna_decode/forward/genome_edit.py:18` (`_CODON`), `dna_decode/forward/inverse.py:45` (`CODON_TABLE`), `dna_decode/typing/codon_map.py:10` (`CODON`), `scripts/forward_inverse_roundtrip.py:67`, and `scripts/fungal_erg11_caller.py:53` each carry their own 64-codon dict. All three in-package copies were compared programmatically and are **byte-identical** standard tables with `TGA/TAA/TAG → '*'` and `AGA → 'R'`. That is NCBI translation table 1, and it is wrong for real clades: Mycoplasma/Spiroplasma read UGA as tryptophan, vertebrate mitochondria read AGA/AGG as stop, ciliates read TAA/TAG as glutamine. Pointed at a Mycoplasma gene today, the forward cell would call a normal tryptophan a nonsense mutation. This is **not** a live defect — every shipped cell (E. coli, Klebsiella, HIV, SARS-CoV-2, Candida, TB) uses standard assignments — it is a latent assumption that becomes a defect precisely when the organism set widens, which is what Phase 5 does.

**Gap 2 — there is no constraint layer, and no discipline that an inert filter must declare itself.** Constraint-shaped logic exists but is scattered and implicit: `genome_edit.py` checks reading frame and REF-base match inline, the FBA cell enforces mass balance via cobrapy, and the HIV/fungal callers self-check reference translation. Nothing collects these as a reusable set, nothing specializes them per clade, and — measured in Step 1 — there is **no shared kill-count primitive** (only ad-hoc `n_excluded` counters in three unrelated scripts). This repo has twice shipped a filter that removed nothing and presented it as a control: the HIV censoring guard filtered `<`/`>` operator-prefixed folds when zero such values exist in any column, and the ResFinder comparison excluded `POINT` rows when zero POINT rows exist across all 1,818 runs. Both had to be retracted. A constraint layer without kill-count reporting would reproduce that failure at scale.

## Codebase Context

**The additive-overlay precedent to mirror.** `dna_decode/data/experimental_drug_rules.py` is a NON-FROZEN, scorer-local rule module created expressly because the frozen `amr_rules.DRUG_RULE` could not express an AND-across-gene-families rule. Its docstring states the contract: "Adding a rule here touches NO frozen file." The new `constraints/` package follows that shape exactly — additive, non-frozen, no frozen-surface edit.

**The universal-plus-specialization API precedent.** `dna_decode/data/mic_tiers.py` already implements this structure: per-drug specifics (`DRUG_BREAKPOINTS:37`, `DRUG_LOCI_BY_MECHANISM:193`) alongside shared cross-drug modifiers (`CO_RESISTANCE_MECHANISMS:280`), reached through accessor functions (`breakpoints_for:71`, `amrfinder_classes_for:169`, `loci_by_mechanism_for:283`). The constraint registry mirrors this accessor-over-table design rather than inventing one.

**The provenance-refusal precedent.** `dna_decode/data/ecoff_catalog.py` requires BOTH `source_url` and `verbatim_quote` per entry and its accessor **raises** `UnsourcedEcoffError` rather than returning a plausible default — verified in Step 1. The constraint registry reuses this exact stance for `sourced` constraints.

**The precedence precedent.** `wiki/decisions-log.md` records a ratified tier-precedence decision from the genome-map virulence overlay ("AMR determinant-phenotype wins tier precedence over virulence"), establishing that in this codebase an ordering between a general and a specialized rule is made explicit and tested, never left implicit.

**Codon consumer surface (regression risk, enumerated by grep).** In-package: `forward/genome_edit.py:35,77`; `forward/inverse.py:155,158,165`; `typing/codon_map.py` (exports `translate`, `subject_aa_by_codon`). Scripts: `fungal_erg11_caller.py:53,116`; `forward_inverse_roundtrip.py:143,149`; `build_holdout_hybrid_manifest.py:88`; `scripts/kaggle/pear_prosst_kernel.py:46`. Two facts shape the plan: `dna_decode/pointfinder/runner.py:18` already imports from `typing.codon_map`, making that module the de-facto canonical home and the public name to preserve; and `scripts/kaggle/pear_prosst_kernel.py` builds its table from a compressed string because it executes standalone on Kaggle with no access to the package, so it cannot import and must instead be test-pinned.

**`scripts/fungal_erg11_caller.py` is the highest-blast-radius consumer.** Per `CLAUDE.md`, its `observed_substitutions` codon-mapper is reused by the HIV RT/PR/IN callers and the SARS-CoV-2 Mpro caller as well as fungal ERG11. Any behavioural change there reaches three validated cells at once.

**Regime already has an owner.** `dna_decode/eval/regime.py` holds eight regimes keyed `(population, endpoint, method)` with a `key` field. The regime axis of this plan **imports** those keys rather than re-declaring them — a second copy is the drift this repo has hit five times.

**Test infrastructure.** `tests/conftest.py` exists; 505 files under `tests/`; `pyproject.toml:241` sets `testpaths = ["tests"]` and `addopts = "-v"`. Convention is one test module per feature with real-data tests guarded by `pytest.skip` when gitignored data is absent.

### Reusable-Code Survey
- `dna_decode/data/experimental_drug_rules.py` — the non-frozen additive-overlay pattern this package copies; proves a new rule layer need not touch a frozen file
- `dna_decode/data/ecoff_catalog.py` — `source_url` + `verbatim_quote` requirement and the raise-on-unsourced accessor, reused verbatim for `sourced` constraints
- `dna_decode/data/mic_tiers.py` — universal-plus-specialization table/accessor shape (`CO_RESISTANCE_MECHANISMS` vs per-drug tables) to mirror for the two bucket axes
- `dna_decode/eval/regime.py` — the eight existing regime keys; imported as the regime axis, never re-declared
- `dna_decode/typing/codon_map.py` — existing `translate` / `subject_aa_by_codon` helpers and the public name `CODON` that must keep resolving for `pointfinder`

## Pre-Change Baseline

Captured 2026-10-01 before any edit; each is re-checkable by the command shown.

- **Codon tables:** three in-package copies, each 64 codons, stops exactly `{TAA, TAG, TGA}`, and all three **byte-identical** — `python -c "from dna_decode.forward import genome_edit as G, inverse as I; from dna_decode.typing import codon_map as C; print(G._CODON==I.CODON_TABLE==C.CODON, len(G._CODON))"` → `True 64`. `TGA → '*'`, `AGA → 'R'` in all three.
- **Codon consumer sites:** 8 call sites across 3 package modules and 4 scripts (enumerated in Codebase Context).
- **Clade-specific translation tables present:** 0.
- **Shared kill-count primitive:** 0 (ad-hoc `n_excluded`/`correction_removed` in `scripts/ar_bank_caur_validate.py`, `scripts/build_phage_receptor_report.py`, `scripts/doubt_layer_per_cell.py` only).
- **Full test suite:** 5157 passed, 11 skipped, 0 failed (`uv run pytest tests/ -q`, measured 2026-10-01).
- **Frozen AMR surface:** `amr_rules.py`, `calibrated_amr_rules.json`, `mic_tiers.py`, `shipped_decoder_surface.py`, `cohort_manifest.py` — must remain byte-identical (`git diff --quiet HEAD -- <paths>`).
- **Prospective lock:** v2 (`lock_date 2026-08-31`) verifies `ok=True`.
- **Evidence surface:** 115 registered cells (`uv run python scripts/project_status.py`).

## Verification Signal

- **Pure-refactor proof:** `constraints.codon_tables.STANDARD` equals the pre-change literal for all three in-package tables, asserted against values captured in the test file, and `typing.codon_map.CODON` still resolves for `pointfinder` (its import is unchanged).
- **The clade table changes behaviour, measurably:** on a CDS containing a mid-frame `TGA`, the standard table yields an internal stop and table 4 yields `W` with no stop — a diff, not a tautology. Run on a real Mycoplasma CDS if the fetch succeeds; synthetic CDS otherwise, labelled as such.
- **Kill counts are reported:** the probe prints `n_evaluated / n_refuted / status` per constraint over real CDS from the committed MG1655 reference and the committed HIV RT/PR/IN references; every constraint with `n_refuted == 0` is labelled `inert`, not silently counted as protection.
- **The conflict detector is non-vacuous:** a deliberately contradictory specialization (general says X, clade says not-X, same key) raises; removing the contradiction passes.
- **Unsourced refusal fires:** a `sourced`-kind constraint without `source_url` + `verbatim_quote` raises `UnsourcedConstraintError`.
- **No regression:** full suite ≥ 5157 passed / 0 failed; frozen AMR surface `git diff --quiet` clean; prospective lock still `ok=True`.

## Implementation Steps

### Step 1: Canonical clade-aware codon tables
Files: dna_decode/constraints/__init__.py, dna_decode/constraints/codon_tables.py, tests/test_constraints_codon_tables.py
Depends on: none

**What changes:**
- New non-frozen `dna_decode/constraints/` package, following the `experimental_drug_rules.py` contract (additive; touches no frozen file).
- `codon_tables.py` defines `STANDARD` (NCBI table 1) with assignments copied verbatim from the existing identical in-package dicts, so consolidation is provably behaviour-preserving.
- `TABLE_IDS` maps a clade key to an NCBI translation-table id; `table_for(clade: str) -> dict[str, str]` returns the table.
- `table_for` **REFUSES** an unrecognised clade (raises `UnknownCladeError`) rather than silently defaulting to standard — a wrong table is worse than a declared unknown.
- Non-standard tables (2, 4, 6) are declared with `source_url` + `verbatim_quote` fields and ship **only after** the Risk-Flag verification step; until then `table_for` raises `UnverifiedTableError` for them. Table 11 (bacterial) ships now because its codon assignments are identical to table 1 (it differs only in permitted start codons, which this layer does not model).

**Test strategy:**
- `tests/test_constraints_codon_tables.py`: assert `STANDARD` has 64 entries and stops exactly `{TAA,TAG,TGA}`; assert `STANDARD['TGA'] == '*'` and `STANDARD['AGA'] == 'R'` (the pre-change values, pinned as literals).
- Assert `table_for('unknown_clade')` raises `UnknownCladeError`, and that an unverified table id raises `UnverifiedTableError` — both proven non-vacuous by also asserting a verified id succeeds.

### Step 2: Constraint registry with two axes, precedence, and conflict detection
Files: dna_decode/constraints/registry.py, tests/test_constraints_registry.py
Depends on: none

**What changes:**
- `Constraint` dataclass: `key`, `scope` (`universal` | `taxon:<bucket>` | `regime:<key>`), `kind` (`axiomatic` | `sourced`), optional `source_url` / `verbatim_quote`, and a `check(ctx) -> ConstraintVerdict` callable. A law (mass balance) is `axiomatic` and needs no citation; a clade fact (NCBI table 4) is `sourced` and does.
- `register()` raises `UnsourcedConstraintError` when `kind == 'sourced'` and either provenance field is missing — the `ecoff_catalog` stance reused.
- `TAXON_BUCKETS`: the taxonomic axis (`virus`, `bacteria`, `archaea`, `fungi`, `plant`, `animal`, and nested keys such as `animal.mammal`), carrying molecular-machinery nuances only.
- The regime axis is **imported** from `dna_decode.eval.regime` (its eight `key` values), never re-declared here.
- `constraints_for(taxon=None, regime=None)` returns universal ∩ taxon-specialization ∩ regime-specialization, with **specialization overriding general on an equal `key`** — the explicit ordering the genome-map tier decision established as this codebase's convention.
- `detect_conflicts()` raises `ConstraintConflictError` when a general and a specialized constraint share a `key` and disagree in a way the precedence rule does not resolve (e.g. two sibling specializations both claiming the same key for the same bucket).

**Test strategy:**
- `tests/test_constraints_registry.py`: a specialized constraint overrides a universal one on the same key; a sibling collision raises `ConstraintConflictError`; removing the collision passes (non-vacuity).
- Assert the regime keys returned match `eval.regime` exactly, so a drifted copy fails.
- Assert `register()` raises on a `sourced` constraint missing either provenance field, and succeeds with both.

### Step 3: Kill-count reporting primitive
Files: dna_decode/constraints/report.py, tests/test_constraint_kill_count.py
Depends on: none

**What changes:**
- `ConstraintReport` accumulates per-constraint `n_evaluated`, `n_refuted`, `n_passed`, `n_inapplicable`.
- `status` is derived, never set: `active` when `n_refuted > 0`, `inert` when `n_evaluated > 0 and n_refuted == 0`, `not_evaluated` when `n_evaluated == 0`. The three states are distinct for the same reason the doubt layer keeps three — collapsing "removed nothing" into "passed" is how a filter poses as a control.
- `as_dict()` emits every constraint's counts and status, so an inert constraint is visible in the artifact rather than absent from it.
- `assert_no_silent_inert()` raises if any constraint is `inert` and the caller has not explicitly acknowledged it via `acknowledge_inert(key, reason)` — making inertness a decision rather than an oversight.

**Test strategy:**
- `tests/test_constraint_kill_count.py`: a constraint that refutes nothing reports `inert` (not `active`, not `passed`); one that was never evaluated reports `not_evaluated`; the two are not conflated.
- `assert_no_silent_inert()` raises on an unacknowledged inert constraint and passes once acknowledged — both directions asserted.
- Encode the two historical retractions as named regression cases: a zero-removal filter must not be reportable as a control.

### Step 4: Universal (axiomatic) constraints
Files: dna_decode/constraints/universal.py, tests/test_constraints_universal.py
Depends on: Step 1, Step 2, Step 3

**What changes:**
- Register the Type A laws, each `kind='axiomatic'`, `scope='universal'`, each returning a `ConstraintVerdict` of `refuted` / `satisfied` / `inapplicable`:
  - `cds_length_multiple_of_three`
  - `no_internal_stop_codon` — takes the clade table from Step 1, so it is universal in form and clade-specialized in data
  - `ref_base_matches_cds` — the check `genome_edit.py` already performs inline, lifted to a named constraint
  - `substitution_reachable_by_single_nt` — the codon-accessibility constraint whose absence previously let a conservativeness finding be over-read
  - `start_codon_present`
- Each constraint returns `inapplicable` rather than `satisfied` when its precondition is absent (no CDS supplied, etc.), so the report can distinguish "checked and clean" from "could not check".

**Test strategy:**
- `tests/test_constraints_universal.py`: each constraint refutes a hand-built violating case and satisfies a hand-built clean case; each returns `inapplicable` on missing input.
- `no_internal_stop_codon` is driven with both the standard table and table 4 over the same mid-frame-`TGA` sequence, asserting opposite verdicts — the clade table is load-bearing, not decorative.

### Step 5: Clade and regime specializations
Files: dna_decode/constraints/specializations.py, tests/test_constraints_specializations.py
Depends on: Step 1, Step 2

**What changes:**
- Taxonomic specializations carrying molecular-machinery nuances, each `sourced` where it asserts a biological fact: genetic-code table per clade (Mycoplasma/Spiroplasma → 4, vertebrate mitochondria → 2, ciliate → 6, bacteria → 11); `splicing_applies` for eukaryotes; `operon_polycistronic` for bacteria; `hgt_acquisition_possible` for bacteria (the property that makes a distributed mobile-element mechanism bacteria-specific); `ploidy` for eukaryotes.
- Regime specializations keyed by the imported `eval.regime` keys, recording which constraint set applies per regime — the axis this project's own measurements found more predictive than taxonomy (natural vs constructed population design; concentrated vs distributed determinant; loss- vs gain-of-function).
- Document inline that "higher life form" is not a modelled axis: clade nuances **swap** rather than accumulate, so no bucket inherits "more constraints" by depth.

**Test strategy:**
- `tests/test_constraints_specializations.py`: assert `constraints_for(taxon='bacteria.mycoplasma')` resolves the genetic-code constraint to table 4 while `constraints_for(taxon='bacteria')` resolves it to table 11 — specialization overrides general.
- Assert every `sourced` specialization carries both provenance fields (fails the unverified ones until the Risk-Flag verification lands, which is the intended gate).
- Assert `detect_conflicts()` is clean across the shipped specialization set.

### Step 6: Point in-package codon consumers at the canonical table
Files: dna_decode/forward/genome_edit.py, dna_decode/forward/inverse.py, dna_decode/typing/codon_map.py, tests/test_codon_consumers_canonical.py
Depends on: Step 1

**What changes:**
- Replace the three literal dicts with imports of `constraints.codon_tables.STANDARD`, preserving each module's existing public name (`_CODON`, `CODON_TABLE`, `CODON`) as an alias so no consumer import changes.
- `dna_decode/pointfinder/runner.py` is deliberately **not** edited — it imports `typing.codon_map`, whose public surface is unchanged.
- Behaviour-preserving by construction: the three dicts were verified byte-identical to `STANDARD` in Step 1.

**Test strategy:**
- `tests/test_codon_consumers_canonical.py`: assert `genome_edit._CODON is constraints.codon_tables.STANDARD` (and the same for the other two) — one object, not three copies.
- Assert each alias still equals the pinned pre-change literal, so the refactor is proven inert behaviourally.
- Re-run the existing forward/inverse and pointfinder suites unchanged; all must pass without edits.

### Step 7: Point script-level codon consumers at the canonical table
Files: scripts/fungal_erg11_caller.py, scripts/forward_inverse_roundtrip.py, scripts/build_holdout_hybrid_manifest.py, tests/test_script_codon_consumers_canonical.py
Depends on: Step 1

**What changes:**
- Same alias-preserving import substitution for the three script-level copies.
- `scripts/fungal_erg11_caller.py` is treated as the highest-risk edit: its codon mapper is reused by the HIV RT/PR/IN and SARS-CoV-2 Mpro callers, so the change is import-only with the alias retained and no call-site edits.
- `scripts/kaggle/pear_prosst_kernel.py` is **left alone** — it executes standalone on Kaggle without the package and cannot import it. It instead gains a drift pin (Step 9).

**Test strategy:**
- Re-run `tests/test_fungal_erg11_caller.py`, `tests/test_hiv_targetsite_caller.py`, `tests/test_sarscov2_caller.py` and `tests/test_hiv_rt_caller.py` unchanged — these carry the reference-integrity self-checks that assert the committed reference translates to the catalog wild-type at every catalogued position, which is the strongest available guard on this edit.
- `tests/test_script_codon_consumers_canonical.py`: assert each script alias equals `STANDARD`.

### Step 8: Verify-in-batch — run the constraints on real data and inspect kill counts
Files: scripts/constraint_kill_count_probe.py, tests/test_constraint_kill_count_probe.py
Depends on: Step 4, Step 5, Step 6, Step 7

**What changes:**
- New probe that loads real CDS — the committed `NC_000913.3` MG1655 reference (via the existing genome/annotation loaders) and the committed HIV RT/PR/IN and SARS-CoV-2 Mpro references under `data/hiv_ref/` and `data/sarscov2_ref/` — evaluates every universal constraint, and emits `wiki/constraint_kill_count_<date>.json` with per-constraint counts and status.
- **Expectation recorded in the artifact before the run:** most constraints are predicted to be `inert` on these references, because they are curated and already satisfy the laws. An inert result is the honest outcome and must be reported as inert, not tuned until something fails.
- Demonstrates the clade case end-to-end: fetch a real Mycoplasma CDS and show the standard table reporting a false internal stop where table 4 reports tryptophan. If the fetch is unavailable the probe falls back to a synthetic CDS and **labels the result `synthetic`** in the artifact, never presenting it as a real-organism measurement.
- Calls `assert_no_silent_inert()` so an unacknowledged inert constraint fails the probe rather than passing quietly.

**Test strategy:**
- `tests/test_constraint_kill_count_probe.py`: offline unit tests over synthetic CDS covering refuted / satisfied / inapplicable paths; real-data assertions guarded by `pytest.skip` when the gitignored references are absent, following repo convention.
- Assert the artifact records the pre-stated inert expectation and the `synthetic` vs real label for the clade demo.

### Step 9: Cross-cutting drift and regression guards
Files: tests/test_codon_table_single_source.py
Depends on: Step 6, Step 7, Step 8

**What changes:**
- A guard that **discovers** codon-table definitions by scanning the package and `scripts/` for 64-entry codon dicts, rather than checking a hand-written list — the hand-enumerated-list failure mode this repo has hit five times.
- Asserts every discovered table is either `constraints.codon_tables.STANDARD` itself or an explicitly allow-listed standalone copy, and that each allow-listed copy **equals** `STANDARD` (this is what pins `scripts/kaggle/pear_prosst_kernel.py`).
- Asserts the allow-list is non-empty and each entry carries a stated reason, so a new unexplained copy fails.
- Asserts the frozen AMR surface paths are untouched by this plan's file set.

**Test strategy:**
- Prove the guard non-vacuous by temporarily introducing a divergent copy in a tmp fixture and asserting it fails, then asserting the live tree passes.
- Assert the discovery step actually finds the known copies (a discovery guard that finds nothing is the inert-filter failure in test form).

## Execution Preview

Wave 0 (3 parallel):  Step 1 — Canonical clade-aware codon tables, Step 2 — Constraint registry with two axes, Step 3 — Kill-count reporting primitive
Wave 1 (4 parallel):  Step 4 — Universal constraints, Step 5 — Clade and regime specializations, Step 6 — In-package codon consumers, Step 7 — Script-level codon consumers
Wave 2 (1 parallel):  Step 8 — Verify-in-batch kill-count probe
Wave 3 (1 parallel):  Step 9 — Cross-cutting drift and regression guards

Critical path: Step 1 → Step 4 → Step 8 → Step 9 (4 waves)
Max parallelism: 4 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — **`scripts/fungal_erg11_caller.py` has the largest blast radius of any file here.** Its codon mapper is reused by the HIV RT/PR/IN callers, the SARS-CoV-2 Mpro caller and the fungal ERG11 caller, so a behavioural change reaches three validated cells simultaneously. Mitigation: Step 7 is import-only with the `_CODON` alias preserved and zero call-site edits; the table was proven byte-identical in Step 1; and those cells' existing reference-integrity self-checks (which assert the committed reference translates to the catalog wild-type at every catalogued position) are re-run unchanged as the guard.
- Severity: high — **NCBI translation tables 2, 4 and 6 are `[unverified]`.** Their codon assignments in this plan come from model knowledge, not a fetched authoritative source. Writing biological reference data from memory is the exact fabrication hazard this project refuses elsewhere (the volunteered ECOFF value was rejected for this reason). Mitigation: `table_for` raises `UnverifiedTableError` for tables 2/4/6 until each ships with a `source_url` + `verbatim_quote` from the NCBI translation-tables page; only tables 1 and 11 (identical codon assignments) are usable before that verification lands. Verification method: fetch the NCBI page and quote the assignment lines verbatim.
- Severity: medium — **`scripts/kaggle/pear_prosst_kernel.py` cannot import the package** (it executes standalone on Kaggle) and so keeps a private codon table, creating a genuine drift surface. Mitigation: Step 9's discovery guard allow-lists it with a stated reason and asserts equality with `STANDARD`, so divergence fails a test rather than going unnoticed.
- Severity: medium — **clade mis-assignment is worse than assuming standard.** Applying table 4 to a non-Mycoplasma genome would corrupt every call. Mitigation: `table_for` refuses unknown clades rather than defaulting, clade selection is explicit opt-in (never inferred from sequence in this plan), and the default path remains the standard table so existing behaviour is unchanged by omission.
- Severity: medium — **the two bucket axes could become a second source of truth for regime.** Mitigation: Step 2 imports the eight regime keys from `dna_decode/eval/regime.py` and a test asserts exact equality, so a drifted copy fails.
- Severity: low — **constraints may all be inert on the curated references**, which would make Step 8 look like it achieved nothing. This is a legitimate and expected outcome, recorded in the artifact before the run; the mitigation is explicitly *not* to tune constraints until something fails, and `assert_no_silent_inert()` forces each inert result to be acknowledged rather than ignored.
- Severity: low — **File overlaps between steps:** none within a wave. Steps 1–5 each create a distinct new file; Step 6 touches only `forward/` and `typing/`; Step 7 touches only `scripts/`; Steps 8–9 each add one new file. Verified against the wave grouping.
- Severity: low — **Transitive imports not captured by `Files:` lists.** Step 6 changes a module that `pointfinder/runner.py` imports, and Step 7 changes a module three viral/fungal callers import, without those files appearing in any `Files:` list. Both are behaviour-preserving alias substitutions and are covered by re-running the dependents' existing suites; noted here because the dependency is real even though no edit is required.
- Severity: low — **External tool integrations.** `cobrapy` (FBA mass balance) is referenced in Problem Statement only and is not touched by any step. The NCBI translation-tables page is a documentation source, not a runtime dependency; its surface is a web page quoted at authoring time, marked `[unverified]` above until fetched. No CLI flags, config keys or exit codes are consumed by this plan. `pytest` surface verified via `pyproject.toml:241-243` (`testpaths`, `addopts`).
- Severity: low — **Restructuring applied.** Per-step tests were folded into each step's Test strategy rather than collected into a single test step, and one cross-cutting guard step (Step 9) was retained for checks that span multiple steps. Steps 6 and 7 were split by directory (package vs scripts) to honour the one-concern-per-step principle and to isolate the high-blast-radius `fungal_erg11_caller.py` edit into its own reviewable step. Each step's new test file was added to its `Files:` line after the toolkit's `test-strategy-leak` gate flagged them as undeclared.
- Severity: low — **No code generators or build tools** write any file in this plan's file set.

## Open Questions

- **Which clades ship in v1?** Recommendation: tables 1 and 11 only (identical codon assignments, no verification debt), with 2/4/6 landing as a follow-on once the NCBI fetch supplies verbatim quotes. This keeps Step 1 shippable without carrying unverified biological data, at the cost of the Mycoplasma demo in Step 8 running against a fetched reference rather than a shipped table.
- **Should the taxonomic axis be assignable from sequence, or always explicit?** This plan makes it explicit-only (refuse on unknown) because inferring a clade from sequence is itself a prediction and would put an unvalidated classifier upstream of every constraint. Worth ratifying, since an explicit-only axis means a caller must know the organism before the constraints can specialize.

## Architecture Impact

**Expected structural changes:**
- Modules added: `dna_decode/constraints/` (`__init__.py`, `codon_tables.py`, `registry.py`, `report.py`, `universal.py`, `specializations.py`), `scripts/constraint_kill_count_probe.py`, and the per-step test modules.
- Modules modified: `dna_decode/forward/genome_edit.py`, `dna_decode/forward/inverse.py`, `dna_decode/typing/codon_map.py`, `scripts/fungal_erg11_caller.py`, `scripts/forward_inverse_roundtrip.py`, `scripts/build_holdout_hybrid_manifest.py` — all import-substitution only.
- Modules deleted: none.
- Dependency direction: `constraints/` is a leaf data/logic layer. It imports only `dna_decode.eval.regime` (for the regime keys). `forward/`, `typing/` and `scripts/` import **from** `constraints/`, never the reverse.
- Boundary changes: one new cross-layer edge, `constraints → eval.regime`. No existing boundary is relocated.

**Risk areas:**
- Shared abstraction touched: the codon table itself, which is the single most widely reused datum in the repo (8 call sites, 7 files, 4+ validated cells downstream).
- Cross-layer calls widened: `scripts/` gains imports from a new package module; previously those scripts were self-contained for codon data.
- Potential cycles: `constraints → eval.regime` is the only inbound edge. `eval/regime.py` must not import `constraints/` or a cycle forms; Step 2's test asserts the regime import is one-directional.

**Sentrux gates:**
- Baseline before Wave 0: no — `n/a — sentrux not installed`
- Compare after each merged wave: no — `n/a — sentrux not installed`
- Acceptable delta: `n/a — sentrux not installed`
- Regression-response stance: structural regression is instead checked by the Step 9 discovery guard plus the frozen-surface `git diff --quiet` assertion; if either fails, investigate and fix before proceeding rather than accept-and-document.

**Anti-drift checks:**
- `eval/regime.py` has zero imports from `dna_decode.constraints` (asserted, prevents the cycle).
- Exactly one 64-entry codon dict exists in `dna_decode/`, and it is `constraints.codon_tables.STANDARD`; every other discovered copy is allow-listed with a reason and equals it (Step 9, discovery-based not hand-listed).
- `dna_decode/pointfinder/runner.py` still imports `typing.codon_map` and that import is unedited.
- The five frozen-surface files show no diff against `HEAD` after the full plan.
- `constraints/` imports nothing from `forward/`, `typing/`, `pointfinder/` or `scripts/`.

## Verification

1. `uv run python -c "from dna_decode.constraints import codon_tables as T; from dna_decode.forward import genome_edit as G, inverse as I; from dna_decode.typing import codon_map as C; print(G._CODON is T.STANDARD, I.CODON_TABLE is T.STANDARD, C.CODON is T.STANDARD)"` → `True True True`.
2. `uv run pytest tests/test_constraints_codon_tables.py tests/test_constraints_registry.py tests/test_constraints_universal.py tests/test_constraints_specializations.py tests/test_constraint_kill_count.py tests/test_codon_table_single_source.py -q` → all pass.
3. `uv run pytest tests/test_fungal_erg11_caller.py tests/test_hiv_rt_caller.py tests/test_hiv_targetsite_caller.py tests/test_sarscov2_caller.py -q` → all pass unedited (the reference-integrity guard on the highest-risk edit).
4. `uv run python scripts/constraint_kill_count_probe.py` → exits 0, writes `wiki/constraint_kill_count_<date>.json`, and prints a per-constraint table in which every `n_refuted == 0` row reads `inert`.
5. Clade behaviour differs: the probe's clade section shows an internal stop under table 1 and `W` under table 4 for the same CDS, labelled `real` or `synthetic`.
6. `git diff --quiet HEAD -- dna_decode/data/amr_rules.py dna_decode/data/calibrated_amr_rules.json dna_decode/data/mic_tiers.py dna_decode/data/shipped_decoder_surface.py dna_decode/eval/cohort_manifest.py` → clean.
7. `uv run python -c "from dna_decode.eval import prospective_lock as P; p,m=P.resolve_active_lock(); print(P.verify_lock(m).ok)"` → `True`.
8. `uv run pytest tests/ -q` → ≥ 5157 passed, 0 failed.

## Save-time amendments

Captured at: 2026-10-01
Source: `/save-plan` arguments

> **Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads
> ONLY `## Implementation Steps` for executable work; amendments are NOT executable instructions.
> If an amendment changes step contracts (file lists, dependencies, structure), re-run
> `/technical-plan` before `/execute-plan`.

- clade-aware codon consolidation
- two-axis constraint registry
- kill-count reporting
- NCBI tables 2/4/6 gated as unverified
- regime axis imported from eval/regime.py

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high,test-strategy-leak -->
