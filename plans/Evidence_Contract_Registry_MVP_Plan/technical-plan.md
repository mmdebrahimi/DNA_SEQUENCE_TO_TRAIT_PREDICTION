# Technical Plan — Evidence-Contract Registry MVP

## Lens status
Probe (repo-grounded, this session) ✅ · Brainstorm (2-round Codex, this session) ✅ · sentrux n/a (not installed). Accepted feedback folded in: cell_id key (probe 2.1), shared-vocab-not-scale (2.2), report-card-reads-from-registry (2.4), trust-layer-theater guardrails (brainstorm).

## Problem Statement
The decoder's trust state is scattered: AMR cells live in `shipped_decoder_surface` (`(organism,drug)`-keyed), PGx/typing/viral cells are *not* represented there, abstention is expressed in **10 heterogeneous vocabularies**, and falsifiers are ad-hoc scripts unlinked to cells. Goal: one **checked-in, test-enforced Evidence-Contract Registry** — a `cell_id`-keyed row per shipped cell declaring its claim, evidence tier, validation slice, label+calibration provenance, abstention vocabulary, falsifier reference, incoming-data gate, and demotion rule — that the validation report card reads from. **Success:** every CLI-routable cell has a contract (test-enforced); the registry's AMR projection equals the frozen `shipped_decoder_surface`; the report card sources its grid+tiers from the registry; no aggregate score / shared confidence scale / fused call (honesty guardrails). **Non-goals:** executing falsifiers (declare only), the monitor (#3/alt A — a later brick), #8 acquisition (parallel non-code track), modifying any frozen surface.

## Codebase Context
- `dna_decode/data/shipped_decoder_surface.py` (73 LOC) — `_FIELDS=(organism,drug,engine,organism_scope,phenotype_source_status,census_group)`, `shipped_decoder_rows()`, `surface_index()`, `all_surface_drugs()`. **AMR-only**; the registry's AMR projection must equal it.
- `scripts/build_validation_report_card.py:47` — `from …shipped_decoder_surface import surface_index`. The Step-6 refactor seam.
- `tests/test_shipped_decoder_surface.py::test_every_cli_drug_in_surface` — the coverage pattern Step 5 generalizes.
- `dna_decode/eval/prospective_lock.py` — `FROZEN_SURFACE_FILES`, `surface_hashes` (the monitor brick will later consume the registry; out of scope here).
- `tests/test_tb_leak_guard.py` — sha256-pins `amr_rules.py` + `calibrated_amr_rules.json`. The registry imports `shipped_decoder_surface` read-only and touches **neither** → leak-guard stays green.
- Abstention terms in-tree (grounded counts): `ABSTAIN`(64) `UNDERPOWERED`(22) `NOT_CENSUSED`(17) `NO_GO`(16) `SUSPEND`(15) `NO_FREE_PHENOTYPE`(12) `LABEL_CONFOUNDED`(11) `phenotype_withheld`(10) `ABSTAINS_BY_DESIGN`(10) `O?`(7).
- Falsifier scripts: `scripts/{amr,cef,cipro_bounded}_falsifier.py`, `dna_decode/eval/cohort_deconfound.py` (declarable refs).

### Reusable-Code Survey
- **`shipped_decoder_surface.py`** — reuse as the AMR projection source (import, don't duplicate).
- **`trust_surface.py` / `build_validation_report_card.py`** — the consumers to refactor onto the registry.
- **The `*_falsifier.py` scripts** — referenced (by path) in the `falsifier_ref` field.
- `None — searched:` graphify-out (absent), `dna_decode/data/`, `dna_decode/eval/`, `scripts/build_*`.

## Pre-Change Baseline
Today: trust state is split across `shipped_decoder_surface` (AMR), `pgx_report_card` (PGx), `build_validation_report_card` (AMR roll-up), and 10 scattered abstention strings; **no single cell-identity registry**; PGx/typing/viral cells absent from the deployed-claim surface. 100 PGx-suite + leak-guard green at HEAD (`8f45730`).

## Verification Signal
- `test-exit-0`: `uv run pytest tests/test_cell_registry.py -q` passes.
- `file-exists`: `dna_decode/data/cell_registry.py` + `wiki/cell_registry_report.{md,json}`.
- Registry AMR projection == `shipped_decoder_surface` rows (consistency test).
- Every CLI-routable cell (AMR drugs + PGx genes + typing schemes) has exactly one contract (coverage test).
- `tests/test_tb_leak_guard.py` green (frozen surface untouched).

## Implementation Steps

### Step 1: Cell-identity scheme + CellContract dataclass
Files: dna_decode/data/cell_registry.py
Depends on: none

**What changes:**
- New module. `cell_id = f"{track}:{organism}:{target}"` where `track ∈ {amr, pgx, typing, viral}`, `target` = drug | gene | scheme. Frozen `CellContract` dataclass with fields: `cell_id, track, organism, target, claim, evidence_tier, validation_slice, label_provenance, calibration_provenance, abstention_vocab, native_abstention, falsifier_ref, incoming_data_gate, demotion_rule`.
- `EVIDENCE_TIERS` enum (independent_measured / near_independent / faithful_to_tool / knowledge_baseline / no_free_source). NO numeric confidence field (guardrail).
- Accessors: `cells()`, `by_cell_id()`, `cli_routable_cell_ids()`.

**Test strategy:**
- dataclass-shape + enum-membership unit tests (pure).

### Step 2: Abstention-vocabulary mapping (shared vocab, NOT a scale)
Files: dna_decode/data/cell_registry_vocab.py
Depends on: Step 1

**What changes:**
- Controlled `ABSTENTION_VOCAB` enum (e.g. `SCORED / ABSTAIN_BY_DESIGN / UNDERPOWERED / WITHHELD_NONCORE / NO_FREE_SOURCE / LABEL_CONFOUNDED`) + `NATIVE_TO_VOCAB` map collapsing the 10 in-tree terms (`SUSPEND`/`phenotype_withheld`/`O?`/`NO_GO`/… → the enum). Each cell keeps its `native_abstention`; the vocab is for cross-cell *grouping*, never a probability.

**Test strategy:**
- every in-tree native term maps to a vocab value; no value implies a numeric scale (assert no float).

### Step 3: Populate AMR + PGx + typing + viral contracts
Files: dna_decode/data/cell_registry.py
Depends on: Step 1, Step 2

**What changes:**
- AMR cells: derive from `shipped_decoder_surface.shipped_decoder_rows()` (import) → one contract each (claim/tier from the existing `phenotype_source_status` + the provdisjoint/lineage report cards).
- PGx: CYP2C19 (`GeT-RM 72/72 ⋈ 1000G`, near_independent calling / faithful_to_cpic phenotype), CYP2C9 (73/73), VKORC1.
- Typing/viral cells: one contract each at their honest tier (faithful_to_tool / knowledge_baseline / no_free_source).

**Test strategy:**
- AMR-projection equals `shipped_decoder_surface`; spot-check PGx tiers vs the report-card JSONs.

### Step 4: Link falsifier refs + incoming-gate + demotion rules
Files: dna_decode/data/cell_registry.py
Depends on: Step 3

**What changes:**
- Populate `falsifier_ref` (path to `scripts/*_falsifier.py` or `none`), `incoming_data_gate` (which of G1–G8 apply, or `n/a`), `demotion_rule` (the SCORED→UNDERPOWERED→… trigger). Declared only (not executed).

**Test strategy:**
- every non-`none` `falsifier_ref` resolves to an existing file; `incoming_data_gate` ⊆ the 8 gate names.

### Step 5: Test-enforced coverage + consistency
Files: tests/test_cell_registry.py
Depends on: Step 2, Step 3, Step 4

**What changes:**
- Coverage: every CLI-routable cell (AMR drug via `all_surface_drugs()` + each PGx `--gene` + typing schemes) has exactly one contract.
- Consistency: registry AMR projection == `shipped_decoder_surface`. Enum validity for tier + vocab. Falsifier-ref existence.

**Test strategy:**
- the tests ARE the deliverable; generalize `test_every_cli_drug_in_surface`.

### Step 6: Report card reads from the registry
Files: scripts/build_validation_report_card.py, scripts/build_pgx_report_card.py
Depends on: Step 5

**What changes:**
- Refactor `build_validation_report_card` to source its grid+tiers from `cell_registry` (the AMR projection preserves current rows — additive, backward-compatible) and emit a unified `wiki/cell_registry_report.{md,json}` (or fold into the existing card). Per-cell tier only; **no aggregate health score** (guardrail). `build_pgx_report_card` reads PGx contracts from the registry.

**Test strategy:**
- report-card states map to registry-declared tiers; existing report-card tests stay green.

## Execution Preview
- **Wave 0:** Step 1 (dataclass + key).
- **Wave 1:** Step 2 (vocab — separate file).
- **Wave 2:** Step 3 (populate — needs 1+2).
- **Wave 3:** Step 4 (falsifier refs — needs 3).
- **Wave 4:** Step 5 (tests — needs 2,3,4).
- **Wave 5:** Step 6 (report card — needs 5).
- Total 6 waves; max parallelism 1 (mostly one module) → effectively sequential; critical path = all 6. Small, single-developer-friendly.

## Risk Flags
- **Frozen surface:** the registry imports `shipped_decoder_surface` **read-only** and touches no frozen file (`amr_rules.py`/`calibrated_amr_rules.json` untouched) → leak-guard stays green. Asserted in Step 5.
- **Theater risk (brainstorm):** the plan forbids an aggregate health score / shared confidence scale / fused call — enforce by *absence* (no such field exists) + a test asserting no numeric confidence field.
- **Tier subjectivity:** `evidence_tier` per cell is a judgment; mitigate by deriving AMR tiers from the existing `phenotype_source_status` + provdisjoint JSONs, not fresh opinion.
- **[unverified]:** the exact typing/viral cell list to enumerate — Step 2 must grep the CLI dispatch (`dna_decode/cli.py` TRAITS) for the authoritative routable set.

## Open Questions
1. Registry as a **new module** (recommended — clean `cell_id` key) vs in-place generalization of `shipped_decoder_surface` (reuses its coverage test). The plan assumes new module + AMR projection.
2. Fold the unified report into the existing `build_validation_report_card` output, or emit a new `cell_registry_report`? Plan assumes a new unified report that the old cards can later be retired into.
3. Should `demotion_rule` be free-text (v0) or a controlled vocabulary now? Plan assumes free-text v0 (the controlled demotion enum is the monitor brick's concern).

## Verification
`uv run pytest tests/test_cell_registry.py tests/test_shipped_decoder_surface.py tests/test_tb_leak_guard.py -q` green + `wiki/cell_registry_report.{md,json}` generated + the AMR-projection-equals-`shipped_decoder_surface` assertion passes.

## Repo grounding

### Captured by: brainstorm @ 2026-06-26
Files read: plans/Evidence_Contract_Registry_MVP_Plan/technical-plan.md, dna_decode/cli.py, dna_decode/amr/cli.py, dna_decode/pgx/cli.py, dna_decode/data/shipped_decoder_surface.py, dna_decode/data/cell_key.py, dna_decode/data/trust_surface.py, scripts/build_validation_report_card.py, scripts/build_pgx_report_card.py, tests/test_shipped_decoder_surface.py, tests/test_trust_surface.py, tests/test_build_validation_report_card.py, pyproject.toml.

Key claims (grounded):
- cli.py TRAITS = {amr, pathotype, plasmid, resfinder, pointfinder, disinfinder, mlst, ktype, salmserovar, pneumoserotype, pgx} — wider than the plan's {amr,pgx,typing,viral} 4-track model; finder/typing tools don't fit track:organism:target. ADD a 5th `finder` track + a `route` field; cell_id is a DISPLAY rendering over structured identity.
- dna-amr --drug routes bacterial+fungal+antimalarial+influenza+HIV+SARS through one option; tests/test_shipped_decoder_surface.py is blind to HIV/SARS. Step 5 coverage must assert an explicit PER-ROUTE cell manifest, not grep-TRAITS.
- dna_decode/data/cell_key.py::canonical_cell_key is the SINGLE canonical (organism,drug) join for the 3 sidecars; the AMR-projection consistency test must join via it, NOT raw-string cell_id equality (organism aliasing: Klebsiella vs Klebsiella_pneumoniae, Escherichia vs Escherichia_coli_Shigella).
- dna_decode/data/trust_surface.py is the RUNTIME badge reading 4 packaged cards force-included in pyproject.toml; Step 6 must rebuild the EXISTING decoder_validation_report_card.json + pgx_report_card.json from the registry under the same filenames/schemas + add a registry->card consistency test + packaging drift test, NOT emit a parallel cell_registry_report.json (avoids two sources of truth). Full trust_surface migration = later brick.
- AMR tier derivation mixes 3 distinct signals (phenotype_source_status claim-status + provdisjoint measured-performance + lineage disclosure); keep claim_status / validation_state / evidence_tier as SEPARATE fields, not one fused tier.

### Captured by: brainstorm @ 2026-06-26 (post-ship v0.1 honesty audit)
Files read: dna_decode/data/cell_registry.py, dna_decode/data/cell_registry_vocab.py, tests/test_cell_registry.py, scripts/build_validation_report_card.py, dna_decode/data/shipped_decoder_surface.py, dna_decode/data/hiv_amr.py, dna_decode/data/sarscov2_amr.py, wiki/hiv_decoder_report_card.json, dna_decode/data/trust_surface.py, dna_decode/pgx/cli.py.

Key claims (grounded) — for a v0.1.1 honesty patch:
- CRITICAL honesty overclaim: _viral_contracts() tags ALL 26 all_supported_hiv_drugs() independent_measured, but delavirdine is CLI-routable yet has NO row in hiv_decoder_report_card.json (25 cells). The registry tags an UNVALIDATED drug as measured + flattens real per-drug caveats (position-based ddI, small CAI n, underpowered subtype transfer). Fix: add EvidenceTier.NOT_CENSUSED; data-drive HIV tiers from the PACKAGED report card (trust_surface.py's loader path) via an explicit adapter; no-card-row -> NOT_CENSUSED; pin the coupling with a test (all_supported_hiv_drugs() - card_drugs == {delavirdine}).
- PGx coverage is FAIL-OPEN: dna_decode/pgx/cli.py:27 hardcodes --gene choices; cli_routable_manifest() + test duplicate the literal -> a 4th gene passes coverage vacuously. Fix: one importable PGX_GENES constant driving CLI choices + manifest + test.
- surface_index() builds all 67 cells (imports hiv/sars/cli) to return the 25 AMR rows -> project from _amr_contracts() directly.
- NOT a gap (verified): drug-level AMR coverage is weak BUT the amr_projection_keys()==frozen-surface consistency test is the real fail-closed organism×drug guard (AMR cells are projected from the surface). No divergent prod consumer of shipped_decoder_surface.surface_index; trust_surface reads the same frozen list the registry projects (two readers, one source).
- Open decision: encode the independent_measured-vs-underpowered threshold as a RULE in code (card row + n>=K + not position-based-only), not prose.
