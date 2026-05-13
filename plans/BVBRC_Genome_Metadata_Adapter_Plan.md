# BVBRC Genome Metadata Adapter — Technical Plan

> Wire `BVBRC_genome.csv` (BV-BRC Genomes-tab export) into the cohort path as a new adapter module, bypassing the wrong-contract `pilot.fetch_ncbi_assembly_quality` scaffold and feeding the existing `--assembly-metadata` wire that `cohort.candidates_from_bvbrc_ast` already accepts.

---

## Problem Statement

Two `NotImplementedError` paths exist in `dna_decode/data/pilot.py`:
- `fetch_bvbrc_drug_counts` (line 104) — resolved via the local-TSV workaround.
- `fetch_ncbi_assembly_quality` (line 157) — still scaffolded; returns only `dict[str, dict[str, int]]` with 2 keys (contig_count + n50).

The user downloaded `BVBRC_genome.csv` (BV-BRC Genomes tab, 2157 rows, 152 E. coli) with columns: `Genome ID`, `MLST`, `Assembly Accession`, `Contigs`, `Contig N50`, `Size`, `CheckM Completeness`, `Collection Year`, `Isolation Country`, `Plasmids`, `CDS`. This is richer than `fetch_ncbi_assembly_quality`'s narrow return schema — it supplies all 6 fields `cohort.candidates_from_bvbrc_ast` expects in its `assembly_metadata` arg (plus signal `fetch_ncbi_assembly_quality` would have dropped).

Without assembly QC metadata, `candidates_from_bvbrc_ast` defaults `contig_count=0` and `n50=0` for every strain, which makes them fail `build_cohort`'s assembly-quality filter (`contig_count_max=500` AND `n50_min=50_000`). Real-data cohort construction is dead-in-the-water until this wire exists.

---

## Design Decisions

### D1: Bypass the scaffold instead of implementing it

**Decision:** Add a new adapter module (`dna_decode/data/bvbrc_genome.py`) that loads the CSV directly into the dict-of-dicts shape `candidates_from_bvbrc_ast` accepts. Leave `pilot.fetch_ncbi_assembly_quality` untouched as a scaffolded stub.

**Rationale:** `fetch_ncbi_assembly_quality`'s 1-arg-signature returns only 2 keys (contig_count + n50). The CSV supplies 7+ richer fields (assembly_accession, mlst, country, year + the 2). Forcing the new code through the legacy scaffold's narrow contract would throw away signal. The scaffold was designed for a live-API path the CSV bypasses entirely.

**Trade-off:** Considered "implement the NotImplementedError in-place" — rejected because the signature is wrong for the data shape. The scaffold can be replaced later when Phase 3's live API integration ships.

### D2: New CLI flag rather than overloading existing `--assembly-metadata`

**Decision:** Add `--assembly-metadata-csv` to `pipeline.py ingest`, mutually exclusive with the existing `--assembly-metadata` (YAML) flag.

**Rationale:** Distinct file formats with distinct loaders; an explicit flag makes the dispatch obvious. The YAML path stays for inline test fixtures and any future hand-authored metadata files; the CSV path is the production real-data wire.

**Trade-off:** Considered auto-detection by file extension — rejected because explicit flags are clearer when both are accepted simultaneously is ambiguous (the mutex group raises argparse 2 with a clear error).

### D3: Coverage-log line surfaces ID-namespace mismatches early

**Decision:** After loading metadata, `cmd_ingest` emits `[ingest] assembly_meta covers M / N AST strain_ids (X%); K will fail QC filter` to stdout.

**Rationale:** The load-bearing assumption is that BV-BRC's AMR Phenotypes export's `genome_id` and the Genomes export's `Genome ID` use the same key namespace (e.g., both `562.12345`). If they don't, the join silently produces `assembly_meta = {}` for every strain. The coverage line catches this on the first real ingest without requiring a separate diagnostic tool.

**Trade-off:** Considered building a separate audit tool to inspect this — rejected as scope creep. A one-line log on the existing path is sufficient.

### D4: `fetch_ncbi_assembly_quality` stays scaffolded

**Decision:** No code change to `pilot.fetch_ncbi_assembly_quality`. `CLAUDE.md` is updated to clarify the scaffold is intentional (live-API path; Phase 3 work).

**Rationale:** The CSV adapter is the supported real-data path. The scaffold's docstring + caller-side `NotImplementedError` remain accurate for the live-API contract. Removing the scaffold would create a stub-versus-no-stub asymmetry with `fetch_bvbrc_drug_counts` (which is also still scaffolded for the live-API path).

**Trade-off:** Considered deleting `fetch_ncbi_assembly_quality` — rejected to preserve the Phase 3 wire location.

---

## Implementation Plan

### Step 1: New adapter module `dna_decode/data/bvbrc_genome.py`
Files: dna_decode/data/bvbrc_genome.py
Depends on: none

**What changes:**
- `dna_decode/data/bvbrc_genome.py` — new module. Exports `load_bvbrc_genome_metadata(path, organism="Escherichia coli") -> dict[str, dict[str, object]]` that:
  - Uses `pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)` for CSV/TSV auto-detect.
  - Tolerant column mapping: `{"Genome ID": "strain_id", "Genome Name": "organism", "Assembly Accession": "assembly_accession", "MLST": "mlst", "Isolation Country": "country", "Collection Year": "year", "Contigs": "contig_count", "Contig N50": "n50", "Species": "species"}`.
  - Optional organism filter: prefers `Species` column when present (more reliable than `Genome Name`), falls back to `Genome Name str.contains(organism)`.
  - Per-row dict construction with int parsing for `contig_count`, `n50`, `year` (returns 0 on parse failure).
  - Returns `dict[str, dict]` keyed by string `strain_id`.

**Key details:**
- Constants: `DEFAULT_ORGANISM = "Escherichia coli"` (matches `ast_data.py:15`).
- Schema: `GENOME_METADATA_KEYS = ("assembly_accession", "mlst", "country", "year", "contig_count", "n50")` — matches the keys `candidates_from_bvbrc_ast` reads at `cohort.py:376-381`.
- Plasmid/chromosome resistance genes intentionally NOT populated — AMRFinder-derived; out of scope.
- Missing-column handling: if a required column is absent, raise `BvBrcGenomeError` with column name.
- Empty-rows handling: returns `{}` (matches `ast_data.py:119-124`).

**Test strategy:**
- Unit: load tiny inline CSV with all columns → verify dict shape.
- Unit: TSV variant (sep auto-detect).
- Unit: extra columns (Plasmids / CDS / Size) — silently dropped.
- Edge: missing critical column → BvBrcGenomeError.
- Edge: organism filter excludes non-target rows.
- Edge: `Contigs="N/A"` / blank → contig_count=0 (no crash).
- Edge: `Genome ID` collisions → last-write-wins, log warning.

### Step 2: `datasources.yaml` registration
Files: config/datasources.yaml
Depends on: none

**What changes:**
- Add `bvbrc_genomes:` block under existing `bvbrc_ast:` block:
  ```yaml
  bvbrc_genomes:
    default_csv_path: ""  # set by user: path to BVBRC_genome.csv
    default_filters:
      organism: "Escherichia coli"
    cache_dir: "data/cache/bvbrc_genomes"
  ```

**Key details:**
- `default_csv_path` empty by default (no auto-discovery; explicit user opt-in via CLI flag).
- `cache_dir` reserved for future (e.g., normalized parquet output); not used by this plan.

**Test strategy:**
- Existing `tests/test_bootstrap.py::test_phase1_drugs_complete` and `::test_foundation_models_complete` continue to pass (new block doesn't affect those).
- New test in `test_bootstrap.py`: assert `bvbrc_genomes` key exists with `default_filters.organism` populated.

### Step 3: `scripts/pipeline.py` accept `--assembly-metadata-csv`
Files: scripts/pipeline.py
Depends on: Step 1

**What changes:**
- `cmd_ingest` (around lines 91-98) — add CSV branch:
  - New CLI flag `--assembly-metadata-csv` (alongside existing `--assembly-metadata` YAML).
  - When CSV flag provided, call `load_bvbrc_genome_metadata(path, organism=cfg["bvbrc_genomes"]["default_filters"]["organism"])` and assign to `assembly_meta`.
  - Mutually exclusive with `--assembly-metadata` (argparse mutex group) — exit 2 with clear error if both given.
- After loading metadata, emit coverage log line: `[ingest] assembly_meta covers M / N AST strain_ids ({fmt%}); K will fail QC filter at contig_count>500 OR n50<50000`.
- argparse: add `--assembly-metadata-csv` to the `ingest` subparser around line 456.

**Key details:**
- Coverage computation: `len(set(ast["strain_id"]) & set(assembly_meta.keys()))` vs `len(set(ast["strain_id"]))`.
- QC-fail estimate: count strains where `contig_count > 500 OR n50 < 50000` (matches `CohortSelectionCriteria` defaults at `cohort.py:78-80`).
- Backward compat: `--assembly-metadata` YAML path remains unchanged.
- Error path: file not found → exit 2 with stderr message (matches existing pattern at `cmd_ingest:85`).

**Test strategy:**
- Unit: argparse accepts `--assembly-metadata-csv path.csv`.
- Unit: argparse rejects both `--assembly-metadata` and `--assembly-metadata-csv` simultaneously (mutex group → SystemExit 2).
- Integration: full `cmd_ingest` with mock AST CSV + mock genome CSV → cohort builds with non-zero contig_count/n50.
- Edge: missing genome CSV file → exit 2 with clear stderr.
- Coverage-log assertion: capture stdout, verify "covers M / N AST strain_ids" line emitted.

### Step 4: Unit tests for the adapter
Files: tests/test_data_bvbrc_genome.py
Depends on: Step 1

**What changes:**
- New file `tests/test_data_bvbrc_genome.py` mirroring the structure of `tests/test_data_ast.py`:
  - Inline mock CSV fixture with header matching BV-BRC's exact column names (Genome ID, MLST, Assembly Accession, Contigs, Contig N50, Isolation Country, Collection Year, Genome Name, Species).
  - 5-7 test rows: mix of E. coli + non-E.coli, with valid + blank metadata fields.
  - 10-12 test functions: load happy path, organism filter, missing-column raise, blank-value tolerance, TSV variant, extra-columns-ignored, key collision warning, empty-file empty-dict.

**Key details:**
- Fixture uses real BV-BRC column names verbatim (case-sensitive: `Genome ID`, not `genome_id`).
- One test pulls a 5-row real-data slice header from the user's `Downloads/BVBRC_genome.csv` (verified columns match): documents schema-of-record.
- No external file dependencies (all fixtures inline strings).

**Test strategy:**
- See above. ~12 tests total.

### Step 5: Pipeline CLI test extension
Files: tests/test_pipeline_cli.py
Depends on: Step 3

**What changes:**
- Append to `tests/test_pipeline_cli.py`:
  - `test_ingest_accepts_assembly_metadata_csv` — happy path with both flags wired.
  - `test_ingest_rejects_both_metadata_flags_simultaneously` — argparse mutex.
  - `test_ingest_csv_path_missing_exits_2` — file-not-found path.
  - Capture stdout, assert coverage-log line present.
- Reuse existing `xgboost = pytest.importorskip(...)` pattern at top of file.

**Key details:**
- Use `tmp_path` for ephemeral CSV + AST TSV.
- Use `monkeypatch.chdir(project_root)` per existing test pattern.
- Coverage-log regex: `r"covers (\d+) / (\d+) AST strain_ids"`.

**Test strategy:**
- See above. 3-4 new tests.

### Step 6: Documentation updates
Files: README.md, CLAUDE.md, TODOS.md
Depends on: Step 3, Step 5

**What changes:**
- **README.md** — Phase 1 quickstart section (around lines 65-90): add note that `pipeline ingest` accepts `--assembly-metadata-csv` as an alternative to `--assembly-metadata <yaml>`; cite the BV-BRC Genomes-tab CSV format.
- **CLAUDE.md** — Gotchas section: update the `fetch_ncbi_assembly_quality` mention to clarify it remains scaffolded *intentionally* — the CSV adapter at `dna_decode/data/bvbrc_genome.py` is the supported real-data path; live API resolution stays Phase 3.
- **TODOS.md** — under "Pre-existing known limitations": narrow the `fetch_ncbi_assembly_quality` bullet to "live BV-BRC + NCBI Datasets REST integration deferred; CSV adapter is the supported path".

**Key details:**
- README: ~5-8 lines of new prose + one code-block example.
- CLAUDE.md: ~3-line gotcha edit.
- TODOS.md: ~2-line bullet edit.

**Test strategy:**
- No unit tests for docs. Visual review.
- Skim existing tests for any that grep for `fetch_ncbi_assembly_quality` text (grep confirms 0 production callers).

---

## Wave Grouping (for /execute-plan)

```
Wave 0 (2 parallel):  Step 1 — adapter module, Step 2 — datasources.yaml
Wave 1 (2 parallel):  Step 3 — pipeline.py flag, Step 4 — adapter tests
Wave 2 (1 step):      Step 5 — pipeline CLI tests
Wave 3 (1 step):      Step 6 — docs
```

Critical path: Step 1 → Step 3 → Step 5 → Step 6 (4 waves)
Max parallelism: 2 agents

---

## Risk Flags

- **File overlaps:** None. Wave 0 is fully disjoint files; Wave 1 splits scripts/pipeline.py and tests/test_data_bvbrc_genome.py cleanly.
- **AST↔genome ID namespace mismatch** [grounded by code; UNVERIFIED on real data]: this plan assumes BV-BRC's AMR Phenotypes export's `genome_id` and the Genomes export's `Genome ID` use the same key namespace. If they don't match, the join silently produces `assembly_meta = {}`. Mitigation: Step 3's coverage-log line surfaces this immediately on the first real ingest. Confirming the assumption requires the AMR file (Option A user action).
- **fetch_ncbi_assembly_quality stays scaffolded:** intentional, documented in CLAUDE.md update (Step 6). Future Phase 3 live-API path can replace it without conflicting with this CSV adapter.
- **`config/datasources.yaml` change is non-breaking:** `bvbrc_genomes` is a new key; existing config consumers unaffected. Verified against `tests/test_bootstrap.py:test_foundation_models_complete` which only asserts a subset.

---

## Verification

End-to-end verification after execution:
1. `uv run pytest tests/ -v` → 297 + N new tests, all green (N ≈ 15-17 new tests).
2. `uv run python -m scripts.pipeline ingest --drugs ciprofloxacin --ast-tsv <BVBRC_amr.csv> --assembly-metadata-csv ~/Downloads/BVBRC_genome.csv --download-genomes` → emits coverage log line + cohort builds with non-zero contig_count/n50.
3. Manual: open generated `cohort.parquet`, confirm `contig_count` and `n50` columns are populated for ≥80% of strains.
4. `git grep "fetch_ncbi_assembly_quality"` → still present as scaffold, no new call sites.
5. Phase 2 entry criteria status in `TODOS.md` reflects updated state.

**Spot-check before executing:** Step 3's coverage-log feature is load-bearing for the namespace-mismatch risk; ensure it lands cleanly. If coverage is 0% on first real ingest, do NOT proceed to Gate B until the namespace mismatch is fixed.
