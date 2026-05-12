# E. coli G2P Platform — Phase 1 Technical Plan

> Build a buildable Phase-1 (months 1-3) MVP for an E. coli genotype-to-phenotype prediction platform with biological interpretability. Antibiotic resistance as the first phenotype: **ciprofloxacin + ceftriaxone + tetracycline**. 17 steps (incl. Step 0.5 pilot gate) grouped into 8 waves. Greenfield project; every file is NEW. **Post-tech-plan brainstorm applied 2026-05-11** — Captum removed from Phase 1 (deferred to Phase 2); drug-first cohort + Mash/ANI phylogeny controls + leaderboard added.

## Problem Statement

Phase 1 deliverable: an end-to-end pipeline that ingests E. coli K-12 reference + a drug-first-validated 450-strain cohort (≥150 strains per drug for cipro/ceftriaxone/tet, broth-microdilution AST only) + CARD/AMRFinder/BV-BRC AST data, computes Evo + DNABERT-2 foundation-model embeddings (cached), trains a frozen-embedding XGBoost binary classifier per drug, evaluates via **leave-one-Mash-clade-out CV** (primary) + LOSO + LOMO (secondary) + clade-only-baseline control (target: embedding-model AUROC ≥0.10 above clade-only on ≥75% of held-out clades), and produces attribution maps via **in-silico mutagenesis only** (gene-level + nucleotide-level saturation on top-K=20 genes) that overlap known resistance loci at top-K=20 precision ≥0.6 for ciprofloxacin (gyrA / parC / qnr / aac(6')-Ib-cr) and ≥0.4 for ceftriaxone (CTX-M / SHV / AmpC). Phase 2-3 extensions (Captum IG with diff MLP head, multi-task across 10-20 drugs, MIC regression, pan-genome graph layer via Panaroo/PyG, interpretability dashboard, comparative genomics) are scoped separately.

Refined-from project ledger at `C:\Users\Farshad\PythonProjects\dna_decode\project_state\dna-decode-2026-05-11.md` — Decisions #1-3 lock the goal + plan; hypotheses H1-H11 cover the empirical questions; observations #1-9 cover the calibration evidence.

## Codebase Context

**Greenfield project.** No existing code; no `.repo-index/`; no prior `wiki/decisions-log.md` for this project (the parent `rca_engine` decisions-log is irrelevant). Step 1 of /technical-plan's standard "research the codebase" workflow is reduced to: verify external dependencies install cleanly + verify data-source accessibility + commit to a directory layout.

**Sibling-project patterns to mirror (NOT lift):**
- `C:\Users\Farshad\PythonProjects\rca_engine\backend\` — clean compute-layer + API-layer split; module-summary discipline; lazy-loaded dependencies pattern. Mirror the structure (one `compute/`-shaped package + one `scripts/`-shaped entry-point dir + one `tests/`-shaped tree).
- `C:\Users\Farshad\PythonProjects\Athena_Development\` — config-file-driven data-source registry. Mirror the pattern for declaring "which strains are in the pan-genome catalog" via a yaml config rather than hardcoding.

**External dependency landscape:**
- Foundation models: `evo-model` (PyPI; loads Evo from HuggingFace Hub via `together/evo-1-131k-base`), `transformers` for DNABERT-2 (`zhihan1996/DNABERT-2-117M`), `sentence-transformers` (already in user's env).
- Sequence data: `biopython` for FASTA/GenBank/GFF3, `pysam` for any VCF later.
- Cache: `h5py` (HDF5) for embeddings; `parquet` via `pyarrow` for tabular data.
- ML: `xgboost`, `scikit-learn`, `pytorch` (CPU-or-GPU build chosen at install time).
- Attribution: `captum` (PyTorch attribution library).
- Viz: `pygenometracks` (CLI; produces PNG/SVG from BED+annotation+attribution tracks).
- Test: `pytest`, `pytest-cov`.
- Project mgmt: `uv` (sibling pattern from rca_engine; faster than pip).

**Data-source landscape:**
- NCBI RefSeq E. coli K-12 reference: `GCF_000005845.2_ASM584v2` (MG1655). Accessible via NCBI Datasets API or FTP.
- CARD database: download from `https://card.mcmaster.ca/latest/data` (free; CC-BY-NC-SA license).
- NCBI AMRFinderPlus database: `ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/`.
- BV-BRC E. coli AST data: download via the BV-BRC FTP (`ftp.bvbrc.org`); strain-level susceptibility test results.
- EnteroBase E. coli MLST + pan-genome strain catalog: REST API at `enterobase.warwick.ac.uk`.

## Implementation Steps

### Step 0.5: Real-data pilot gate
Files: scripts/pilot_gate.py, dna_decode/data/pilot.py, tests/test_pilot_gate.py
Depends on: Step 1

**What changes:**
- `dna_decode/data/pilot.py` — `run_pilot_gate(drugs: list[str], filters: CohortSelectionCriteria) -> PilotReport` performs metadata-only HTTP calls (no genome downloads, no embedding compute): fetches BV-BRC AST manifest, NCBI assembly summaries (TSV), CARD + AMRFinder reference catalogs. Computes per-drug counts at each filter stage: raw → broth-microdilution filter → assembly-quality filter → 3-drug-intersection. Estimates download volume + embedding compute time projections for the full pipeline.
- `scripts/pilot_gate.py` — `python -m scripts.pilot_gate [--drugs ciprofloxacin,ceftriaxone,tetracycline] [--target-per-drug 150]`: writes `data/processed/pilot_report.md` summarizing per-stage counts + GO/NO-GO recommendation. Returns non-zero exit if any target-per-drug count is unmet.
- `tests/test_pilot_gate.py` — mocked HTTP responses; verify per-stage count math + report format + exit-code semantics

**Key details:**
- **HARD pre-execution gate.** Step 6 (and every step that ingests cohort genomes) depends on Step 0.5. If per-drug counts are unmet, `pilot_gate.py` exits non-zero → `/execute-plan` halts at this step and reports the failure. The human MUST rescope (relax a filter / pick alternate drug agent / drop a drug from Phase 1) and re-run Step 0.5 before downstream steps can proceed. This is intentional — fail-fast on data insufficiency is cheaper than discovering it after multi-hour embedding runs.
- Output report format: per-drug counts at each filter stage in a markdown table; ANI/MLST diversity estimate; expected download volume in GB; expected embedding compute time on RTX 4090 (calibrated from a 10-strain micro-pilot run within `pilot_gate.py` itself).
- Release condition (exit code 0): all three Phase 1 drugs (cipro / ceftriaxone / tet) have ≥150 strains after broth-microdilution + assembly-quality + per-drug-AST filters; 3-drug intersection ≥75 strains; estimated embedding-compute fits the user's compute ceiling.
- Why it's a hard gate: the original plan assumed sufficient labeled data without verifying. Post-tech-plan brainstorm surfaced that broth-microdilution + assembly-quality filters could collapse the cohort. Hard-gating Phase 1 at the metadata-only stage costs <1 hour; a soft-gate that lets execution continue would burn 4-8 hours of ingestion + cache compute before failure became visible.

**Test strategy:**
- Unit: with mocked manifests, verify filter math; verify HARD-gate exit-code semantics at the 150-per-drug threshold (exit 0 on pass, non-zero on fail)
- Edge cases: zero strains after filtering (report flags failure + suggests filter relaxation; exits non-zero); network failure mid-fetch (report flags incomplete; exits non-zero with diagnostic)

### Step 1: Project bootstrap
Files: pyproject.toml, .gitignore, README.md, dna_decode/__init__.py, dna_decode/data/__init__.py, dna_decode/models/__init__.py, dna_decode/interp/__init__.py, dna_decode/eval/__init__.py, dna_decode/viz/__init__.py, tests/__init__.py, tests/conftest.py, config/datasources.yaml
Depends on: none

**What changes:**
- `pyproject.toml` — declare project name `dna_decode`, Python ≥3.11, dependencies (biopython, transformers, torch, xgboost, scikit-learn, h5py, pyarrow, **bitsandbytes** for 4-bit Evo quantization on RTX 4090 24GB VRAM, pygenometracks-via-system-binary-or-pip-equivalent, pytest, pyyaml), tooling (ruff, mypy optional). NOTE: captum REMOVED from Phase 1 deps per post-tech-plan brainstorm C1 (Step 12 deferred to Phase 2); add when the Phase 2 differentiable-MLP-head work begins.
- `.gitignore` — ignore `data/`, `*.h5`, `*.parquet`, model checkpoints, `__pycache__/`, `.venv/`
- `README.md` — one-paragraph project description + setup instructions (uv sync; uv run pytest)
- `dna_decode/__init__.py` and module-level `__init__.py` files — empty placeholders defining the package structure (`data`, `models`, `interp`, `eval`, `viz` subpackages)
- `tests/__init__.py` + `tests/conftest.py` — empty package marker + a conftest.py with one fixture that returns the path to `tests/fixtures/` (created in Step 15)
- `config/datasources.yaml` — declarative config listing the canonical data sources (RefSeq accession, CARD URL, AMRFinder FTP path, BV-BRC FTP host, EnteroBase API base). Used by Steps 2-6 to avoid hardcoding URLs across modules.

**Key details:**
- Python ≥3.11 (matches user's other projects)
- Use `uv` for env management (sibling pattern from rca_engine `backend/`)
- Pin specific versions for `evo-model`, `DNABERT-2`, `bitsandbytes`, `mash` (system binary; document install) to avoid foundation-model + binary-tool API churn that broke prior projects
- 4-bit Evo quantization via bitsandbytes is the Phase 1 compute default (RTX 4090 24GB VRAM); Phase 2 may switch to full-precision on rented A100 if 4-bit quantization is shown to hurt attribution precision (validate this empirically before committing to A100 spend)

**Test strategy:**
- Unit: `tests/test_bootstrap.py` — verifies the package imports cleanly (`import dna_decode; import dna_decode.data; ...` for all 5 subpackages); verifies `config/datasources.yaml` loads as valid YAML with the expected top-level keys

### Step 2: NCBI RefSeq downloader + cache
Files: dna_decode/data/refseq.py, tests/test_data_refseq.py
Depends on: Step 1

**What changes:**
- `dna_decode/data/refseq.py` — `download_genome(accession: str, cache_dir: Path) -> Path` that fetches the FASTA + GFF3 + GenBank trio from NCBI Datasets API or FTP, caches under `cache_dir/<accession>/{genome.fna, annotations.gff3, annotations.gbk}`, returns the cache directory. Idempotent: short-circuits if cache exists and `--force` not passed. `list_cached() -> list[str]` enumerates cached accessions. `default_ecoli_k12_accession()` returns `"GCF_000005845.2"` (MG1655 reference).
- `tests/test_data_refseq.py` — unit tests with mocked HTTP (no live network in CI); fixture: a tiny FASTA + GFF3 pair in `tests/fixtures/ecoli_mini/` (created in Step 15)

**Key details:**
- Network calls go through `requests` with explicit timeout (30s) + retry-with-backoff (max 3)
- Cache contents are write-once: a `.complete` sentinel file is the last thing written, atomically; partial caches re-trigger download
- Error handling: 404 → `RefSeqAccessionNotFound`; 5xx + timeout → `RefSeqDownloadError` after retries
- No CLI here — that's Step 14

**Test strategy:**
- Unit: download a small mocked genome, verify cache structure, verify idempotence (second call no-ops), verify `.complete` sentinel
- Edge cases: 404 propagates correct exception; partial cache (no sentinel) re-downloads

### Step 3: Genome annotation parser
Files: dna_decode/data/annotations.py, tests/test_data_annotations.py
Depends on: Step 1

**What changes:**
- `dna_decode/data/annotations.py` — `parse_gff3(path: Path) -> AnnotationTable` produces a typed table (pandas DataFrame with stable column schema: `seqid, source, type, start, end, strand, gene_id, locus_tag, product`); `parse_genbank(path: Path) -> AnnotationTable` same shape from a GenBank file. `extract_cds_sequences(genome_fasta: Path, annotations: AnnotationTable) -> dict[str, str]` returns gene-id → nucleotide-sequence mapping for CDS features. `extract_intergenic_regions(...)` returns gene-id-pair → intergenic-sequence mapping.
- `tests/test_data_annotations.py` — uses `tests/fixtures/example.gff3` + `example.fna` (a 10-gene synthetic E. coli-like fixture; created in Step 15)

**Key details:**
- Use `biopython.SeqIO` for FASTA + GenBank, raw line parsing for GFF3 (biopython's GFF3 support is in `bcbio-gff` which is a separate dep — avoid)
- AnnotationTable is a typed DataFrame; type hints via `pandas-stubs` if added later
- Strand handling: revcomp for `-` strand CDS extraction
- Error handling: malformed GFF3 line → `AnnotationParseError` with line number context

**Test strategy:**
- Unit: parse 10-gene fixture, verify row count, verify a known gene's start/end/strand, verify CDS extraction produces correct revcomp for a `-`-strand gene
- Edge cases: empty file (returns empty table); single-gene file; gene with no `locus_tag` (allowed)

### Step 4: Resistance database loaders (CARD + AMRFinder)
Files: dna_decode/data/resistance_db.py, tests/test_data_resistance_db.py
Depends on: Step 1

**What changes:**
- `dna_decode/data/resistance_db.py` — `load_card(cache_dir: Path) -> ResistanceCatalog` downloads + parses the CARD JSON model (the `card.json` artifact from `card.mcmaster.ca/latest/data`); `load_amrfinder(cache_dir: Path) -> ResistanceCatalog` parses the NCBI AMRFinderPlus reference TSV. Unified `ResistanceCatalog` dataclass: `entries: list[ResistanceEntry]` where each has `gene_symbol, gene_family, drug_class, resistance_mechanism, source_db, source_id`. `map_gene_to_resistance(gene_symbol: str) -> list[ResistanceEntry]` does case-insensitive lookup with aliases.
- `tests/test_data_resistance_db.py` — fixture-driven; uses `tests/fixtures/card_mini.json` + `amrfinder_mini.tsv` (5-entry mock corpora, created in Step 15)

**Key details:**
- Re-download cadence: cache for 30 days unless `--force-refresh`
- CARD's data model is rich; v1 extracts only the AMR Gene Family hierarchy + drug class associations
- Gene-symbol normalization: CARD uses `gyrA` while AMRFinder uses `GyrA`; normalize to title-case + maintain an alias map

**Test strategy:**
- Unit: load mock CARD JSON + mock AMRFinder TSV, verify entries parsed, verify case-insensitive lookup, verify alias mapping
- Edge cases: missing `drug_class` field (defaults to "unknown"); gene with multiple resistance mechanisms (returns all entries)

### Step 5: AST phenotype data loader (BV-BRC)
Files: dna_decode/data/ast_data.py, tests/test_data_ast.py
Depends on: Step 1

**What changes:**
- `dna_decode/data/ast_data.py` — `load_bvbrc_ast(cache_dir: Path, organism: str = "Escherichia coli") -> ASTTable` downloads the BV-BRC AST TSV from FTP, filters to the organism, parses into a typed DataFrame with columns: `strain_id, antibiotic, susceptibility_label (S/I/R), mic_value, mic_units, measurement_method, source`. `get_drug_list(ast: ASTTable, min_strains: int = 50) -> list[str]` returns antibiotics with sufficient labeled strains for training. `binarize_susceptibility(label: str) -> int` returns 0 for S/I, 1 for R (the v1 binary classification target).
- `tests/test_data_ast.py` — fixture: `tests/fixtures/bvbrc_ast_mini.tsv` (20 rows; created in Step 15)

**Key details:**
- BV-BRC AST file is large (~1GB); cache with `.complete` sentinel
- Missing MIC values are common; v1 keeps S/I/R labels and treats missing MIC as NaN (used only when MIC regression is implemented in Phase 2)
- Categorical→binary mapping: S=0, I=0 (intermediate treated as susceptible for v1 binary task), R=1

**Test strategy:**
- Unit: load mock TSV, verify row filtering to E. coli, verify `get_drug_list` thresholding, verify binarize_susceptibility mapping
- Edge cases: empty AST file; rows with NaN MIC; rows with non-standard susceptibility labels (logged + dropped)

### Step 6: Strain/AST cohort catalog (drug-first)
Files: dna_decode/data/cohort.py, tests/test_data_cohort.py
Depends on: Step 1, Step 2, Step 4, Step 5, Step 0.5

**Renamed from "Pan-genome strain catalog" per post-tech-plan brainstorm M1.** Pan-genome *clustering* (Panaroo/Roary) is Phase 2; Phase 1 only builds a strain catalog with AST labels.

**What changes:**
- `dna_decode/data/cohort.py` —
  - `build_cohort(drugs: list[str], target_per_drug: int = 150, criteria: CohortSelectionCriteria) -> StrainCohort` selects strains drug-first: for each target drug (default Phase 1: `ciprofloxacin`, `ceftriaxone`, `tetracycline`), ensure ≥`target_per_drug` strains with non-NaN AST labels for that drug; maximize the joint intersection (strains with all 3 drugs labeled) up to `target_per_drug // 2`; maximize MLST/clade diversity within each per-drug pool.
  - `CohortSelectionCriteria` — controls: `assembly_quality_threshold(contig_count_max=500, n50_min=50_000)` (REPLACES the prior "complete-circle only" gate per M1), `plasmid_localization_required: bool` (when True, prefers strains where AMRFinder reports plasmid-contig vs chromosome-contig localization for resistance genes; informational, not hard filter), `measurement_method_filter: list[str]` (default `["broth_microdilution"]` per failure-mode #4).
  - `download_cohort_genomes(cohort: StrainCohort) -> dict[str, Path]` uses `refseq.download_genome` (Step 2) to materialize each strain.

- `tests/test_data_cohort.py` — mocked NCBI + BV-BRC responses; verifies drug-first selection produces ≥target_per_drug per drug; verifies assembly-quality filter; verifies broth-microdilution filter

**Key details:**
- Cohort target: 150 strains/drug × 3 drugs Phase 1, with ≥75 strains in the 3-drug intersection (parameterized)
- Drug-first selection means a strain can be in the cohort for *one* drug even if not labeled for the other two; per-drug training uses only the drug-labeled subset
- Assembly-quality threshold (contig_count ≤ 500, N50 ≥ 50K) preserves cohort size vs the rejected "complete-circle only" gate; AMRFinder plasmid/chromosome metadata is informational (helps interpret attribution downstream, doesn't gate selection)
- Output: `cohort_v1.parquet` with strain metadata (id, accession, MLST, country, year, contig_count, n50, AST_label_<drug> for each Phase 1 drug, plasmid_resistance_genes, chromosome_resistance_genes)
- Materializing genomes parallelizable via `ThreadPoolExecutor` (max_workers=4)

**Test strategy:**
- Unit: with mocked candidate strains containing per-drug AST coverage, verify the cohort meets ≥150 per drug; verify joint-intersection target; verify assembly-quality filter excludes contig_count > 500
- Edge cases: fewer candidate strains than target (returns all + warning); zero strains for one drug (raise `CohortConstructionError` with diagnostic per-drug counts); broth-microdilution filter eliminates all strains (raise)

### Step 7: Foundation model wrappers (Evo + DNABERT-2)
Files: dna_decode/models/foundation.py, tests/test_models_foundation.py
Depends on: Step 1

**What changes:**
- `dna_decode/models/foundation.py` — abstract base `FoundationModel` with `embed(sequence: str) -> np.ndarray` returning a fixed-dimension sequence embedding. Concrete: `EvoModel` (loads `together/evo-1-131k-base` via HuggingFace), `DNABERT2Model` (loads `zhihan1996/DNABERT-2-117M`). `embed_batch(sequences: list[str]) -> np.ndarray` batches efficiently with auto-truncation/sliding-window for sequences exceeding model context. `model_factory(name: str) -> FoundationModel` dispatches.
- `tests/test_models_foundation.py` — uses model `MockFoundationModel` (returns hash-based deterministic embeddings) to test the wrapper interface contract; the real Evo/DNABERT-2 loaders are tested in Step 15's smoke test (too heavy for unit tests)

**Key details:**
- Lazy model load: the model is only loaded when `embed` or `embed_batch` is called the first time (avoid loading on import)
- Device selection: `torch.device("cuda" if available else "cpu")`; configurable via env var `DNA_DECODE_DEVICE`
- Sliding-window for sequences > model context: stride = context-length / 2, mean-pool overlapping windows
- The two models have different embedding dimensions (Evo 4096, DNABERT-2 768); downstream code must respect this — embed return shape is `(n_windows_or_1, embedding_dim)`

**Test strategy:**
- Unit: instantiate MockFoundationModel, verify embed returns correct shape, verify embed_batch is faster than serial embed (basic perf sanity), verify sliding-window over a 10K-bp synthetic sequence with 1K context produces correct number of windows
- Edge cases: empty sequence (raises `ValueError`); single-window-fitting sequence; sequence shorter than context

### Step 8: Embedding cache (HDF5)
Files: dna_decode/models/cache.py, tests/test_models_cache.py
Depends on: Step 1, Step 7

**What changes:**
- `dna_decode/models/cache.py` — `EmbeddingCache` class wrapping an HDF5 file at `cache_dir/embeddings_<model_name>.h5`. Methods: `put(strain_id: str, gene_id: str, embedding: np.ndarray) -> None`, `get(strain_id: str, gene_id: str) -> np.ndarray | None`, `has(strain_id, gene_id) -> bool`, `list_strains() -> list[str]`, `bulk_get(pairs: list[tuple[str, str]]) -> np.ndarray` for fast retrieval during training. `populate(model: FoundationModel, strain_genomes: dict[str, Path], annotations: dict[str, AnnotationTable]) -> None` runs end-to-end: for each strain × each gene, compute embedding if not cached.
- `tests/test_models_cache.py` — uses MockFoundationModel + in-memory HDF5 file

**Key details:**
- HDF5 layout: `/strains/<strain_id>/<gene_id>` datasets containing 1D embedding arrays; per-strain `attrs` with metadata (model name, model version, embedding dim, timestamp)
- Concurrent write: HDF5 SWMR mode for parallel population (or batch-then-write pattern, simpler — populate in single thread, parallelism via embed_batch in Step 7)
- File-format upgrade: if existing cache has different `model_version` attr, refuse to write (force user to delete + repopulate); v0.2 may add migration
- Bulk get returns a stacked `(n, embedding_dim)` array, NaN-fills missing pairs

**Test strategy:**
- Unit: put → get round-trip; bulk_get with mix of cached and missing pairs; verify version-mismatch refusal
- Edge cases: empty cache; corrupted HDF5 file (graceful exception); concurrent open + read

### Step 9: Baseline classifiers
Files: dna_decode/models/classifiers.py, tests/test_models_classifiers.py
Depends on: Step 7, Step 8

**What changes:**
- `dna_decode/models/classifiers.py` — `train_xgboost_classifier(X: np.ndarray, y: np.ndarray, drug_name: str, params: XGBParams | None = None) -> TrainedClassifier` fits XGBoost on (n_strains, embedding_dim) features → binary R/S label. `predict_proba(model: TrainedClassifier, X: np.ndarray) -> np.ndarray` returns calibrated probabilities. `feature_importance(model: TrainedClassifier) -> np.ndarray` returns embedding-dim-aligned importances.
- `tests/test_models_classifiers.py` — synthetic feature/label data (sklearn `make_classification`); verifies training + prediction + importance shapes

**Key details:**
- Default XGBParams: `n_estimators=200, max_depth=6, learning_rate=0.1, subsample=0.8, eval_metric='auc'`. These are pre-Step-10 placeholders; Step 10's CV harness will tune.
- Probability calibration via `sklearn.calibration.CalibratedClassifierCV` (sigmoid method)
- Per-gene embedding aggregation: v1 uses mean-pooling of all gene embeddings within a strain into a single fixed-dim feature vector. Alternative gene-specific feature engineering deferred to Step 11 (mutagenesis) / Phase 2.

**Test strategy:**
- Unit: train on synthetic 100-sample 256-dim data, verify model predicts above-random on training set, verify importance vector shape matches feature dim
- Edge cases: all-one-class labels (raises clear error); ≤10 samples (raises minimum-samples error)

### Step 10: Evaluation harness (LOSO + LOMO + Mash/ANI + clade-only baseline + per-clade reporting)
Files: dna_decode/eval/cv.py, dna_decode/eval/metrics.py, dna_decode/eval/phylogeny.py, dna_decode/eval/clade_baseline.py, tests/test_eval.py, tests/test_eval_phylogeny.py
Depends on: Step 1

**Strengthened per post-tech-plan brainstorm M2.** MLST-only phylogeny control is low-resolution and can give false comfort when most MLSTs are label-pure. Phase 1 ships with 3 phylogeny controls.

**What changes:**

- `dna_decode/eval/cv.py` — 3 CV strategies:
  - `leave_one_strain_out_cv(features, labels, strain_ids, train_fn, predict_fn) -> CVResult` — baseline (LOSO).
  - `leave_one_mlst_out_cv(features, labels, strain_ids, mlst_assignments, train_fn, predict_fn) -> CVResult` — holds out entire MLST sequence-types.
  - `leave_one_clade_out_cv(features, labels, strain_ids, clade_assignments, train_fn, predict_fn) -> CVResult` — holds out entire Mash/ANI-distance phylogenetic clusters (higher resolution than MLST). **Phase 1 primary control.**

- `dna_decode/eval/phylogeny.py` — `compute_mash_distances(strain_genomes: dict[str, Path]) -> DistanceMatrix` runs Mash on each strain genome (`mash sketch` + `mash dist`) to produce pairwise ANI-like distances. `cluster_by_ani(distance_matrix: DistanceMatrix, threshold: float = 0.02) -> dict[str, int]` runs hierarchical clustering (scipy `linkage` + `fcluster`) at the configured ANI threshold; returns strain → cluster-id mapping. Default 0.02 ANI ≈ 98% genome identity (sub-species-level groupings).

- `dna_decode/eval/clade_baseline.py` — `train_clade_only_classifier(strain_ids, clade_assignments, labels, drug) -> CladeOnlyModel` trains a classifier on JUST one-hot clade-membership features (no sequence embeddings). `predict_clade_only(model, strain_ids, clade_assignments) -> np.ndarray`. The clade-only model is the **null baseline**: if the embedding model's per-clade-holdout AUROC ≤ clade-only baseline AUROC + 0.05, the embedding model is learning clade signature, NOT mechanistic resistance signal.

- `dna_decode/eval/metrics.py` — `compute_metrics(y_true, y_score) -> Metrics` (AUROC, AUPRC, accuracy at 0.5 threshold, F1, Brier score, calibration ECE). `compute_attribution_precision(predicted_loci, known_loci, k) -> float`. `compute_per_clade_metrics(cv_result: CVResult, clade_assignments) -> dict[clade_id, Metrics]` aggregates per-held-out-clade. `compute_within_mlst_permutation_control(cv_result, labels, mlst_assignments) -> ControlReport` — within-MLST label-shuffle control per H8.

- `tests/test_eval.py` + `tests/test_eval_phylogeny.py` — synthetic data with known labels + known clades + known importance loci

**Key details:**
- Phase 1 primary CV strategy is **leave-one-Mash-clade-out**, not LOSO or LOMO. LOSO + LOMO ship as secondary metrics for reference.
- Mash CLI (`mash`) dependency: add to `pyproject.toml` install instructions; pinned version. Fallback: pure-Python `pyani` if Mash binary unavailable.
- Clade-only baseline output should be reported alongside the embedding-model output in the same `CVResult`. The platform's "validation gate" is: embedding-model per-clade AUROC ≥ clade-only AUROC + 0.10 on ≥75% of held-out clades.
- Per-clade reporting catches the "0.90 on 8 clades / 0.45 on 2 clades hides clade-specific failure" pattern Codex flagged.
- AUROC, AUPRC, ECE use sklearn implementations.
- Attribution-precision: predicted loci ranked by absolute ISM importance (NOT Captum, per C1); known loci per drug: ciprofloxacin = {gyrA, gyrB, parC, parE, qnr-family, aac(6')-Ib-cr}; ceftriaxone = {CTX-M family, SHV family, AmpC overproduction loci}; tetracycline = {tetA, tetB, tetM, tetW family}; precision = (predicted top-K ∩ known) / K.

**Test strategy:**
- Unit: synthetic 30-strain dataset across 6 clades; verify LOSO produces 30 folds, LOMO produces some number, leave-one-clade-out produces 6 folds; clade-only baseline trains + scores; per-clade metrics aggregate correctly.
- Mash subprocess wrapped + mocked in tests (CI doesn't run real Mash binary on synthetic fixtures).
- Edge cases: single-class folds (warn + NaN); single-clade dataset (leave-one-clade-out raises clear error); zero-distance clade (handled by clustering threshold).

### Step 11: In-silico mutagenesis (gene-level + saturation)
Files: dna_decode/interp/mutagenesis.py, tests/test_interp_mutagenesis.py
Depends on: Step 7, Step 8

**What changes:**
- `dna_decode/interp/mutagenesis.py` — Phase 1's **sole attribution mechanism** (Captum/IG deferred to Phase 2 per post-tech-plan brainstorm critique C1: XGBoost on mean-pooled embeddings is non-differentiable end-to-end, and discrete tokenizer input doesn't yield clean position-level IG).
  - `gene_level_mutagenesis(model: FoundationModel, classifier: TrainedClassifier, strain_genome: Path, annotations: AnnotationTable, drug: str) -> GeneEffectTable` performs in-silico knockout of each gene (replace gene region with N's or shuffled background), recomputes embedding, recomputes prediction, records prediction-delta per gene.
  - `saturation_mutagenesis(model, classifier, gene_id: str, sequence: str, alt_bases: tuple = ("A","C","G","T")) -> PositionEffectTable` does single-base substitutions across the gene length (length × 3 alternative bases × prediction-delta). Always run on the top-K (default K=20) genes ranked by gene-level mutagenesis — full-genome saturation is too expensive.
  - `motif_recovery(saturation_table: PositionEffectTable, known_motifs: list[Motif]) -> MotifRecoveryReport` aligns saturation-mutagenesis high-impact windows against known motifs from RegulonDB / JASPAR-like sources and reports overlap.

- `tests/test_interp_mutagenesis.py` — MockFoundationModel + synthetic classifier with known importance loci

**Key details:**
- Gene-level is cheaper (one prediction per gene); saturation runs on top-K genes only.
- Use parallel batched inference via `embed_batch` from Step 7.
- Output `GeneEffectTable` columns: `gene_id, locus_tag, prediction_delta, baseline_probability, knockout_probability` sorted by absolute delta desc.
- `PositionEffectTable` columns: `gene_id, position, ref_base, alt_base, prediction_delta`; pivot for visualization downstream.
- ISM is theoretically sound for any predictor (no differentiability required), works for XGBoost + frozen-embedding pipeline, and is the published consensus for sequence-bio interpretability — defensible Phase 1 ground truth.

**Test strategy:**
- Unit: synthetic mock model where gene X is known to drive prediction; gene-level mutagenesis ranks gene X first; saturation-mutagenesis within gene X identifies the seeded motif positions.
- Edge cases: empty annotations (returns empty table); gene with zero coverage (skip + warn); saturation-mutagenesis on gene longer than model context (sliding-window aggregation).

### Step 12: [REMOVED] Captum attribution wrapper — deferred to Phase 2

**Why removed (post-tech-plan brainstorm C1):**
- Step 9's XGBoost classifier on mean-pooled embeddings is non-differentiable → Captum IG cannot backprop through it.
- Foundation-model attribution via IG requires a differentiable end-to-end path from tokens → embeddings → classifier → score. The plan as originally written did not provide one.
- Discrete BPE/k-mer tokenizer inputs don't give clean per-nucleotide IG even with a differentiable head — outputs are token-level and need back-resolution.

**Phase 2 backlog entry:**
- Add a small differentiable MLP head (e.g., 2-layer MLP) trained on the same frozen embeddings as Step 9's XGBoost — runs as an alternate attribution head, NOT a replacement for the XGBoost benchmark classifier.
- Captum IG runs through MLP-head + foundation-model embedding path with end-to-end gradient.
- Phase 2 cross-validation: IG attribution vs ISM attribution on the same top-K genes — if they disagree systematically, the model is exposing different signal under different attribution probes.

### Step 13: Genome-browser visualization
Files: dna_decode/viz/browser.py, tests/test_viz_browser.py
Depends on: Step 3, Step 11

**What changes:**
- `dna_decode/viz/browser.py` — `render_ism_browser(annotations: AnnotationTable, gene_effects: GeneEffectTable, position_effects: PositionEffectTable, drug_name: str, output_path: Path) -> Path` generates a `pygenometracks`-compatible ini config + BED tracks (annotation track, gene-level-delta heatmap track, position-level-delta bedgraph track for the top-K saturation-mutagenesis genes), invokes pygenometracks CLI, returns PNG path. `render_strain_comparison(annotations, gene_effects_per_strain: dict[str, GeneEffectTable], region: GenomicRegion, output_path: Path) -> Path` stacks multiple-strain gene-level-effect tracks for visual comparison.
- `tests/test_viz_browser.py` — verifies ini + BED files are syntactically valid (pygenometracks parses them); skips actual PNG render in CI (heavy)

**Key details:**
- pygenometracks is a CLI tool; call via `subprocess.run` with explicit env (the `make_tracks_file` and `pyGenomeTracks` binaries)
- Track config is a Python-generated `.ini` written to a tempfile
- BED-format conversion: ISM position-level prediction-delta → bedgraph (chrom, start, end, value); gene-level prediction-delta → BED + score column for heatmap rendering
- Annotation track uses the GFF3 directly via pygenometracks' built-in GFF support
- Inputs come from Step 11's `GeneEffectTable` + `PositionEffectTable` (NOT Captum, which is Phase 2 only per C1)

**Test strategy:**
- Unit: generate an ini config + bed file from a synthetic 1Kbp ISM saturation map + 5-gene annotation, verify the ini is parseable, verify the BED is sorted + non-overlapping
- Edge cases: empty ISM tables (single zero-valued track); pygenometracks-not-installed (raises clear setup error)

### Step 14: CLI entry points
Files: scripts/ingest.py, scripts/train.py, scripts/predict.py, scripts/attribute.py
Depends on: Step 2, Step 4, Step 5, Step 6, Step 7, Step 8, Step 9, Step 10, Step 11

**What changes:**
- `scripts/ingest.py` — `python -m scripts.ingest --target-strains 200 [--force-refresh]`: orchestrates Steps 2, 4, 5, 6 (download reference + CARD + AMRFinder + AST + pan-genome strain genomes); writes a manifest at `data/processed/ingest_manifest.json`.
- `scripts/train.py` — `python -m scripts.train --drug fluoroquinolone --model evo [--cv-strategy loso|lomo]`: loads ingested data, computes embeddings (Steps 7 + 8), trains classifier (Step 9), runs CV (Step 10), prints metrics + saves trained model to `data/processed/models/<drug>_<model>.pkl`.
- `scripts/predict.py` — `python -m scripts.predict --fasta <path> --model <path> --output <report.json>`: loads a trained model, embeds the input FASTA, predicts, writes a JSON report (`{drug, probability, calibration_method, model_metadata}`).
- `scripts/attribute.py` — `python -m scripts.attribute --strain-id <id> --drug <drug> [--saturation-top-k 20]`: loads trained model, runs gene-level + saturation ISM (Step 11), prints top-20 loci + writes a position-level effect table for the top-K genes.

**Key details:**
- All CLIs use `argparse` (stdlib, no click/typer dep for v1)
- All CLIs write to `data/processed/` (gitignored)
- All CLIs read `config/datasources.yaml` for paths + cache locations
- Error handling: each CLI catches `dna_decode` package exceptions and prints user-friendly messages

**Test strategy:**
- Unit: import each script's main function, call with mocked deps, verify exit code + output structure
- Smoke: full end-to-end run via Step 15

### Step 15: Smoke pipeline + fixtures
Files: scripts/smoke_pipeline.py, tests/fixtures/ecoli_mini/genome.fna, tests/fixtures/ecoli_mini/annotations.gff3, tests/fixtures/card_mini.json, tests/fixtures/amrfinder_mini.tsv, tests/fixtures/bvbrc_ast_mini.tsv, tests/fixtures/example.gff3, tests/fixtures/example.fna, tests/test_smoke.py
Depends on: Step 14

**What changes:**
- `scripts/smoke_pipeline.py` — end-to-end run on the bundled mini-fixture: ingest synthetic 5-strain dataset, embed via MockFoundationModel (no real Evo download), train a logistic-regression placeholder for fluoroquinolones, run CV, run attribution, render a browser PNG. Completes in <60s on CPU. Used both as a smoke test AND as user-facing documentation of the full pipeline.
- `tests/fixtures/ecoli_mini/` — 5 synthetic 10-Kbp E. coli-like sequences with 10 genes each, designed so that gene "synth_gyrA" carries an artificial "resistance signal" (a specific motif present in resistant strains)
- `tests/fixtures/{card_mini.json, amrfinder_mini.tsv, bvbrc_ast_mini.tsv, example.gff3, example.fna}` — small fixture files referenced by Steps 2-5 unit tests
- `tests/test_smoke.py` — invokes smoke_pipeline.py, asserts exit code 0 + all output files exist + AUROC on the synthetic resistance signal is ≥0.85 (model SHOULD ace this constructed signal)

**Key details:**
- Smoke pipeline uses `MockFoundationModel` (deterministic hash-based embedding) to avoid the Evo download cost (~10GB) in CI
- The synthetic resistance signal is a 50-bp motif inserted into "synth_gyrA" of resistant strains and absent from susceptible strains — a learnable signal that validates the end-to-end harness without requiring real biology
- Fixture file sizes <100KB each (acceptable to commit)

**Test strategy:**
- Smoke: full pipeline runs in <60s on CPU; AUROC ≥0.85 on synthetic signal; attribution top-1 locus = synth_gyrA
- This is the canonical regression test for any future change

### Step 16: Documentation
Files: docs/ARCHITECTURE.md, docs/HOW_TO_ADD_ORGANISM.md, docs/HOW_TO_RUN.md, README.md
Depends on: Step 15

**What changes:**
- `docs/ARCHITECTURE.md` — 1-2 page architecture overview with the diagram from the conversation, the 4-layer pipeline (ingest → embed → predict → interpret), and pointers to each module
- `docs/HOW_TO_ADD_ORGANISM.md` — playbook for extending beyond E. coli (which steps need new code vs config-only changes, what data sources to find)
- `docs/HOW_TO_RUN.md` — end-user setup: install + ingest + train + predict + attribute commands
- `README.md` (UPDATE Step 1's stub) — replace placeholder with: 1-paragraph project description + quickstart pointing at HOW_TO_RUN.md + status banner ("v0.1: Phase 1 MVP — antibiotic resistance prediction on E. coli")

**Key details:**
- README.md is the only file overlap with Step 1 — Step 1 wrote a stub; Step 16 finalizes
- ARCHITECTURE.md mirrors the rca_engine/CLAUDE.md "what's here" + "architecture" sections pattern
- No external docs site for v1; everything lives in `docs/` markdown

**Test strategy:**
- Unit: `tests/test_docs.py` (optional) — verify each docs file exists + has expected headings (smoke test against doc drift)

### Step 17: Comparative model benchmarking leaderboard
Files: scripts/leaderboard.py, dna_decode/eval/leaderboard.py, tests/test_leaderboard.py
Depends on: Step 7, Step 8, Step 9, Step 10, Step 14

**Added per Adjustment C — scope-limited per Codex's partial-agree pushback.** NOT an open-ended model bakeoff; a fixed leaderboard against the SAME cohort + cache + CV protocol.

**What changes:**
- `dna_decode/eval/leaderboard.py` — `run_leaderboard(cohort: StrainCohort, drugs: list[str], models: list[str] = ["evo", "dnabert2", "nucleotide-transformer", "gena-lm"], cv_strategy: str = "leave_one_clade_out") -> LeaderboardReport` orchestrates: for each model in `models`, ensure embeddings are cached for the cohort (compute if missing via Step 8's `populate`), train Step-9 XGBoost per drug, run Step-10 CV, compute attribution-precision via Step-11 ISM, aggregate into a single ranked table.
- `scripts/leaderboard.py` — `python -m scripts.leaderboard [--drugs ...] [--models ...]`: writes `data/processed/leaderboard.md` with per-model × per-drug AUROC + AUPRC + clade-only baseline gap + attribution-precision + embedding compute time. Phase 1 success: at least one model achieves ≥0.85 AUROC on ciprofloxacin AND ≥0.10 gap vs clade-only AND ≥0.6 attribution-precision @ K=20.
- `tests/test_leaderboard.py` — mocked MockFoundationModel × N models; verifies report structure + ranking logic

**Key details:**
- **Scope-limited per Codex:** runs against the SAME pilot-validated cohort (Step 0.5) + same embedding-cache protocol (Step 8) + same CV strategy (Step 10). Adding a new model = adding to the `models` list; new cohorts require explicit rerun.
- Embedding compute is the bottleneck: Evo embeddings on 450 strains × ~5K genes per strain ≈ 2.25M sequences. Pre-compute once, cache (Step 8), reuse across all CV folds.
- Real-data leaderboard runs are expensive (4 models × ~6 hours embedding each); Phase 1 ships with leaderboard initially run for Evo + DNABERT-2 only; NT + GENA-LM added incrementally
- Compare-and-rank logic considers: primary metric (AUROC per drug), secondary (attribution-precision), tertiary (embedding compute cost). Tie-breaking documented in the leaderboard report.

**Test strategy:**
- Unit: with 3 MockFoundationModels of differing simulated performance, verify ranking produces correct order; verify report has all required columns
- Integration (in Step 15 smoke): leaderboard.py on synthetic 5-strain fixture across 2 mock models completes <30s

## Execution Preview

```
Wave 0 (1 step):     Step 1 — Project bootstrap
Wave 1 (7 parallel): Step 0.5 — Real-data pilot gate
                       Step 2 — RefSeq downloader
                       Step 3 — Annotation parser
                       Step 4 — Resistance database loaders
                       Step 5 — AST phenotype data loader
                       Step 7 — Foundation model wrappers
                       Step 10 — Evaluation harness
Wave 2 (2 parallel): Step 6 — Strain/AST cohort catalog (drug-first)
                       Step 8 — Embedding cache
Wave 3 (2 parallel): Step 9 — Baseline classifiers
                       Step 11 — In-silico mutagenesis (gene-level + saturation)
Wave 4 (1 step):     Step 13 — Genome-browser visualization
Wave 5 (1 step):     Step 14 — CLI entry points
Wave 6 (2 parallel): Step 15 — Smoke pipeline + fixtures
                       Step 17 — Comparative model benchmarking leaderboard
Wave 7 (1 step):     Step 16 — Documentation
```

Step 12 (Captum attribution wrapper) **removed from Phase 1** per post-tech-plan brainstorm C1 — deferred to Phase 2 with a dedicated differentiable MLP head. ISM (Step 11) is Phase 1's sole attribution mechanism.

Critical path: Step 1 → Step 7 → Step 8 → Step 9 → Step 14 → Step 15/17 → Step 16 (7 waves)
Max parallelism: 7 agents (Wave 1)

Note: dna_decode is a fresh project — parallel execution requires `git init` + a configured remote. If not set up, /execute-plan falls back to sequential mode.

## Risk Flags

- **No git remote yet for `dna_decode/`.** /execute-plan parallel mode requires a remote. Either run `git init` + create an empty GitHub/local-bare repo + `git remote add origin ...` before /execute-plan, or accept sequential execution.
- **Foundation model download cost.** Evo is ~10GB+ on first load; first-time `train.py` invocation will be slow. Mitigation: smoke test uses MockFoundationModel; first real Evo load is a one-time cost.
- **GPU memory.** Evo 7B at full precision needs ~28GB VRAM. Use 4-bit quantization (bitsandbytes) on a single RTX 4090 (24GB), or rent an A100 (40GB) on-demand for Phase 1 calibration runs. Already covered in the Pending Decisions on the ledger.
- **BV-BRC AST data labeling inconsistency.** Real-world AST data has measurement-method variance (disk-diffusion vs broth-microdilution; CLSI vs EUCAST breakpoints) → label noise. Phase 1 filters to broth-microdilution only (Step 5) per failure-mode #4; Phase 2 should add method-aware reweighting + EUCAST-vs-CLSI breakpoint reconciliation.
- **Attribution-precision threshold (0.6 at K=20 for ciprofloxacin; 0.4 for ceftriaxone) is calibrated against literature, not vetted experimentally.** This is an honest ML eval, not a biological claim. β-lactam attribution-precision is set lower because resistance signal is distributed across plasmid β-lactamases + porin mutations + efflux pumps. Reflected in project ledger H7 + H9.
- **Step 8 file overlap with Step 7.** Both touch `dna_decode/models/`. Wave 2 has Step 8 depending on Step 7 already — no overlap conflict at file level (different files within the same package).
- **Phylogeny-control sufficiency.** Phase 1 ships with 3 controls (LOSO, LOMO, leave-one-Mash-clade-out) + a clade-only baseline + within-MLST permutation. A "passes all 4" model is still not proven causal — only correlational. Frame outputs as "associated with" not "causes" per the spec's causality-stance.
- **Mash CLI external binary dependency.** Step 10's `phylogeny.py` shells out to `mash sketch` + `mash dist`. Add Mash install instructions to README. Fallback: pure-Python `pyani` (slower but no binary needed).
- **Plasmid vs chromosomal encoding masking.** Step 6's cohort filter uses AMRFinder's plasmid-vs-chromosome annotation only when present; assemblies with no annotation get included with `plasmid_localization=unknown`. Phase 1 attribution interpretation must be careful: a high gene-level attribution on a plasmid-borne CTX-M means "this strain has the plasmid" — not "this DNA sequence causes resistance independently."
- **Annotation drift across Bakta versions.** Bakta version pinned in `pyproject.toml` (Step 1). Reannotation with a new Bakta version requires recomputing the embedding cache. Phase 2 should add a Bakta-version field to cache metadata and refuse-on-mismatch (matches Step 8's existing model-version refuse-on-mismatch pattern).

**Restructuring applied (cumulative):**
1. Step 5 (AST data) moved earlier in the dependency chain so Step 6 (cohort) can prioritize strains with AST labels.
2. Step 12 (Captum attribution) REMOVED from Phase 1 per post-tech-plan brainstorm C1; deferred to Phase 2 with a dedicated differentiable MLP head.
3. Step 11 (ISM) extended with nucleotide-level saturation-mutagenesis on top-K genes (was gene-level only).
4. Step 6 renamed from "Pan-genome strain catalog" → "Strain/AST cohort catalog (drug-first)" per M1; selection logic rewritten as drug-first per-drug coverage; "complete-circle only" filter replaced with assembly-quality threshold + AMRFinder plasmid/chromosome metadata.
5. Step 10 strengthened with Mash/ANI clustering + clade-only baseline + per-held-out-clade reporting per M2.
6. Step 0.5 (Real-data pilot gate) added at Wave 1 per T2 — soft gate that estimates per-drug isolate counts after all filters before the full pipeline runs.
7. Step 17 (Comparative model benchmarking leaderboard) added at Wave 6 per Adjustment C — scope-limited to the pilot-validated cohort + same CV protocol; not an open-ended bakeoff.
8. Phase 1 drug list specified: ciprofloxacin + ceftriaxone (specific β-lactam agent per Codex Adjustment-A refinement; CTX-M family attribution validation target) + tetracycline.

## Verification

After all 16 steps (Step 0.5 + Steps 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17 — Step 12 removed) merge:

1. `uv run pytest tests/ -v` — all unit tests pass.
2. `uv run python -m scripts.pilot_gate` — real-data pilot gate completes; per-drug counts ≥150 for cipro + ceftriaxone + tetracycline after broth-microdilution + assembly-quality filters. GO recommendation.
3. `uv run python scripts/smoke_pipeline.py` — synthetic-fixture end-to-end completes in <60s with AUROC ≥0.85 on the constructed signal.
4. `uv run python -m scripts.ingest --drugs ciprofloxacin,ceftriaxone,tetracycline` — real ingestion of NCBI + CARD + AMRFinder + BV-BRC + cohort-validated strains, completes in <4 hours on home internet.
5. `uv run python -m scripts.train --drug ciprofloxacin --model evo` — full real-data training on cipro resistance, completes overnight on a single GPU, prints **leave-one-Mash-clade-out CV AUROC ≥0.80** (Phase 1 SLO; ≥0.85 is the target).
6. `uv run python -m scripts.train --drug ciprofloxacin --model evo --include-clade-baseline` — also reports clade-only baseline AUROC; embedding-model AUROC must be ≥0.10 above clade-only on ≥75% of held-out clades (per-clade reporting Step 10).
7. `uv run python -m scripts.train --drug ceftriaxone --model evo` — same as #5 for ceftriaxone; attribution-precision target is ≥0.4 (lower than cipro per H9).
8. `uv run python -m scripts.train --drug tetracycline --model evo` — same as #5 for tetracycline.
9. `uv run python -m scripts.attribute --strain-id GCF_000005845.2 --drug ciprofloxacin` — produces top-20 loci with gyrA + parC in top 5 (validates H7 cipro target).
10. `uv run python -m scripts.attribute --strain-id <CTX-M-positive-strain-id> --drug ceftriaxone` — produces top-20 loci with CTX-M family in top 10 (validates H7 ceftriaxone target).
11. `uv run python -m scripts.leaderboard --drugs ciprofloxacin,ceftriaxone,tetracycline --models evo,dnabert2` — runs the Phase 1 leaderboard; writes `data/processed/leaderboard.md` ranking models by primary metric (AUROC per drug).
12. Visual: render ISM-attribution browser PNGs for each of the 3 drugs on the reference genome; eyeball that high-attribution regions overlap with known resistance-gene neighborhoods.

**Phase 1 ships** when verifications 1, 2, 3, 4, 5, 6, 9 pass. Verifications 7, 8, 10, 11, 12 are nice-to-have for the Phase 1 launch artifact; gaps are documented in Phase 2 backlog. Verification 6's clade-only-baseline gap is the **mechanistic-signal validation gate** — failing this means the model is learning phylogeny, not resistance biology.
