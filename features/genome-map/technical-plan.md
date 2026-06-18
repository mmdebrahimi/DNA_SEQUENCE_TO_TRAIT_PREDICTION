# Genome Map v1 — Bakta Honesty Report — Technical Plan (v2, brainstorm-hardened)

> Single-genome evidence-tiered function/QC map (Bakta annotation re-tiered for honesty + AMRFinder determinant call-outs behind a hard join-quality gate + a DB-labelled unknown rate), bacterial-AMR-only v1, with a prevent-wrong-inference gate.

## Lens status
Inputs: `features/genome-map/feature.md` (status spec-ready) + its `## Repo grounding` (intake/design/spec/probe/brainstorm) + the 2026-06-17 pre-exec /brainstorm (this conversation, 2 Codex rounds — C1-C5 + M1-M2 folded here).
Degradations: repo-index unconfigured (direct reading); Bakta/AMRFinder GFF+main.tsv field sets [unverified] until the Step-2 tool-surface manifest; no spec.md (spec in feature.md); sentrux/project-rules/DESIGN absent.

## Problem Statement
The north star ("DNA → what its parts do") is achievable only as an evidence-tiered MOLECULAR-FUNCTION / QC map, not a learned phenotype predictor (closed negative). v1 = the "Bakta honesty report": a personal exploration/QC tool that RE-TIERS Bakta's own annotation for honesty + overlays the existing AMRFinder determinant cells, on the ALREADY-INSTALLED Bakta(db-light) + AMRFinder + the existing GFF parser/runner. Phenotype claims appear ONLY behind the validated-determinant wall; the unknown rate carries its db-light tooling caveat.

The pre-exec /brainstorm folded in (all [grounded]): **C1** AMRFinder must be RUN (reuse `_run_amrfinder`) + an explicit organism; **C2** Bakta's embedded `##FASTA` block breaks `parse_gff3`/`load_annotation_table` → a SHARED `##FASTA`-stripping loader for both the Bakta-run and the offline provided-GFF path; **C3** the determinant→feature join needs a `DeterminantHit` from raw `main.tsv` (coords/protein) + a HARD join-quality gate (symbol-fallback joins are visible but EXCLUDED from `determinant-phenotype` + G1; spike NO-GOs if all-fallback); **C4** the map retains raw fields + `classification_reason` so G1 computes from the map alone; **C5** a Wave-0/1 tool-surface manifest GATES the tier + overlay steps (not parallel/guessed); **M1** TB is DROPPED from the v1 spike (VCF-vs-GFF contract mismatch); **M2** v1 accepts Bakta-compatible GFF only.

Non-goals: TB/fungal overlay in the v1 spike (schema future-compat only); hmmer/Pfam/eggNOG tiers; the pathway-module tier; arbitrary GFF/GBK normalization; a visual browser; multi-genome comparison; ANY edit to the frozen AMR surface; a learned/embedding model.

## Codebase Context
- `dna_decode/data/annotations.py::parse_gff3` (`:63`) + `load_annotation_table` (`:163`, dispatches `.gff3`→`parse_gff3` at `:173`) — the existing parser, but it does NOT strip Bakta's embedded `##FASTA` block → choked output. `scripts/pathotype_laptop_pipeline.py::parse_bakta_gff3` (`:70`) strips `##FASTA` then calls `parse_gff3`; `bakta_annotate` (`:46`) is the existing Bakta runner.
- `scripts/drug_mechanism_audit.py::_run_amrfinder` (`:83`, sig `(fasta, out_dir, organism="Escherichia") -> (main_tsv, mutations_tsv)`; `-O` is organism-specific for QRDR) + `dna_decode/amr/cli.py::_run_amrfinder_for_genome` (`:215`) — the existing AMRFinder runner.
- `dna_decode/eval/amr_rules.py::call_resistance` / `cipro_determinants_from_main` (`:68`) — return determinant dicts with element symbol/name/class, **NOT coordinates**. AMRFinder `main.tsv` uses the `Gene symbol` column (`concordance/core.py:71`, `resistance_db.py:133`).
- `dna_decode/organism_rules/tb_amr.py` + `tb_vcf.py` (`CHROM=NC_000962.3`) — TB scores VCF variants vs H37Rv COORDINATES (a different overlay contract; out of v1 spike).
- `tools/docker_runner.py::run` — docker primitive. Test infra: pytest under `tests/`; offline-safe (BLAST/vf_diff) degradation pattern.

### Reusable-Code Survey
- `scripts/pathotype_laptop_pipeline.py` (`parse_bakta_gff3` ##FASTA-strip + `bakta_annotate`) — Step 1 ingest + annotate.
- `scripts/drug_mechanism_audit.py::_run_amrfinder` + `dna_decode/amr/cli.py::_run_amrfinder_for_genome` — Step 1 AMRFinder run.
- `dna_decode/data/annotations.py::parse_gff3` — the inner GFF parse (wrapped by the ##FASTA strip).
- `dna_decode/eval/amr_rules.py` + `dna_decode/organism_rules/` — determinant catalogs (read-only; the overlay parses raw main.tsv itself for coords).
- None — searched: graphify-out/GRAPH_REPORT.md (absent), src/lib/utils dirs (absent). Grounding entries consulted: probe @ + brainstorm @ in `features/genome-map/feature.md`.

## Pre-Change Baseline
- No genome-map exists. The determinant cells give a genome-level R/S verdict; nothing tiers a genome's features or surfaces the unknown rate.
- **Bakta + AMRFinder annotation are installed-but-UNPROVEN end-to-end on this host** (Bakta smoke deferred 2026-05-15; Docker-mount-corruption history) — the Step-2 tool-surface manifest is the feasibility gate.
- Differentiation baseline (what the map must BEAT, per G1): the raw Bakta GFF3 `product` wording read directly.
- Frozen AMR surface (`dna_decode/eval/amr_rules.py` + `calibrated_amr_rules.json`) byte-frozen — leak-guard target.

## Verification Signal
- **Tool surface (Step 2):** `wiki/genome_map_tool_surface_<date>.json` lists the real Bakta GFF fields + product-vocabulary examples + the AMRFinder `main.tsv` headers (incl. whether protein-id/contig/coords are present), OR a documented `BAKTA_ANNOTATION_BLOCKED` / `AMRFINDER_BLOCKED`. The tier + overlay steps consume it.
- **Map (Steps 3-5):** a per-feature JSON map + flat table; one primary tier ∈ {determinant-phenotype, curated-molecular-function, homology-only-hypothesis, unknown}; raw fields (`raw_product`/`raw_gene_symbol`/`raw_locus_tag`/`raw_feature_type`/`source_tool`/`classification_reason`/`secondary_evidence[]`) retained; phenotype field non-empty ONLY on determinant-phenotype features that cleared a HIGH-confidence join (symbol-fallback excluded); top-level `unknown_under_bakta_db_light` + per-tier counts + join-quality counts (`n_main_rows/n_high_confidence_join/n_symbol_fallback/n_unjoined`).
- **Gate (Steps 6-7, make-or-break):** `wiki/genome_map_spike_verdict_<date>.md` reports, per 3 bacterial prototype genomes, the DB-labelled unknown rate + ≥3 prevent-wrong-inference G1 features (≥1 of type demote-homology/surface-determinant; "unknown rate exists" alone fails) + G2 spot-check + GO/NO-GO. **NO-GO if a genome's determinant joins are all symbol-fallback** (the gene-symbol-trap guard). Honest NO-GO/BLOCKED allowed.
- Frozen AMR surface byte-unchanged (leak-guard); offline path (provided Bakta-compatible GFF) emits via the SAME ##FASTA-stripping loader.

## Implementation Steps

### Step 1: Tool plumbing — shared ##FASTA-safe GFF loader + Bakta + AMRFinder runners
Files: dna_decode/genome_map/__init__.py, dna_decode/genome_map/ingest.py, dna_decode/genome_map/annotate.py, dna_decode/genome_map/amrfinder.py, tests/test_genome_map_ingest.py
Depends on: none

**What changes:**
- New `genome_map` package. **C2:** `ingest.py::load_genome_gff(path)` — the SHARED loader that strips Bakta's embedded `##FASTA` block (cf. `parse_bakta_gff3`) BEFORE `parse_gff3`; used by both the Bakta-run path AND the offline provided-GFF path (v1 = Bakta-compatible GFF only, **M2**).
- `annotate.py`: thin Bakta runner (adapt `bakta_annotate`; pinned `oschwengers/bakta:v1.11.4`, db-light, entrypoint quirk, D:-cache, skip-existing).
- **C1:** `amrfinder.py::run_amrfinder(fasta, out_dir, organism)` — reuse `drug_mechanism_audit._run_amrfinder`; organism is an EXPLICIT required input (no Bakta auto-detect in v1).

**Test strategy:**
- `load_genome_gff` strips a synthetic `##FASTA`-bearing GFF correctly + parses the records; a plain GFF is unaffected. Bakta/AMRFinder wrappers build correct docker args without running (mock `docker_runner.run`).

### Step 2: Tool-surface manifest (feasibility smoke + field/vocab/header inventory)
Files: scripts/genome_map_tool_surface.py, tests/test_genome_map_tool_surface.py
Depends on: Step 1

**What changes:**
- **C5:** run Bakta + AMRFinder on ONE genome (via Step 1) and INVENTORY: Bakta GFF fields + product-vocabulary examples (the real db-light wording) + AMRFinder `main.tsv` headers (presence of protein-id/contig/start/end). Emit `wiki/genome_map_tool_surface_<date>.json` — the artifact that GATES Steps 3 (tier vocab) + 4 (overlay join keys). A wedge → `BAKTA_ANNOTATION_BLOCKED` / `AMRFINDER_BLOCKED`, never a fake manifest.

**Test strategy:**
- Manifest builder unit-tested on cached/synthetic Bakta GFF + AMRFinder main.tsv → correct field/header inventory; a missing-tool path → BLOCKED status, not a crash.

### Step 3: Function/confidence tier classifier (vocab seeded from the manifest)
Files: dna_decode/genome_map/tiers.py, dna_decode/genome_map/tier_vocab.py, tests/test_genome_map_tiers.py
Depends on: Step 2

**What changes:**
- `tier_vocab.py`: the low-confidence wording patterns (`putative`/`probable`/`by similarity`/`domain-containing`/`uncharacterized`) + the `hypothetical`/empty unknown patterns — **seeded/validated against the Step-2 manifest's real product vocabulary** (not guessed).
- `tiers.py::classify_feature_tier(product, gene_symbol) -> (tier, classification_reason)` over the 4 v1 tiers (pathway deferred). Returns the REASON (the matched pattern) — feeds `classification_reason` (C4).

**Test strategy:**
- Table-driven over the manifest's real product examples + synthetic: named-gene+specific→curated; putative/DUF-domain→homology-hypothesis; hypothetical/empty→unknown; each returns its classification_reason.

### Step 4: Determinant overlay — DeterminantHit + hard join-quality gate
Files: dna_decode/genome_map/phenotype_overlay.py, tests/test_genome_map_overlay.py
Depends on: Step 2

**What changes:**
- **C3:** parse the raw AMRFinder `main.tsv` into a `DeterminantHit` retaining ALL columns (Gene symbol, class/subclass, + any protein-id/contig/start/end). Reuse `amr_rules.cipro_determinants_from_main` for the drug/class semantics + provenance, but the JOIN uses the raw coords/protein.
- Join each hit to a parsed feature by an explicit hierarchy: protein-id → coordinate-overlap (contig+start/end) → symbol-fallback, carrying `join_confidence` ∈ {high, coord, symbol_fallback}. **Symbol-fallback joins are VISIBLE but do NOT earn the `determinant-phenotype` primary tier and do NOT count for G1.** Report `n_main_rows/n_high_confidence_join/n_symbol_fallback/n_unjoined`.
- A determinant cell ABSTAIN/SUSPEND propagates to the feature (no forced call). Read-only; no frozen edit.

**Test strategy:**
- Synthetic main.tsv (with + without coords) → protein/coord joins = high-confidence determinant-phenotype; a coords-absent hit → symbol_fallback, NOT determinant-phenotype, surfaced with join_confidence; ABSTAIN propagates; the per-genome join-quality counts are emitted.

### Step 5: Map assembler + raw-field JSON/table + DB-labelled unknown rate
Files: dna_decode/genome_map/build_map.py, tests/test_genome_map_build.py
Depends on: Step 3, Step 4

**What changes:**
- Assemble: `load_genome_gff` features → primary tier = `determinant-phenotype` ONLY if a HIGH-confidence/coord determinant join hit it, else `classify_feature_tier`. **C4:** retain `raw_product`/`raw_gene_symbol`/`raw_locus_tag`/`raw_feature_type`/`source_tool`/`classification_reason`/`secondary_evidence[]` per feature. Emit JSON map + flat table.
- Top-level metrics: per-tier counts + `unknown_under_bakta_db_light` + the determinant-phenotype list + the join-quality counts (from Step 4).
- Phenotype WALL structural: `phenotype` field populated ONLY on high-confidence `determinant-phenotype` features. Offline-safe via the Step-1 shared loader.

**Test strategy:**
- Synthetic features + overlay → correct primary tiers (high-confidence determinant beats function beats homology beats unknown); symbol-fallback determinant does NOT get the phenotype field; raw fields + classification_reason retained; unknown rate DB-labelled; offline (provided GFF) path works.

### Step 6: Differentiation audit + G1/G2 gate (make-or-break)
Files: dna_decode/genome_map/gate.py, tests/test_genome_map_gate.py
Depends on: Step 5

**What changes:**
- `gate.py`: raw-vs-tiered audit using the retained `raw_product` + `classification_reason`. Classify prevent-wrong-inference cases: (a) a raw product taken as fact (`putative`/`by similarity`) DEMOTED to `homology-only-hypothesis`; (b) a HIGH-confidence determinant SURFACED as `determinant-phenotype` that raw lists as a plain gene. `evaluate_gate(map) -> {g1_features[], g1_pass (>=3, >=1 of type a/b), g2_spotcheck, all_joins_symbol_fallback, verdict}`.
- **Gate guard (C3):** `verdict = NO_GO` if `all_joins_symbol_fallback` (the gene-symbol-trap guard) OR g1 fails OR g2 fails.

**Test strategy:**
- Synthetic maps: 3 demoted-homology + 1 surfaced-determinant → GO; unknown-rate-only "win" → NO-GO (gameability guard); a map whose determinant joins are all symbol-fallback → NO-GO; G2 (no non-determinant phenotype field) checked.

### Step 7: Prototype spike (3 bacterial genomes) + GO/NO-GO verdict
Files: scripts/genome_map_spike.py, tests/test_genome_map_spike.py, wiki/genome_map_spike_verdict_2026-06-17.md
Depends on: Step 2, Step 6

**What changes:**
- **M1:** run the pipeline on 3 BACTERIAL genomes (Open Question A: E. coli cohort strain + a 2nd bacterium [Klebsiella/Pseudomonas] + a homology-heavy bacterium) — annotate+amrfinder (Step 1) → manifest (Step 2) → tier (3) → overlay (4) → assemble (5) → gate (6). Emit `wiki/genome_map_spike_verdict_<date>.md`: per-genome tier counts + `unknown_under_bakta_db_light` + join-quality counts + the G1 evidence + G2 + GO/NO-GO. BLOCKED-gates honestly if Step 2 is BLOCKED.

**Test strategy:**
- Orchestration test on synthetic/cached inputs (no live tools): the verdict aggregates per-genome results + computes GO/NO-GO from the gate (incl. the all-symbol-fallback NO-GO); a BLOCKED manifest → BLOCKED verdict. Live 3-genome run is the manual deliverable.

## Execution Preview

Wave 0 (1):  Step 1 — tool plumbing (shared loader + Bakta + AMRFinder runners)
Wave 1 (1):  Step 2 — tool-surface manifest (feasibility + inventory)
Wave 2 (2):  Step 3 — tier classifier, Step 4 — determinant overlay + join gate
Wave 3 (1):  Step 5 — map assembler (raw-field schema)
Wave 4 (1):  Step 6 — G1/G2 gate
Wave 5 (1):  Step 7 — prototype spike + verdict

Critical path: Step 1 → 2 → 3 → 5 → 6 → 7 (6 waves)
Max parallelism: 2 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags
- Severity: high — **Bakta + AMRFinder unproven end-to-end on this host** (smoke deferred; Docker-mount-corruption history). The Step-2 manifest is the feasibility gate; a wedge → BLOCKED, and Step 7 BLOCKED-gates rather than faking a map. [grounded]
- Severity: high — **Determinant join quality is the integrity crux** (the gene-symbol trap). Mitigated by the `DeterminantHit` coord/protein join + the hard gate (symbol-fallback excluded from determinant-phenotype + G1; all-fallback → NO_GO). [grounded — brainstorm]
- Severity: medium — **External tool surfaces [unverified] until Step 2:** the AMRFinder `main.tsv` columns ACTUALLY present (does this version emit protein-id/contig/coords?) determine whether high-confidence joins are even possible; if coords are absent, MOST joins fall to symbol-fallback → likely NO-GO (an honest outcome, surfaced). Bakta db-light product vocabulary similarly drives the tier boundary. The manifest pins both. [grounded]
- Severity: medium — **db-light unknown-rate is tooling-coverage** → the field is DB-labelled `unknown_under_bakta_db_light`. [grounded]
- Severity: low — **TB/fungal out of the v1 spike** (VCF-vs-GFF contract) — schema future-compat only; not a defect.
- Severity: low — **File overlap:** none within a wave (`genome_map/__init__.py` created once in Step 1; Steps 3/4 write disjoint modules). Cross-wave import edges captured by `Depends on:`.

## Open Questions
- A. **The 3 bacterial prototype genomes** — E. coli cohort strain (rich AMR) + a 2nd bacterium (Klebsiella/Pseudomonas) + a homology-heavy/hypothetical-protein-rich bacterium (the honesty stress test; must be genuinely UNFAMILIAR for the G1 test). Confirm exact accessions.
- B. (Resolved by the brainstorm: symbol-fallback determinant joins are EXCLUDED from `determinant-phenotype` + G1 credit; TB dropped from v1; provided-GFF = Bakta-compatible only. The spec record was reconciled to match.)

## Verification
1. `uv run pytest tests/test_genome_map_ingest.py tests/test_genome_map_tool_surface.py tests/test_genome_map_tiers.py tests/test_genome_map_overlay.py tests/test_genome_map_build.py tests/test_genome_map_gate.py tests/test_genome_map_spike.py -q` — all green.
2. `uv run pytest tests/ -q` (excluding `tests/test_models_foundation.py`) — 0 regressions; frozen AMR surface byte-unchanged.
3. Tool surface: `wiki/genome_map_tool_surface_<date>.json` exists (Bakta fields + vocab + AMRFinder headers) or a documented BLOCKED.
4. Map: a per-feature JSON map with the phenotype wall holding (only high-confidence determinant features carry phenotype), `unknown_under_bakta_db_light` present, join-quality counts emitted.
5. Gate: `wiki/genome_map_spike_verdict_<date>.md` reports per-genome DB-labelled unknown rate + join-quality + ≥3 G1 prevent-wrong-inference (≥1 type a/b) + G2 + GO/NO-GO; NO-GO if all determinant joins are symbol-fallback — or an honest NO-GO/BLOCKED.

## Save-time amendments

Captured at: 2026-06-17
Source: `/save-plan` arguments

> Audit-notes-only: `/execute-plan` reads ONLY `## Implementation Steps`. These amendments are provenance — the fixes below are ALREADY folded into the Steps above (this is the v2 post-/brainstorm regeneration, not a v1+amendments stack).

- genome-map v2: brainstorm-hardened — AMRFinder-run step, shared ##FASTA loader, DeterminantHit + hard join-quality gate, raw-field schema, tool-surface manifest gating, bacterial-only v1

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high mode=feature:genome-map -->
