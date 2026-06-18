# Session handoff — 2026-06-17: EXECUTE the genome-map plan (fresh session)

Everything a cold session needs to implement the genome-map v1 spike. Repo: `C:\Users\Farshad\PythonProjects\dna_decode`
(origin `mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`, branch `main`, all work pushed). HEAD ≈ `1d4d1c6`.

## ⇒ START HERE
1. `cd C:\Users\Farshad\PythonProjects\dna_decode`
2. Read the plan: **`features/genome-map/technical-plan.md`** (status: candidate; 7 steps / 6 waves; v2 brainstorm-hardened — the executable contract).
3. Read the spec/decisions: `features/genome-map/feature.md` (AC1–AC13 + the `## Repo grounding` log = intake/design/spec/probe/brainstorm evidence).
4. **Settle the one open decision** (Open Question A — the 3 prototype genomes; see below) — a one-line "use these accessions" unblocks Step 7.
5. Run:
```
/execute-plan features/genome-map/technical-plan.md
```
   It will run **sequential mode** by convention (commit-to-main, no PRs). The whole pre-exec chain
   (idea-anchor→intake→design→spec→technical-plan→probe→design↺→spec↺→brainstorm→technical-plan↺→save-plan)
   is COMPLETE — no more planning is owed; this is pure execution.

## What this builds (1 paragraph)
A **personal exploration/QC tool**: point at ONE microbial genome → an honest, evidence-tiered per-feature
map. It RE-TIERS Bakta's own annotation (4 tiers: determinant-phenotype > curated-molecular-function >
homology-only-hypothesis > unknown) + overlays the existing AMRFinder determinant cells, with phenotype
claims ONLY behind a validated-determinant wall and a DB-labelled unknown rate. It is the honest,
label-free, achievable form of the project north star — NOT a learned phenotype predictor (that arm is a
closed negative; see `wiki/north_star_distance_brainstorm_2026-06-17.md`). v1 is a **spike** that ends in a
GO/NO-GO verdict on whether to invest further.

## The 7 steps (waves) — the contract is in the plan; this is the map
- **Step 1 (W0)** `dna_decode/genome_map/{__init__,ingest,annotate,amrfinder}.py` — shared `##FASTA`-stripping GFF loader + Bakta runner + AMRFinder runner.
- **Step 2 (W1)** `scripts/genome_map_tool_surface.py` — run Bakta+AMRFinder on 1 genome, emit `wiki/genome_map_tool_surface_<date>.json` (Bakta fields/vocab + AMRFinder headers). **This is the feasibility GATE** + it seeds Steps 3/4.
- **Step 3 (W2)** `dna_decode/genome_map/{tiers,tier_vocab}.py` — tier classifier (vocab seeded from the Step-2 manifest).
- **Step 4 (W2)** `dna_decode/genome_map/phenotype_overlay.py` — `DeterminantHit` + the **hard join-quality gate**.
- **Step 5 (W3)** `dna_decode/genome_map/build_map.py` — assembler + raw-field JSON/table + `unknown_under_bakta_db_light`.
- **Step 6 (W4)** `dna_decode/genome_map/gate.py` — G1 (prevent-wrong-inference) + G2 + verdict.
- **Step 7 (W5)** `scripts/genome_map_spike.py` — 3-genome run → `wiki/genome_map_spike_verdict_<date>.md` (GO/NO-GO).

## ⚠️ The 5 brainstorm catches — DO NOT re-introduce (they're already in the plan)
1. **Bakta GFF has an embedded `##FASTA` block** that `parse_gff3` / `load_annotation_table` choke on → ALWAYS go through the shared `ingest.load_genome_gff` (strips `##FASTA` first; pattern = `scripts/pathotype_laptop_pipeline.py::parse_bakta_gff3:70`). The offline provided-GFF path uses it too.
2. **AMRFinder must be RUN** (the overlay needs `main.tsv`) — reuse `scripts/drug_mechanism_audit.py::_run_amrfinder(fasta, out_dir, organism)` (`-O` is organism-specific for QRDR). **Organism is an explicit input — do NOT auto-detect from Bakta in v1.**
3. **Determinant→feature join = the integrity crux.** `cipro_determinants_from_main`/`call_resistance` return symbol/class **NOT coordinates**; joining by `gene_symbol` is the documented 0%-overlap trap. Parse the raw `main.tsv` into a `DeterminantHit` (retain protein-id/contig/start/end), join protein-id → coord-overlap → symbol-fallback with `join_confidence`. **Symbol-fallback joins are VISIBLE but EXCLUDED from `determinant-phenotype` + G1; the spike NO-GOs if ALL determinant joins are symbol-fallback.**
4. **Map schema retains raw fields** (`raw_product`/`raw_gene_symbol`/`raw_locus_tag`/`raw_feature_type`/`source_tool`/`classification_reason`/`secondary_evidence[]`) so Step 6's G1 computes from the map alone.
5. **Step 2's manifest GATES Steps 3+4** — the tier wording rules + the AMRFinder join keys are seeded from the REAL inventory, never guessed.

## ⚠️ Open Question A (settle before Step 7) — the 3 prototype genomes
Must be **bacterial** (TB is OUT of v1 — VCF-vs-GFF contract mismatch) and the 3rd genuinely UNFAMILIAR (so G1 isn't confirmation bias). Recommended:
- (i) a **cipro-R E. coli** from the existing N=147 cohort (`data/processed/stage2_n150_cipro_cohort.parquet` — pick an ST131 R strain; rich AMR overlay). Organism `Escherichia`.
- (ii) a **Klebsiella pneumoniae** (carbapenem-R; the project has Kleb cells). Organism `Klebsiella_pneumoniae`.
- (iii) a **homology-heavy / hypothetical-protein-rich bacterium** you don't know (e.g. a less-studied environmental isolate) to stress the homology + unknown tiers. Organism per AMRFinder's `-O` list (or run AMR with no `-O` / a generic — the manifest will show coverage).
Resolve by confirming the 3 accessions; genomes download via NCBI Datasets (the project's `download_cohort_genomes` / refseq path).

## Tooling reality (feasibility)
- **INSTALLED** (CLAUDE.md "Stage 2 toolchain"): Bakta `oschwengers/bakta:v1.11.4` + **db-LIGHT** at `C:/Users/Farshad/dna_decode_stage2/bakta_db/db-light`; AMRFinderPlus `ncbi/amr:4.2.7` + DB at `C:/Users/Farshad/dna_decode_stage2/amrfinder_db`; BLAST+ native at `C:/Users/Farshad/ncbi-blast/bin`. NOT installed: hmmer/Pfam/eggNOG (out of v1 scope by design).
- **Bakta annotation was NEVER smoke-tested on this host** (deferred 2026-05-15, CPU-heavy) → **Step 2 is where it's first proven**. Expect it to be SLOW; a wedge (the documented WSL2/Docker-mount-corruption gotcha — `wsl --shutdown` to recover) → emit `BAKTA_ANNOTATION_BLOCKED`, never a fake map.
- Bakta `bakta_db` needs the `--entrypoint /bin/bash -c "..."` wrapper; direct docker from Git Bash needs `MSYS_NO_PATHCONV=1` (see CLAUDE.md gotchas + `wiki/stage2_install_artifact_2026-05-15.md`). The reuse helpers already handle these.
- AMRFinder `main.tsv` uses the `Gene symbol` column; whether THIS version emits protein-id/contig/coords is `[unverified]` → **Step 2's manifest decides it.** If coords are absent, most joins fall to symbol-fallback → likely an honest NO-GO (surfaced, not hidden).

## FROZEN — never edit (leak guard asserts byte-equality)
`dna_decode/eval/amr_rules.py` · `dna_decode/data/calibrated_amr_rules.json` (reproducibility freeze 2026-06-13).
The overlay READS these; it must not modify them.

## Verification (the execute-plan acceptance)
1. `uv run pytest tests/test_genome_map_*.py -q` — all green.
2. `uv run pytest tests/ -q` (exclude `tests/test_models_foundation.py` — host torch-paging limit) — 0 regressions; frozen AMR surface byte-unchanged.
3. `wiki/genome_map_tool_surface_<date>.json` exists (or documented BLOCKED).
4. A per-feature JSON map with the phenotype wall holding + `unknown_under_bakta_db_light` + join-quality counts.
5. `wiki/genome_map_spike_verdict_<date>.md` = GO/NO-GO (NO-GO if all-symbol-fallback) — an honest NO-GO/BLOCKED is a valid outcome.

## Working conventions (this project)
- `python` / `uv run` (NOT `python3`). Commit straight to `main` (= the cross-machine sync channel; user syncs ~weekly). Commit footer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- Stage only your own files. Pre-existing dirty files to LEAVE: `uv.lock`, `wiki/ciprofloxacin_mechanism_audit_2026-06-05.*`, `bash.exe.stackdump`, `research_outputs/eukaryotic-...unsupported.md`, `wiki/3 idea anchors...rtf`, the untracked stale `plans/TB_AMR_Decoder_CRyPTIC_Technical_Plan.md` (superseded; safe to delete).
- Codex `/brainstorm` JSON-parse can hit a Windows cp1252 crash on `→` chars — force `PYTHONIOENCODING=utf-8` + `errors="replace"` and DON'T `rm` the codex_out before parsing.
- Project ledger: `project_state/dna-decode-2026-05-11.md` (action log through row 121). Log a `--append-action` row per executed batch.

## Parked / not-this-thread
- **TB decoder** is at a clean waypoint (RIF lineage-collapsed sens **0.739** after the MNV-matcher fix, `TB_SUBSET_PLUMBING`; full BASELINE needs the ~1.6 TB regeno fetch → D: via `scripts/populate_tb_regeno_detached.bat`; independent number needs a hand-curated post-2023 gold set). Plan: `plans/TB_AMR_Decoder_RIF_INH_On_CRyPTIC_Plan/`. Do NOT mix TB into the genome-map spike.
- The learned/embedding decoder is a CLOSED NEGATIVE on free data — do not reopen without a new label-clean substrate.
