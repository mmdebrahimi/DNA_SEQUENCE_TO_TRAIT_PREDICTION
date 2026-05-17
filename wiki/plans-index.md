# Plans Index
<!-- Auto-maintained by /save-plan. Do not edit manually. -->

## [plan_file: Cipro_Decision_Bundle_Technical_Plan.md] 2026-05-17
**Summary:** Implementation blueprint for the post-SUSPEND_CONDITION_4 decision bundle, scoped down per /review reductions: collapse Tier 0 Steps 1-2 into one script, drop Bakta + circular variants, defer mean+max preflight v3.
**Key decisions:**
- 12 steps across 5 waves: Wave 0 (3 parallel foundation refactors) → Wave 1 (5 parallel test+consumer) → Wave 2 (2 parallel tests) → Wave 3 (manual census+manifest run) → Wave 4 (conditional label-sensitivity LOSO)
- Critical path: Step 1 → Step 7 → Step 9 → Step 11 → Step 12 (5 waves; 3 code + 2 runtime); max parallelism = 5 agents (Wave 1)
- /review reductions absorbed: collapse census + manifest into one script (Step 7), drop Bakta smoke + mechanism completeness + 3 of 5 manifest variants + --ignore-gate flag, defer mean+max preflight v3 to a conditional follow-up
- Pre-conditions PC1 (D9 framing lock) + PC2 (D4 numeric threshold) are user-locked BEFORE Wave 3 runs (not enforced by code)
- Source plan drift flagged: `plans/Cipro_Decision_Bundle_Plan.md` still contains un-reduced spec; should be edited in-place before /execute-plan per 2026-05-14 HIGH-salience lesson

---

## [plan_file: Cipro_Decision_Bundle_Plan.md] 2026-05-17
**Summary:** Tier-0 cheap decision bundle on N=38 + BV-BRC MIC census, fired BEFORE any Databricks burst or per-gene NT diagnostic. Replaces the rejected binary Path A vs Path B framing after SUSPEND_CONDITION_4 verdict landed 2026-05-17.
**Key decisions:**
- D1: Replace binary A/B with 4-tier decision bundle (Tier 0 cheap → Tier 3 large spend gated)
- D2: BV-BRC census must cover multiple phenotype policies (HIGH_R/HIGH_S + CLSI-strict + EUCAST-strict)
- D3: Frozen label-policy manifest before any relabeling experiment (single source of truth, no hand-coded overrides)
- D4: Primary estimand for relabeling = per-strain error concentration + rank-order stability, NOT max AUROC
- D5: Mean+max attribution preflight v3 is a closeout falsifier (Tier 1), not a fork-decider
- D6: Curated baseline informational run needs 3rd verdict field `given_suspended_gate: INFORMATIONAL_ONLY`
- D7: Bakta 4-strain smoke selection requires a negative-control strain (borderline_S no-mech)
- D8: Mechanism completeness — AMRFinder differential test first, manual blastn only on discordant rows
- D9 (open tradeoff): Phase 1 EP1 deliverable framing ("publish" vs "ship working classifier") affects Tier 2/3 weights — user must lock before Tier 3 fires

---

## [plan_file: Return_Decision_Tree_Patch_Plan.md] 2026-05-16
**Summary:** Apply the /review synthesis's 3 surgical correctness fixes + structural restructure to `wiki/return_decision_tree_2026-05-16.md`, plus the one-line `-u` flag fix to `run_stage1b_detached.bat`. Implementation-ready under 4 steps; max parallelism 2 (doc + bat are independent).
**Key decisions:**
- D1: Eng correctness patches come FIRST, restructure SECOND (sequenced Step 1 → Step 2 because both modify the same doc file)
- D2: `-u` flag (not `PYTHONUNBUFFERED=1` env var) for the .bat fix — explicit + tighter blast radius
- D3: Merge enumerated 8 sub-steps into 4 actual steps; adjacent same-file ops serialize anyway
- D4: Don't touch the running Stage 1b process; Step 3 affects FUTURE relaunches only

---

## [plan_file: Sidework_Sequence_Ship_Path_Plan.md] 2026-05-13 (decisions locked 2026-05-14)
**Summary:** Delta from `Sidework_Sequence_Plan.md` after `/review` (2026-05-13) + post-save `/brainstorm` (2026-05-14). Resolves the load-bearing B-scope problem (B-B locked: drop clade-only from smoke; 12 unique MLST → singleton clades make clade-only degenerate). Fixes deterministic-hashing reproducibility bug. Narrows ARCHITECTURE.md to per-line judgment. Locks --per-class=20. Excludes wiki/GATE_A_REPORT.md from C scope.
**Key decisions (all locked 2026-05-14):**
- B-resolution = B-B (drop clade-only from smoke; 4 variants) (D1)
- Step C scope = 13 files (11 originals + LESSONS_LEARNED + docs/ARCHITECTURE.md with per-line judgment); wiki/GATE_A_REPORT.md NOT in scope (D2)
- Step C is edit-then-stage, not stage-only (D3)
- `mlst_to_clade_id` helper lives in `dna_decode/data/cohort.py` (D4)
- Test scope = 1 parametrized in tests/test_data_cohort.py + 1 integration in tests/test_pipeline_cli.py (both files exist; extend) (D5)
- Drop `fallback_counter: dict`; use deterministic hashing (D6)
- Multi-scheme MLST + scheme-collision: scheme-aware tuple hash via `zlib.crc32` (NOT Python `hash()` — process-salted) (D7)
- Step M demoted out of wave graph (D8)
- Time-box C to 30 min (D9)
- B commit message includes numerical before/after of clade-only AUROC (D10)
- Post-populate slow tests gate the smoke gate (D11)
- Step D `--per-class 20` locked (N=40 with 5R margin) (D12)
- Register `slow` pytest marker (D13)

---

## [plan_file: Sidework_Sequence_Plan.md] 2026-05-13
**Summary:** Ordered work to do while NT v2 100M populate runs in background (~45 min remaining). Brainstorm-revised after Codex critique surfaced under-scoped C, mixed-scope TODOS hunks, and B's deeper-than-15-min reality.
**Key decisions:**
- Sequence A → C → E → B; D optional last (D1)
- C scope = current-state files only; skip archived plans + project_state snapshots (D2)
- TODOS.md hunk staging via `git add -p` across A, C, E (D3)
- B = narrow scope (helper extraction + unit test); defer per_clade_baseline strain-keying (D4)
- pytest discipline during populate = CPU-only target tests via `-m "not slow"` (D5)
- Auto-memory user_environment.md update is separate from C (D6)

---

## [plan_file: Phase2_Decision_Gate_Plan.md] 2026-05-13 (Step 1/2 updated 2026-05-14 for B-B lock)
**Summary:** Split the originally-conflated "12-strain decision gate" into a 12-strain smoke/falsification gate (4 variants, clade-only dropped per B-B lock 2026-05-14) + a tiered N=50 → N=150 staged decision gate (local screen → Databricks burst). Stage 2 acceptance: NT ≥ best classical + 5 pp AUROC AND top-10 NT attribution includes gyrA / parC / parE (biological-plausibility check).
**Key decisions:**
- 12-strain = smoke gate, tiered N=50 → N=150 = staged decision gate (D1)
- Fix clade-only `hash(mlst) % 10` placeholder BEFORE running smoke gate (D2)
- Keep 12-strain smoke at existing 5 variants — no RF/TabPFN/SNP-table additions yet (D3)
- Document the smoke result as smoke, not decision (D4)
- Tiered Option-C threshold: N=50 local screen "any positive lift" → N=150 Databricks "5 pp + biology check" (D5)
- TODOS additions for 4 deferred /research items; SNP-table scope corrected to "parse AMRFinderPlus POINT* rows" per Codex 2026-05-13 (D6)

---

## [plan_file: NT_Deferral_Docs_Cleanup_Ship_Path_Plan.md] 2026-05-13
**Summary:** Scope-tightened delta from `NT_Deferral_Docs_Cleanup_Plan.md` after `/review` synthesis — drops D1's `/save-plan` hedge (already disproven on disk), trims the over-prose `[BLOCKED]` bullet, sharpens the deferral annotation, attempts cheap NT revision retrieval before recording the gap, and adds an explicit untracked-file staging note.
**Key decisions:**
- Drop the `/save-plan` hedge — direct manual edit only (D1)
- Trim `[BLOCKED]` bullet to ≤3 sentences + separate Environment line (D2)
- Attempt NT revision retrieval before recording the gap (D3)
- Sharpen "gate failed" → "equivalence test failed at model load" (D4)
- Step 3 explicit `git add` reminder for untracked plan file (D5)

---

## [plan_file: NT_Deferral_Docs_Cleanup_Plan.md] 2026-05-13
**Summary:** Docs-only follow-up to commit `d4a4652` — apply the three issues surfaced by `/brainstorm` against the just-shipped NT AutoModel refactor deferral: stale plans-index, conflated TODOS scope, missing reproducibility metadata.
**Key decisions:**
- Re-run `/save-plan` before manual index edit (D1)
- Split the TODOS entry into specific (BLOCKED) + general (OPEN) (D2)
- Lean reproducibility metadata, one line (D3)
- Diagnostic spike deferred, NOT killed (D4)

---

## [plan_file: Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md] 2026-05-13
**Summary:** Scope-reduced delta from `Audit_Calibration_NT_AutoModel_Plan.md` after `/review` synthesis — 5 steps / 4 waves → 2 commits, 2 waves. Drops dual-verdict columns (institutionalization risk), drops wiki update (no text to replace), splits NT refactor into a separate gated commit (equivalence test required).
**Key decisions:**
- Asymmetric warning banner replaces dual-verdict columns (D1)
- Drop wiki/GATE_B_REPORT.md update entirely (D2)
- Split NT refactor into a separate, gated commit (D3)
- Keep default-semantics test as the regression lock (D4)
- thresholds_block(rules) helper, not inline string list (D5)
**Status:** Commit 1 shipped (473b8eb); Commit 2 deferred 2026-05-13 — equivalence test failed at model load (AutoModel.from_pretrained state_dict mismatch on NT v2 100M trust_remote_code checkpoint). See plans/NT_Deferral_Docs_Cleanup_Plan.md and TODOS.md [BLOCKED] NT AutoModel swap.

---

## [plan_file: Audit_Calibration_NT_AutoModel_Plan.md] 2026-05-13
**Summary:** Fix a credibility bug in the just-shipped audit cohort generator (`scripts/audit_cohort.py`) — the "GO" verdict was emitted under silently-relaxed thresholds; defaults produce "WARN" — AND simultaneously replace `NucleotideTransformerModel`'s `AutoModelForMaskedLM` with `AutoModel` to eliminate the `output_hidden_states=True` workaround.
**Key decisions:**
- Audit report header MUST surface threshold values (D1)
- Two verdict columns — Phase 1 production + Gate B infra-only (D2)
- Pin default semantics in tests (D3)
- NT switches to `AutoModel`, not `AutoModelForMaskedLM` (D4)
- Pooling-strategy tag stays "single_seq_mean" (D5)
**Status:** Commit 1 shipped (473b8eb); Commit 2 deferred 2026-05-13 — equivalence test failed at model load (AutoModel.from_pretrained state_dict mismatch on NT v2 100M trust_remote_code checkpoint). See plans/NT_Deferral_Docs_Cleanup_Plan.md and TODOS.md [BLOCKED] NT AutoModel swap.

---

## [plan_file: BVBRC_Genome_Metadata_Adapter_Plan.md] 2026-05-12
**Summary:** Wire `BVBRC_genome.csv` (BV-BRC Genomes-tab export) into the cohort path as a new adapter module, bypassing the wrong-contract `pilot.fetch_ncbi_assembly_quality` scaffold and feeding the existing `--assembly-metadata` wire that `cohort.candidates_from_bvbrc_ast` already accepts.
**Key decisions:**
- Bypass the scaffold instead of implementing it (D1)
- New CLI flag rather than overloading existing `--assembly-metadata` (D2)
- Coverage-log line surfaces ID-namespace mismatches early (D3)
- `fetch_ncbi_assembly_quality` stays scaffolded (D4)

---

## [plan_file: Ecoli_G2P_Phase1_Closeout_Plan.md] 2026-05-12
**Summary:** Wrap up the stalled `/execute-plan` epilogue for `Ecoli_G2P_Phase1_Ship_Path_Plan.md` — toolchain restore, doc reconciliation, first authoritative test pass, archive, state cleanup, push, final report.
**Key decisions:**
- Selective expansion over hold scope (D1)
- Real-data validation = Phase 2 entry criterion, not Phase 1 closeout (D2 — pending)
- Toolchain restore approach: uv vs pip (D3 — pending)
- Archival convention: status-header + git tag recommended (D4 — pending)
- `/documentation` before commit (D5)
- Retrospective re-derivation, not skip (D6)
- Test outcome recorded, not gated (D7)
- Delete both state files at end (D8)
- `.claude/execute-plan-state/` added to `.gitignore` (D9)

---

## [plan_file: Ecoli_G2P_Phase1_Ship_Path_Plan.md] 2026-05-12
**Summary:** Contracted path to ship Phase 1 of `Ecoli_G2P_Platform_Technical_Plan.md`. Captures the `/review` synthesis verdict (HOLD scope + selective contraction within remaining steps) plus the deferred Wave 3.5 hardening fixes from the post-Wave-3 `/brainstorm`. Estimated remaining work: ~700 LOC across 5 implementation steps + 4 hardening edits.
**Key decisions:**
- HOLD scope, do not expand (D1)
- Reorder — Step 15 (smoke + fixtures) BEFORE Step 14 (CLI) (D2)
- Step 14 collapses to one `scripts/pipeline.py` with subcommands (D3)
- Step 13 visualization uses matplotlib + TSV export, NOT pygenometracks (D4)
- Step 17 leaderboard collapses to a shell loop over `pipeline.py train` (D5)
- Step 16 docs trimmed to README + ARCHITECTURE.md only (D6)
- Apply Wave 3.5 hardening BEFORE Step 14 wiring fires (D7)
- Add quantization-fidelity micro-step (selective addition) (D8)

---

## [plan_file: Gene_Presence_AUROC_Bug_Fix_Plan.md] 2026-05-14
**Summary:** Strengthen the diagnostic, confirm the strain-unique-identifier-domination hypothesis on real data, then add a `gene_symbol` column to `AnnotationTable` so the gene-presence smoke variant returns a non-degenerate AUROC at N=12.
**Key decisions:**
- Add `gene_symbol` column to AnnotationTable; do NOT rewrite `gene_id` (preserves embedding cache key + fixes parse_gff3/parse_genbank asymmetry) (D1)
- Strengthen diagnostic before mounting F: drive (absolute counts + per-prefix namespace breakdown + side-by-side dual-extractor AUROC) (D2)
- Strengthen synthetic falsifier with strain-unique-blocks + shared-core LOSO + all-zero held-out row case (D3)
- Add `INDETERMINATE_IDENTIFIER_OOV` smoke verdict as defense-in-depth guardrail (D4)

---

## [plan_file: Stage1_N40_Cipro_Engineering_Screen_Plan.md] 2026-05-14
**Summary:** Run a 4-experiment matrix (NT-XGBoost gate + NT-logreg sanity + k-mer-XGB classical + NT+k-mer-fusion-logreg diagnostic) under LOSO on the N=40 cipro cohort (effective N=38) with paired bootstrap CI, MLST diagnostic appendix, and a 3-bucket verdict to decide whether to spend Stage 2 N=150 Databricks burst budget.
**Key decisions:**
- Restore NT-XGBoost as primary gate-bearing head; add NT-logreg as sanity-check baseline; fusion is diagnostic-only NOT gate-bearing (D1)
- All variants run with `calibrate=False` for primary AUROC (uniform calibration discipline matching smoke-gate; calibration is small-N footgun) (D2)
- Add diagnostic appendix: MLST distribution + per-strain LOSO predictions + paired bootstrap CI (B=1000) + 3-bucket verdict (≥5 pp CLEAN / 3-5 pp NOISY / <3 pp FAIL) (D3)
- Gene-presence + AMRFinderPlus POINT* baselines explicitly out of scope; result packet notes 'best classical' is bounded (D4)

---

## [plan_file: Stage1_Refactor_And_Test_Hardening_Plan.md] 2026-05-14
**Summary:** Convert the /review synthesis into a 3-step refactor: pre-commit decision rules in the Stage 1 plan, reduce `scripts/stage1_n40_cipro.py` to thin orchestration over existing infrastructure, and pin the two critical untested behaviors (fusion-exclusion from gate + `calibrate=False` discipline).
**Key decisions:**
- Treat Stage 2 burst as atomic; pre-commit deterministic per-bucket actions (CI-lower-bound rule converts borderline NOISY PASS → FAIL) (D1)
- Refactor runner to reuse existing infrastructure: `leave_one_strain_out_cv` for NT variants, factored `dna_decode/eval/loso_kmer.py` for k-mer + fusion, `_train_baseline_logreg(calibrate=False)` for logreg path (D2)
- Replace silent mean-fallback on `ClassifierTrainingError` with re-raise; eliminate the failure-masking pattern (D3)
- Add fusion-exclusion + `calibrate=False` discipline regression tests; pin /brainstorm-flagged failure modes (D4)
- Loud MLST handling (raise on None instead of "unknown" fallback) + bootstrap-skip-count reporting (D5, D6)

---

## [plan_file: Stage2_N150_Prep_Plan.md] 2026-05-14
**Summary:** Resolve the three deferred Stage 2 decisions (annotation source, AMRFinderPlus integration, Databricks vs local) and ship the infrastructure needed for a Stage 2 N=150 cipro decision-gate run, so the gate runs cleanly once Stage 1 PASSes.
**Key decisions:**
- Annotation source = Bakta re-annotation for cross-strain stable gene symbols (defer Roary; accept-degenerate rejected) (D1)
- AMRFinderPlus POINT* SNP-table baseline IS in scope for Stage 2 (gyrA/parC/parE textbook signal; load-bearing comparator) (D2)
- Compute = Databricks burst for N=150 NT populate (~3-5 hr A100); local CPU for everything else (Bakta + AMRFinder + analysis) (D3)
- Stage 2 cohort = N=150 expanded from gate_b_cohort.parquet (67 strains) via audit-cohort pipeline with relaxed assembly-quality thresholds (D4)

---

## [plan_file: Stage2_Docker_Tools_Install_Plan.md] 2026-05-14
**Summary:** After user starts Docker Desktop, install Mash + Bakta + AMRFinderPlus via pinned Docker images, write a single tools/docker_runner.py Python orchestration module (NOT .sh wrappers), and smoke-validate on K-12 + one cipro-R strain. Resolves the Phase A.1 / A.2 / A.5 install steps from Stage2_N150_Prep_Plan.md.
**Key decisions:**
- Docker (containers) — NOT WSL2 Ubuntu (D1)
- One Python tools/docker_runner.py module — NOT three .sh wrappers (D2)
- Pin Docker image tags — NOT :latest (D3)
- Correct AMRFinderPlus invocation: amrfinder_update for DB, --database (not --database_path), --mutation_all <file> takes a path (D4)
- Stage Bakta DB on C: drive if room; verify before flagging install complete; record versions+SHA digests (D5)
- Smoke-test on TWO strains: K-12 (binary works) + one cipro-R (POINT-row parsing actually exercised) (D6)

---
