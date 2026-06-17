# TB AMR Decoder (RIF + INH) on CRyPTIC — Technical Plan

> First M. tuberculosis decoder cell: deterministic WHO-catalogue determinant rule for RIF (rpoB) + INH (katG + inhA promoter), lineage-collapsed on CRyPTIC as a labelled knowledge-baseline (lineage via a pinned barcode caller, NOT de-novo SNP distance), with an independent post-2023 gold-set arm in-scope.

## Lens status
Inputs: conversation (handoff settled-design + 2026-06-16 authority decisions + 2026-06-17 pre-exec /brainstorm fixes C1-C3 + M1-M2). Regenerates the v1 plan; supersedes the untracked draft `plans/TB_AMR_Decoder_CRyPTIC_Technical_Plan.md`.
Degradations: repo-index unconfigured — direct file reading only; lineage-barcode SNP set + REGENOTYPED_VCF callability schema [unverified] until built; TB Portals independent cohort access-gated [unverified].

## Problem Statement
The deterministic AMR decoder counts curated resistance determinants → R/S and is strongest on point-mutation mechanisms (cipro QRDR 0.925, Klebsiella 1.0). The public-label AMR track on existing organisms is saturated; the binding constraint is labels. CRyPTIC provides a free, deep TB substrate (12,287 *M. tuberculosis* isolates × 13 drugs, reference broth-microdilution MIC + binary phenotype + per-isolate VCF vs H37Rv NC_000962.3). The feasibility probe (2026-06-16) returned GREEN; coordinate alignment to the WHO catalogue is VERIFIED (rpoB S450L → 761155 C>T; katG S315T → 2155168 C>G match cached VCFs exactly).

This plan ships the FIRST TB decoder cell — deterministic, two mechanisms (RIF = `rpoB` RRDR; INH = `katG` + `inhA` promoter), rule from the WHO *M. tuberculosis* mutation catalogue v2 (2023, pinned commit 0bb39143), scored over a clonality-corrected CRyPTIC cohort.

Two ratified authority decisions define "done":
- **Deliverable (a):** ship the cell as a labelled in-distribution **knowledge-baseline** (WHO catalogue built partly from CRyPTIC → recalls its own training phenotypes). All results tagged `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`.
- **Deliverable (b):** the post-2023 **independent gold set is in-scope now** — build the independent-cohort scoring harness + acquisition runbook, scored separately, labelled `INDEPENDENT_VALIDATION`, never conflated with the baseline.

Pre-exec /brainstorm (2026-06-17) folded in: **C1** the VCF call rule is `FILTER==PASS` + `GT` non-reference (no `GCP` FORMAT key exists); **C2** lineage collapse is a pinned barcode caller feeding the existing `cluster_weighted_confusion` (NOT de-novo SNP distance + representative-dedup); **C3** the baseline label gates on a full per-drug prevalence-preserving cohort + a callability denominator; **M1** report within-lineage mixed-prediction counts; **M2** pin masked-vs-REGENOTYPED VCF per purpose.

Non-goals: no learned/embedding model; no edit to the frozen E. coli surface; no drugs beyond RIF + INH; no de-novo SNP-distance transmission clustering (deferred to v1c).

## Codebase Context
- TB has **no code yet**: `dna_decode/organism_rules/` is absent (new package). Only TB assets are `scripts/cryptic_feasibility_probe.py` + `scripts/tb_coordinate_alignment_probe.py` + `scripts/stage_tb_vcf_subset.py` + the GREEN probe artifacts.
- **VCF schema (verified):** cached CRyPTIC masked VCFs have FORMAT `GT:DP:DPF:COV:FRS:GT_CONF:GT_CONF_PERCENTILE` and FILTERs `MIN_DP`/`MIN_FRS`/`MIN_GCP`/`mask`. There is **no `GCP` FORMAT field** — the quality floor IS the FILTER set, so `FILTER==PASS` encodes it. Cached files are DECOMPRESSED text despite the `.vcf.gz` name. The reuse table has a `VCF` (masked) and a `REGENOTYPED_VCF` (col 44, joint-genotyped with explicit reference calls) column — **no lineage column**.
- **Catalogue (pinned + verified):** `data/raw/who_tb_catalogue/` has the v2 master (`drug/gene/mutation/variant/tier/genomic position` + `INITIAL CONFIDENCE GRADING` = `1) Assoc w R` / `2) Assoc w R - Interim` / …) + the genomic-coordinate file (`variant → chromosome/position/ref/alt`), pinned at commit 0bb39143 in `CHECKSUMS`. 438 distinct grade-1/2 variants. Join key = `variant` (e.g. `rpoB_p.Ser450Leu`).
- **Lineage reuse:** `eval/clonality.py::cluster_weighted_confusion(preds, labels, clusters)` collapses each cluster to ONE vote — `cluster_class` (R/S/**DISCORDANT**, never label-majority-votes) + `_cluster_prediction` (member-prediction majority) + `wilson_ci` + `effective_lineage_n`. The barcode caller supplies `clusters` (strain → lineage_id); the existing fn does the rest. `greedy_representative_clusters_from_matrix` is NOT used in v1b (it needs FASTA-Mash distances TB doesn't have).
- Test infra: pytest under `tests/`; `uv run pytest tests/ -q` (exclude `tests/test_models_foundation.py`). Frozen E. coli suite must stay green.

### Reusable-Code Survey
- `dna_decode/eval/clonality.py` — `cluster_weighted_confusion` (the v1b metric, reused verbatim), `cluster_class`, `wilson_ci`, `effective_lineage_n`.
- `dna_decode/data/experimental_drug_rules.py` — non-frozen overlay pattern (branding constants, scorer-local rule, frozen-helper reuse) the TB modules mirror.
- `scripts/tb_coordinate_alignment_probe.py` — the verified catalogue↔VCF join + grade-1/2 filter logic Step 2/3 build on.
- `scripts/cryptic_feasibility_probe.py` — VCF fetch/cache (`fetch_vcf`) + reuse-table loader reused by acquisition.
- None — searched: dna_decode/organism_rules (absent), graphify-out/GRAPH_REPORT.md (absent), dna_decode/eval, dna_decode/data, scripts.

## Pre-Change Baseline
- No TB decoder cell exists. The only TB number is the probe PoC (RIF rpoB-RRDR window rule, sens/spec 1.0 on N=30) — explicitly NOT a validated metric.
- Frozen E. coli/Klebsiella/S. aureus/C. auris surface = 6 SCORED cells; `dna_decode/eval/amr_rules.py` + `dna_decode/data/calibrated_amr_rules.json` byte-frozen at commit b3761c8 (2026-06-13) — the byte-equality target for the leak guard.
- Cache state: 300 masked RIF VCFs staged (150R/150S — convenience sampling, NOT a v1b cohort); WHO catalogue pinned; coordinate alignment verified; no lineage caller, no scoring orchestrator.

## Verification Signal
- **v1a (plumbing):** `dna_decode/organism_rules/tb_amr.py` exists; cell tests pass; the coordinate-alignment fixture (`rpoB S450L` + `katG S315T` → real VCF records via the verified probe path) passes; per-sentinel R/S/ABSTAIN calls + `coverage_scope` + `excluded_grade12_loci` reported. NO sens/spec beyond fixtures.
- **v1b (baseline metric, deliverable a):** results artifact reports lineage-collapsed sens/spec (via `cluster_weighted_confusion`) + raw sens/spec + raw→lineage shrinkage + `n_discordant` + `n_clusters_mixed_prediction` (M1) + Wilson CI + `effective_lineage_n` + per-drug `n_uncallable` (callability denominator), tagged `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE` — emitted ONLY on the full per-drug prevalence-preserving HIGH-quality cohort; else a `TB_SUBSET_PLUMBING` or `LINEAGE_COLLAPSE_BLOCKED_NO_LINEAGE_CALL` status (never the baseline label on convenience data).
- **Leak guard:** `tests/test_tb_leak_guard.py` green — no CRyPTIC phenotype column read during rule construction, WHO files checksum-pinned, frozen E. coli surface byte-untouched.
- **Independent arm (deliverable b):** `scripts/score_tb_independent_goldset.py` + acquisition runbook exist; an independent run scores separately labelled `INDEPENDENT_VALIDATION`, else `INDEPENDENT_VALIDATION_BLOCKED_NO_GOLDSET`.
- Frozen E. coli suite: 0 regressions.

## Implementation Steps

### Step 1: TB VCF acquisition + genotype-aware parser + callability spec
Files: dna_decode/organism_rules/__init__.py, dna_decode/organism_rules/tb_vcf.py, tests/test_tb_vcf.py
Depends on: none

**What changes:**
- New `organism_rules` package. `tb_vcf.py`: acquisition (reuse `cryptic_feasibility_probe.fetch_vcf`; record provenance). **Ratified B:** the 300-VCF cache is plumbing-only; v1b fetches the FULL per-drug prevalence-preserving HIGH-quality eligible cohort (~10-12k VCFs, ~2.6 GB) into a **D:-drive cache** (C: too full) — restartable/skip-existing background fetch.
- **Genotype-aware parser (C1 — corrected):** a non-reference CALL requires `FILTER==PASS` AND `GT` allele index ≥1. **There is NO `GCP` FORMAT field** — the quality floor is already the `MIN_DP`/`MIN_FRS`/`MIN_GCP` FILTERs, so PASS subsumes it. Expose `DP`/`FRS`/`GT_CONF_PERCENTILE` as parsed provenance fields ONLY (never a second hard floor). Normalize MNVs/indels (needed for `inhA` promoter + `katG`). "position+ALT present" is NOT a call.
- **Callability spec (M2):** parse BOTH the masked `VCF` (canonical for determinant calls — matches the cache + the verified alignment) and, where present, the `REGENOTYPED_VCF` (explicit ref/no-call at union sites — the source of truth for whether a determinant window is callable). Define `is_callable(isolate, window)`: callable iff the regenotyped VCF has an explicit ref-or-alt call across the window; absent/`./.`/masked → uncallable. Document that masked VCFs are variant-only (a determinant-window absence is ambiguous without the regenotyped track).

**Test strategy:**
- Unit tests over ≥1 cached VCF: correct handling of `0/0`-with-ALT, `./.`, multi-allelic, FILTER-failed records; assert no code path reads a `GCP` FORMAT key.
- Sentinel: rpoB/katG variant parses as a non-reference call; a reference site does not; an uncallable window returns `is_callable=False`.

### Step 2: WHO catalogue load + pin-verify + RIF/INH determinant join
Files: dna_decode/data/tb_who_catalogue.py, tests/test_tb_who_catalogue.py
Depends on: none

**What changes:**
- `tb_who_catalogue.py`: load the already-pinned catalogue (`data/raw/who_tb_catalogue/`), verify against `CHECKSUMS` (commit 0bb39143 + per-file sha256), and JOIN master grade rows ↔ the genomic-coordinate file on `variant`. Build the RIF + INH grade-1/2 determinant table: RIF = `rpoB` RRDR grade-1/2; INH = **all WHO grade-1/2 INH loci** (ratified A — `katG` + `inhA`/`fabG1` promoter are the dominant contributors but `inhA` coding + `ahpC` grade-1/2 are included; report per-locus contribution + `coverage_scope`). Reuse the verified `tb_coordinate_alignment_probe` join/grade logic.

**Test strategy:**
- Join correctness on sentinels (`rpoB_p.Ser450Leu`, `katG_p.Ser315Thr` → grade-1, correct coords).
- Checksum-pin test fails if catalogue files change unpinned; grade-1/2 count == 438 (regression pin).

### Step 3: TB decoder cell (RIF + INH determinant scoring) — v1a plumbing
Files: dna_decode/organism_rules/tb_amr.py, tests/test_tb_amr.py
Depends on: Step 1, Step 2

**What changes:**
- `tb_amr.py` (non-frozen organism-routed cell, TMP-SMX overlay pattern): `score_rif(calls, callable_fn)` + `score_inh(calls, callable_fn)`; branding (`RULE_STATUS`, `RULE_SCOPE="organism_routed"`, `input_type="vcf_h37rv"`, catalogue commit/checksum). Returns R / S / **ABSTAIN** per drug — **ABSTAIN when the determinant window is uncallable** (C3; never susceptible-by-absence). Report `coverage_scope` + `excluded_grade12_loci`.
- Coordinate-alignment fixture (load-bearing, already verified): `rpoB S450L` + `katG S315T` real records resolve through parse→join→call.

**Test strategy:**
- Fixture-driven (20–50 sentinel VCFs: known determinants + negatives + an uncallable-window case). Each scores as expected; uncallable → ABSTAIN; out-of-scope variant → abstain/partial. NO sens/spec beyond fixtures.

### Step 4: TB lineage-barcode caller (pinned SNP barcode, VCF-native)
Files: dna_decode/organism_rules/tb_lineage.py, dna_decode/data/tb_lineage_barcode.py, tests/test_tb_lineage.py
Depends on: Step 1

**What changes:**
- **C2b — replaces the de-novo SNP-distance builder.** `tb_lineage_barcode.py`: the **Napier-2020 lineage barcode** (ratified F) sourced from the TBProfiler `tbdb` barcode file, **pinned by commit SHA + checksum** (same discipline as the WHO catalogue), as genomic positions + lineage alleles vs H37Rv NC_000962.3. Use the barcode DATA applied directly to our VCFs — NOT the TBProfiler caller (avoids the wrapper-trap). Coll-2014 (62-SNP) is the coarser fallback if the tbdb file can't be cleanly pinned.
- `tb_lineage.py`: assign each isolate a lineage/sublineage by matching barcode positions in its parsed Step-1 calls. Output `strain_id → lineage_id` (the dict shape `cluster_weighted_confusion`'s `clusters` consumes). Deterministic; no distance matrix, no threshold, no assembly. Isolates with no barcode hit → an explicit `UNASSIGNED` bucket (reported, never silently merged into one cluster).

**Test strategy:**
- Barcode-pin test; assign known-lineage sentinels (H37Rv → lineage 4; a lineage-2/Beijing barcode pattern → lineage 2); an isolate with no barcode hits → `UNASSIGNED`.

### Step 5: v1b lineage-collapsed scoring orchestrator + cohort/callability gate (deliverable a)
Files: scripts/score_tb_cryptic.py, tests/test_score_tb_cryptic.py
Depends on: Step 3, Step 4

**What changes:**
- **Cohort gate (C3):** build the scored set as the FULL per-drug HIGH-quality eligible cohort (prevalence-preserving — RIF cohort ≠ INH cohort; do NOT artificially balance). The 150R/150S RIF cache is plumbing substrate only. If the eligible cohort isn't fully fetched → emit `TB_SUBSET_PLUMBING` (never the baseline label).
- **Callability denominator (C3):** per drug, partition isolates into scored (callable determinant window) vs `n_uncallable` (ABSTAIN — excluded from sens/spec, reported). A susceptible call requires a callable window with no grade-1/2 determinant — never absence-by-mask.
- **Aggregation (C2a):** build `{preds, labels, clusters}` where `preds` = Step-3 R/S/ABSTAIN, `labels` = CRyPTIC binary phenotype, `clusters` = Step-4 lineage map; call `clonality.cluster_weighted_confusion` → lineage-collapsed sens/spec + `n_discordant` + Wilson CI + `effective_lineage_n`. **M1:** also report `n_clusters_mixed_prediction` (same-label clusters whose member predictions disagree — within-lineage determinant heterogeneity).
- Emit `wiki/tb_rif_inh_cryptic_results_<date>.{md,json}`: lineage-collapsed + raw sens/spec + raw→lineage shrinkage + discordant + mixed-prediction + uncallable counts, tagged `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`. If lineage assignment is unavailable → `LINEAGE_COLLAPSE_BLOCKED_NO_LINEAGE_CALL` (never raw-only as the headline).

**Test strategy:**
- Cohort-gate test (convenience subset → `TB_SUBSET_PLUMBING`, not baseline label); callability test (uncallable window → ABSTAIN/excluded, not S); aggregation reuses `cluster_weighted_confusion` (pin the {preds,labels,clusters} wiring); BLOCKED-status test; honesty label present; raw + lineage both reported.

### Step 6: Frozen-surface leak guard
Files: tests/test_tb_leak_guard.py
Depends on: Step 2, Step 3

**What changes:**
- Assert no CRyPTIC phenotype column (`*_BINARY_PHENOTYPE`/`*_MIC`) is read during rule construction (rules come from the WHO catalogue, not fitted on CRyPTIC labels).
- Assert WHO catalogue + lineage-barcode files are checksum/version-pinned.
- Assert the frozen E. coli surface (`dna_decode/eval/amr_rules.py` + `dna_decode/data/calibrated_amr_rules.json`) is byte-untouched.

**Test strategy:**
- The three assertions; honest-scope docstring (proves no-label-fitting + pinning, NOT biological independence).

### Step 7: Independent post-2023 gold-set scoring harness + acquisition runbook (deliverable b)
Files: dna_decode/organism_rules/tb_goldset.py, scripts/score_tb_independent_goldset.py, tests/test_tb_goldset.py, wiki/tb_independent_goldset_acquisition_2026-06-17.md
Depends on: Step 3

**What changes:**
- `tb_goldset.py`: ingest an independent (non-CRyPTIC) TB cohort (per-isolate VCF/assembly→VCF + WGS-paired DST) into the same `calls` shape Step 3 scores. Independence = from the WHO-catalogue BUILD, i.e. post-2023 isolates (temporal hold-out), since WHO v2 swept most public pre-2023 TB WGS+pDST.
- `score_tb_independent_goldset.py`: reuse the Step 3 cell + Step 4 lineage caller + Step 5 collapse; emit `wiki/tb_rif_inh_independent_results_<date>.md` labelled `INDEPENDENT_VALIDATION` — scored SEPARATELY, never merged with the baseline. BLOCKED-gate → `INDEPENDENT_VALIDATION_BLOCKED_NO_GOLDSET` when absent.
- Acquisition runbook (ratified E): hand-curate a ~30-isolate post-2023 gold set FIRST (public WGS + reference BMD DST, confirmed outside the WHO v2 build) for the first independent number; submit the TB Portals/NIAID DAR in PARALLEL as the larger confirmatory follow-up. Never block deliverable (b)'s first number on the DAR. If hand-curation can't reach N≈30, the DAR is the fallback.

**Test strategy:**
- 2–3 isolate synthetic fixture in the gold-set shape → scores via the Step 3 cell; BLOCKED-status test (no gold set → BLOCKED label, not metrics); independence-label present + distinct from the baseline label.

## Execution Preview

Wave 0 (2 parallel):  Step 1 — TB VCF parser + callability, Step 2 — WHO catalogue join+pin-verify
Wave 1 (2 parallel):  Step 3 — TB decoder cell (v1a), Step 4 — lineage-barcode caller
Wave 2 (3 parallel):  Step 5 — v1b lineage-collapsed score + cohort/callability gate, Step 6 — leak guard, Step 7 — independent gold-set arm

Critical path: Step 1 → Step 3 → Step 5 (3 waves)
Max parallelism: 3 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — **Circularity (not leakage):** WHO catalogue built partly from CRyPTIC → the v1b number is in-distribution. Mitigated by the mandatory `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE` label (Step 5) + the separate independent arm (Step 7). A bare "0.99 on CRyPTIC" headline is forbidden.
- Severity: high — **Lineage-collapse validity:** the v1b number is only as good as the barcode lineage call. The barcode is coarse (lineage/sublineage, not transmission clusters) — honest for lineage disclosure but it does NOT dedup near-identical transmitted isolates within a sublineage. Report `effective_lineage_n`; if it is tiny, Step 5 demotes the cell to "smoke test." De-novo SNP-distance transmission clustering is deferred (v1c) precisely because masked-VCF union-site distance is untrustworthy.
- Severity: high — **Callability denominator [unverified]:** the C3 "ABSTAIN on uncallable window" gate depends on the `REGENOTYPED_VCF` carrying explicit ref/no-call at the determinant windows. The regenotyped-VCF schema is NOT yet inspected — Step 1 must verify it before the callability rule is trusted; if it lacks explicit ref calls, fall back to a conservative ABSTAIN-on-determinant-window-absence + a mask-track check.
- Severity: medium — **External tool / data surfaces:** (1) WHO catalogue — pinned at 0bb39143 + sha256 (verified). (2) lineage-barcode SNP set (Coll-2014 / Napier-2020) — `[unverified]` until the exact positions/alleles are sourced + pinned in Step 4. (3) TB Portals/NIAID independent cohort — access-gated DAR, mixed-method DST `[unverified]` — Step 7 BLOCKED-gates honestly.
- Severity: low — **File overlap:** none within any wave (disjoint files; `organism_rules/__init__.py` created once in Step 1). New-package import edges captured by `Depends on:`.
- Severity: low — **Bedaquiline/clofazimine determinant-invisible** — scoped OUT; first cell is RIF + INH only.

## Open Questions
All ratified 2026-06-17 (Soraya best-judgment per user delegation — USER-OVERRIDABLE before /execute-plan). Folded into the Steps above.
- A. **INH scope — RESOLVED: all WHO grade-1/2 INH loci** (NOT just katG + inhA-promoter). A hand-narrowed subset isn't "the catalogue rule" and manufactures scope-artifact false-negatives; score the full grade-1/2 INH set, report per-locus contribution (katG + inhA-promoter named as dominant).
- B. **v1b cohort cut — RESOLVED: full per-drug prevalence-preserving HIGH-quality cohort**, cache routed to D: (~2.6 GB; C: too full). Bounded subset can't earn the baseline label (C3 gate) and risks too few effective lineages after collapse. The 300-VCF RIF cache stays plumbing-only.
- E. **Independent gold-set route — RESOLVED: hand-curate ~30 post-2023 isolates FIRST** (clean, controllable, no access gate), submit the TB Portals/NIAID DAR in PARALLEL as the larger confirmatory follow-up; never block deliverable (b)'s first number on the DAR. Caveat: hand-curation feasibility (30 clean post-2023 not-in-WHO-build isolates) is [unverified] until attempted; DAR is the fallback if N can't be reached.
- F. **Lineage barcode — RESOLVED: Napier-2020 barcode** from the TBProfiler `tbdb` barcode file, pinned by commit SHA + checksum; apply the SNP positions DIRECTLY to our VCFs (use the DATA, not the TBProfiler caller — avoids the wrapper-trap). Finer sublineage than Coll-2014 → more honest (less aggressive) clonality correction. Coll-2014 is the coarser fallback.
- (Resolved earlier: C → Step 1 `FILTER==PASS`+`GT` rule, no GCP floor. D → M2 masked VCF canonical for calls + REGENOTYPED_VCF for callability.)

## Verification
1. `uv run pytest tests/test_tb_vcf.py tests/test_tb_who_catalogue.py tests/test_tb_amr.py tests/test_tb_lineage.py tests/test_score_tb_cryptic.py tests/test_tb_leak_guard.py tests/test_tb_goldset.py -q` — all green.
2. `uv run pytest tests/ -q` (excluding `tests/test_models_foundation.py`) — frozen E. coli suite 0 regressions; leak guard confirms `amr_rules.py` + `calibrated_amr_rules.json` byte-untouched.
3. v1a: `dna_decode/organism_rules/tb_amr.py` exists; coordinate-alignment fixture passes; uncallable window → ABSTAIN.
4. v1b: `wiki/tb_rif_inh_cryptic_results_<date>.md` reports lineage-collapsed (via `cluster_weighted_confusion`) + raw sens/spec + shrinkage + `n_discordant` + `n_clusters_mixed_prediction` + `n_uncallable` + Wilson CI + `effective_lineage_n`, tagged `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE` — OR a `TB_SUBSET_PLUMBING` / `LINEAGE_COLLAPSE_BLOCKED_NO_LINEAGE_CALL` status.
5. Deliverable (b): `scripts/score_tb_independent_goldset.py` + acquisition runbook exist; independent run labelled `INDEPENDENT_VALIDATION` (or BLOCKED), scored separately from the baseline.

## Save-time amendments

Captured at: 2026-06-17
Source: `/save-plan` arguments

> Audit-notes-only: `/execute-plan` reads ONLY `## Implementation Steps` for executable work. These amendments are provenance, NOT executable instructions. The fixes below are ALREADY folded into the Steps above (this plan is the post-/brainstorm regeneration, not a v1+amendments stack).

- TB CRyPTIC v2: brainstorm fixes folded — FILTER==PASS rule (no GCP), lineage-barcode collapse + cluster_weighted_confusion reuse, cohort+callability gate, in-scope gold-set arm

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high -->
