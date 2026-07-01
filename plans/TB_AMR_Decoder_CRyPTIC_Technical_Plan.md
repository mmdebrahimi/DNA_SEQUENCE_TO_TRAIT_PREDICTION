# TB AMR Decoder (RIF + INH) on CRyPTIC — Technical Plan

**Status:** candidate (pre-execution; pre-exec /brainstorm incorporated 2026-06-16)
**Date:** 2026-06-16
**Owner:** Soraya-driven (attended); author dna_decode
**Pipeline:** /idea-anchor + /probe (done) → /technical-plan (this) → pre-exec /brainstorm (done — findings folded in) → /save-plan → /execute-plan (or `/soraya --until-mvp`)

## Goal

Ship the FIRST M. tuberculosis AMR decoder cell for `dna_decode` — deterministic, interpretable, two
mechanisms: **RIF** (`rpoB` RRDR) and **INH** (`katG` + `inhA` promoter) — scored over the **CRyPTIC reuse
cohort** (~12,287 isolates) using the **WHO M. tuberculosis mutation catalogue v2 (2023)** as the
determinant source. Lineage-collapsed metrics. Honestly labeled a knowledge-baseline.

## Acceptance bar (ratified 2026-06-16)

- Knowledge-baseline first; **independent validation DEFERRED** (post-2023 temporal hold-out / gold set is a
  named follow-up, never conflated with this metric).
- Honesty label (load-bearing): all results tagged **`WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`** — the
  WHO catalogue was partly built from CRyPTIC, so scoring on CRyPTIC is a knowledge-baseline, NOT independent
  validation.

## Hard boundaries

- **FROZEN surface — never edited:** `dna_decode/eval/amr_rules.py`, `dna_decode/data/calibrated_amr_rules.json`
  (reproducibility freeze 2026-06-13). TB rules live in their OWN new module + registry.
- New code: `dna_decode/organism_rules/tb_amr.py` (+ a TB rule registry, + tests). Mirrors the TMP-SMX
  experimental-overlay pattern (`dna_decode/data/experimental_drug_rules.py`).
- Deterministic rules only — NO learned/embedding model.

## Why this plan is STAGED (pre-exec /brainstorm findings)

The naive framing ("write determinant rules + score the local CSV") is **wrong**: the CRyPTIC reuse CSV is
*only* `*_BINARY_PHENOTYPE`/`*_MIC`/`*_PHENOTYPE_QUALITY` + IDs (12,288 rows) — **no variant data**; and
`data/raw/cryptic/vcf_cache/` holds only **30** masked VCFs (the RIF PoC), not 12,287. Lineage collapse has
**no existing substrate** (`eval/clonality.py` builds its matrix via Mash-on-FASTA through `eval/phylogeny.py`;
TB has only masked VCFs, no FASTAs). The real work is a data-engineering pipeline with several correctness
traps. Hence four gated stages, each with its own checkable artifact, and a `BLOCKED` status rather than fake
metrics when substrate is missing.

---

## Stage 0 — VCF acquisition + genotype/filter/callability spec (the foundation)

**Goal:** turn CRyPTIC VCF *paths* into genotype-aware, normalized variant calls + a documented callability spec.

**Approach:**
1. **Acquisition:** stream masked VCFs from the EBI/CRyPTIC FTP into `data/raw/cryptic/vcf_cache/` (gitignored).
   v1 uses a **bounded subset** (see Open Tradeoff B) — NOT all 12,287 up front. Record provenance (URL, fetch
   date, file count).
2. **Genotype-aware parser** (CRITICAL — C2): parse `GT` allele indices + FILTER policy. A cached VCF contains
   `0/0`-with-ALT, `./.`, multi-allelic records, `mask-compass`/`MIN_FRS` filters. **"position+ALT present" is
   NOT a resistance call** (the existing `rrdr_variant_present` shape would mis-score). Declared default call
   policy (Open Tradeoff C): a site is a non-reference call iff `FILTER==PASS` AND `GT` is non-reference (allele
   index ≥1) AND a minimum quality floor (FRS/DP/GCP) holds. Normalize MNVs/indels (essential for `inhA`
   promoter + `katG`).
3. **Callability spec** (C-M2): define which records count as "callable" for downstream SNP distance — masked
   VCFs are union-genotyped variant sites, not whole-genome; state explicitly that pairwise distance is over
   polymorphic union sites (usable for clustering, ≠ full callable-genome SNP distance).

**Files:** `dna_decode/organism_rules/tb_vcf.py` (acquisition + genotype-aware parse + callability), tests.
**Tests:** unit tests over ≥1 cached VCF asserting correct handling of `0/0`-with-ALT, `./.`, multi-allelic,
filtered records; a sentinel test that a known rpoB/katG variant parses as a non-reference call and a known
reference site does not.
**Acceptance:** parser passes the genotype/filter unit suite; callability spec documented in-module.

## Stage 1 — WHO catalogue join + pin (the rule source)

**Goal:** a pinned, joined WHO determinant table for RIF + INH.

**Approach (CRITICAL — C4):** the WHO repo `GTB-tbsequencing/mutation-catalogue-2023` has **no releases**;
output files changed Feb/May 2024 while the master didn't. So:
1. Fetch into `data/raw/who_tb_catalogue/`; **pin a commit SHA + per-file checksums** for every file used
   (record in the artifact + a checksum manifest). NOT "download latest v2".
2. **JOIN** (not a file read): grade-bearing Excel/master rows ↔ the genomic-coordinate VCF (which "maps
   coordinates to variant names only"). Unit-test the join on RIF + INH sentinel determinants.
3. Build the TB rule registry from the join: RIF = `rpoB` RRDR grade-1/2 determinants; INH = `katG` +
   `inhA`-promoter grade-1/2 determinants (scope: Open Tradeoff A).

**Files:** `dna_decode/data/tb_who_catalogue.py` (loader + checksum-pin + join), `data/raw/who_tb_catalogue/CHECKSUMS`.
**Tests:** join correctness on sentinels; checksum-pin test (fails if catalogue files change unpinned).
**Acceptance:** joined determinant table built + pinned; sentinel join tests pass.

## Stage 2 — v1a "plumbing cell" (sentinel fixture, NO metric claims)

**Goal:** end-to-end determinant scoring on a curated fixture — proves the pipeline, NOT cohort performance.

**Approach (M1):** assemble **20–50 handpicked sentinel VCFs** (known RIF/INH determinants + known negatives;
NOT a first-row phenotype sample). Run: Stage-0 parse → Stage-1 determinant match → R/S call per drug. The
artifact reports per-sentinel calls + the honesty label. **Explicitly NOT sens/spec beyond fixtures.**

**Files:** `dna_decode/organism_rules/tb_amr.py` (the decoder cell: `score_rif`, `score_inh`), `tests/test_tb_amr.py`.
**Tests:** fixture-driven — each sentinel scores as expected; abstain/partial-call behavior on out-of-scope
variants; `coverage_scope` + `excluded_grade12_loci` reported (M5).
**Acceptance (v1a checkable):** `file-exists dna_decode/organism_rules/tb_amr.py` · `test-exit-0` the cell
tests · `file-exists` a v1a plumbing artifact (sentinel calls + coverage scope + honesty label).

## Stage 3 — v1b lineage-collapsed score (the real metric, BLOCKED-gated)

**Goal:** lineage-collapsed sens/spec for RIF + INH over the bounded cohort subset.

**Approach:**
1. **SNP-distance matrix from VCFs** (C3 — chosen: direct VCF-derived, NOT consensus-FASTA+Mash): build pairwise
   SNP distances over callable polymorphic union sites (per Stage-0 spec) into a provenance-stamped artifact.
   Prefer **sparse threshold-neighbor clustering** over a dense 12k×12k matrix (~1.2 GB float64) — M3.
2. **Lineage collapse:** feed the matrix to the EXISTING `clonality.greedy_representative_clusters_from_matrix`.
3. **Intra-cluster aggregation** (Open Tradeoff D — chosen: **representative-isolate dedup**, phenotype-blind,
   by lowest missingness / ID order). Do **NOT** majority-vote R/S for TB (intra-cluster discordance may be real
   acquired resistance). Also report raw metrics + discordant-cluster counts.
4. **BLOCKED-gating (C3):** if VCF acquisition or a provenance-stamped SNP-distance artifact is missing, emit
   **`LINEAGE_COLLAPSE_BLOCKED_NO_DISTANCE_MATRIX`** — NEVER raw/fake metrics.

**Files:** `dna_decode/organism_rules/tb_snp_distance.py` (matrix builder + provenance), `scripts/score_tb_cryptic.py`
(orchestrator), `tests/test_tb_snp_distance.py`.
**Tests:** matrix builder unit tests (missing/het/masked handling); dedup determinism; BLOCKED-status test
(no matrix → BLOCKED, not metrics).
**Acceptance (v1b checkable):** `file-exists` a results artifact (e.g. `wiki/tb_rif_inh_cryptic_results_<date>.md`)
reporting lineage-collapsed sens/spec + raw metrics + discordant-cluster counts + the honesty label, OR a
documented `BLOCKED` status.

---

## Leak test (frozen regression — strengthened per brainstorm)

`tests/test_tb_leak_guard.py`:
1. Assert **no CRyPTIC phenotype columns** (`*_BINARY_PHENOTYPE`/`*_MIC`) are read during rule construction
   (rules come from the WHO catalogue, not fitted on CRyPTIC labels).
2. Assert the WHO catalogue files are **pinned by checksum**.
3. Assert the **frozen E. coli surface is byte-untouched** (`amr_rules.py` + `calibrated_amr_rules.json`).
> Honest scope: this proves no-label-fitting + pinning; it CANNOT prove biological independence (that's the
> deferred independent-validation follow-up).

## Open tradeoffs — drafted defaults for RATIFICATION (domain-science calls)

- **A. INH scope:** *default* = `katG` + `inhA`-promoter WHO grade-1/2 only, with `coverage_scope` +
  `excluded_grade12_loci` reported (honest partial coverage). *Alt:* all WHO grade-1/2 INH loci. **Ratify.**
- **B. v1 cohort cut:** *default* = bounded subset for v1b (not all 12,287 up front), scaled later. *Alt:*
  full cohort in v1b. **Ratify.**
- **C. VCF call policy:** *default* = `FILTER==PASS` + `GT` non-reference + a min-quality floor; MNV/indel
  normalized. **Ratify the floor values** (FRS/DP/GCP).
- **D. SNP-distance threshold + canonical VCF:** *default* = masked VCF (matches cache) as canonical, tested
  vs rpoB/katG/inhA sentinels; SNP threshold predeclared (~5–12 SNPs, TB transmission range) + a sensitivity
  check. **Ratify** the canonical VCF (masked vs REGENOTYPED) + the threshold.

## Deferred (named, not silently dropped)

- Independent validation (post-2023 temporal hold-out / hand-curated gold set) — the next milestone after this.
- Full 12,287-cohort scale beyond the v1b subset.
- CRyPTIC lineage/sublineage metadata as a fallback stratification (if SNP-distance clustering underperforms).
- Additional TB drugs/mechanisms beyond RIF + INH.

## Test plan summary

Stage 0 genotype/filter/callability unit suite · Stage 1 catalogue-join + checksum-pin tests · Stage 2
sentinel-fixture cell tests · Stage 3 matrix-builder + dedup-determinism + BLOCKED-status tests · the leak guard.
All new tests under `tests/`; frozen E. coli suite must stay green (0 regressions).

## Execution note

Soraya drives this attended via `/soraya --until-mvp` against a tb-decoder family ledger whose MVP bar = the
v1a + v1b acceptance predicates above. The money/destructive gates apply unchanged; VCF/catalogue fetches are
`auto` (data download, not money). No frozen-surface edits.
