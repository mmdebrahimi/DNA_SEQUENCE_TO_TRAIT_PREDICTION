## Lens status
- **probe:** applied — grounded the 3 unknowns against code (string-match leakage gate; assembly-only ingestion; PROBAC reads-only + BioSample bridge).
- **brainstorm:** applied (pre-exec, 2 rounds on the Fix-C draft) — accepted 4 grounded contract fixes (C1 organism rationale; C2 ensure_run-by-accession; C3 classify_tier tiers not binary; C4 bidirectional Entrez/ENA resolver). All folded into the Steps below.
- **review:** not applied.

## Problem Statement
Re-validate the FROZEN v0.5.0 deterministic AMR decoder (`call_resistance` for cipro/cef/gent) on an INDEPENDENT, non-US, clinically-measured-MIC E. coli cohort (pilot = Oxford ENA/NCBI `PRJNA604975`; fallback = Spain PROBAC `PRJEB62601`), to test whether the shipped tool generalizes outside its US-NCBI-PD tuning provenance. Results land in a SEPARATE `wiki/external_validation_*.json` namespace with `evidence_tier=external_clinical` and a dedicated roll-up — the FROZEN report card + `compute_lineage_metrics` are NOT modified. The decoder is NOT modified — this validates it.

## Codebase Context
- `scripts/provenance_disjoint_validate.py` — THE scoring core to mirror. Its per-strain loop is `ensure_run(acc, own_runs, gcache, amrfinder_organism, reuse_glob) → call_resistance(main_tsv, drug, organism=registry_organism) → (prediction, y)` collected into `independent_cohort_validate._conf`, with leakage exclusion via `cohort_manifest`. The external arm reuses this exact loop; it differs only in (a) cohort SOURCE (external project + measured MIC vs NCBI-PD metadata) and (b) an ADDED BioSample-level leakage check.
- `scripts/independent_cohort_validate.py::_conf` — takes `(prediction_str, label_int)` pairs (prediction ∈ {"R","S","ABSTAIN"}) → `n_scored/tp/fp/tn/fn/abstain/acc/sens/spec`. REUSE THIS `_conf` (NOT `organism_drug_validate._conf`, wrong shape). ABSTAIN is excluded from `n_scored`.
- `scripts/organism_drug_validate.py::ensure_run(acc, own_runs, gcache, amr_org, reuse_glob)` — downloads BY ACCESSION (`refseq.download_genome(acc)` + `fasta_path(acc)`), runs AMRFinder with `-O amr_org`, returns `main.tsv`. It does NOT accept a FASTA path — the key passed in MUST be a downloadable GCA/GCF accession.
- `dna_decode/eval/amr_rules.py::call_resistance(main_tsv, drug, organism=...)` — frozen decoder, reused UNCHANGED. For E. coli, `calibrated_rule_for` returns None (`Escherichia` has no registry entry, amr_rules.py:165) → the validated `DRUG_RULE` defaults apply (cipro thr 2 QRDR; cef thr 1 + extended-spectrum subclass; gent thr 1 + gent subclass). `organism=` is provenance + (registry lookup, a no-op for E. coli) — NOT a registry activator here.
- `dna_decode/data/mic_tiers.py::{breakpoints_for, classify_tier}` — `breakpoints_for(drug)` returns CLSI/EUCAST cutoffs only; `classify_tier(mics: list[float], distinct_calls: set[str], breakpoints) → {HIGH_R,DECISIVE_R,HIGH_S,DECISIVE_S,BORDERLINE,AMBIGUOUS,CONFLICT,NO_MIC}` is the binary-label authority. Strict pass = HIGH_R/HIGH_S; relaxed pass = + DECISIVE_R/DECISIVE_S.
- `dna_decode/data/refseq.py::download_genome` — GCA/GCF assembly fetch via NCBI Datasets. ASSEMBLY-ONLY. Called internally by `ensure_run`; the external arm does NOT call it directly.
- `dna_decode/eval/cohort_manifest.py::{build_manifest, prior_accessions}` — leakage gate; EXACT accession-STRING match (scans every `data/raw/*/selected.tsv` + `data/processed/*.parquet`). Left UNCHANGED; the new preflight resolves these accessions to BioSamples for an ADDED same-isolate cross-check.
- `dna_decode/eval/clonality.py::{greedy_representative_clusters, cluster_weighted_confusion, wilson_ci, effective_lineage_n}` — clonality MATH, reused DIRECTLY in the new roll-up (NOT the `compute_lineage_metrics` orchestrator, which globs provenance_disjoint_validation_* + requires `_provdisjoint_` dir names).
- `scripts/build_validation_report_card.py` + `scripts/compute_lineage_metrics.py` — FROZEN consumers; key cells by `canonical_cell_key(organism,drug)` → collision risk. NOT modified; external arm uses its own namespace + roll-up (Fix C).
- ENA portal API `https://www.ebi.ac.uk/ena/portal/api/filereport` — `result=read_run` gives `run_accession` + `sample_accession` (BioSample) but NOT `assembly_accession`. NCBI Entrez `esearch/elink` (BioSample↔Assembly) is the primary BioSample→GCA resolver; ENA portal `result=assembly&query=sample_accession=` is the fallback.

### Reusable-Code Survey
- `provenance_disjoint_validate.py` scoring loop — mirrored by the new external scorer (ensure_run → call_resistance → _conf + manifest leakage).
- `_conf` (independent_cohort_validate.py) — reused (correct shape, ABSTAIN-aware).
- `ensure_run` / `_run_amrfinder` (organism_drug_validate.py) — reused; downloads by accession.
- `call_resistance` (amr_rules.py) — frozen decoder, reused unchanged.
- `breakpoints_for` + `classify_tier` (mic_tiers.py) — reused for tiered MIC labeling.
- `cohort_manifest.build_manifest/prior_accessions` — reused to enumerate tuning accessions for the BioSample leakage check.
- `clonality.py` math fns — reused inline in the roll-up.
- None new — searched: scripts/, dna_decode/eval/, dna_decode/data/, graphify-out/ (absent).

## Pre-Change Baseline
- Frozen v0.5.0 decoder validated ONLY on NCBI-PD provenance-disjoint cohorts (10 SCORED cells; report card). No external paper-cohort re-validation exists.
- Full suite: 1056 passed (excluding tests/test_models_foundation.py host-torch limit).
- Leakage gate matches accession strings only → blind to GCA-assembly vs ERR-run same-isolate identity (closed here by a BioSample preflight, beside the frozen gate).

## Verification Signal
- Gate-0 preflight emits for the pilot: assembly-availability (BioSamples with a GCA = FREE vs no-assembly = ASSEMBLY-REQUIRED), persisted `assembly_to_biosample` + `biosample_to_assemblies` maps, a BioSample-level leakage report (FAIL-CLOSED if >5% of tuning accessions are unresolved-to-BioSample, OR any cohort/tuning BioSample overlap, OR Entrez/ENA disagree on a tuning accession → counted unresolved), and a MIC-openness flag.
- A new `wiki/external_validation_<cohort>_<drug>_<date>.json` (schema `external-validation-v1`; `evidence_tier=external_clinical`) per cipro/cef/gent with strict + relaxed `_conf` blocks, per-tier bucket counts, leakage_control, and a cross-country/method caveat.
- A new `wiki/external_validation_report_card.{md,json}` rendering external cells + per-cell cluster-weighted strict sens/spec + Wilson CI (clonality math reused inline).
- New tests green; existing 1056 unchanged; FROZEN files byte-unchanged: amr_rules.py, mic_tiers.py, cohort_manifest.py, build_validation_report_card.py, compute_lineage_metrics.py.

## Implementation Steps

### Step 1: Gate-0 preflight + bidirectional BioSample resolver
Files: scripts/external_cohort_preflight.py, dna_decode/eval/biosample_resolver.py, tests/test_external_cohort_preflight.py, tests/test_biosample_resolver.py
Depends on: none

**What changes:**
- New `biosample_resolver.py` with BOTH directions + a persisted cache `data/raw/_biosample_cache.json` (idempotent; network only on miss; cache records source + per-entry):
  - `runs_for_project(project_acc) → [(run_accession, biosample)]` via ENA portal `filereport result=read_run` (run→BioSample only; NO assembly field).
  - `biosample_to_assemblies(biosample) → [gca,...]` — PRIMARY: NCBI Entrez `esearch/elink` (BioSample→Assembly, works for SAMEA + SAMN); FALLBACK: ENA portal `result=assembly&query=sample_accession=`. Empty list (no linked GCA/GCF) is a VALID outcome, not an error.
  - `assembly_to_biosample(gca) → biosample|None` — Entrez primary, ENA fallback (used by the leakage check).
  - Pure parse/intersection helpers split from the HTTP calls for offline unit tests.
- New `external_cohort_preflight.py`: for a project accession →
  - (a) assembly-availability: per cohort BioSample, FREE if `biosample_to_assemblies` non-empty else ASSEMBLY-REQUIRED; emit both counts.
  - (b) BioSample leakage cross-check: enumerate tuning accessions via `cohort_manifest.build_manifest` → `assembly_to_biosample` each → intersect with cohort BioSamples. FAIL-CLOSED if ANY overlap, OR >5% of tuning accessions unresolved-to-BioSample, OR Entrez/ENA disagree on a tuning accession's BioSample (disagreement counts as unresolved).
  - (c) MIC-openness flag (manual human-confirmed input — the one fact code can't settle).
  - Emits `wiki/external_preflight_<cohort>_<date>.json` with an overall PASS/FAIL verdict.

**Test strategy:**
- Unit (no network): synthetic ENA-TSV + Entrez/ENA mapping fixtures → FREE/ASSEMBLY-REQUIRED counts; leakage intersection (overlap→FAIL; disjoint+resolved→PASS); >5%-unresolved→FAIL-CLOSED; Entrez/ENA disagreement→unresolved; no-assembly BioSample→ASSEMBLY-REQUIRED (not error); cache hit/miss idempotency.

### Step 2: External-cohort MIC → tiered R/S labels
Files: dna_decode/data/external_mic_labels.py, tests/test_external_mic_labels.py
Depends on: none

**What changes:**
- New module. `CANONICAL_DRUG` alias map (cipro/CIP→ciprofloxacin; ceftriaxone/CRO→ceftriaxone; gentamicin/CN/GEN→gentamicin); reject drugs outside {ciprofloxacin, ceftriaxone, gentamicin}.
- Parse a per-isolate MIC table → `{isolate: {canonical_drug: (mics:list[float], distinct_calls:set[str])}}`; tolerant of censored (`>`,`<=`) + unit normalization (record the raw token).
- Label via `mic_tiers.classify_tier(mics, distinct_calls, breakpoints_for(drug))` — NOT raw thresholding. Score only decisive tiers: STRICT = HIGH_R/HIGH_S; RELAXED = + DECISIVE_R/DECISIVE_S. Exclude BORDERLINE/AMBIGUOUS/CONFLICT/NO_MIC.
- Emit per drug: `selected_strict.tsv` + `selected_relaxed.tsv` (accession TAB R/S) AND a `buckets_<drug>.json` with the count per tier (so the excluded fraction is visible before sens/spec). Record the breakpoint version used.

**Test strategy:**
- Unit: alias normalization (incl. reject unknown drug); censored-value handling; tier assignment via classify_tier (HIGH_R/DECISIVE_R/HIGH_S/DECISIVE_S/BORDERLINE/AMBIGUOUS/CONFLICT/NO_MIC fixtures); strict-vs-relaxed selected.tsv contents; bucket counts sum to N.

### Step 3: External-cohort genome resolver (accession list, not FASTA)
Files: dna_decode/data/external_cohort_genomes.py, tests/test_external_cohort_genomes.py
Depends on: Step 1

**What changes:**
- New module. Per scored cohort isolate (from Step 2): BioSample → `biosample_to_assemblies` (Step 1). Returns a `{isolate: gca_accession}` map for the FREE subset (one chosen GCA per BioSample; deterministic pick when multiple) + an `ASSEMBLY_REQUIRED` list (BioSamples with no GCA) — counted + reported, never silently dropped.
- Returns ACCESSIONS, not fasta paths: the downstream scorer feeds each GCA accession to `ensure_run`, which performs the download. No download happens in this module.

**Test strategy:**
- Unit (no network): BioSample→GCA mapping (mocked resolver); multiple-GCA deterministic pick; no-GCA → ASSEMBLY_REQUIRED bucket (not dropped); output keys are GCA accessions consumable by `ensure_run`.

### Step 4: External-cohort scorer (mirror provdisjoint loop, separate artifact)
Files: scripts/external_cohort_revalidate.py, tests/test_external_cohort_revalidate.py
Depends on: Step 1, Step 2, Step 3

**What changes:**
- New script: require preflight PASS (Step 1) → tiered labels (Step 2) → GCA accessions (Step 3) → mirror the `provenance_disjoint_validate` per-strain loop: for each GCA accession `ensure_run(acc, own_runs, gcache, AMRFINDER_ORGANISM, reuse_glob)` → `call_resistance(main_tsv, drug, organism=REGISTRY_ORGANISM)` → collect `(prediction, y)` → `independent_cohort_validate._conf`.
- `AMRFINDER_ORGANISM` + `REGISTRY_ORGANISM` are taken VERBATIM from the frozen E. coli cipro/cef/gent provdisjoint runs (read from the committed `wiki/provenance_disjoint_validation_*ecoli*` artifacts / the provdisjoint invocation), NOT invented — so external numbers are comparable to the frozen cells.
- Score BOTH strict and relaxed label sets; write `wiki/external_validation_<cohort>_<drug>_<date>.json` (`_schema=external-validation-v1`, `evidence_tier=external_clinical`, organism/drug, `strict`+`relaxed` `_conf` blocks, `buckets`, `leakage_control`, `independence_tier` with cross-country/method caveat) + `data/raw/<cohort>_extval_<drug>/selected_{strict,relaxed}.tsv`.
- Does NOT write `provenance_disjoint_validation_*.json` (no collision). Fail-closed if preflight FAILED or MIC-openness=controlled, unless `--allow-degraded` (stamps `independence_degraded=true`).

**Test strategy:**
- Unit: end-to-end on a tiny fixture cohort (mocked ensure_run + call_resistance) → artifact carries external-validation-v1 + evidence_tier + strict/relaxed blocks + buckets; uses `_conf`'s `n_scored/acc/sens/spec` shape with ABSTAIN handling; organism strings passed verbatim; fail-closed paths; selected_{strict,relaxed}.tsv written.

### Step 5: External-validation roll-up + inline clonality
Files: scripts/build_external_validation_report.py, tests/test_build_external_validation_report.py
Depends on: Step 4

**What changes:**
- New roll-up: glob `wiki/external_validation_*.json` → render `wiki/external_validation_report_card.{md,json}` (separate from the frozen decoder report card). Per cell: raw strict + relaxed acc/sens/spec + cluster-weighted STRICT sens/spec + Wilson CI computed INLINE via `clonality.py::{greedy_representative_clusters, cluster_weighted_confusion, wilson_ci}` on the cohort genomes (Mash via the existing Docker path). NOT via compute_lineage_metrics.
- Honest tier string: external clinical cohort, different country/lab/method than tuning; lineage-disclosed; strict-tier primary, relaxed secondary.

**Test strategy:**
- Unit: external_validation_*.json fixtures → roll-up MD/JSON carries cells + strict/relaxed + CI + evidence_tier; clonality math invoked on a synthetic distance matrix; assert frozen `build_validation_report_card` output unchanged when run alongside.

### Step 6: Docs + ledger
Files: CLAUDE.md, README.md
Depends on: Step 4, Step 5

**What changes:**
- CLAUDE.md gotcha: the external-cohort re-validation arm (separate external_validation namespace; mirrors provdisjoint loop; organism strings verbatim from frozen cells; bidirectional Entrez/ENA BioSample resolver + fail-closed >5%; classify_tier strict/relaxed scoring; frozen decoder + report card + lineage UNCHANGED).
- README: command entries for `external_cohort_preflight.py`, `external_cohort_revalidate.py`, `build_external_validation_report.py`.

**Test strategy:**
- Docs only. Final full-suite run captured in Verification.

## Execution Preview
- Wave 0: Step 1, Step 2.
- Wave 1: Step 3 (deps Step 1).
- Wave 2: Step 4 (deps 1,2,3).
- Wave 3: Step 5 (deps 4).
- Wave 4: Step 6 (deps 4,5).
- Total waves: 5. Max parallelism: 2 (Wave 0). Critical path: 1→3→4→5→6 (length 5).

## Risk Flags
- **GCA availability unconfirmed for the pilot** [unverified]: if Oxford `PRJNA604975` BioSamples mostly lack linked GCAs, Step 3 routes them to ASSEMBLY-REQUIRED and the free pilot N collapses. Step 1 preflight is the wave-0 go/no-go; if it fails on Oxford, swap to Spain/Denmark BEFORE building further. EXECUTE STEP 1 FIRST.
- **Exact organism triple must be read, not invented** [grounded]: Step 4's `AMRFINDER_ORGANISM` + `REGISTRY_ORGANISM` MUST come from the frozen E. coli provdisjoint artifacts. Hardcoding a guessed `-O` string would make the external numbers non-comparable to the frozen cells. Verified two distinct strings exist (group/registry vs amrfinder_organism) in `provenance_disjoint_validate.py`.
- **MIC table may be MTA-gated** [grounded]: per-isolate MIC could need a data-access agreement even with open genomes → Step 4 fail-closes; cohort drops to user-gated. Surfaced by Step 1's MIC-openness flag.
- **Breakpoint version drift** [inferred]: cohort MIC may use a different EUCAST/CLSI year than mic_tiers (CLSI 2024/EUCAST 14.0); Step 2 records the version so a disagreement is attributable to breakpoint vs decoder.
- **Mash/Docker dependency for Step 5 clonality** [grounded]: inline clonality needs the Docker Mash path; on a non-Docker host Step 5 degrades to raw (uncorrected) sens/spec with a flag.

## Open Questions
- Whether to cross-link the external_validation_report_card from the frozen decoder_validation_report_card (cosmetic; deferred).
- Eventual Fix A (key the frozen consumers by (organism,drug,evidence_source)) for the multi-cohort breadth state — deferred until the pilot proves the path.

## Verification
- `uv run pytest tests/ -q --ignore=tests/test_models_foundation.py` → 0 regressions + new tests green.
- Manual (Gate-0): run Step 1 preflight on PRJNA604975 → assembly-availability + BioSample leakage + MIC-openness verdict (the empirical go/no-go).
- Manual (full): run the scorer on the FREE GCA subset → external_validation_*.json (strict + relaxed); build the external roll-up → external cells + CI render.
- Frozen invariant: git diff shows amr_rules.py + mic_tiers.py + cohort_manifest.py + build_validation_report_card.py + compute_lineage_metrics.py byte-unchanged.

## Save-time amendments

Audit-notes-only: this block is provenance for human readers. `/execute-plan` reads ONLY `## Implementation Steps`. If an amendment changes a Step's contract, re-run `/technical-plan` before `/execute-plan`.

Captured at: 2026-06-15
Source: `/save-plan` arguments (re-emit superseding the first save; Steps now carry the fixes literally)

- re-emit: C1-C4 folded into literal Steps
- mirror provdisjoint loop
- ABSTAIN-aware _conf
- organism triple read-not-invented

Captured at: 2026-06-14 (superseded — recorded for history; these were the first-save intent, now baked into the Steps)
Source: `/save-plan` arguments

- Fix C namespace separation
- real BioSample resolver + fail-closed >5%
- pinned _conf + canonical drug map + exact organism key
- evidence_tier=external_clinical

## Repo grounding

### Captured by: brainstorm @ 2026-06-14

Files read: dna_decode/eval/amr_rules.py, scripts/organism_drug_validate.py, dna_decode/data/mic_tiers.py, dna_decode/data/refseq.py, scripts/independent_cohort_validate.py, scripts/provenance_disjoint_validate.py, plans/Oxford_Cohort_External_Revalidation_Plan/technical-plan.md

Key claims:
- [grounded] `call_resistance` does NOT activate a calibrated registry for E. coli — `calibrated_rule_for` returns None for `Escherichia` (amr_rules.py:165) and falls through to `DRUG_RULE`. The DRUG_RULE defaults (cipro thr 2 QRDR; cef thr 1 + extended-spectrum subclass; gent thr 1 + gent subclass) ARE the validated E. coli path. `organism=` is for AMRFinder `-O` / provenance comparability only — must match the frozen 10 cells' `-O` string.
- [grounded] `ensure_run(acc, own_runs, gcache, amr_org, reuse_glob)` (organism_drug_validate.py:95) downloads by accession via `refseq.download_genome(acc)` + `fasta_path(acc)`; it does NOT accept a FASTA path. Step 3 returns GCA accessions (FREE subset) + an ASSEMBLY-REQUIRED bucket, not `{accession: fasta_path}`.
- [grounded] `breakpoints_for` returns cutoffs only (mic_tiers.py:71); `classify_tier` (mic_tiers.py:92) already models HIGH_R/DECISIVE_R/HIGH_S/DECISIVE_S/BORDERLINE/AMBIGUOUS/CONFLICT/NO_MIC. Step 2 tier-labels + scores only decisive tiers (strict=HIGH_R/HIGH_S primary; relaxed=+DECISIVE secondary) + reports excluded bucket counts — NOT naive binary thresholding.
- [grounded] ENA `read_run` filereport gives run→sample(BioSample) only; Step 1 needs BOTH directions — `assembly_to_biosample` (leakage check) AND `biosample_to_assemblies` (Step 3 genome resolution). Primary resolver = NCBI Entrez esearch/elink BioSample→Assembly; ENA portal `result=assembly&query=sample_accession=` as fallback. A BioSample with no linked GCA/GCF → ASSEMBLY-REQUIRED (cohort outcome), not a resolver error.
- [grounded] `provenance_disjoint_validate.py` is the canonical scoring core (ensure_run → call_resistance(organism=reg_org) → independent_cohort_validate._conf + cohort_manifest leakage); the external scorer mirrors it. `_conf` is ABSTAIN-aware on `(prediction_str, label_int)` pairs.
<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified -->
