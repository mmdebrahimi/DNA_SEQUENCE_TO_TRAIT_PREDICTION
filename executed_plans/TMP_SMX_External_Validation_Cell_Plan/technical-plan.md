## Lens status
- Mode: conversation-driven (post-probe + post-brainstorm synthesis). **Scope divergence (intentional):** the invocation arg said "trimethoprim-sulfa AND-rule + cefepime; non-frozen overlay catalog" — the `/brainstorm` revised this to **TMP-SMX ONLY** (cefepime dropped: 0.986 ceftriaxone-concordant + the reused {CEPH,CARB} rule mis-calls the 30 distinguishing cef-R/fep-S isolates + CLSI SDD 4-8 gray zone) and **scorer-local in the EXISTING external-validation arm** (not a new overlay catalog: the `(sul AND dfr)` rule is a new shape frozen `DRUG_RULE` cannot express, and in-hand-cohort validation already lands in the non-frozen `external_validation_*` namespace).
- Repo index: unconfigured — direct file reading.
- sentrux: n/a (not installed). gh + codex present.
- project-rules.md / DESIGN.md: absent — Project Rules Check omitted; no UI.

## Problem Statement
The frozen deterministic AMR decoder covers 6 drugs. Two in-hand independent measured-MIC cohorts (Oxford ~2900, Sci234 234) carry co-trimoxazole MIC + complete genotype, and trimethoprim-sulfamethoxazole (TMP-SMX) is a genuinely-new resistance MECHANISM (acquired folate-pathway genes: `sul*` + `dfrA/B*`) distinct from all 6 frozen cells. Add TMP-SMX as an honestly-validated EXPERIMENTAL cell, scored by a scorer-local `(>=1 sul) AND (>=1 dfr)` rule, validated on BOTH cohorts, branded distinct from deployed claims, with NO edit to any frozen file. Success: a strict-tier TMP-SMX cell on both cohorts whose mechanistic strata reproduce across cohorts; frozen files byte-unchanged; the cell renders in the external report card marked EXPERIMENTAL_SCORED. Non-goals: cefepime (deferred — needs a cefepime-specific rule + SDD abstain), promotion to the frozen deployed surface, editing the frozen 6-drug catalogs/engine.

## Codebase Context
- `dna_decode/eval/amr_rules.py` (FROZEN): `DRUG_RULE` = `threshold + subclass_any` (count/OR semantics) — CANNOT express AND-across-two-gene-families; do not edit.
- `dna_decode/data/mic_tiers.py` (FROZEN): `classify_tier(mics, distinct_calls, breakpoints)` is a PURE function taking breakpoints as a param (reusable with non-frozen cotrimoxazole breakpoints — no frozen edit); `STRICT_MIC_TIERS={HIGH_R,HIGH_S}`, `RELAXED` adds DECISIVE; `breakpoints_for`/`supported_drugs` cover only the 6 frozen drugs (TMP-SMX absent — confirms scorer-local).
- `scripts/build_external_validation_report.py` (NON-frozen, Fix-C arm): `build_cell` reads TOP-LEVEL `drug/strict/relaxed`; run-scoped roll-up via `--run-id`/`--artifacts`; lineage block degrades to "unavailable" for reads-only cohorts (Sci234/Oxford have no GCA FASTAs on disk). Editable.
- `scripts/sci234_score.py` (mine, 2026-06-16): `load_cohort` already extracts per-isolate acquired genes (incl. sul/dfr) + MIC by study Key. Reuse its loader.
- `scripts/oxford_score.py` (mine): `group_amrfinder` groups `amrfinder.tsv` rows by guuid + header-normalizes `Gene symbol -> Element symbol`. Reuse for Oxford genotype.
- `scripts/independent_cohort_validate.py`: `_conf(pairs)` → n/tp/fp/tn/fn/acc/sens/spec. Reuse.
- `dna_decode/data/external_mic_labels.py` (non-frozen): `build_drug_labels` raises for non-PILOT drugs (cipro/cef/gent) — do NOT reuse for TMP-SMX; the scorer tiers via `classify_tier` + explicit cotrimoxazole breakpoints instead.
- Data (gitignored): Oxford `amrfinder.tsv` (3539 sul/dfr rows) + `main_data.csv` (`Cotrim_upper`, MIC=2^upper, TMP-component units {1,2,4,32}); Sci234 supplement (`sul*`/`dfrA-B*` presence columns + `TRISUL` MIC).
- Breakpoint: co-trimoxazole Enterobacterales (EUCAST v14.0 / CLSI M100 2024) is expressed as the TRIMETHOPRIM COMPONENT: S<=2 / R>=4. Both cohort columns report TMP-component units (verified: Oxford Cotrim {1,2,4,32}).

### Reusable-Code Survey
- `scripts/oxford_score.py::group_amrfinder` — Oxford genotype grouping + header normalization. REUSE.
- `scripts/sci234_score.py::load_cohort` — Sci234 acquired-gene + MIC loader. REUSE.
- `scripts/independent_cohort_validate.py::_conf` — confusion/metrics. REUSE.
- `dna_decode/data/mic_tiers.py::classify_tier` (frozen, breakpoints-param) — strict/relaxed tiering. REUSE without edit.
- `scripts/build_external_validation_report.py` — external roll-up. EXTEND (non-frozen) for branding passthrough.
- `scripts/fam_subclass_resolver.py` — NOT needed (sul->SULFONAMIDE, dfr->TRIMETHOPRIM are direct classes; no subclass refinement).
- None — searched: graphify-out/GRAPH_REPORT.md (absent), src/lib/utils dirs (n/a — flat package).

## Pre-Change Baseline
- TMP-SMX currently has ZERO decoder cells (net-new). The invariant to PRESERVE: the 6 frozen-drug catalogs + report card + lineage modules byte-unchanged (reproducibility freeze, commit b3761c8). Measured prior signal (brainstorm, Sci234, TMP-component R>=4/S<=2): sul+dfr 69/70 R; sul-only 0/40 R; dfr-only 1/6 R; neither 1/117 R — the AND rule's expected separation, to be reproduced on Oxford. Test suite green at session start (126 amr/mic + 20 resolver + external-arm tests).

## Verification Signal
- A strict-tier TMP-SMX `external_validation_*` cell on BOTH Sci234 and Oxford with sens/spec + n.
- THE GATE (de-risk): the 4 genotype strata (sul+dfr / sul-only / dfr-only / neither) reproduce the Sci234 pattern on Oxford — sul+dfr high R-rate, the other three low — computed + asserted by the scorer. If Oxford contradicts (e.g. sul-only often R), the cell is emitted as INDETERMINATE, not SCORED.
- Frozen files byte-unchanged: `git diff --stat` shows no change to `dna_decode/eval/amr_rules.py`, `dna_decode/data/mic_tiers.py`, `dna_decode/eval/cohort_manifest.py`, `scripts/build_validation_report_card.py`, `scripts/compute_lineage_metrics.py`.
- No frozen leak: `trimethoprim-sulfamethoxazole` absent from `mic_tiers.supported_drugs()` AND `amr_rules.DRUG_RULE` (asserted by test).
- Artifact carries `rule_status=EXPERIMENTAL_SCORED`, `rule_scope=scorer_local`, `not_in_shipped_surface=true`, exact rule text + strata; renders in `external_validation_report_card` visibly distinct from frozen-decoder external cells.
- Full `uv run pytest tests/ -q` no new failures.

## Implementation Steps

### Step 1: Non-frozen TMP-SMX rule + cotrimoxazole breakpoint module
Files: dna_decode/data/experimental_drug_rules.py
Depends on: none

**What changes:**
- New non-frozen module. `COTRIMOXAZOLE_BREAKPOINTS = {"clsi_r":4.0,"clsi_s":2.0,"eucast_r":4.0,"eucast_s":2.0}` (TMP-component; docstring cites EUCAST v14.0 + CLSI M100 2024, flagged for authoritative pin).
- Explicit family regexes (M1): `SUL_RE = re.compile(r"^sul[1-4]$", re.I)`; `DFR_RE = re.compile(r"^dfr[AB]\d+", re.I)` applied AFTER stripping allele/suffix decorations; both EXCLUDE regulators/look-alikes (e.g. `sulR`, bare `dfrA` w/o digit) by default.
- `tmp_smx_call(gene_symbols)`: R iff `(any SUL_RE) AND (any DFR_RE)` else S; returns `{prediction, matched_sul, matched_dfr, rule_text}`. This is the scorer-local AND rule (the new shape).
- `cotrimoxazole_tier(mic_tokens, distinct_calls)`: wraps frozen `mic_tiers.classify_tier(..., COTRIMOXAZOLE_BREAKPOINTS)`.
- Branding constants: `RULE_STATUS="EXPERIMENTAL_SCORED"`, `RULE_SCOPE="scorer_local"`, `DRUG="trimethoprim-sulfamethoxazole"`.

**Test strategy:**
- Unit (Step 4): AND truth table; regex excludes `sulR`/bare-`dfrA`; tiering maps MIC 32->HIGH_R, 1->HIGH_S; breakpoint values pinned.

### Step 2: TMP-SMX external scorer (both cohorts, single-cell top-level artifacts)
Files: scripts/tmp_smx_external_validate.py
Depends on: Step 1

**What changes:**
- New scorer. For Sci234: reuse `sci234_score.load_cohort` for per-isolate acquired genes + TRISUL MIC. For Oxford: reuse `oxford_score.group_amrfinder` for per-guuid Element symbols + `Cotrim_upper` MIC (2^upper).
- Per isolate: `tmp_smx_call(gene_symbols)` -> prediction; `cotrimoxazole_tier` -> strict/relaxed tier labels; binary R/S (MIC>=4 / <=2). Compute `_conf` for strict + relaxed + binary.
- Compute the 4-strata breakdown (sul+dfr / sul-only / dfr-only / neither: n, R, S, R-rate) per cohort; set `strata_reproduced` True iff sul+dfr R-rate is the max stratum AND sul-only R-rate < 0.5 (the gate).
- Emit ONE TOP-LEVEL `external-validation-v1` artifact PER cohort: `wiki/external_validation_<cohort>tmpsmx_<run_id>_<date>.json` with top-level `drug="trimethoprim-sulfamethoxazole"`, `strict`/`relaxed`/`binary` cells, `run_id`, `cohort`, `organism`, `rule_status`/`rule_scope`/`not_in_shipped_surface`/`rule_text`, `strata`, `strata_reproduced`, `independence_tier`, `leakage_status`, `fidelity_caveat` (folP/folA point-mutation TMP-R is a curated-determinant blind-spot). If `strata_reproduced` is False -> stamp headline INDETERMINATE (not a SCORED claim).
- `--run-id` arg; print per-cohort cells + strata; exit 0.

**Test strategy:**
- Unit (Step 4): synthetic 2-cohort fixture -> strata math + artifact schema (top-level drug/strict + branding + run_id) + `_conf` correctness + INDETERMINATE stamping when strata fail.

### Step 3: External report-card branding passthrough (non-frozen builder)
Files: scripts/build_external_validation_report.py
Depends on: none

**What changes:**
- `build_cell`: pass through `rule_status`, `rule_scope`, `not_in_shipped_surface`, `strata_reproduced` when present (default None -> frozen-decoder external cell, unchanged behavior).
- `render_md`: add a `scope` column (or trailing marker) so an `EXPERIMENTAL_SCORED`/`scorer_local` cell is visibly distinct from frozen-decoder external cells; existing cells render unchanged (marker blank).
- No change to globbing/run-scoping/lineage; additive only.

**Test strategy:**
- Covered in Step 4 (the builder test is listed in Step 4's Files): `build_cell` carries branding fields; `render_md` shows the experimental marker for a branded cell and leaves a frozen-decoder cell unmarked.

### Step 4: Unit tests for the rule module + scorer + builder passthrough + frozen-leak guard
Files: tests/test_experimental_drug_rules.py, tests/test_tmp_smx_external_validate.py, tests/test_build_external_validation_report.py
Depends on: Step 1, Step 2, Step 3

**What changes:**
- `test_experimental_drug_rules.py`: AND truth table (sul+dfr->R; sul-only/dfr-only/neither->S); regex exclusions (sulR, dfrB w/o digit, intI1); cotrimoxazole tiering + breakpoint pins.
- `test_tmp_smx_external_validate.py`: synthetic-fixture strata + artifact schema + branding + run_id + INDETERMINATE-on-strata-fail; FROZEN-LEAK GUARD: assert `"trimethoprim-sulfamethoxazole" not in mic_tiers.supported_drugs()` and `not in amr_rules.DRUG_RULE`.
- `test_build_external_validation_report.py` (extend): `build_cell` carries the Step-3 branding fields; `render_md` marks a branded cell experimental and leaves a frozen-decoder cell unmarked.

**Test strategy:**
- `uv run pytest tests/test_experimental_drug_rules.py tests/test_tmp_smx_external_validate.py tests/test_build_external_validation_report.py -q` green; full suite no new failures.

### Step 5: Run validation on both cohorts + run-scoped roll-up + result memo
Files: wiki/external_validation_tmpsmx_result_2026-06-16.md
Depends on: Step 2, Step 3

**What changes:**
- Run `python -m scripts.tmp_smx_external_validate --run-id tmpsmx<date>` -> emit the two cohort artifacts.
- Run `python scripts/build_external_validation_report.py --run-id tmpsmx<date> --no-clonality` -> render the external report card with the two TMP-SMX cells (lineage unavailable for reads-only cohorts).
- Write the result memo (companion to the Oxford/Sci234 memos): per-cohort strict/binary cells, the 4-strata table for BOTH cohorts, strata-reproduction verdict, frozen-untouched confirmation, fidelity caveat, EXPERIMENTAL_SCORED scope statement.

**Test strategy:**
- Manual: verify both artifacts emitted + report card renders 2 distinct EXPERIMENTAL cells; `git diff --stat` shows zero frozen-file changes; strata reproduce across cohorts.

## Execution Preview
- Wave 0 (parallel x2): Step 1 (rule module) + Step 3 (builder passthrough) — independent files.
- Wave 1: Step 2 (scorer) — needs Step 1.
- Wave 2 (parallel x2): Step 4 (tests) — needs 1+2+3; Step 5 (run+memo) — needs 2+3.
- Total waves: 3. Max parallelism: 2. Critical path: Step 1 -> Step 2 -> Step 4 (toolkit); Step 5 (run+memo) is the parallel wave-2 leaf.

## Risk Flags
- **Breakpoint authority [unverified-in-plan]:** co-trimoxazole TMP-component S<=2/R>=4 is well-established but must be pinned against EUCAST v14.0 / CLSI M100 2024 before the cell is treated as a claim (Step 1 docstring + this flag). Oxford Cotrim units verified TMP-component.
- **folP/folA blind-spot:** target point-mutation TMP-R is invisible to the acquired-gene AND rule (the sul+dfr-negative R isolates: Sci234 neither 1/117, dfr-only 1/6). Bounded FN, documented in `fidelity_caveat`.
- **dfr-only stratum n=6 (Sci234) is tiny:** Oxford strata reproduction is the de-risk; the scorer emits INDETERMINATE (not SCORED) if Oxford contradicts.
- **Shared genotype-caller vocabulary:** both cohorts use curated acquired-gene callers — independence is partial (independent of MY pipeline, not of the caller class); documented.
- **cefepime DEFERRED (arg divergence):** dropped per brainstorm; needs a cefepime-specific rule + CLSI-SDD-band abstain — out of scope here.

## Open Questions
- T1 surface: this plan emits single-cell top-level TMP-SMX artifacts (render directly). Separately patching `build_cell` to FLATTEN the existing `drugs`-map Oxford/Sci234 artifacts (so the cipro/cef/gent cells also render) is a deferred follow-up, not scoped here.
- Promotion: graduating TMP-SMX from EXPERIMENTAL to the frozen deployed surface (DRUG_RULE/DRUG_BREAKPOINTS/shipped_decoder_surface + freeze re-stamp) is a separate deliberate decision, out of scope.

## Verification
- `uv run pytest tests/ -q` — no new failures; new rule/scorer/builder tests green.
- `git diff --stat dna_decode/eval/amr_rules.py dna_decode/data/mic_tiers.py dna_decode/eval/cohort_manifest.py scripts/build_validation_report_card.py scripts/compute_lineage_metrics.py` — empty (frozen invariant).
- Frozen-leak test asserts TMP-SMX absent from `supported_drugs()` + `DRUG_RULE`.
- Both cohort artifacts emitted with `_schema=external-validation-v1` + branding; external report card renders 2 EXPERIMENTAL cells distinct from frozen cells.
- Strata reproduce across cohorts (sul+dfr high-R, others low) — else INDETERMINATE.

## Save-time amendments

Captured at: 2026-06-16
Source: `/save-plan` arguments

Audit-notes only — `/execute-plan` reads executable work from `## Implementation Steps`, not this block.

- TMP-SMX-only (cefepime deferred)
- scorer-local sul-AND-dfr rule in the non-frozen external-validation arm
- single-cell top-level external_validation_* artifact branded EXPERIMENTAL_SCORED
- reuse frozen classify_tier with non-frozen cotrimoxazole breakpoints
- strata-reproduction gate
- frozen files untouched
<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
