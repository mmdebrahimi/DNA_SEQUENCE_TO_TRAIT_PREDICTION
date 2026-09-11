# ECOFF-Anchored Tiering — Evaluation Only — Technical Plan

> Measure whether an ECOFF (wild-type) anchor agrees with our genotype calls better than the clinical breakpoint, without touching the frozen surface.

## Lens status
Inputs: conversation (candidate row 5 of `project_state/evidence-surface-2026-08-31.md`)
Degradations: sentrux unavailable (not installed) — Architecture Impact gates n/a; repo-index unconfigured — Step 1 research done by direct file reading; EUCAST ECOFF values `[unverified]` until Step 1 sources them over the network.

## Problem Statement
`dna_decode/data/mic_tiers.py` anchors every tier decision to **clinical breakpoints** (CLSI 2024 + EUCAST 14.0) and versions nothing per entry. The EUCAST WGS subcommittee's position, surfaced by the 2026-09-03 prior-art scan, is that a genotype-based decoder should be anchored to the **ECOFF** instead, because the two answer different questions:

- a **clinical breakpoint** answers *will treatment succeed* — a function of dose, site and PK/PD as well as the organism;
- an **ECOFF** answers *does this isolate carry an acquired mechanism* — wild-type vs non-wild-type.

A determinant-based decoder predicts the second and is scored against the first. The same scan recorded that ResFinder 4.0's residual errors cluster **one dilution above the ECOFF**, which is the signature this mismatch would produce.

Scope is **evaluation only**. Adopting ECOFF anchoring would edit `mic_tiers.py`, which is sha256-pinned in `wiki/prospective_lock_manifest_2026-08-31.json`, and would retire the v2 lock — a user-authority call of the same class as the gentamicin v2 deployment. This plan stops at the measurement that informs that call.

Non-goals: changing any shipped rule; re-tiering the 10 frozen SCORED cells (they carry categorical R/S from NCBI-PD, not numeric MICs, so they cannot support an ECOFF re-tiering at all); deriving an ECOFF ourselves.

## Codebase Context
- `dna_decode/data/mic_tiers.py` — FROZEN. `classify_tier(mics, distinct_calls, breakpoints)` takes the breakpoint dict as a **parameter**, so an alternative anchor can be evaluated without editing the module. `DRUG_BREAKPOINTS` covers cipro/cef/tet/gent/meropenem/oxacillin; its own comment already notes tetracycline "uses ECOFF only" and carries `eucast_r`/`eucast_s` as `None` there.
- `classify_tier` returns R/S tiers with 4x safety margins and hardcodes the keys `clsi_r`, `clsi_s`, `eucast_r`, `eucast_s`. ECOFF semantics are **binary WT/NWT with no intermediate zone**, so passing an ECOFF in as `clsi_s` would be silently wrong — the same failure class as the log2-unit trap below. A separate classifier is required.
- `data/raw/oxford/` — the only substrate on disk with numeric per-isolate MICs plus a genotype call: `main_data.csv` (2,897 isolates, per-drug log2 dilution intervals) and `amrfinder.tsv` (4,979 scanned isolates). MIC is `2**upper`, per `scripts/oxford_score.py:14`.
- `dna_decode/eval/error_rates.py::vme_me` — VME (R called S) / ME, shipped 2026-09-10; the clinical-microbiology convention and the right headline metric here.
- `tests/test_mic_tiers.py` — 46 tests pinning current tier semantics; they must stay green and untouched.

### Reusable-Code Survey
- `dna_decode/data/mic_tiers.py` — `classify_tier` takes breakpoints as a parameter, so the clinical arm needs no new code.
- `dna_decode/data/experimental_drug_rules.py` — precedent for a NON-frozen scorer-local overlay; TMP-SMX reused the frozen `classify_tier` with a non-frozen breakpoint dict passed in and touched no frozen file.
- `scripts/gentamicin_rmt_candidate.py` — precedent for evaluating a candidate rule offline, scorer-local, frozen surface untouched and asserted by test.
- `scripts/oxford_score.py` — `load_mic_labels` already parses the log2 encoding and joins MIC to AMRFinder by `guuid`.
- `dna_decode/eval/error_rates.py` — VME/ME computation, already tested.

## Pre-Change Baseline
Measured 2026-09-11 from the committed deposit, before any code is written:

- Oxford's MIC grid is **collapsed, not a full dilution series** — gentamicin takes only {1, 2, 4, 32} (no 8, no 16); ciprofloxacin {0.125, 0.25, 0.5, 8}; ceftriaxone {0.5, 1, 2, 4, 32}.
- Clinical-breakpoint class counts (CLSI): gentamicin 192 R / 2,681 S; ciprofloxacin 370 R / 2,471 S; ceftriaxone 238 R / 2,630 S / 6 I.
- Size of the stratum a lower anchor could move, by candidate cutoff: gentamicin 25 isolates at cutoff 2 and 431 at cutoff 1; ciprofloxacin 89 at 0.25 and 1,836 at 0.125; ceftriaxone 448 at 0.5.
- Current shipped tiering is clinical-breakpoint-only; `mic_tiers.py` sha256 is `47d6c70103673736580c510cf1484b09917665981c6f8594c3b8521e850bc283` per the v2 lock manifest.

## Verification Signal
- A committed `wiki/ecoff_tiering_result_2026-09-11.json` whose verdict is re-derivable by a test from its own numbers, against a bar frozen in a separate file before any result existed.
- Per drug, a feasibility verdict from Step 2 (`SCOREABLE`, `DEGENERATE_ECOFF_BELOW_PANEL`, or `UNDERPOWERED_STRATUM`) and, for scoreable drugs, the determinant-carriage rate among clinically-S isolates the ECOFF calls non-wild-type versus those it calls wild-type, with VME/ME under both anchors.
- `mic_tiers.py` sha256 unchanged and `verify_lock` still passing — asserted by test, not by inspection.
- Full suite green.

## Implementation Steps

### Step 1: Source the EUCAST ECOFFs and pin them with provenance
Files: dna_decode/data/ecoff_catalog.py, wiki/ecoff_sources_2026-09-11.json
Depends on: none

**What changes:**
- New NON-frozen module holding per-drug E. coli ECOFFs, each entry carrying `value`, `organism`, `source_url`, `eucast_version`, `retrieved` and `verbatim_quote`.
- Values are **fetched from the EUCAST MIC-distribution website, never recalled**. Writing a biological reference value from memory is the standing fabrication hazard in this repo, and a wrong ECOFF would produce a fully self-consistent, entirely wrong evaluation with nothing downstream to catch it.
- Any drug whose ECOFF cannot be sourced is recorded with `value: null` and `status: "unsourced"` and is excluded downstream; `ecoff_for(drug)` raises rather than defaulting.
- No import from, and no edit to, `mic_tiers.py`.

**Test strategy:**
- Every entry with a value must carry a non-empty `source_url` and `verbatim_quote`; a test fails if either is missing.
- `ecoff_for` raises on an unsourced or unknown drug rather than returning a default.

### Step 2: Substrate feasibility gate, with the degenerate outcomes named in advance
Files: scripts/ecoff_tiering_feasibility.py, wiki/ecoff_tiering_feasibility_2026-09-11.json
Depends on: Step 1

**What changes:**
- For each drug, compare the sourced ECOFF against the cohort's **measured** MIC grid and emit one of three pre-named verdicts: `SCOREABLE`; `DEGENERATE_ECOFF_BELOW_PANEL` (the ECOFF is at or below the lowest dilution the panel reports, so every isolate is non-wild-type and the anchor cannot discriminate); `UNDERPOWERED_STRATUM` (fewer than 20 isolates separate the two anchors, matching the repo's existing per-class floor).
- Reports the discriminating stratum size per drug, and refuses to emit a verdict for a drug whose ECOFF is unsourced.
- Runs **before** the acceptance bar is frozen, deliberately: two bars this session were mis-specified because a threshold was frozen before the mechanism analysis that says what a correct result can look like.

**Test strategy:**
- Synthetic grids exercise all three verdicts.
- Non-vacuity: a grid whose ECOFF sits below the lowest reported dilution must return `DEGENERATE_ECOFF_BELOW_PANEL`, never `SCOREABLE` — the live risk for ciprofloxacin, whose panel bottoms out at 0.125.

### Step 3: Freeze the acceptance bar, naming every outcome
Files: wiki/ecoff_tiering_acceptance_bar.json
Depends on: Step 2

**What changes:**
- Records the substrate as measured in Step 2 (per-drug stratum sizes), then freezes the verdict rule with all outcomes named: `SUPPORTED` (determinant carriage is higher among clinically-S/ECOFF-NWT isolates than among clinically-S/ECOFF-WT, and the gap exceeds its permutation null); `WEAK_DIRECTIONAL` (right direction, inside the null); `FALSIFIED` (carriage equal or lower); `INDETERMINATE` (no drug reaches `SCOREABLE`, or strata fall below the floor).
- States in advance that a `SUPPORTED` verdict argues the ECOFF is the better anchor **for this decoder on this cohort** and is not by itself grounds to edit the frozen surface.
- Written before any carriage number is computed, and not edited afterwards.

**Test strategy:**
- A test asserts the bar names all four outcomes and carries a `FROZEN BEFORE` status string, so a bar naming only the outcome that happened cannot pass as a pre-registration.

### Step 4: Wild-type classifier, separate from the clinical tiering rather than a reuse of it
Files: dna_decode/eval/wildtype_tiering.py
Depends on: Step 1

**What changes:**
- `classify_wildtype(mics, ecoff)` returning `WT`, `NWT` or `NO_MIC`, in a NON-frozen module, using the same median-of-MICs convention as `classify_tier` so the two arms differ only in the anchor.
- Deliberately **not** implemented by passing an ECOFF into `classify_tier` as `clsi_s`: that function returns R/S tiers with 4x safety margins and an intermediate zone, none of which exist in ECOFF semantics, and the abuse would typecheck while silently misleading.
- Carries an explicit note that WT/NWT is not S/R and must never be rendered as a clinical call.

**Test strategy:**
- Boundary cases at exactly the ECOFF (at-ECOFF is wild-type; one dilution above is not).
- A test pins that this module does not import `mic_tiers`, and that `mic_tiers` is unchanged.

### Step 5: Run the evaluation on the Oxford cohort
Files: scripts/ecoff_tiering_evaluate.py, wiki/ecoff_tiering_result_2026-09-11.json
Depends on: Step 3, Step 4

**What changes:**
- For each `SCOREABLE` drug, join MIC to AMRFinder by `guuid` reusing the log2 convention from `oxford_score.load_mic_labels` (MIC is `2**upper`), classify each isolate under both anchors, and call determinant presence with the deployed rule.
- Headline comparison on the discriminating stratum: determinant-carriage rate among clinically-S/ECOFF-NWT versus clinically-S/ECOFF-WT isolates, against a permutation null on the anchor labels holding stratum size fixed.
- Reports **VME/ME under both anchors** via `error_rates.vme_me`, since the dangerous direction is the one that matters clinically.
- Verdict computed mechanically by a `verdict_from_bar` function applying Step 3's frozen rule, never authored in prose.
- Refuses with a non-zero exit and writes no artifact if no drug is `SCOREABLE`, rather than emitting an empty artifact that would read as a null result.

**Test strategy:**
- The committed verdict is re-derivable from the artifact's own numbers.
- A run where the two anchors agree on every isolate must return `INDETERMINATE`, not `FALSIFIED` — an empty stratum is not evidence against the claim.

### Step 6: Tests, non-vacuity proofs, and the frozen-surface assertion
Files: tests/test_ecoff_catalog.py, tests/test_wildtype_tiering.py, tests/test_ecoff_tiering_evaluate.py
Depends on: Step 4, Step 5

**What changes:**
- Guards for each decision point above, each proven non-vacuous by re-injecting the defect it guards and confirming the guard fails.
- An explicit assertion that the five frozen-surface files are byte-unchanged and that `prospective_lock.verify_lock` still passes, mirroring `tests/test_tb_leak_guard.py`.
- A guard that the unit convention is honoured: feeding raw log2 indices instead of `2**upper` must fail loudly, since doing so silently returned zero resistant isolates on all three drugs during Step 1 research.

**Test strategy:**
- Full suite green, with the 46 existing `test_mic_tiers.py` tests untouched and passing.

### Step 7: Memo, artifact, and ledger row
Files: wiki/ecoff_tiering_result_2026-09-11.md, project_state/evidence-surface-2026-08-31.md
Depends on: Step 5, Step 6

**What changes:**
- Memo leading with the per-drug feasibility verdict, then the carriage comparison and VME/ME under both anchors, then the limits.
- States explicitly that the cohort is **single-source** — Oxford tripped G7 and is not-applicable on G2 for exactly that reason (`wiki/oxford_gate_screen_2026-09-11.md`) — so the result is scoped to one UK region, 2008-2018, on a collapsed MIC panel.
- Ledger row appended with a `[done:...]` marker retiring candidate row 5.
- Ends with the adoption question stated as a user-authority decision, not as a recommendation to execute.

**Test strategy:**
- A test asserts the memo's headline numbers match the committed JSON, so the prose cannot drift from the artifact.

## Execution Preview

Wave 0 (1 parallel):  Step 1 — Source the EUCAST ECOFFs
Wave 1 (2 parallel):  Step 2 — Substrate feasibility gate, Step 4 — Wild-type classifier
Wave 2 (1 parallel):  Step 3 — Freeze the acceptance bar
Wave 3 (1 parallel):  Step 5 — Run the evaluation
Wave 4 (1 parallel):  Step 6 — Tests and the frozen-surface assertion
Wave 5 (1 parallel):  Step 7 — Memo, artifact, ledger row

Critical path: Step 1 → Step 2 → Step 3 → Step 5 → Step 6 → Step 7 (6 waves)
Max parallelism: 2 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags
- Severity: critical — Fabrication hazard on the ECOFF values. A recalled ECOFF would produce a fully self-consistent, entirely wrong evaluation with nothing downstream to catch it. Step 1 requires `source_url` plus `verbatim_quote` per entry and excludes anything unsourced; the values stay `[unverified]` until that step runs against EUCAST.
- Severity: critical — The frozen surface must not change. `mic_tiers.py` is sha256-pinned in the v2 prospective lock; editing it retires the lock and spends the prospective evidence. Every arm here passes breakpoints in rather than editing, and Step 6 asserts the five files byte-unchanged.
- Severity: high — The substrate may be degenerate for some drugs. Oxford's ciprofloxacin panel bottoms out at 0.125; if the sourced ECOFF is at or below that, every isolate is non-wild-type and the anchor cannot discriminate. Step 2 names that outcome in advance rather than discovering it mid-analysis.
- Severity: high — Semantic mismatch between the two anchors. ECOFF is binary WT/NWT; `classify_tier` returns R/S tiers with safety margins. Passing an ECOFF in as `clsi_s` would typecheck and silently mislead, which is why Step 4 is a separate classifier and not a reuse.
- Severity: medium — Unit convention. Oxford's MIC columns are log2 dilution indices, not mg/L; applying breakpoints directly returns zero resistant isolates on all three drugs, which reads as a clean all-susceptible cohort rather than an error. MIC is `2**upper` per `scripts/oxford_score.py:14`, guarded in Step 6.
- Severity: medium — Single-source cohort. Oxford is one study, one UK region, 2008-2018, and holds zero `rmt` carriers in 4,979 isolates. Any result is scoped accordingly and must not be stated as a population-general property.
- Severity: medium — Network dependency. Step 1 needs the EUCAST site. If it is unreachable the plan stalls at Wave 0 rather than proceeding on recalled values.
- Severity: low — File overlaps. None within a wave; Step 2 and Step 4 share only their dependency on Step 1 and touch disjoint files.
- Severity: low — External tool integrations. No CLI or service is invoked beyond the EUCAST fetch in Step 1. `classify_tier`'s signature was verified by reading `dna_decode/data/mic_tiers.py:92` directly rather than inherited from a reference.
- Severity: low — Restructuring applied. Step 2 (feasibility) was ordered before Step 3 (freeze the bar) rather than the reverse, so the threshold is set with the mechanism analysis in hand.

## Open Questions
- Adoption is out of scope and is yours to decide. Even a `SUPPORTED` verdict only argues the ECOFF is the better anchor for this decoder on this cohort; acting on it edits `mic_tiers.py`, retires the v2 lock and restarts the prospective clock, exactly as the gentamicin v2 deployment did.
- If Oxford proves degenerate for every drug, should the evaluation move to a BV-BRC measured-MIC cohort (finer dilution grid, but a different and less controlled provenance), or stop and report the substrate wall?
- Which EUCAST release should the ECOFFs be pinned to — the current one, or 14.0 to match the clinical breakpoints already in `mic_tiers.py`? Mixing releases across the two arms would confound the comparison.

## Verification
1. `uv run python scripts/ecoff_tiering_feasibility.py` prints a per-drug verdict and writes its artifact; a drug with an unsourced ECOFF is refused rather than defaulted.
2. `uv run python scripts/ecoff_tiering_evaluate.py` exits non-zero and writes nothing when no drug is `SCOREABLE`.
3. `uv run pytest tests/test_ecoff_catalog.py tests/test_wildtype_tiering.py tests/test_ecoff_tiering_evaluate.py tests/test_mic_tiers.py -q` is green.
4. `git diff --stat HEAD -- dna_decode/data/mic_tiers.py dna_decode/data/amr_rules.py dna_decode/data/calibrated_amr_rules.json dna_decode/data/shipped_decoder_surface.py dna_decode/eval/cohort_manifest.py` is empty.
5. `uv run python -m scripts.prospective_lock_validate` still verifies the v2 lock.
6. Full suite `uv run pytest tests/ -q` is green.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high -->
