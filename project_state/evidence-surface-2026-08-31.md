# F-C · Evidence Surface (L3)
<!-- project-schema: 0.1 -->

> Initialized 2026-08-31. Project ID: evidence-surface-2026-08-31. Originating goal (verbatim): "Expose the evaluation machinery — de-confounding, nulls, denominators, leakage gates, provenance and source concentration — as a first-class product surface rather than internal scaffolding."

## Project Context
- **Project ID:** evidence-surface-2026-08-31
- **Project root:** C:/Users/Farshad/PythonProjects/dna_decode
- **Captured:** 2026-08-31
- **Originating goal:** Expose the evaluation machinery — de-confounding, nulls, denominators, leakage gates, provenance and source concentration — as a first-class product surface rather than internal scaffolding
- **Refined goal:** Make the L3 evidence machinery reachable and legible from the shipped CLI, so a user can ask "how was this validated, and what could it not have seen?" and get a machine-readable answer without reading `wiki/`.
- **Horizon (months):** 6
- **Schema:** project-schema 0.1
- **Family of:** dna-decode-2026-05-11 · **blocked_by:** (none — independent)

## Empirical Concerns
- **Verdict:** N-A
- **Check status:** not-applicable
- **Provisional:** NO
- **Findings:** (Imperative goal; no factual-shape claim about the external world. Its premise — that this machinery exists and is under-exposed — is verified against committed modules: `cohort_deconfound.py`, `clonality.py`, `cohort_manifest.py`, `prospective_lock.py`, `trust_surface.py`, plus the three namespace-separate disclosure layers on the report card.)

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Bounded — exposes machinery that already exists. Builds no new evaluation method; the deliverable is reach and legibility.

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:**
  - **C1 (selected)** — Every disclosure layer on the report card is reachable from the CLI, not only from `wiki/`. *Falsifier:* a layer exists on the card with no CLI route.
  - **C2 (selected)** — A cell's answer to "what could this validation NOT have seen?" is machine-readable (source concentration, lineage collapse, prospective status). *Falsifier:* the answer requires prose interpretation.
  - **C3** — The four gate families (G1-G10 in the negative-results map) are runnable as a screen against a candidate dataset, not just readable as a memo.
  - **C4 (rejected)** — A single aggregate "trustworthiness score". REJECTED: the report card deliberately has NO aggregate headline, and adding one would destroy the per-cell honesty the surface exists for.

## Goal Hierarchy
### Long-term (12+ months tier)
The evaluation discipline is the product's differentiator, and a user can see it without reading the repository.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Every card disclosure layer has a CLI route | a test enumerates layers and asserts each is reachable | 2 months |
| 2 | "What could this not have seen?" is machine-readable per cell | a structured limits block per cell | 3 months |
| 3 | The rejection gates run as a screen, not just a memo | `screen_candidate_dataset(...)` returns per-gate verdicts | 4 months |

### Short-term (≤1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Enumerate every disclosure layer and its current reachability | research | Soraya | days |
| 2 | Add the missing CLI routes, augment-only | edit-local-code | Soraya | weeks |
| 3 | Guard test: no layer exists without a route | run-tests | Soraya | days |

## State Snapshot
### Assumptions
- The evaluation discipline is the real differentiator — **medium-high** (three critique papers describe its absence; that is inference from omission, not measurement).
- Exposure changes nothing about correctness, only about reach — **high**.
- Augment-only is the right discipline here too — **high** (three prior layers all held it).

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| E1 | Three disclosure layers exist, all namespace-separate and augment-only | wiki/decoder_validation_report_card.json | high | 2026-08-29 |
| E2 | Source concentration measurably explains a real blind spot | wiki/provdisjoint_source_concentration_2026-08-28.md | high | 2026-08-28 |
| E3 | Ten rejection gates are written but are prose, not code | wiki/negative_results_map_2026-06-13.md | high | 2026-08-26 |
| E4 | `trust_block()` is the existing always-safe record accessor to mirror | dna_decode/data/trust_surface.py:293 | high | 2026-08-31 |
<!-- project-state:end:evidence -->

### Unknowns
- Whether exposing limits reduces user trust more than it earns — untested, and plausibly the reason nobody ships this.
- Whether the gate screen can be coded at all, or is irreducibly a judgment call per dataset.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | Every disclosure layer can be routed without changing any cell's tier | open | (untested) |
| H2 | The rejection gates are codeable as a screen | open | (untested) |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| No aggregate trustworthiness headline, ever | 2026-08-31 | the card's per-cell honesty is the point; an aggregate destroys it |
| Augment-only: exposure never changes a tier | 2026-08-31 | same discipline as lineage / prospective / source-concentration layers |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Does a single-source cell warrant demotion, or only disclosure? | Soraya | **user authority** | current answer: disclose; demoting is a scope call |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Machinery built and load-bearing; reach is the gap — most of it is visible only by reading `wiki/`.

### Target state / terminal condition
Every disclosure layer is CLI-reachable and every cell answers "what could this not have seen?" in machine-readable form.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` + `gates-passed`. 2026-08-31: 2 / 2 MVP criteria met (reachability test suite + layer inventory artifact). C1 + C2 met at record AND human level. **2026-09-01: C3 MET** — the ten gates run
(`scripts/screen_candidate_gates.py`), refuse on the 2 judgment gates, and reproduce both committed hand
verdicts; the reproduction check caught a G2 applicability-ordering bug and an overstated PEAR headline.
All three refinement candidates are now closed; residual work is n=2 schema breadth, not the mechanism.

### Candidate next actions
Actions 1-3 completed 2026-08-31 (Action Log rows 2-4). C1 + C2 are met at the record level AND the
human-readable level; the follow-on set is below.

| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | DONE 2026-09-01 — C3 shipped: `eval/rejection_gates.py` + `scripts/screen_candidate_gates.py --verify` | edit-local-code | high | high | resolved | — |
| 1b | Screen a THIRD candidate (n=2 worked examples is the schema's stated limit) | research | med | high | med | 1-2hr |
| 2 | Check reachability of the NON-report-card surfaces (hiv/tb/pgx cards) | research | med | med | med | 1-2hr |
| 3 | Decide whether the doubt/lineage lines belong in `dna-decode` routes too | propose | low | med | med | 1-2hr |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- **Default:** after any action class fires.
- **Family-specific:** whenever a NEW disclosure layer lands on the report card — that is exactly when reachability drifts.

## MVP Criteria
- `test-exit-0 uv run pytest tests/test_evidence_surface_reachable.py -q`
- `file-exists wiki/evidence_surface_layer_inventory.json`

Attempt budget: 3.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` / `research` / `write-plan` / `run-tests` / `ask-user` / `stop` — auto; `edit-local-code` — per-action approval.

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-08-31 | propose | ledger created (project-init protocol applied by hand) | eligible; independent of F-A |
| 2 | 2026-08-31 | research | enumerate the 4 card disclosure layers vs CLI reachability | 2 of 4 card-only; prospective surfaced only when it CONTRADICTED |
| 3 | 2026-08-31 | edit-local-code | attach lineage + source_concentration + full prospective to trust_block | all 4 now record-reachable; augment-only verified by diff |
| 4 | 2026-08-31 | edit-local-code | human-readable renderers (doubt / lineage / source concentration) | all 4 now print; a JSON-only disclosure is not one |
| 5 | 2026-09-01 | edit-local-code | C3: rejection gates as a runnable screen (eval/rejection_gates.py + scripts/screen_candidate_gates.py) | 8 mechanical / 2 judgment; G1+G3 REFUSE without a human reading. Reproduces both hand verdicts; caught a G2 applicability-ordering bug + an overstated PEAR headline |
<!-- project-state:end:action-log -->

## Open Questions for User
- Whether a single-source SCORED cell warrants more than disclosure (3 of 10 rest on one BioProject).

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-08-31
- **Progress signal:** (none yet — init only)
