# F-E · Label Acquisition
<!-- project-schema: 0.1 -->

> Initialized 2026-08-31. Project ID: label-acquisition-2026-08-31. Originating goal (verbatim): "Address the binding constraint — labels, not models — by identifying and, if authorized, acquiring a label source that clears the eight rejection gates by construction."

## Project Context
- **Project ID:** label-acquisition-2026-08-31
- **Project root:** C:/Users/Farshad/PythonProjects/dna_decode
- **Captured:** 2026-08-31
- **Originating goal:** Address the binding constraint — labels, not models — by identifying and, if authorized, acquiring a label source that clears the eight rejection gates by construction
- **Refined goal:** Maintain a screened, ranked list of candidate label sources scored against the ten rejection gates, so that IF the user authorizes an acquisition the target is already chosen and its gate profile already known. Executor scope is screening only; acquisition is authority.
- **Horizon (months):** 12
- **Schema:** project-schema 0.1
- **Family of:** dna-decode-2026-05-11 · **blocked_by:** (none — but its terminal move is external/authority)

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** skipped-by-flag
- **Provisional:** NO
- **Findings:** "Labels, not models" is the reproducibility freeze's own recorded conclusion after every public-label expansion closed (E1). The named candidate (PEAR) is described from a published source; its accession-level reachability is **UNVERIFIED** and is recorded as U1 rather than assumed — a paper's cited accessions have resolved empty or to the wrong organism before.

## Project vs Research-Program
- **Verdict:** FAIL
- **Provisional:** NO
- **Classification:** hybrid
- **Rationale:** The SCREENING half is a bounded project (a ranked list, gate-scored, with a verification step). The ACQUISITION half is unbounded from the executor's side — it depends on an external party, possibly money, and a user decision. Split accordingly: this ledger owns screening only.

## Refinement Candidates
- **Verdict:** FAIL
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:** (FAIL — the originating goal mixes an executor task with an authority action.)
  - **C1 (selected)** — A ranked candidate list exists, each entry scored against all ten rejection gates. *Falsifier:* an entry ships with a gate unscored.
  - **C2 (selected)** — Every candidate's accessions are RESOLVED before it is recommended (right organism, links to data, on both ENA and NCBI). *Falsifier:* a recommendation rests on a paper's citation alone.
  - **C3 (authority)** — Acquire a non-public wet-lab or clinical label source. NOT an executor task.
  - **C4 (rejected)** — Reopen public-label expansion. REJECTED: closed with recorded negatives; the provdisjoint grid is saturated.

## Goal Hierarchy
### Long-term (12+ months tier)
When a label source is worth acquiring, the target is already screened, verified, and ranked — the decision is the only thing left.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Candidate list scored against all ten gates | no entry has an unscored gate | 3 months |
| 2 | PEAR accessions resolved or the candidate demoted | resolution recorded per accession, both archives | 4 months |
| 3 | The free forward path (prospective accrual) kept alive | periodic sweep re-run; accrual status current | ongoing |

### Short-term (≤1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Score PEAR against all ten gates, explicitly | research | Soraya | days |
| 2 | Resolve 2-3 PEAR accessions on ENA + NCBI before recommending | research | Soraya | days |
| 3 | Re-run the prospective accrual sweep (free, no authority) | run-tests | Soraya | days |

## State Snapshot
### Assumptions
- Labels, not models, are the binding constraint — **high** (the freeze's own conclusion).
- PEAR is constructed variation at scale on an AMR target — **medium** (from its publication; unverified in-hand).
- The prospective-lock path remains free and un-gated — **high** (it fired 2026-08-24 and produced real scores).

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| E1 | Every public-label AMR expansion is closed; constraint is LABELS | wiki/reproducibility_freeze_2026-06-13.md | high | 2026-06-13 |
| E2 | Ten reusable rejection gates exist and screen candidates | wiki/negative_results_map_2026-06-13.md | high | 2026-08-26 |
| E3 | The prospective lock produced real scores on 63 post-lock isolates | wiki/prospective_lock_first_accrual_2026-08-24.md | high | 2026-08-24 |
| E4 | Constructed variation is the regime with measured positives | wiki/organism_gp_regime_correction_2026-08-29.md | high | 2026-08-29 |
<!-- project-state:end:evidence -->

### Unknowns
- **U1** — Whether PEAR's accessions actually resolve, to the right organism, with data attached. UNVERIFIED. A paper's cited accessions have been empty or wrong-organism before; verify 2-3 on both ENA and NCBI BEFORE any multi-hour fetch.
- Whether PEAR's growth measurements survive gate G6 (MIC censoring) and G10 (rule scoreability against genotype).
- Whether any acquisition target is reachable without money.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | PEAR clears all ten gates | open | (unscored) |
| H2 | A free path exists that clears the label wall (prospective accrual) | confirmed | 2026-08-24 |
| H3 | Some acquisition target is reachable without money | open | (untested) |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| Screening is executor scope; acquisition is authority | 2026-08-31 | the split that makes this family workable at all |
| Verify accessions before recommending, never after | 2026-08-31 | cited accessions have resolved empty or wrong-organism before |
| Nothing downstream may be planned as if acquisition will land | 2026-08-31 | it is external and may never happen |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Authorize a PEAR acquisition? | Soraya | **user authority** + possible money | ~23,000 E. coli, single-copy blaCTX-M-14 variants, measured growth |
| Authorize any non-public wet-lab/clinical label source? | Soraya | **user authority** + money | clears the label gates by construction |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
One named candidate, unverified; the free forward path (prospective accrual) is alive and produced real scores; acquisition itself is authority-gated.

### Target state / terminal condition
A gate-scored, accession-verified, ranked candidate list exists — so an acquisition decision has nothing left to research.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` + `gates-passed`. At init: 0 retired, 0 of 2 MVP criteria met.

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | Score PEAR against all ten gates | research | med | high | med | 1-2hr |
| 2 | Resolve PEAR accessions on ENA + NCBI (retires U1) | research | high | high | med | hours |
| 3 | Re-run the prospective accrual sweep | run-tests | med | med | low | hours |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- **Default:** after any action class fires.
- **Family-specific:** re-run the accrual sweep periodically — the free path accrues on NCBI-PD ingestion lag, not on anything we control.

## MVP Criteria
- `file-exists wiki/label_source_candidates_gate_scored.json`
- `test-exit-0 uv run pytest tests/test_label_source_gate_screen.py -q`

Attempt budget: 3.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` / `research` / `write-plan` / `run-tests` / `ask-user` / `stop` — auto; `edit-local-code` — per-action approval.
- **Note:** any acquisition step is money/authority — hard pause, no flag disables it.

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-08-31 | propose | ledger created (project-init protocol applied by hand) | screening scope only; acquisition is authority |
<!-- project-state:end:action-log -->

## Open Questions for User
- Whether to authorize a PEAR acquisition, or any non-public label source. Both are external and may involve money.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-08-31
- **Progress signal:** (none yet — init only)
