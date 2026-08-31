# F-D · Learned, Narrow (L4 restraint)
<!-- project-schema: 0.1 -->

> Initialized 2026-08-31. Project ID: learned-narrow-2026-08-31. Originating goal (verbatim): "Hold the learned layer to its measured regime — molecular endpoints and constructed variation — and keep the boundary enforced rather than merely written down."

## Project Context
- **Project ID:** learned-narrow-2026-08-31
- **Project root:** C:/Users/Farshad/PythonProjects/dna_decode
- **Captured:** 2026-08-31
- **Originating goal:** Hold the learned layer to its measured regime — molecular endpoints and constructed variation — and keep the boundary enforced rather than merely written down
- **Refined goal:** Encode the measured regime boundary as a checkable artifact, so a proposal to extend L4 outside it is REFUSED by a screen rather than by whoever happens to remember the negative.
- **Horizon (months):** 6
- **Schema:** project-schema 0.1
- **Family of:** dna-decode-2026-05-11 · **blocked_by:** (none — independent)
- **Family kind:** RESTRAINT. Its deliverable is a boundary that stays enforced, not a build.

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** skipped-by-flag
- **Provisional:** NO
- **Findings:** The goal's premises are measured and each cites a committed artifact (E1-E4). The regime split is not asserted from memory — it was CORRECTED from memory once (see D1), and the corrected version is what this ledger encodes.

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Bounded and unusual in direction: the success criterion is that a specific class of work does NOT get built without clearing a screen. Deliverable is a screen plus its test, not a capability.

## Refinement Candidates
- **Verdict:** FAIL
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:** (FAIL because "hold the boundary" is not falsifiable as stated and had to be rewritten as an artifact.)
  - **C1 (selected)** — A regime-classification screen exists: given a proposed learned-decoder cell, it returns the regime and whether that regime has a measured positive. *Falsifier:* the screen passes a natural-population zero-shot proposal.
  - **C2 (selected)** — The screen is wired into the same place a new cell registers, so it cannot be skipped by not remembering it.
  - **C3** — The regime map is derived from the artifacts at read time, not written in prose. Prose went stale three times.
  - **C4 (rejected)** — "Never build learned decoders." REJECTED and factually wrong: three regimes have measured positives (segregant cross 12/12 r 0.46-0.80; TEM-1 genome-edit 0.761; FBA per-condition MCC 0.70-0.74).

## Goal Hierarchy
### Long-term (12+ months tier)
No learned-decoder proposal reaches build without its regime being named and its regime's evidence stated.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Regime map derived from artifacts, not prose | a script prints the map; prose cites the script | 2 months |
| 2 | Regime screen returns a verdict per proposal | screen refuses a natural-population zero-shot proposal | 3 months |
| 3 | Screen wired where new cells register | a new cell cannot register without a regime field | 5 months |

### Short-term (≤1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Derive the regime map from committed artifacts | edit-local-code | Soraya | days |
| 2 | Pin the corrected regime statement by test | run-tests | Soraya | days |
| 3 | Draft the screen's refusal criteria | propose | Soraya | days |

## State Snapshot
### Assumptions
- The discriminating variable is POPULATION DESIGN, not organism complexity — **high** (corrected 2026-08-29; the earlier "organism complexity" reading was wrong).
- Scale is dead in this regime; modality is live — **high** (650M > 3B > 15B, reproduced on our own full-benchmark run).
- A boundary written in prose will be violated — **high** (violated three times by the same compression error).

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| E1 | Natural-population + zero-shot: 0-for-5 de-confounded | wiki/organism_gp_regime_correction_2026-08-29.md | high | 2026-08-29 |
| E2 | Constructed variation works: segregant cross 12/12, r 0.46-0.80 | same | high | 2026-08-29 |
| E3 | Orthogonal modality lifts; scale does not (ESM2+GEMME+ProSST 90.5% paired) | wiki/forward_modality_hybrid_2026-07-17.json | high | 2026-07-17 |
| E4 | Antagonistic endpoints INVERT: ESM 0.454 below chance vs catalog 0.926 | wiki/hiv_esm_vs_catalog_2026-07-09.md | high | 2026-07-09 |
<!-- project-state:end:evidence -->

### Unknowns
- Where the constructed/natural boundary actually sits for an intermediate design (a structured pedigree, a mapping population).
- Whether the condition-SWITCH cell (FBA) is a fourth regime or a coverage problem inside an existing one.
- Whether a screen can classify a proposal's regime from its description, or needs the dataset in hand.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | Population design, not organism complexity, is the discriminator | confirmed | 2026-08-29 |
| H2 | A proposal's regime is classifiable from its description alone | open | (untested) |
| H3 | The FBA condition-switch cell is a coverage problem, not a new regime | under-investigation | 2026-08-29 |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| **D1** — "organism-level g2p is a closed negative" is WRONG and must not be repeated | 2026-08-29 | a 12/12 positive at r 0.46-0.80 exists; the compression hid a live direction three separate times |
| Do NOT extend L4 to organism-level natural populations | 2026-08-31 | 0-for-5, de-confounded, independently confirmed at 24,000 genomes where more data does not rescue it |
| Do NOT pretrain; scale is measured dead here and the field is crowded | 2026-08-31 | gLM2-650M is a download, not a research programme |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Score off-the-shelf gLM2-650M against the curated baseline? | Soraya | none (executor work) | the cheapest decisive test of the gene-LLM idea; no training run |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Boundary correct and measured, but enforced only by memory — and memory has failed it three times.

### Target state / terminal condition
A regime screen exists, is derived from artifacts, and is wired where new cells register.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` + `gates-passed`. At init: 0 retired, 0 of 2 MVP criteria met.

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | Derive the regime map from artifacts | edit-local-code | med | med | low | 1-2hr |
| 2 | Pin the corrected regime statement by test | run-tests | high | low | low | 30min |
| 3 | Score gLM2-650M vs the curated baseline | research | med | high | high | days |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- **Default:** after any action class fires.
- **Family-specific:** whenever a learned-decoder proposal appears — that is precisely when the boundary is load-bearing.

## MVP Criteria
- `test-exit-0 uv run pytest tests/test_regime_boundary.py -q`
- `file-exists wiki/learned_regime_map.json`

Attempt budget: 3.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` / `research` / `write-plan` / `run-tests` / `ask-user` / `stop` — auto; `edit-local-code` — per-action approval.

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-08-31 | propose | ledger created (project-init protocol applied by hand) | restraint family; eligible |
<!-- project-state:end:action-log -->

## Open Questions for User
- None requiring authority. This family's whole content is a boundary the user already set; the work is making it checkable.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-08-31
- **Progress signal:** (none yet — init only)
