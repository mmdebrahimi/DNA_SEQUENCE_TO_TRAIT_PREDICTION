# F-B · Catalog Curation
<!-- project-schema: 0.1 -->

> Initialized 2026-08-31. Project ID: catalog-curation-2026-08-31. Originating goal (verbatim): "Close the measured catalog-completeness gaps, and define the procedure by which a curated biological fact is allowed to enter a shipped catalog."

## Project Context
- **Project ID:** catalog-curation-2026-08-31
- **Project root:** C:/Users/Farshad/PythonProjects/dna_decode
- **Captured:** 2026-08-31
- **Originating goal:** Close the measured catalog-completeness gaps, and define the procedure by which a curated biological fact is allowed to enter a shipped catalog
- **Refined goal:** For each measured gap, land the curation ONLY when (a) every entry is sourced to a named external authority, (b) the recovery is measured against the F-A doubt-layer baseline, and (c) the lock/freeze consequence is stated explicitly and authorized.
- **Horizon (months):** 6
- **Schema:** project-schema 0.1
- **Family of:** dna-decode-2026-05-11 · **blocked_by:** doubt-layer-2026-08-31

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** skipped-by-flag
- **Provisional:** NO
- **Findings:** The goal's factual premises are all measured in-repo and each cites a committed artifact (E1-E3). The specific curation targets — six HIV NNRTI driver substitutions and the `rmt*` methyltransferase family — are named, counted, and attributable to external authorities (Stanford HIVDB; AMRFinder Subclass assignment across 158 cached rows). No premise here depends on an unresolved external fact.

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Bounded: two named gaps, a countable set of entries each, and a measurable recovery. Not a program — it does not propose curating the catalog in general.

## Refinement Candidates
- **Verdict:** FAIL
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:** (FAIL because the originating goal conflates two moves with DIFFERENT authority consequences, and had to be split.)
  - **C1** — HIV NNRTI curation. `hiv_amr.py` is **NOT** in the prospective-lock manifest, so curating it does **not** invalidate the lock or the freeze. Drivers named and counted: V179D x12, A98G x10, H221Y x7, F227C x5, V108I x4, V179E x3. *Falsifier:* an entry lands with no external source citation.
  - **C2** — Gentamicin `rmt` curation. This DOES invalidate the prospective lock AND the reproducibility freeze; it requires an unfrozen revision, its own validation, and a NEW lock. Strictly separate from C1.
  - **C3** — The curation PROCEDURE itself: a written, tested rule for what may enter a shipped catalog. This is the durable deliverable; C1/C2 are its first two instances.
  - **C4 (rejected)** — Curating the 40 unrecorded colour loci. REJECTED for now: fabrication hazard unless every locus is OMIA/literature-sourced, and the colour family is frozen at 19 cells pending a user scope call.

## Goal Hierarchy
### Long-term (12+ months tier)
No curated biological fact enters a shipped catalog without a named source, a measured recovery, and a stated lock consequence.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Curation procedure written and testable | a test refuses a catalog entry lacking a source citation | 1 month |
| 2 | HIV NNRTI gap closed (C1) | the 53 missed isolates measured before/after against the F-A baseline | 2 months |
| 3 | Gentamicin `rmt` revision + NEW lock (C2) | v2 lock manifest exists; old lock explicitly retired, not silently broken | authority-gated |
| 4 | Over-call risk bounded, not assumed | S-labelled carriers sought in a NEW source; absence stops being arithmetic | 4 months |

### Short-term (≤1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Write the curation procedure + its refusal test | edit-local-code | Soraya | days |
| 2 | Draft the HIV NNRTI entries with per-entry Stanford citations | propose | Soraya | days |
| 3 | Surface the C1-vs-C2 authority split to the user | ask-user | Soraya | now |

## State Snapshot
### Assumptions
- Curation is the winning framing for the blind spot — **high** (framing sweep: incumbent framing 0 survivors / 2 executed kills; curation framing 3 survivors).
- The HIV drivers are real DRMs, not label noise — **medium** (named in the measured set; not independently re-derived here).
- Specificity cost of the `rmt` fix is unknown, not zero — **high** (see U1).

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| E1 | `hiv_amr.py` is NOT pinned by the prospective lock | wiki/innovate_blindspot_framing_sweep_2026-08-31.md | high | 2026-08-31 |
| E2 | `rmt` fix recovers +0.369 sens on 131 disjoint isolates | wiki/gentamicin_rmt_disjoint_validation_2026-08-28.md | high | 2026-08-28 |
| E3 | armA is ALREADY counted (24/24 GENTAMICIN); the gap is `rmt*` only | CLAUDE.md accrual correction 2026-08-29 | high | 2026-08-29 |
| E4 | Zero S-labelled `rmt` carriers exist across three datasets | wiki/gentamicin_rmt_label_hunt (63 public carriers, PPV 62/63) | high | 2026-08-28 |
<!-- project-state:end:evidence -->

### Unknowns
- **U1** — Over-calling risk of the `rmt` rule. "Specificity unchanged" is ARITHMETIC over an absence, not evidence: no S-labelled `rmt` carrier exists in any dataset checked, so the denominator is empty.
- Whether curating HIV NNRTI positions recovers the 53 or merely relabels them.
- Whether a curation procedure can be enforced by test at all, or is irreducibly a review discipline.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | The blind spot is a curation gap, not a computation gap | under-investigation | 2026-08-31 |
| H2 | Curation recovery exceeds the F-A doubt-layer baseline | open | blocked on F-A |
| H3 | An empty S denominator can be filled from a new source | open | (untested) |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| C1 and C2 are separate decisions with different authority consequences | 2026-08-31 | conflating them would smuggle a lock-invalidating change past a lock-safe one |
| Curation must be measured against F-A, not asserted | 2026-08-31 | this is why F-B is blocked_by F-A |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Edit a shipped catalog at all (C1)? | Soraya | **user authority** | lock-safe, but a scope call about what ships |
| Gentamicin v2 lock (C2)? | Soraya | **user authority** | invalidates the prospective lock AND the reproducibility freeze |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Blocked on F-A by construction; both curation targets named and counted; both terminal moves are user-authority calls.

### Target state / terminal condition
A written, tested curation procedure exists, and each measured gap is either closed under it or explicitly declined.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` + `gates-passed`. At init: 0 retired, 0 of 2 MVP criteria met.

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | Write curation procedure + refusal test | edit-local-code | med | med | low | 1-2hr |
| 2 | Draft HIV NNRTI entries with citations | propose | med | med | med | 1-2hr |
| 3 | Seek S-labelled `rmt` carriers in a new source (retires U1) | research | high | high | high | days |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- **Default:** after any action class fires.
- **Family-specific:** re-evaluate the moment F-A produces a measured baseline — that is what unblocks this family.

## MVP Criteria
- `file-exists wiki/catalog_curation_procedure.md`
- `test-exit-0 uv run pytest tests/test_catalog_curation_procedure.py -q`

Attempt budget: 3.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` / `research` / `write-plan` / `run-tests` / `ask-user` / `stop` — auto; `edit-local-code` — per-action approval.

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-08-31 | propose | ledger created (project-init protocol applied by hand) | blocked_by doubt-layer-2026-08-31 |
<!-- project-state:end:action-log -->

## Open Questions for User
- **C1** — may a curated biological fact enter the shipped `hiv_amr.py`? Lock-safe, but it changes what the tool asserts.
- **C2** — is the gentamicin v2 lock authorized? It invalidates the prospective lock and the reproducibility freeze; a fix is an unfrozen revision needing its own validation and a NEW lock.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-08-31
- **Progress signal:** (none yet — blocked on F-A by construction)
