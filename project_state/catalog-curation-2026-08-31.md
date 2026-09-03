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
- **U1 (SHARPENED 2026-09-02, still open)** — Over-calling risk of the deployed `rmt` rule. The denominator is no longer empty *by absence*: a cache-independent sweep of 20,816 gentamicin-labelled PD isolates found **60** S-labelled carriers — the first ever. But **all 60 come from a single BioProject (PRJNA1322038)**, and a pre-registered `aac(3)` control shows that project calls the undisputed gentamicin determinant R only **2%** of the time vs **97%** elsewhere (and calls no-gene isolates R 86%) → `LABEL_ARTIFACT`; they cannot test the rule. Outside it: **146/146 carriers are R** across 23 BioProjects (was 62/63, cache-bounded). So the rule's evidence is stronger AND its over-call risk is still UNTESTED. Next lever is an INDEPENDENT archive, not PD. `wiki/gentamicin_rmt_specificity_hunt_2026-09-02.md`
- Whether curating HIV NNRTI positions recovers the 53 or merely relabels them.
- RESOLVED 2026-09-01: PARTLY. 1 of 4 conditions is enforced by test (per-MODULE external authority; `tests/test_catalog_provenance.py`, all 8 shipped catalogs already pass). The other 3 -- per-ENTRY sourcing, measurement against the doubt-layer baseline, and reviewing a derivation as biology -- are review discipline. Per-entry is not representable without restructuring the bare `set[str]` catalogs -- **but AMRrules (2026-09-03 prior-art scan) is the existence proof that it IS representable at production scale: 24 fields per rule incl. PMID, ECO `evidence code`, `evidence grade`, `evidence limitations`. The constraint was OUR schema, not the problem.** Stating the split is load-bearing: a green suite checking only the weak condition reads as more assurance than it is.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | The blind spot is a curation gap, not a computation gap | falsified | 2026-09-01 |
| H2 | Curation recovery exceeds the F-A doubt-layer baseline | falsified | 2026-09-01 |
| H4 | A literature-anchored (not data-derived) curation could still be worth it | open | (untested) |
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
| RESOLVED 2026-09-01: Edit a shipped catalog at all (C1)? | Soraya | none -- measured | Answered by measurement, not authority: data-derived curation LOSES to the free doubt layer on every variant tested (recovery 0.000-0.500 vs 0.604) and the 3x variant would drop canonical Y181C. No edit made. `wiki/hiv_nnrti_curation_verdict_2026-09-01.md` |
| RESOLVED 2026-08-31: Gentamicin v2 lock (C2)? | Soraya | none -- authorized + shipped | User-authorized; deployed with symbol_rescue + a NEW lock manifest (prospective_lock_manifest_2026-08-31.json). E. coli N=131 sens 0.523 -> 0.892. `wiki/gentamicin_v2_lock_2026-08-31.md` |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Both named curation targets are RESOLVED: gentamicin `rmt` shipped under a v2 lock (2026-08-31); HIV NNRTI measured and DECLINED (2026-09-01, loses to the free doubt layer). Only C3, the curation procedure itself, remains open -- and is lower-value now that the first candidate was rejected on measurement.

### Target state / terminal condition
A written, tested curation procedure exists, and each measured gap is either closed under it or explicitly declined.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` + `gates-passed`. 2026-09-01: the family's central question is ANSWERED (curation does not earn its place; measured 3 ways). 2 hypotheses falsified. The MVP criteria as written (a curation PROCEDURE + its refusal test) are now lower-value than when drafted, since the first candidate curation was measured and rejected.

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | DONE 2026-09-01 -- wiki/catalog_curation_procedure.md + tests/test_catalog_provenance.py (18). Resolves U3: per-ENTRY citation is NOT representable (bare set[str] catalogs); per-MODULE authority IS, and all 8 already pass -> a regression guard, not a migration | edit-local-code | med | med | resolved | -- |
| 2 | Draft HIV NNRTI entries with citations | propose | med | med | med | 1-2hr |
| 3 | DONE 2026-09-02 -- 60 S-labelled carriers FOUND (first ever) over 20,816 labelled isolates; ALL 60 from ONE BioProject and killed by an aac(3) control -> LABEL_ARTIFACT. U1 SHARPENED, not retired | research | high | high | resolved | -- |
| 4 | DONE 2026-09-03 -- BV-BRC (independent on BOTH axes; 162/169 carriers new). 67 S carriers found; control verdict SPECIFIC_TO_RMT. Klebsiella PPV 0.475, E. coli 12/12 clean | research | high | high | resolved | -- |
| 5 | Re-call the 67 BV-BRC S-carrier genomes with AMRFinder -- the deployed rule consumes AMRFinder, not CARD | run-tests | high | high | med | 1-2hr |
| 6 | Investigate WHY a full-length rmtB sits at MIC<=1 in Klebsiella (silencing / expression / plasmid context) | research | med | high | high | days || 5 | Adopt an AMRrules-shaped per-entry schema for ONE catalog as a pilot (PMID + ECO evidence code + grade + limitations); costed against touching the frozen surface | propose | med | high | med | days |
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
| 2 | 2026-09-01 | edit-local-code | scripts/hiv_nnrti_mutant_catalog.py -- deconfounded NNRTI derivation | built; NNRTI was the only HIV class without one |
| 3 | 2026-09-01 | run-tests | measure curation vs the F-A doubt-layer baseline, 3 ways | NO-SHIP: every variant recovers < the free flag's 0.604 |
| 4 | 2026-09-01 | edit-local-code | fix self-to-self CompMutList parsing defect (L234L/K238K/M230M/R72R) | headline was flattering by 8->5 additions before the fix |
| 5 | 2026-09-01 | edit-local-code | C3: wiki/catalog_curation_procedure.md + tests/test_catalog_provenance.py | 4 conditions from 2 executed instances (gentamicin shipped / NNRTI declined); 1 enforced, 3 review discipline -- split stated explicitly |
| 6 | 2026-09-02 | run-tests | inverted rmt specificity hunt over PD AMR_genotypes + aac(3) project control | 60 counter-examples found and killed by the control (LABEL_ARTIFACT); PPV outside 146/146; v2 validation verified uncontaminated (0 overlap) |
<!-- project-state:end:action-log -->

## Open Questions for User
- **C1** — may a curated biological fact enter the shipped `hiv_amr.py`? Lock-safe, but it changes what the tool asserts.
- **C2** — is the gentamicin v2 lock authorized? It invalidates the prospective lock and the reproducibility freeze; a fix is an unfrozen revision needing its own validation and a NEW lock.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-08-31
- **Progress signal:** (none yet — blocked on F-A by construction)
