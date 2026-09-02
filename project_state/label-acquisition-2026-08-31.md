# F-E · Label Acquisition
<!-- project-schema: 0.1 -->

> Initialized 2026-08-31. Project ID: label-acquisition-2026-08-31. Originating goal (verbatim): "Address the binding constraint — labels, not models — by identifying and, if authorized, acquiring a label source that clears the eight rejection gates by construction."

## Project Context
- **Project ID:** label-acquisition-2026-08-31
- **Project root:** C:/Users/Farshad/PythonProjects/dna_decode
- **Captured:** 2026-08-31
- **Originating goal:** Address the binding constraint — labels, not models — by identifying and, if authorized, acquiring a label source that clears the eight rejection gates by construction
- **Refined goal:** Maintain a screened, ranked list of candidate label sources scored against the ten rejection gates, so that IF the user authorizes an acquisition the target is already chosen and its gate profile already known. Executor scope is screening only; acquisition is authority.
- **SCOPE CORRECTION 2026-08-31:** screening PEAR reclassified it OUT of this family's premise. PEAR is
  constructed variation with a continuous molecular endpoint (`constructed_molecular` -> regime `WORKS`),
  it is PUBLIC and FREE, and it therefore neither needs acquisition authority nor addresses the AMR label
  wall. It is a candidate **L4 forward-cell replication** substrate. This family retains the screening
  METHOD; PEAR itself belongs to the forward/inverse line. See `wiki/pear_substrate_screen_2026-08-31.md`.
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
| 1 | DONE 2026-09-01 -- PEAR scored mechanically by `scripts/screen_candidate_gates.py`; verdict INCOMPLETE (G6 unscreened), not `clears` | research | Soraya | done |
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
| E5 | PEAR BioProject PRJNA687219 resolves: E. coli K-12 MG1655, 45 SRA experiments, 478 Gbases raw reads | NCBI BioProject (fetched) | high | 2026-08-31 |
| E6 | PEAR is constructed single-gene variants w/ continuous relative-growth readout -> constructed_molecular regime | wiki/pear_substrate_screen_2026-08-31.md | high | 2026-08-31 |
| E7 | PEAR clears every applicable gate EXCEPT G6, which is UNSCREENED -> status INCOMPLETE, not `clears` (corrected 2026-09-01 by the mechanical screen; the memo headline overstated it, its own table did not). Blocker remains ARTIFACT FORMAT | wiki/pear_substrate_screen_2026-08-31.md + wiki/rejection_gate_screen_2026-09-01.md | high | 2026-09-01 |
| E8 | PEAR G6 PASSES on real values: mode-share 0.0019 / 2,106 distinct levels (n=2,114); ten-gate verdict CLEARS | wiki/pear_g6_screen_2026-09-01.md | high | 2026-09-01 |
| E9 | The "C: at 99%" blocker was a wrong-drive problem: D: has 4.1 TB; R installed there in minutes | wiki/pear_g6_screen_2026-09-01.md | high | 2026-09-01 |
| E10 | Forward cell EXTERNALLY replicated on a 2nd beta-lactamase: ESM2 CTX-M-14/cefotaxime rho 0.352 vs TEM-1 0.761 -- magnitude is protein-specific | wiki/pear_forward_replication_2026-09-02.md | high | 2026-09-02 |
| E11 | ESM2 beats BLOSUM62 on PEAR (0.352 vs 0.198): the learned model earns its keep, so the regime's DIRECTION holds | wiki/pear_forward_replication_2026-09-02.md | high | 2026-09-02 |
| E12 | Ceftazidime rho 0.078 is STRUCTURAL: CTX-M-14 is a cefotaximase and CAZ resistance is gain-of-function (IQR 0.080 vs CTX 0.180, max +3.68); a damage predictor cannot score gain | wiki/pear_forward_replication_2026-09-02.md | high | 2026-09-02 |
| E13 | No noise ceiling is derivable from PEAR's published tables -- Figure2B vs Figure3A correlate at EXACTLY 1.0 (same numbers, monotone transform, not replicates) | wiki/pear_forward_replication_2026-09-02.md | high | 2026-09-02 |
<!-- project-state:end:evidence -->

### Unknowns
- ~~**U1** — Whether PEAR's accessions resolve~~ **RETIRED 2026-08-31.** BioProject `PRJNA687219`
  resolves: *E. coli* K-12 MG1655 (taxid 511145), correct title, 45 SRA experiments, Sun Yat-sen
  University. GitHub `woson2020/CTXM-14` exists. The discipline paid, but not as written — the
  accessions are fine; the ARTIFACT FORMAT is the blocker (see U4).
- ~~**U4** — Whether the per-variant fitness values are extractable~~ **RETIRED 2026-09-01: YES.**
  R 4.6.1 installed via micromamba to `D:/tools/r_env` (D: has 4.1 TB free — the "C: at 99%" blocker was
  a WRONG-DRIVE problem, not a real wall). The ggplot objects carry their source data in `$data`;
  extracted to `D:/dna_decode_cache/pear/extracted/*.tsv` by `scripts/pear_extract_fitness.R`.
  Journal supplementary tables never needed.
- ~~**U5**~~ **ANSWERED 2026-09-01, partly unfavourably.** The plot objects carry the AGGREGATED
  per-variant effect sizes for the full scanned range — Figure.2A/2B nucleotide-level (792 positions x
  A/C/G/T = 3,168 substitutions per drug) and Figure3.A per-variant (2,114 in `C648T` notation with CTX
  + CAZ). They do **NOT** carry the ~23,000 raw barcoded strains. Quoting "23,000 variants" off this
  extraction would be wrong.
- ~~**U6** — G6 censoring~~ **RETIRED 2026-09-01 by measurement.** `assay_degeneracy` run on the real
  extracted values: mode-share **0.0019** over **2,106 distinct levels** (Figure3.A cefotaxime, n=2,114);
  all four variant tables pass. **G6 PASSES**; the ten-gate screen now returns CLEARS.
  TRAP RECORDED: including the wild-type baseline rows (`gt == isWt`, relative growth 1.0 by
  construction) alone puts mode-share at 0.2678 and TRIPS the bar — a censored-assay verdict manufactured
  entirely by the normalizer. Same defect class as the NNRTI `L234L` self-to-self entries, opposite
  direction. `wiki/pear_g6_screen_2026-09-01.md`.
- Whether any acquisition target is reachable without money.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | PEAR clears all ten gates | confirmed | 2026-09-01 |
| H4 | PEAR is an L1 label source that helps the AMR label wall | falsified | 2026-08-31 |
| H5 | PEAR's processed per-variant fitness table is directly downloadable | falsified | 2026-08-31 |
| H2 | A free path exists that clears the label wall (prospective accrual) | confirmed | 2026-08-24 |
| H7 | The forward cell's 0.761 is a general property of the genome-edit path | falsified | 2026-09-02 |
| H3 | Some acquisition target is reachable without money | open | (untested) |
| H6 | HBV has a free measured-phenotype source usable as a decoder label | falsified | 2026-09-01 |
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
| RESOLVED 2026-08-31: Authorize a PEAR acquisition? | Soraya | none -- moot | PEAR is PUBLIC and FREE (SRA + GitHub, no DUA). The authority gate was a consequence of misclassifying it as an acquisition. |
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
| 1 | DONE 2026-09-01 -- G6 screened on real extracted values; PASSES (mode-share 0.0019 / 2,106 levels). Ten-gate verdict CLEARS | research | med | high | resolved | -- |
| 1b | DONE 2026-09-02 -- ESM2 CTX rho 0.352 (BLOSUM62 0.198) vs the TEM-1 benchmark 0.761; regime DIRECTION holds, MAGNITUDE is protein-specific | run-tests | high | high | resolved | -- |
| 1c | DONE 2026-09-02 (Kaggle T4) -- ESM2+ProSST hybrid LOSES here (0.204 vs ESM2 0.352) because ProSST alone is at chance (-0.040); orthogonality premise fails on this protein. GEMME (MSA) not run | run-tests | med | high | resolved | -- |
| 2 | DONE 2026-08-31 -- U1 retired; PRJNA687219 resolves (E. coli K-12 MG1655, 45 SRA experiments) | research | high | high | resolved | -- |
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
| 2 | 2026-08-31 | research | screen PEAR: resolve accessions + classify its regime | RECLASSIFIED -- not an L1 label source and NOT an acquisition; see wiki/pear_substrate_screen_2026-08-31.md |
| 3 | 2026-08-31 | research | score PEAR against all ten rejection gates | clears every gate that applies; G6 (censoring) is the one open risk |
| 4 | 2026-08-31 | research | verify the processed data is reachable | BLOCKED on format: the GitHub .RData is a serialized ggplot object, not a table; R not installed |
| 5 | 2026-09-01 | research | screen HBV as a 6th viral cell against G1 (label source) | NO-BUILD: every free HBV resource is a rule or a prevalence table; no measured-phenotype source |
<!-- project-state:end:action-log -->

## Open Questions for User
- Whether to authorize a PEAR acquisition, or any non-public label source. Both are external and may involve money.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-08-31
- **Progress signal:** (none yet — init only)
