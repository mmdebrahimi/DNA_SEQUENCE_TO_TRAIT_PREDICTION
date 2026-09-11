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
| 1 | Screen a THIRD candidate through the gate screen (n=2 worked examples is the schema's stated limit) | research | med | high | med | 1-2hr |
| 2 | Check reachability of the NON-report-card surfaces (hiv/tb/pgx cards) | research | med | med | med | 1-2hr |
| 3 | Decide whether the doubt/lineage lines belong in `dna-decode` routes too | propose | low | med | med | 1-2hr |
| 4 | Report VME/ME (very major = R called S) beside sens/spec on the report card -- the clinical-microbiology convention that puts the dangerous error in the headline | edit-local-code | med | med | low | 1-2hr |
| 5 | Evaluate ECOFF-anchored tiering vs our clinical-breakpoint-only `mic_tiers.py` (EUCAST WGS subcommittee principle; ResFinder 4.0 errors cluster one dilution above the ECOFF) | research | med | high | med | days |
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
| 6 | 2026-09-09 | edit-local-code | salmserovar trust surface: corrected the 'resolves uniquely' contract falsehood on 4 surfaces (module docstring, printed CLI caveat, router validation string, report card) + the stale 'coverage 80' prose the argparse-default drift guard is structurally blind to; added 4 non-vacuity-proven guards; stratified the published 0.8400 by resolution route; pinned the SeqSero2 comparator's result headers; extended the prose guard to all 7 typing cells | 168 hits are 165 unique / 3 arbitrary-tie-break / 0 fallback -> headline SURVIVES, tie-break-independent accuracy 0.8250. Recorded against it: the arbitrary-winner policy buys 3 hits by converting 5 abstentions into confidently-wrong calls (3/8), the losing side of this cell's own pre-registered asymmetry -- left OPEN as an acceptance-bar question. Fallback zero is a measurement (reached 6x, failed each) so C3 disclosure work has zero measured headroom. Largest no-call cluster = 4 isolates on 4:i:- (monophasic Typhimurium), named not taken. Suite 4495 pass / 0 fail; 3 commits pushed; frozen AMR surface byte-unchanged. |
| 7 | 2026-09-10 | edit-local-code | Derived data inventory shipped (scripts/data_inventory.py + dna_decode/data/inventory_notes.py + wiki/data_inventory.{md,json}); pointers in CLAUDE.md orientation + LESSONS_LEARNED.md. Logged here as the machinery-legibility family. CORRECTION (same day, before any other run): this row first claimed dna-decode-2026-05-11 was schema-drifted with TWO end:action-log markers. FALSE -- I grepped the bare marker NAME, which also matched Action Log row 668, whose TEXT describes inserting that marker. Counted with the delimiters, as /project-state Step 3 actually counts, all six markers there are unique and that ledger is fine. This was the SIXTH instance in one day of prose-about-a-thing being read as the thing, and the only one I published as fact before checking -- the pre-publish falsification rail exists for exactly this and I skipped it. | 248 datasets across 4 roots; consumed_by/reported_in DERIVED per dataset; 59 referenced NOWHERE (largest 13.2GB) = the forgotten-asset bucket; 28 curated, 220 UNDOCUMENTED and reported as such. Self-contamination bug hit FIVE times (emitted artifacts -> orphans collapsed 59 to 8 while printing success; then the notes module; then the scanner's own comment; then the guard test; then the sentence describing the limitation) -- hand-listed exclusion replaced by a DERIVED filename pattern. Suite 4522 pass / 0 fail. |
| 8 | 2026-09-10 | edit-local-code | Soraya --advance: acted on the /innovate survivor DOUBT-bacterial-gap. RETRACTED its headline (the kill-test called target_site_doubt directly; no shipped surface makes that call for a bacterial drug, and bacterial cells DO carry a scored doubt_layer via trust_block) and fixed the two real defects the retraction uncovered. Then measured the ARBANK-source-conc survivor. Also relocated my own previous commit's inverted-polarity kill-tests out of tests/. | (1) lenacapavir (mutant-level CAI) was told "this catalog is position-based" two lines above its own MUTANT-LEVEL caveat; (2) NOT-MEASURED rendered as SILENCE on 4 cells (sarscov2-mpro + 2 fungal since 2026-09-02, plus lenacapavir) -- doubt_one_line knew applicable/assessed but not measured. Both fixed, 15 guards, non-vacuous (4 fail on revert). Registering a cell in the doubt map ALONE crashed a real call with a bare KeyError -> parity test. AR Bank source concentration MEASURED: 7 of 13 fail the 0.60 bar, NOT the predicted 13 of 13 -- gono arm single-panel (1.000), Klebsiella/E. coli PASS; grouped by CDC panel_id, a stated proxy for BioProject. 4539 pass / 0 fail after moving the 6 by-design-failing probes to probes/. |
| 9 | 2026-09-10 | edit-local-code | Soraya --advance: executed candidate-actions row 4 -- VME/ME (very major error = R called S) rendered beside sens/spec on the report card, ordered worst-VME-first. New pure module dna_decode/eval/error_rates.py reusing clonality.wilson_ci; augment-only under its own `error_rates` key on SCORED cells. | Adds NO information by construction (VME = 1-sens, ME = 1-spec exactly) -- the value is FRAMING, and it changes what the card leads with: klebsiella x meropenem reads "sens 0.467" in the old table and VME 0.533 [0.361-0.698] on 30 resistant isolates at the TOP of the new one, i.e. more than half of carbapenem-resistant isolates would be reported susceptible. E. coli cipro shows the inverse profile (VME 0.067 / ME 0.300). AUGMENT-ONLY VERIFIED BY DIFF: 27 cells before/after, identical key set, ZERO non-error_rates field changes, 10 cells gained the block. NO acceptance bar asserted (regulatory ceilings exist but are not verifiable from this repo) and a test forbids one appearing. 9 tests + report-card suite 38 green. [done:cc2777c1] |
| 10 | 2026-09-10 | edit-local-code | Soraya --advance: executed candidate row 2 (reachability of the NON-report-card surfaces) and went past it -- new dna_decode/data/cell_evidence_line.py renders a cell's MEASURED registry evidence as one line, wired into 11 CLI routes beside (never replacing) their hand-written caveats. Also RECOVERED an 11-file NUL-fill corruption + a corrupt .git/index found mid-run. | Two hypotheses I checked were WRONG and checking beat publishing: HCMV IS reached (via a dynamic card pattern) and the 7 unread cards belong to cells with a different call shape. The real gap was one level down -- non-AMR CLIs read neither trust_surface nor cell_registry, so measured evidence reached no user. 11 of 44 routes wired; the 33-route shortfall is asserted BY TEST so a green suite over a partial subset cannot read as completeness. Guard caught my own silent no-op (console script dna-pneumo-serotype vs package pneumoserotype -> None -> printed nothing). Verified on real FASTA runs, not just tests. Corruption: 11 tracked files 100% NUL + corrupt index, all recovered from HEAD (zero loss), disk-full ruled out at 35.5GB free, root cause NOT established. Suite 4582 pass / 0 fail. [done:5ee2565a] |
| 11 | 2026-09-10 | run-tests | Soraya --advance: executed the deployable half of the /innovate survivor POOLING-dilution on cached embeddings (no new compute), against a bar frozen BEFORE any number existed (wiki/pooling_vs_locus_acceptance_bar.json). | Verdict REFUTED_LOCUS_IDENTITY, code-computed. N=123 cipro strains / 71 lineages, leave-one-lineage-out: whole-genome pooled 0.8029, QRDR-restricted 0.8386, random-4-gene null max 0.9020 -- A RANDOM DRAW OF 4 ARBITRARY GENES BEAT THE 4 THAT CAUSE THE RESISTANCE. Narrowing helps; locus IDENTITY does not. The random-gene CONTROL was the load-bearing arm and it killed the claim. Does NOT contradict the 0-for-5: leave-one-lineage-OUT is weaker than the within-lineage question, and the strict test is barely answerable here (5 of 71 lineages mixed-label, largest 25R/1S) -- infeasibility computed into the artifact, not asserted. TWO defects in my own code, both found by RUNNING: the null drew on strain-unique gene_id so it silently ran ZERO draws (0%-overlap trap hitting the control, not the model), and two annotation vocabularies name gyrA differently with one colliding with parC (gyrase IS a type-II topoisomerase) -- fixing both took usable strains 78 -> 123. 16 tests, both key guards proven non-vacuous. Sweep artifact corrected: 2nd of 6 survivors to fall once acted on. Suite 4598 pass / 0 fail. |
| 12 | 2026-09-10 | run-tests | Soraya --advance: executed COUPLING-not-popdesign, the last surviving mechanism claim of the /innovate sweep, as a DIRECT measurement replacing the proxy its kill-test used -- per-determinant lineage concentration over 621 E. coli genomes / 104 lineages, both arms from the SAME genomes and SAME lineage partition, against a bar frozen before any score existed. | SUPPORTED on all four pre-registered bands, every gap exceeding the MAXIMUM of 1000 arm-label permutations (headline point 0.4708 n=24 vs acquired 0.0788 n=92, gap 0.3920 vs perm max 0.1762). FIRST survivor of this sweep to survive a STRONGER test rather than fall (2 of the other 3 acted on were refuted). TWO defects found, and the first was one step from BEING the headline: raw coupling falls with prevalence in BOTH arms, which reads as 'the determinants driving most resistance are near-chance' -- a METRIC ARTIFACT. The null holds carrier count fixed, which controls the EXPECTATION but NOT the MAXIMUM; the largest lineage is 53 of 621, so a 400-carrier determinant cannot exceed 53/400 and has ~0.048 of coupling available TOTAL. Normalized by achievable, gyrA_S83L is 0.895 not 0.041 -- near-SATURATED, not decoupled -- while acquired at the same prevalence sit BELOW zero (blaTEM-1 -0.134). THE ARM GAP IS LARGER UNDER NORMALIZATION (0.5227 vs perm max 0.2184), and the unrestricted normalized mean EXCEEDS the prevalence-matched one, showing the matching was only ever compensating for the ceiling. Frozen bar's rule is on the RAW score and was applied unchanged; normalization ships as a SECONDARY band, never a retro-fit. Defect 2 found by reconciling my row count (156,471) against the BAR's (137,288): they differ EXACTLY by 19,183 unparseable symbols, and 12 of those are determinants AMRFinder CALLS (indels/frameshifts/nonsense/negative-coordinate promoter) -- ~a fifth of the chromosomal arm. Corrected filter ships as a SENSITIVITY band (30 families vs 24, verdict unchanged), not a silent reshaping of a frozen substrate. Establishes the mechanism's PREMISE only; population design and coupling stay perfectly correlated so regime.py's attribution is still under-determined. 21 tests, both fixes proven non-vacuous by re-injection. Suite 4619 pass / 0 fail. Frozen AMR surface byte-unchanged. |
| 13 | 2026-09-11 | edit-local-code | Soraya --advance: executed candidate row 1 -- screened a THIRD candidate (Oxford E. coli bacteraemia cohort, L1) through the rejection-gate screen, the schema's own stated n=2 limit. Chose the empty diagonal deliberately: PEAR is L4/CLEARS and HBV is L1/REJECTED-at-G1, so nothing had ever exercised the L1 gates PAST G1. | FOUND A REAL DEFECT IN THE GATES. The memo defines G2 as a class x SOURCE CONTINGENCY ('one BioProject supplies most of ONE CLASS') but the implementation tested largest_source_share of the WHOLE COHORT -- a different quantity that diverges hardest at one source: Oxford scores share 1.00 and TRIPPED, though with no source variation the source cannot explain any label variance, so the named confound is STRUCTURALLY IMPOSSIBLE. It rejected a cohort we had already validated against (gent acc 0.990, 2026-06-15) for a confound it cannot have, and made G2 unsatisfiable by ANY single-source cohort. Fixed: n_sources<=1 -> not_applicable, memo's contingency preferred, cohort share demoted to a labelled weaker fallback. Single-source is NOT declared harmless -- it is a GENERALIZABILITY limit owned by the source_concentration disclosure layer + source_diverse_validate, and Oxford proves the concern real (0 rmt carriers in 4,979). Both committed screens still REPRODUCE (--verify). RESULT: Oxford REJECTED on G7 -- a THIRD failure mode. HBV dies because no measured phenotype EXISTS; Oxford's labels are as good as this project has (wet-lab BMD MIC, ascertained on bacteraemia = independent of the phenotype, 192R/2681S gent, ZERO breakpoint-censored) and it dies on THIN METADATA: the deposit carries only guuid, so no provenance-disjoint split can be built from it. G7 REPORTED not explained away -- the distinction that keeps the G2 fix from being special pleading is structurally-impossible (bug) vs genuinely-absent (finding). UNIT TRAP: the MIC columns are log2 dilution indices (negative values), and applying CLSI breakpoints directly returned R=0 on all three drugs -- a clean powered all-susceptible cohort rather than an error; the convention (MIC=2**upper) was READ OUT of scripts/oxford_score.py rather than guessed. Limits: screens the DEPOSIT not the cohort (ENA may hold the submitter fields, not fetched), G8 unmeasured, n=3 still small. 9 tests, fix proven non-vacuous by re-injection (4 fail). Frozen AMR surface byte-unchanged. [done:cd09424e] |
<!-- project-state:end:action-log -->

## Open Questions for User
- Whether a single-source SCORED cell warrants more than disclosure (3 of 10 rest on one BioProject).

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-08-31
- **Progress signal:** (none yet — init only)
