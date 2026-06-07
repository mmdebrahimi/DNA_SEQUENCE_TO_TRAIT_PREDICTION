# Eukaryotic Trait-Decoding Cycle
<!-- project-schema: 0.1 -->

> Initialized 2026-06-07. Project ID: eukaryotic-trait-decoding-cycle-2026-06-07. Originating goal (verbatim user input): "Eukaryotic trait-decoding cycle for dna_decode — extend the validated deterministic genome→phenotype decoder beyond bacterial AMR into the eukaryotic kingdom via two parallel substrates across two machines with phase gates. Path A (laptop, no GPU): fungal AMR — Candida auris azole resistance ... Path B (workhorse GPU): Arabidopsis thaliana flowering-time embedding test ... phase gates G0/G1/G2 ... money gate on paid compute; do not route personal code through the Bombardier/DLP machine."

## Project Context
- **Project ID:** eukaryotic-trait-decoding-cycle-2026-06-07
- **Project root:** C:\Users\Farshad\PythonProjects\dna_decode
- **Captured:** 2026-06-07
- **Originating goal:** Extend the validated deterministic genome→phenotype decoder beyond bacterial AMR into the eukaryotic kingdom via two parallel substrates across two machines with phase gates (Path A fungal AMR / Path B Arabidopsis embedding test).
- **Refined goal (from 3c top-ranked candidate):** Across two machines with no duplication + share-and-decide phase gates: (A, laptop, no compute) validate a deterministic fungal-AMR decoder on C. auris azole resistance (G1: acc≥0.80/sens≥0.80 on a de-confounded WGS+MIC cohort OR documented efflux/aneuploidy failure mode); (B, GPU workhorse, gated on compute) run the first methodologically-clean embedding-niche test on Arabidopsis flowering-time (G2: embedding R² beats SNP-PRS + kinship baselines under clade-stratified CV, PASS/FAIL).
- **Horizon (months):** 12
- **Schema:** project-schema 0.1

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** attempted (web-verified by the 2026-06-07 substrate-feasibility survey, `research_outputs/eukaryotic-multimodal-substrate-feasibility-2026-06-07.md`)
- **Provisional:** NO
- **Findings:** (none detected) — C. auris ERG11↔MIC linkage (181/188), India/S.Africa WGS+MIC depth (188-350), Arabidopsis flowering-time depth (1,003 acc), plant-DNA-FM VRAM class (24-80GB GPU), and the project's internal embedding 0-for-3 record are all sourced/verified.

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Bounded — TWO named substrates (C. auris azole; Arabidopsis flowering-time) with measurable, dated phase gates (G0 caller works / G1 fungal-AMR acc≥0.80 / G2 embedding PASS-FAIL), not the unbounded "decode all eukaryotes." The broader kingdom ambition lives only in the long-term tier; the actionable scope is two falsifiable tests.

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:**
  1. **Path A — fungal AMR decoder (C. auris azole)** validated at G1 (no compute; determinant-scan transfer). **Top-ranked / primary win.**
  2. **Path B — Arabidopsis flowering-time embedding test** at G2 (GPU-gated; the embedding thesis's decisive 4th test).
  3. **Fungal-AMR infra** (catalog ✅ + BLAST caller G0) as the enabling sub-goal for #1.

## Goal Hierarchy
### Long-term (12+ months tier)
Establish whether the dna_decode decoder generalizes across the eukaryotic kingdom — both via the proven deterministic determinant-scan (fungi) and, decisively, whether frozen DNA-foundation-model embeddings ever earn their keep (plant quantitative phenotype) — feeding the broader "decode any trait, any organism" north star.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | G0 — fungal BLAST caller works (laptop) | `scripts/fungal_erg11_caller.py` calls a known ERG11 mutation on a C. auris reference genome | weeks |
| 2 | G1 — fungal AMR validated (laptop) | deterministic azole decoder acc≥0.80 / sens≥0.80 on a clade-de-confounded C. auris WGS+MIC cohort, OR documented efflux/aneuploidy failure mode | 1-2 months |
| 3 | G2 — Arabidopsis embedding PASS/FAIL (workhorse, GPU-gated) | embedding R² beats SNP-PRS + kinship baselines under clade-stratified CV (PASS), or 0-for-4 + close embedding frontier (FAIL) | gated on compute |

### Short-term (≤1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Verify C. auris WGS+MIC cohort is extractable (supplementaries → accessions) — iron-law gate | research | laptop | days |
| 2 | Build `scripts/fungal_erg11_caller.py` (BLAST ERG11/FKS1 → translate → call vs catalog) → G0 | edit-local-code | laptop | days |
| 3 | Pre-stage Path-B AraGWAS download manifest + baseline spec for the workhorse | write-plan | laptop | days |
| 4 | Resolve workhorse identity (Precision 7780 vs Bombardier/DLP) — safety gate for Path B | ask-user | user | now |

## State Snapshot
### Assumptions
- Fungal azole-R is determinant-detectable enough for a BLAST scan (ERG11↔MIC near-perfect, but multi-locus efflux/aneuploidy caps sensitivity). confidence: medium
- A usable C. auris WGS+MIC cohort with downloadable accessions is extractable from the S.Africa/India supplementaries. confidence: medium
- A plant DNA-FM fits (or nearly) the workhorse's ~12GB GPU; else Path B → cloud/money. confidence: low
- git origin is a safe + sufficient cross-machine channel. confidence: high
### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| 1 | C. auris S.Africa 181/188 MIC>32 had ERG11 mutations | research_outputs/eukaryotic-multimodal-substrate-feasibility-2026-06-07.md | high | 2026-06-07 |
| 2 | Arabidopsis flowering-time GWAS: 1,003 accessions, public (AraPheno/AraGWAS) | same | high | 2026-06-07 |
| 3 | plant DNA-FM (PlantCaduceus) needs 24-80GB GPU class | same | high | 2026-06-07 |
<!-- project-state:end:evidence -->
### Unknowns
- Is the C. auris WGS+MIC cohort actually downloadable at depth (accessions in supplementaries)?
- Does a plant DNA-FM fit 12GB, or is cloud (money) required for Path B?
- What fraction of C. auris azole-R is ERG11-only vs efflux/aneuploidy (sets fungal sensitivity ceiling)?
### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | The deterministic determinant-scan method transfers to fungi (C. auris azole) at acc≥0.80 | open | never |
| H2 | Frozen DNA-FM embeddings beat PRS+kinship baselines on Arabidopsis flowering-time (the embedding niche exists) | open | never |
<!-- project-state:end:hypotheses -->
### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| Path C ratified (A now / B queued on compute) | 2026-06-07 | user-ratified via /soraya |
| Two-machine split, no duplication, git-only sync | 2026-06-07 | plans/Eukaryotic_DualMachine_Coordination.md |
<!-- project-state:end:decisions-made -->
### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Workhorse identity (Precision 7780 vs Bombardier/DLP) | Soraya | user | SAFETY: Path B handoff assumes personal Precision 7780; do NOT route personal code through Bombardier |
| Path B compute: 12GB GPU sufficient OR cloud budget | Soraya | user | money gate — no paid compute without explicit OK |
| Pre-commit: G2-FAIL closes the embedding frontier permanently (no 5th attempt) | Soraya | user | guards against the diminishing-returns trap on a 4th embedding bet |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Fungal determinant catalog shipped (`dna_decode/data/fungal_amr.py`, 7 tests); EP7/EP8 + dual-machine coordination plan + workhorse handoff written + pushed. Path A ready to build the BLAST caller; Path B queued on compute.

### Target state / terminal condition
G1 (fungal AMR decoder validated on C. auris, or documented failure) AND G2 (Arabidopsis embedding PASS/FAIL) both resolved, each shared + decided at its gate. Either substrate may independently succeed/fail/decline — both resolved = cycle complete.

### Progress proxy
- **v0.1 metric:** unknowns-retired (3 open) + gates-passed (0 of G0/G1/G2).

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | C. auris cohort feasibility check | research | med | high (gates G1 viability) | med | low |
| 2 | Build fungal_erg11_caller.py → G0 | edit-local-code | high | med | low | low |
| 3 | Pre-stage Path-B manifest | write-plan | low | low | low | low |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- Re-run `/project-state` after each gate (G0/G1/G2) + share result packet.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` — auto · `research` — auto · `write-plan` — auto · `edit-local-code` — REQUIRES per-action human approval · `run-tests` — auto if local · `ask-user` — auto · `stop` — auto

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-06-07 | propose | /project-init invoked (eukaryotic cycle) | ledger created; 3a PASS / 3b PASS(project) / 3c PASS |
| 2 | 2026-06-07 | research | Path A iron-law gate: C. auris WGS+MIC cohort feasibility | PASS — 841 C. auris assemblies on NCBI (taxon 498019, downloadable via refseq path) + MIC/ERG11 labels from S.Africa(188)/India(350) supplementaries. G1 substrate viable; the cohort-extractable unknown is RETIRED. Next: build scripts/fungal_erg11_caller.py (G0). |
| 3 | 2026-06-07 | edit-local-code | Gate G0 (MACHINERY) REACHED: fungal ERG11 BLAST caller works | scripts/fungal_erg11_caller.py: blastn(CDS-vs-genome) + gap-aware codon-translate + catalog match (tblastn absent so blastn used; C. auris ERG11 intronless = colinear). Validated against a KNOWN planted mutation via REAL makeblastdb+blastn: resistant genome (codon132 TAC->TTC) -> R, ERG11:Y132F detected; wild-type -> S with efflux/aneuploidy blind-spots; offline -> INDETERMINATE. 2 caller tests + 7 catalog tests green. G0-COMPLETION (next) = swap synthetic ref for real C. auris ERG11 allele + validate on a real resistant genome (confirms catalog numbering vs reference). Then G1 = cohort validation. |
| 4 | 2026-06-07 | edit-local-code | Gate G0-COMPLETION reached: fungal ERG11 caller validated on REAL data | Fetched real C. auris ERG11 reference (RefSeq XM_029033208.2, 525aa) + confirmed catalog numbering matches reference (Y@132, K@143, V@125). Validated call_erg11 end-to-end (real makeblastdb+blastn) on real GenBank isolate alleles: PV630306 WT->S (efflux/aneuploidy blind spots surfaced), PV630305 Y132F->R, PV630302 K143R->R = 3/3 correct. Numbering-mismatch risk RETIRED. Committed reference + 3 public allele fixtures (data/fungal_ref/) + 3 real-data regression tests (5 caller tests green). Next: G1 = build C.auris WGS+MIC cohort + validate decoder (acc>=0.80/sens>=0.80 OR documented efflux/aneuploidy failure). |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Workhorse identity (SAFETY):** is the GPU "workhorse" the personal Precision 7780, or the Bombardier/DLP machine? Path B handoff assumes the former; if the latter, Path B is blocked on principle (no personal code on a DLP machine).
- **Path B compute (MONEY):** is the ~12GB RTX 3500 Ada acceptable to attempt, or do you want a cloud A100 (paid — requires explicit budget approval)?
- **Embedding pre-commitment:** confirm that a clean G2-FAIL closes the embedding frontier permanently (this is the thesis's 4th test after 0-for-3; pre-committing guards against the diminishing-returns trap).

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-06-07
- **Progress signal:** (init only — catalog shipped, caller next toward G0)
