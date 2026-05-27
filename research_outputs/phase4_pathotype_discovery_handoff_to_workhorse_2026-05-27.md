# Phase 4 Pathotype Discovery — Handoff Back to Workhorse Machine (2026-05-27)

## Purpose

Closes the loop on the handoff your machine sent at `dna_decoder_phase4_handoff_to_other_machine_2026-05-27.md`. This is the discovery-lane deliverable in the form the workhorse asked for: dataset access feasibility + label recoverability/class map + cohort-shape proposal + first heavy execution task.

**Role split confirmed:**
- Discovery lane (GTX 860M laptop) = this machine — produced this memo + the supporting artifacts
- Execution workhorse (Precision laptop) = your machine — receives this memo + executes the named first heavy task

## TL;DR

E. coli pathotype prediction IS the right first non-AMR phenotype to promote — with one important scoping discipline: **do NOT commit 3 months of execution time on the v0 architecture without two cheap pre-build gates passing first**. The substrate is real and accessible; the deterministic-resolver architecture is shippable on top of pinned VirulenceFinder + ETECFinder; the wedge is narrow but defensible.

**Read-order safety note for the workhorse:** Section 4 names the first heavy execution task (Horesh cohort build) but that task is **gated** on Section 5's Gate A + Gate B both passing. If you skim top-down and miss Section 5, you will start building a 100-200 GB cohort that may be wasted work if the gates fail. Read Section 5 before acting on Section 4.

## 1. Dataset access feasibility — VALIDATED with caveats

| Source | Accessibility | Per-record label-provenance status | H1 (≥70% independent) verdict | Tier |
|---|---|---|---|---|
| **Horesh 2021 curated 10,146 E. coli collection** | Public Figshare DOI 10.6084/m9.figshare.13270073 (4.8 MB CSV File F1 + assemblies via NCBI `Assembly_name`) | BUILT-IN: Pathotype column has `(predicted)` suffix marking gene-rule-derived labels. Independent labels = 2,077 records (20.47% of full collection); 97.6% publication-attributed | **PASS** on the 2,077-record independent-label subset (by construction). FAIL on full collection. | **Tier-1a** (promoted from Tier-2 post-audit) |
| **Whittam STEC Center DECA collection** | Direct contact required (stec@cvm.msu.edu); MTA + non-export terms unknown | MLST clonal-lineage labels predate modern VF DB; expected H1-passing | **PASS** (expected) | Tier-1b — load-bearing for STEC/EHEC/EAEC supplement |
| **von Mentzer 2021 ETEC reference** | Public NCBI BioProject + 7 long-read reference genomes | Clonal-lineage labels + CF/toxin annotation | **PASS** | Tier-1b |
| **N=5 prototype reference strains** | Public type-strain repositories (EAEC 042 / ETEC H10407 / EPEC E2348-69 / STEC EDL933 / EIEC E11) | Curated literature labels | **PASS** | Tier-1b — positive controls only (N too small to drive cohort) |
| **NCBI Pathogen Detection** | Public; SOLR-indexed isolate browser | No discrete pathotype field; `host_disease` + `isolation_source` + AMRFinderPlus virulence call only | Indirect (needs facet query) | Tier-3 — pending Action 4 on this machine |
| **EnteroBase** | Public; large E. coli database with BlastFrost pathovar assignment | STRUCTURALLY CIRCULAR (BlastFrost uses same gene markers as v0 resolver would) | **FAIL** | **Tier-4** — usable for clonal-lineage + cohort assembly ONLY, NOT pathotype labels |

**Headline:** the Horesh File F1 audit converted "Horesh is hybrid-labeled and partially circular" (substrate-survey memo Row 17) into a clean per-record verdict. The 2,077-record independent-label subset is H1-passing by construction. Detailed audit at `research_outputs/horesh-2021-file-f1-audit-2026-05-27.md`.

## 2. Label recoverability / class map — RECOMMENDED 11-CLASS HONEST SURFACE

Post-/brainstorm round 2 + T2 resolution, the v0 output contract is:

| Output class | Required cluster combination |
|---|---|
| `EHEC_COMPATIBLE` | `STX1` ∨ `STX2` ∧ `LEE` |
| `STEC_NON_LEE` | `STX1` ∨ `STX2`, no `LEE` |
| `tEPEC_COMPATIBLE` | `LEE` ∧ `BFP_EAF`, no `STX` |
| `aEPEC_COMPATIBLE` | `LEE`, no `STX`, no `BFP_EAF` |
| `ETEC_COMPATIBLE` | `LT` ∨ `ST` |
| `EAEC_COMPATIBLE` | `aggR` + (`AAF_*` ∨ `aatA` ∨ `aaiC`) OR `AAF_*` + (`aatA` ∨ `aaiC`) |
| `UPEC_COMPATIBLE` | ≥2 strong ExPEC markers (papC/G, sfa, afa, hlyA, cnf1) + isolation context recommended |
| `HYBRID` | ≥2 primary DEC modules confidently present |
| `AMBIGUOUS` | partial markers / weak ID-cov / contradictory / single weak support |
| `UNCLASSIFIED` | EIEC/DAEC/AIEC-like markers OR unmapped |
| `COMMENSAL_LOW_MARKER_BURDEN` | no DEC module + sub-threshold ExPEC + assembly QC PASS |

`--legacy-6class` flag collapses to the originating goal's {EPEC/EHEC/ETEC/UPEC/EAEC/commensal} surface for backward compatibility.

Full marker catalog (≥23 markers) + decision table + abstention rules + provenance JSON schema is in `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md` under `## v0 Output Contract`. **Treat that section as a RECOMMENDATION from the discovery machine, not authoritative spec** — your machine owns the final implementation contract.

### Smallest credible first class map (if anything in Tier-1b doesn't materialize)

Drop to 7-class on Tier-1a-only: `ExPEC_COMPATIBLE / EPEC_COMPATIBLE / ETEC_COMPATIBLE / AMBIGUOUS / UNCLASSIFIED / COMMENSAL_LOW_MARKER_BURDEN / OUT_OF_SCOPE_PATHOTYPE` (where STEC/EHEC/EAEC roll into OUT_OF_SCOPE_PATHOTYPE). v0 ships honest at 7-class; widens to 11-class when Whittam DECA arrives.

## 3. Cohort-shape proposal — N≈350-500 strains, 5+ clade-balanced folds

| Pathotype class | Tier-1a (Horesh independent) | Tier-1b (DECA + von Mentzer + prototypes) | Tier-3 (NCBI Pathogen Detection) | Target | H2 floor |
|---|---:|---:|---:|---:|---:|
| ExPEC / UPEC_COMPATIBLE | 1574 (downsample to ~150) | — | ~50 (commensal-context urinary isolates) | ~150-200 | ≥50 ✓ |
| EPEC (tEPEC + aEPEC) | 269 | ~10 prototype + curated | — | ~50-100 | ≥50 ✓ |
| ETEC | 183 | ~50 (von Mentzer 7 lineages × representatives) | — | ~50-100 | ≥50 ✓ |
| EAEC | 2 | ~10 from DECA EAEC clones | ~30 (if available) | ≥50 (pending DECA) | ≥50 (PENDING) |
| EHEC | 6 | ~30 from DECA EHEC clones (O157:H7 + O26 + O111) | — | ~50 (pending DECA) | ≥50 (PENDING) |
| STEC (non-LEE) | 31 | ~20 from DECA STEC clones | — | ~50 (pending DECA) | ≥50 (PENDING) |
| Commensal | 2 | ~10 reference commensals (K-12 etc.) | ~30-75 from NCBI host_disease=healthy | ≥75 | ≥75 (PENDING) |
| **Total** | ~430 from Horesh | ~130 from DECA/von Mentzer/prototypes | ~80-150 from NCBI | **~350-500** | — |

### Exclusion rules
- Drop records where `Pathotype` is `(predicted)` from Horesh (gene-rule-derived → H1 circular).
- Drop records where `Isolation` is `Unknown` from Horesh (provenance-weak).
- Drop records with assembly QC fail (N50 < 50kbp OR contigs > 500).
- Drop EnteroBase records entirely from the label-substrate pool (Tier-4: clonal-lineage context only, not labels).

### Split constraints
- 5+ Mash-clade-balanced folds (reuse Mash-cluster-N=147 infrastructure pattern from Phase 1).
- No single fold may be missing >2 target classes.
- ExPEC must be ST-stratified across at least ST73 + ST131 + ST95 + others (top-3 ST share = 43.2% on Horesh independent subset — manageable).
- STEC must be ST-diversified via Tier-1b (Horesh STEC is 90% single-ST = severe leakage risk without DECA supplement).
- Within-fold class balance allowed to vary; across-fold marginal pathotype distribution must match overall cohort.

## 4. First heavy execution task for the workhorse — RUNS ONLY AFTER GATE A PASSES (see §5)

> **PRECONDITION (load-bearing):** **DO NOT START THIS TASK UNTIL GATE A IN §5 PASSES.** Gate A is a ~4-day no-code sanity check (manual VF 3.2.0 + ETECFinder on 5 known-pathotype strains; hand-write the 11-class decision-table verdict). If Gate A fails — VF/ETECFinder install friction, DB schema drift, decision table not computable from raw outputs, or Horesh per-strain labels not usable — **the heavy lift below is wasted work** and the project needs scope re-evaluation, not cohort assembly. Gate B (cold-email user-pain validation, §5) runs in parallel with Gate A. Heavy lift fires only when BOTH gates PASS.

> **Build the Horesh 2021 H1-passing cohort substrate.** Download FASTA assemblies for the 2,077 independent-label records from Horesh File F1 using the `Assembly_name` column (NCBI Datasets API). Run pinned VirulenceFinder 3.2.0 + ETECFinder DB calls on each. Materialize per-strain cluster profile + provenance JSON in the v0 schema drafted in the project ledger. Output: `data/processed/pathotype_horesh_h1_cohort/` directory with N≤2,077 genomes + annotations + VF/ETECFinder outputs + provenance JSONs.

**Estimated effort:** 1-2 days at workhorse bandwidth/storage. **Disk budget:** ~100-200 GB (assembled E. coli genomes average 5 MB compressed; FASTA + Bakta annotations + VF outputs roughly 30-50 MB per strain).

**Why this first AMONG post-gate tasks:** unblocks every downstream falsifier (H2 / H3 / H4 / H5). Without the cohort built + characterized, no decision-table tuning, no abstention-rule calibration, no clade-balanced validation, no VF-concordance check. This is the single biggest workhorse-shaped task — but it is NOT the first thing to do; Gate A is.

## 5. Two pre-build gates (RUN BEFORE COMMITTING TO HEAVY LIFT)

A parallel `/idea-validation-council` 6-lens review on this discovery machine produced verdict **"Pursue test"** — meaning don't commit 3 months of v0 implementation before two cheap gates pass:

### Gate A — Week-1 manual sanity (~4 days, no code, workhorse-doable)
1. Manually run VF 3.2.0 + ETECFinder on 5 known-pathotype strains from Horesh 2021 (1 each: EHEC / ETEC / EAEC / UPEC / commensal). Save raw outputs.
2. Hand-write the 11-class decision-table verdict for those 5 strains using ONLY the raw outputs.
3. **PASS:** tools install + run on workhorse + per-strain decision table is computable + Horesh per-strain labels are usable.
4. **FAIL any step → STOP, reassess.** Possible failure: VF DB schema may have drifted; ETECFinder integration friction; decision table may not be computable from raw outputs alone.

### Gate B — Cold-email user-pain validation (~1 week, no code, discovery-machine-doable but workhorse may have better contacts)
1. Email 5 target users (PulseNet, GenomeTrakr, ECDC-adjacent lab, nf-core E. coli pipeline maintainer, clinical-micro researcher).
2. Send the v0 one-pager + ask: "Would you install + use this in your pipeline within 60 days?"
3. **PASS:** ≥2/5 say yes.
4. **FAIL:** ≤1/5 say yes → KILL or PIVOT (consider PR-against-ECTyper as alternative distribution path).

### Why the gates matter

- Competitor research surfaced a serious **PathotypeR canary**: somebody shipped a deterministic E. coli pathotype resolver in R, got 0 stars / 2 commits / dormant. This is a strong negative adoption signal.
- ECTyper v2.0 (May 2025) covers DEC pathotype + per-gene %ID/contig/position with public-health adoption + Galaxy/IRIDA integration. The wedge is narrower than the brief framed — the audit-packet/abstention/honest-surface layer is empty in 2026, but DEC pathotype calling itself is solved.
- Without Gate B's user-pain signal, the project risks shipping a thin wrapper nobody installs.

Full lens-by-lens evidence + sources are in the discovery-machine conversation log; the verdict ("Pursue test") + the 5-bullet rationale are the load-bearing outputs.

## 6. What NOT to do on the workhorse (mirroring your earlier list back)

Do **not**:
- Skip the two pre-build gates and commit 3 months of v0 implementation directly.
- Treat the discovery machine's v0 output contract as authoritative spec — the workhorse owns the final implementation contract. Refine + adjust as the cohort reveals real edge cases.
- Build the trained-classifier sidecar (Candidate 2 / Candidate 3) before H1 PASS is confirmed on at least one cohort subset.
- Spend time on diff layer + side-by-side ambiguity-mapping in v0 — defer to v0.1+ per Executor lens scope-cut.
- Commit to 4 falsifiers (H2 + H3 + H4 + H5) in 3 months — cut to H2 only for v0; defer the other three to v0.1.

## 7. Pending open items on the discovery machine (this machine continues to handle)

| # | Open item | Owner | ETA |
|---|---|---|---|
| 1 | **NCBI Pathogen Detection `host_disease` facet query** for E. coli commensal substrate density | this machine (claude) | hours |
| 2 | **Whittam STEC Center direct-contact follow-up** (user-driven email send) | user | 1-2 weeks response window |
| 3 | EP-1 SUSPEND-gate audit-gate reuse + pathotype-specific opacity definition | this machine (claude) | hours |
| 4 | Cold-email Gate B (5 target users) | user (with discovery-machine support drafting one-pager) | 1 week response window |

Items 1 + 3 fire on this machine without further user input. Item 2 (Whittam) email draft is at `research_outputs/whittam-stec-center-contact-draft-2026-05-27.md` and ready to send. Item 4 needs user to identify the 5 target contacts.

## 8. Cross-machine state pointers

- Project ledger (v3): `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md`
- Substrate survey (Mission Control L1 run): `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` + raw + unsupported sidecars
- Horesh F1 per-record audit: `research_outputs/horesh-2021-file-f1-audit-2026-05-27.md`
- Whittam STEC Center email draft: `research_outputs/whittam-stec-center-contact-draft-2026-05-27.md`
- Followup queue (4 source memos, 20 active candidates): `research_outputs/_followup_queue.md`
- Run audit trail: `mission-control-runs/2026-05-27-0140-research-ecoli-pathotype-substrate/{intent-contract,audit-trail}.md`
- This memo: `research_outputs/phase4_pathotype_discovery_handoff_to_workhorse_2026-05-27.md`

Git commit `8971d22` on `main` carries all of the above (except this memo — separate commit incoming after writing).

## One-paragraph transfer summary

E. coli pathotype prediction IS the right first non-AMR phenotype: dataset access is validated (Horesh 2021 publication-extracted independent-label subset = 2,077 records, audit-confirmed accessible from public Figshare; Whittam DECA in flight; NCBI Pathogen Detection facet pending); class map is locked at 11-class honest surface with `--legacy-6class` flag (smallest-credible-fallback = 7-class on Tier-1a alone if Tier-1b doesn't materialize); cohort target is ~350-500 strains balanced across 6-7 classes with 5+ Mash-clade folds (Mid-term M4); first heavy execution task for the workhorse is **build the Horesh 2021 H1-passing cohort substrate** (N≤2,077 genomes + Bakta annotations + pinned VF 3.2.0 + ETECFinder calls + provenance JSONs, ~1-2 days at workhorse bandwidth/storage). BUT: a parallel 6-lens validation council returned verdict "Pursue test" — DO NOT commit 3 months of v0 implementation before two cheap pre-build gates pass (Gate A = manual sanity on 5 strains; Gate B = cold-email user-pain validation on 5 target contacts). Run the gates first; the heavy lift fires only on both-PASS.
