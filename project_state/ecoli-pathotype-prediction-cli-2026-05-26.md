# E. coli Pathotype Prediction CLI
<!-- project-schema: 0.1 -->

> Initialized 2026-05-26. Project ID: ecoli-pathotype-prediction-cli-2026-05-26. Originating goal (verbatim user input): "E. coli pathotype prediction from a user-supplied genome assembly: an open-source CLI tool that takes FASTA + optional GFF3 and emits a multiclass pathotype call (EPEC / EHEC / ETEC / UPEC / EAEC / commensal) with audit-grade provenance, including which acquired virulence-gene clusters drove the call (stx for EHEC; LEE for EPEC; afa/papC for UPEC; ETEC enterotoxins) and a side-by-side comparison against CGE VirulenceFinder gene-call output."

## Project Context
- **Project ID:** ecoli-pathotype-prediction-cli-2026-05-26
- **Project root:** C:\Users\Farshad\PythonProjects\dna_decode
- **Captured:** 2026-05-26
- **Originating goal:** E. coli pathotype prediction from a user-supplied genome assembly: an open-source CLI tool that takes FASTA + optional GFF3 and emits a multiclass pathotype call (EPEC / EHEC / ETEC / UPEC / EAEC / commensal) with audit-grade provenance, including which acquired virulence-gene clusters drove the call (stx for EHEC; LEE for EPEC; afa/papC for UPEC; ETEC enterotoxins) and a side-by-side comparison against CGE VirulenceFinder gene-call output.
- **Refined goal (if 3c produced one):** (no rewrite — 3c PASS; candidate refinements recorded in `## Refinement Candidates`). **Updated 2026-05-27 post-/research substrate-survey:** v0 architecture FINAL LOCK = deterministic multilabel cluster resolver + abstention (Candidate 5); substrate-tier strategy LOCKED (Tier-1 DECA + von Mentzer + prototypes / Tier-2 Horesh 2021 pending File F1 audit / Tier-3 NCBI host_disease / Tier-4 EnteroBase clonal-lineage only); classifier track (Candidates 2+3) deferred to v0.1+ conditional on H1 PASS. Concrete v0 output contract at `## v0 Output Contract` (23 markers + 11-class decision table + abstention rules + provenance JSON).
- **Horizon (months):** 12
- **Schema:** project-schema 0.1
- **Continuation context:** continuation of session "DNA Decode AI Project"; this project is EP-4 candidate per `plans/EP_4_Non_AMR_Phenotype_Candidates.md` + the long-horizon roadmap at `plans/Trait_Decoding_Roadmap.md`.

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** attempted
- **Provisional:** NO
- **Findings:**
  - VirulenceFinder is current and actively maintained (CGE service at cge.food.dtu.dk/services/VirulenceFinder/; **PyPI Virulencefinder 3.2.0 released 2025-06-04**; ETECFinder DB revision 2024; PulseNet 2.0 STEC validation 2025). Output schema = BLAST-table with %ID, query length, contig position, predicted phenotype, accession. **Side-by-side diff is mechanically feasible.** Version-pin discipline: provenance schema MUST capture VF software version + DB version + DB commit/checksum (rule-table drift across DB revisions is real per 2024 ETECFinder revision: 455 new alleles, 50 replaced/renamed, 2 removed).
  - Precision concern (not a goal-malforming error): the parenthetical "LEE for EPEC; stx for EHEC" is a simplification. LEE is shared between EPEC and EHEC; the actual discriminator is `LEE+ stx- (+ bfp+ = typical, bfp- = atypical) → EPEC` vs `LEE+ stx+ → EHEC`. Goal's "drove the call" framing implies cluster-combination logic, so this is recoverable at design time.
  - Precision concern: EAEC is listed as a pathotype but no defining gene cluster is named in the parenthetical. EAEC is defined by AAF (aggregative adherence fimbriae) variants + aggR regulator + aaiC + aatA. The CLI's gene-cluster catalog will need to include EAEC markers explicitly.
  - Precision concern: ETEC "enterotoxins" is under-specified. The discriminating markers are LT (heat-labile) and ST (heat-stable) toxin genes. ETECFinder DB (2024 revision) is the authoritative gene catalog.
  - **Biology-vs-framing tension:** the multiclass framing assumes mutually-exclusive pathotypes, but hybrid pathotypes exist in nature (e.g., the 2011 German outbreak strain O104:H4 was STEC/EAEC hybrid). This is a real factual-shape issue but is more naturally addressed as an architecture fork in `## Refinement Candidates` than as a goal rewrite, because the goal as written is actionable; the question is whether the implementation should output multiclass (with hybrid as a separate class) or multilabel (gene-cluster profile + decision rule).

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** bounded scope (single CLI surface, FASTA-in → pathotype-out + provenance + VirulenceFinder diff), measurable success criteria implied (multiclass label + audit), 12-month horizon, decomposable into milestones. No "understand all" / "decode any" framing. Substrate is well-defined (E. coli pathotype-labeled genomes, EnteroBase candidate). Output contract is concrete (CLI tool, not a research thesis).

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:** the originating goal is well-formed and actionable; refinement here surfaces an architecture fork. **Updated 2026-05-26 post-/brainstorm round 2:** the v0 default is Candidate 5 (deterministic multilabel cluster resolver + abstention); Candidates 3+4 are demoted to v0.1+ falsifier-track work; Candidate 2 (classifier) is deferred until label-provenance audit confirms non-circular labels. All candidates ship the same external CLI surface; they differ in the internal pathotype-decision layer.

  1. **Rule-based pathotype CLI (lowest-novelty, highest-determinism) — DEMOTED (subsumed by Candidate 5).** Pipeline: Bakta annotate → virulence-gene catalog match (VirulenceFinder DB + ETECFinder DB + curated additions) → hand-coded cluster-combination rules (`LEE+ ∧ stx+ → EHEC`; `LEE+ ∧ stx- ∧ bfp+ → tEPEC`; `LT+ ∨ ST+ → ETEC`; `AAF+ ∧ aggR+ → EAEC`; `papC+ ∨ afa+ ∨ hlyA+ → UPEC`; else → commensal). VirulenceFinder side-by-side is a parallel run + diff. **Contribution:** engineering polish (clean CLI, audit packet, hybrid handling). **Risk:** low scientific novelty; a 2010-era tool with nicer output. **Status:** subsumed by Candidate 5, which uses VirulenceFinder + ETECFinder as the pinned caller layer rather than reimplementing the gene-call layer with Bakta-anchored BLAST.

  2. **Classifier-based pathotype decoder (highest-novelty, highest-risk) — DEFERRED to v0.1+ pending label-provenance audit.** Train multiclass classifier on labeled E. coli genome corpus (EnteroBase pathotype labels OR ECOR + curated subset) using NT embeddings or gene-presence vectors. SUSPEND-style audit gate from EP1 surfaces model-vs-VirulenceFinder disagreements. **Contribution:** decoder-pattern continuation; discovers atypical / hybrid cases rules miss. **Risk:** EP-2 tet finding (NT fails on distributed mobile-element mechanisms) is a direct warning shot — EAEC virulence is distributed across mobile elements; classifier may underperform on it. Label substrate quality is the cohort-construction risk. **Status:** does not ship v0; if labels are gene-rule-derived (likely for EnteroBase / NCBI Pathogen Detection), classifier uplift is circular and meaningless. Re-evaluate after H1 passes.

  3. **Hybrid disagreement-detector — DEMOTED to v0.1+ (longer-term output architecture, not v0).** Rules are PRIMARY (always emit a rule-based call + provenance). Classifier runs in PARALLEL. Tool's loud output is the disagreement set: cases where rules + classifier conflict, or where VirulenceFinder marker presence is ambiguous (e.g., LEE+ stx+ bfp+ → genuine ambiguity). **Contribution:** the tool's value-add is finding edge cases, not re-doing pathotype calls. **Fit with dna_decode north star:** matches "AI DNA decoder, not papers" — the AI surfaces what the rules miss. **Status:** correct as v0.1+ direction, but requires Candidate 2's classifier to exist first; depends on H1 passing.

  4. **Multilabel virulence-cluster profile + commensal-residual — PROMOTED into Candidate 5.** Drop multiclass entirely; emit a multilabel vector over {LEE, stx, bfp, LT, ST, afa/papC, AAF/aggR, hlyA, cnf1, ...} + a "commensal-residual" flag if no virulence cluster fires above threshold. **Status:** the multilabel-internal core insight is preserved inside Candidate 5; the "no derived primary" version is rejected because clinicians/epi want a primary label.

  5. **Deterministic multilabel cluster resolver + abstention (RECOMMENDED v0 DEFAULT).** Pipeline: pinned VirulenceFinder + ETECFinder gene-calls → multilabel cluster profile (≥23 marker clusters; see `## v0 Output Contract` section below) → deterministic decision-table → derived pathotype label with explicit abstention rules (AMBIGUOUS / HYBRID / UNCLASSIFIED / COMMENSAL_LOW_MARKER_BURDEN). Audit packet + side-by-side VF diff. **No trained classifier required for v0.** **Contribution:** honest, version-pinnable, clinically legible CLI that ships in 3 months rather than 12. Classifier work (Candidate 2) and disagreement-detection (Candidate 3) layer on top in v0.1+ as falsifier tracks. **Fit with north star:** "AI DNA decoder tool, not papers; failure-tolerant iteration" — v0 is the engineering scaffold; AI augmentation is the next iteration, gated on labels being non-circular. **Risk:** the "side-by-side VF comparison" framing in the originating goal becomes weaker (VF is both internal caller AND baseline — see Pending Decision T1).

  **Recommended v0 default:** Candidate 5. Recommended v0.1+ track: Candidate 5 + Candidate 2 + Candidate 3 layered as falsifier sidecar. Recommended out-of-scope: Candidate 1 (subsumed) + Candidate 4 (subsumed).

## Goal Hierarchy
### Long-term (12+ months tier)
Ship an open-source E. coli pathotype prediction CLI that ingests a genome assembly (FASTA + optional GFF3) and emits an auditable pathotype call with virulence-cluster provenance + a side-by-side comparison against CGE VirulenceFinder. **v0 (3-6 mo) = deterministic multilabel cluster resolver + abstention (no trained classifier).** **v0.1+ (6-12 mo) = classifier sidecar + disagreement-detector layer, conditional on H1 PASS (labels are not gene-rule-derived).**

### Mid-term (3-12 months) — REORDERED 2026-05-27 (substrate-first; 3-month v0 ship target; Status column added)
| # | Milestone | Success Criterion | Horizon | Status |
|---|---|---|---|---|
| 1 | **Label-provenance audit + substrate selection** | H1 PASS verdict on chosen substrate (≥70% independent-label fraction); substrate source decision recorded in `## Decisions Made`; per-class isolate count availability documented at H2 floors (N≥50 each pathotype + N≥75 commensal) | ≤1 month | **PARTIAL 2026-05-27** — 4-tier substrate strategy LOCKED via /research substrate-survey; H1 PASS verdict pending Horesh File F1 audit + Whittam direct-contact follow-up |
| 2 | Architecture-fork final lock | Decision recorded in `## Decisions Made` confirming Candidate 5 (deterministic resolver + abstention) or revised candidate; H1 + T1 + T2 closed | ≤1 month (post-M1) | **COMPLETED 2026-05-27** — Candidate 5 FINAL LOCK; T1+T2 closed; substrate strategy provides path even with H1 partial |
| 3 | v0 marker catalog + decision rules + parser | ≥23 markers per `## v0 Output Contract` codified; 11-row decision table + abstention rules implemented; pinned VF 3.2.0 + ETECFinder DB + parser tests passing on frozen fixtures; H5 PASS | 1-3 months | open |
| 4 | Labeled-substrate cohort built + clade-balanced folds | N≥350 E. coli genomes (per-class minimums from H2); MLST + Mash distances populated for all isolates; ≥5 clade-balanced folds constructed; H3 substrate prerequisite satisfied | 2-4 months | open |
| 5 | Per-class precision/recall + abstention quality evaluation | H2 + H3 + H4 ship-gate verdicts produced; SUSPEND-style audit gate adapted with pathotype-specific opacity definitions | 3-5 months | open |
| 6 | v0 CLI ship + falsifier | Public-facing CLI w/ docs; falsifier verdict on held-out test set + external dataset (GenomeTrakr / NCBI Pathogen Detection labeled subset); 11-class output surface + optional `--legacy-6class` flag | **3-6 months** (down from 6-9; classifier track deferred to v0.1+) | open |
| 7 | (v0.1+ — DEFERRED) Classifier sidecar + disagreement-detector layer | Candidate 2 + Candidate 3 layered as falsifier track; only fires if H1 PASSES (otherwise classifier uplift is circular) | 6-12 months | deferred (conditional on H1 PASS) |

### Short-term (≤1 month) — REORDERED 2026-05-27 post-/research substrate-survey completion
| # | Action | Class | Owner | Horizon | Status |
|---|---|---|---|---|---|
| ~~1~~ | ~~Substrate / label-provenance survey~~ | research | claude | ≤2 weeks | **COMPLETED 2026-05-27** — `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` (4-tier strategy locked) |
| ~~2~~ | ~~Resolve T1 + T2~~ | ask-user | user | ≤1 week | **COMPLETED 2026-05-26** — T1 = option C; T2 = 11-class honest default |
| ~~4~~ | ~~Architecture-fork final lock~~ | ask-user | user | ≤1 week (post-Action-1) | **COMPLETED 2026-05-27** — Candidate 5 locked as v0 default |
| 1 (new) | **Horesh 2021 File F1 supplementary audit** — download File F1 from PMC8208696; parse per-record metadata to estimate independent-label fraction (target: ≥70% for Tier-2 status); count per-pathotype × per-ST distribution for H3 fold-construction feasibility | research | claude | ≤1 week | open |
| 2 (new) | **Whittam STEC Center direct-contact follow-up** — email stec@cvm.msu.edu requesting (a) DECA per-pathotype strain counts, (b) per-clone roster + curation provenance, (c) MTA terms + non-export restrictions | ask-user | user | ≤2 weeks | open |
| 3 | Confirm VirulenceFinder integration pathway (CLI wrap vs PyPI 3.2.0 import vs Bitbucket DB direct match); pin VF software + DB versions + DB commit/checksum | research | claude | ≤2 weeks | open |
| ~~4 (new)~~ | ~~NCBI Pathogen Detection host_disease facet query~~ | research | claude | ≤1 week | **PARTIAL 2026-05-27** — HONEST-GAP verdict per `research_outputs/ncbi-pathogen-detection-host-disease-facet-2026-05-27.md`. Programmatic facet access blocked. Three manual paths documented (Isolates Browser UI / BigQuery / pd-help email); user owns next step. |
| ~~5~~ | ~~Cross-check audit-gate reuse~~ | research | claude | ≤2 weeks | **COMPLETED 2026-05-27** — verdict ADAPTED_REUSE per `research_outputs/ep1-suspend-gate-pathotype-reuse-feasibility-2026-05-27.md`. Workhorse-side reuse recommendation provided. |

## State Snapshot
### Assumptions
- (confidence: high) VirulenceFinder is the right comparison baseline (active, well-known in field, gene-call output well-defined per 2024-2025 sources).
- (confidence: high) EnteroBase has labeled pathotype substrate sufficient for cohort construction at N≥150.
- (confidence: medium) The EP-1 mechanism × phenotype audit-gate pattern (SUSPEND verdict, opacity gating) generalizes from AMR-mechanism-call to pathotype-cluster-call.
- (confidence: medium) The chosen pathotypes (EPEC/EHEC/ETEC/UPEC/EAEC + commensal) cover ≥90% of real-world clinical E. coli isolates. (DAEC and AIEC are missing; intentional Phase-0 scope cap.)
- (confidence: low) NT embeddings will transfer to pathotype prediction (EP-2 tet finding suggests NT struggles on distributed mobile-element mechanisms, and EAEC is mobile-element-distributed).
- (confidence: low) Multiclass framing is OK as the EXTERNAL contract even if the internal representation is multilabel.
- (confidence: medium) The dna_decode Phase 1 infrastructure (Bakta + AMRFinder + cohort builder + audit-gate) covers ≥60% of what this project needs to ship.
- (confidence: high) Users typing the CLI command already have an assembled FASTA — reads-in is out of scope.

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| 1 | VirulenceFinder is active 2025, output schema = BLAST-table format with %ID/length/position/phenotype/accession | WebSearch 2026-05-26: PulseNet 2.0 STEC paper (PubMed 40572198, June 2025); ETECFinder revision (PMC11237473, 2024); PyPI Virulencefinder 3.1.0 | high | 2026-05-26 |
| 2 | ETECFinder DB (CGE-hosted, 2024 revision) is the authoritative ETEC gene catalog; added 455 new alleles | ResearchGate / PMC11237473 (2024) | high | 2026-05-26 |
| 3 | EP-2 finding: NT-frozen-whole-genome-pooling fails on distributed mobile-element mechanisms (tet 0.400 anti-predictive) | `wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md` | high | 2026-05-26 (transferred from prior project context) |
| 4 | EnteroBase E. coli pathotype labels are derived via BlastFrost gene-presence calls (same markers v0 deterministic resolver would track) → STRUCTURALLY CIRCULAR with v0 | `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Row 1 (PMC9393565) | high | 2026-05-27 |
| 5 | NCBI Pathogen Detection has NO discrete pathotype field; virulence calls come from AMRFinderPlus (gene-rule-derived) | `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Rows 5+6 (NCBI Pathogen Detection help) | high | 2026-05-27 |
| 6 | Horesh 2021 curated 10,146-genome E. coli collection: pathotype assignment uses HYBRID (ariba + VirulenceFinder DB + isolation-source refinement) | `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Row 17 (PMC8208696) | high | 2026-05-27 |
| 7 | Horesh 2021 pathotype representation: EPEC ≈ 3% + ETEC ≈ 2% of 10,146 = ~304 EPEC + ~203 ETEC strains | `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Row 16 (PMC8208696) | high | 2026-05-27 |
| 8 | Horesh 2021 lineage concentration: six STs (11, 131, 73, 10, 95, 21) = 50% of 10,146 isolates → H3 lineage-leakage risk | `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Row 18 (PMC8208696) | high | 2026-05-27 |
| 9 | Whittam STEC Center repository: 16,000 strains; DECA collection = 15 representative diarrheagenic clones (MLST clonal-lineage labeled, INDEPENDENT of v0 marker rules); per-pathotype counts NOT publicly documented | `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Rows 10+11 (Whittam STEC Center + Hazen 2012 DECA) | high | 2026-05-27 |
| 10 | von Mentzer 2021: 7 major ETEC lineages with long-read reference genomes; clonal-lineage + CF/toxin labels (INDEPENDENT of v0 marker rules) | `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Row 13 (Nature Sci Rep) | high | 2026-05-27 |
| 11 | Prototype reference strains with curated literature labels (N=5): EAEC 042, ETEC H10407, EPEC E2348-69, STEC EDL933, EIEC E11 | `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Row 12 (PMC7593697) | high | 2026-05-27 |
| 12 | Phase 4 pathotype Gate A PASS on workhorse (Precision 7780) 2026-05-28: real VirulenceFinder runtime (`python -m virulencefinder` + cloned VF DB + local blastn); 5-strain frozen panel all called correctly (ehec_edl933→EHEC_COMPATIBLE; etec_h10407→ETEC_COMPATIBLE; eaec_042→EAEC_COMPATIBLE; upec_cft073→UPEC_COMPATIBLE; commensal_mg1655→COMMENSAL_LOW_MARKER_BURDEN); caller_runtime_ok + raw_outputs_usable + decision_table_computable all true | `research_outputs/pathotype_gate_a_result_2026-05-28.json` + `dna_decoder_phase4_gate_a_pass_2026-05-28.md` | high | 2026-05-29 |
| 13 | Gate A runtime finding: working VF invocation is `python -m virulencefinder` + cloned VF DB + external blastn — NOT separate virulencefinder/etecfinder executables or Docker images; current etecfinder branch is an alias to the VF runtime scoped to virulence_ecoli (sufficient for the 5-strain decision table; stricter standalone ETECFinder DB is a future v0-spec decision) | `research_outputs/dna_decoder_phase4_gate_a_pass_2026-05-28.md` | high | 2026-05-29 |
<!-- project-state:end:evidence -->

### Unknowns
- ~~Whether EnteroBase has per-pathotype labels at the resolution needed~~ **RESOLVED 2026-05-27 via substrate survey:** EnteroBase pathotype labels are derived via BlastFrost gene-presence calls → structurally CIRCULAR with v0 deterministic resolver. EnteroBase remains usable for clonal-lineage + cohort-assembly only, NOT as pathotype-label source.
- Whether VirulenceFinder's gene-call output is stable enough across releases to diff reliably (DB versioning discipline) — partially mitigated by H5 + provenance schema; full validation deferred to v0 implementation.
- Hybrid pathotype prevalence in available labeled substrate (would there be ≥5 for evaluation?) — partial answer from substrate survey: hybrid pathotypes (STEC/EAEC, UPEC/EAEC) ARE documented in the literature (e.g., O104:H4; ST131 UPEC/EAEC hybrid in Brazilian Amazon study); curated panel feasibility requires direct contact with Whittam STEC Center.
- Whether the chosen pathotype-decision layer's failure modes can be detected by the EP-1 audit-gate pattern (mechanism-class-bounded SUSPEND verdict).
- DAEC and AIEC inclusion: out-of-scope v0 per T2 resolution — emit UNCLASSIFIED.
- **NEW unknown (post-substrate-survey):** Horesh 2021 supplementary File F1 per-record label-provenance distribution — specifically, what fraction of the 10,146 records have isolation-source-based pathotype labels vs gene-rule-only labels. Required to estimate H1 floor on the Horesh substrate.
- **NEW unknown:** Whittam STEC Center DECA per-pathotype strain counts (EPEC / EHEC / ETEC / EAEC / EIEC breakdown of the 15 clones × N strains each). Direct-contact follow-up via stec@cvm.msu.edu required (MTA + non-export restrictions documented).
- **NEW unknown:** NCBI Pathogen Detection `host_disease` field population density for E. coli — whether it's dense enough to serve as independent-label proxy for ExPEC/UPEC + commensal. Requires SOLR-index facet query on NCBI Pathogen Detection isolate browser.

### Hypotheses (Active) — REWRITTEN 2026-05-26 post-/brainstorm round 2 (sharper falsifiers per Codex review)
| ID | Statement | Pass/fail threshold | Substrate requirement | Confound controlled | Status | Last-tested |
|---|---|---|---|---|---|---|
| H1 | Label substrate is NOT mostly circular marker-derived labels | PASS if ≥70% of evaluation records have label provenance independent of v0 marker rules (clinical dx / outbreak investigation / lab pathotype assay / curated literature). FAIL if provenance missing for >40% OR labels explicitly gene-rule-derived | Per-record provenance field: label source, label method, curator, whether marker genes were used in label derivation | Prevents "rules validate against labels made by the same rules" — Phase 1 EP-1 SUSPEND-gate analog | **under-investigation (substrate-tier strategy locked 2026-05-27)** | 2026-05-27 |
| H2 | Deterministic resolver has acceptable per-class performance on independent labels | Ship gate: N≥50 each for EHEC/STEC-LEE, EPEC, ETEC, EAEC, UPEC-compatible/ExPEC + N≥75 commensal/low-marker. Precision floor: ≥0.90 EHEC/STEC + ETEC; ≥0.85 EPEC + EAEC; ≥0.80 UPEC-compatible; ≥0.85 commensal-low-marker. Recall floor: ≥0.80 non-UPEC; ≥0.70 UPEC + commensal-low-marker. Smoke gate: N≥20 acceptable but cannot ship the claim | Independent labels (H1-passing) + raw FASTA assemblies passing QC | Avoids macro-F1 hiding weak EAEC/UPEC/commensal performance | open | 2026-05-26 |
| H3 | Performance survives clade-balanced validation | PASS if clade-out precision floors from H2 are met, OR degrade by ≤10 pp vs random stratified split while retaining ≥0.80 precision on confident calls. FAIL if one dominant Mash/MLST/serotype cluster explains class success | Mash distances OR MLST for all isolates; folds built so close genomic neighbors do not cross train/test; ≥5 clade folds; no fold missing >2 target classes | Lineage leakage + O157:H7-style serotype shortcutting (reuses Mash-cluster N=147 infra) | open | 2026-05-26 |
| H4 | Abstention improves reliability rather than hides failures | PASS if confident-call precision ≥0.95 overall + abstention recall ≥0.80 on known hybrid/ambiguous/low-QC cases + abstention rate ≤15% on curated unambiguous. FAIL if abstention >25% overall without enrichment for known ambiguity | ≥30 known ambiguous/hybrid/out-of-scope isolates + ≥200 curated unambiguous | Prevents AMBIGUOUS from becoming an unmeasured escape hatch | open | 2026-05-26 |
| H5 | VirulenceFinder/ETECFinder → resolver concordance is stable under pinned versions | PASS if parser reproduces 100% of frozen VF/ETECFinder marker calls from fixture outputs; derived pathotype rules agree with synthetic marker-profile fixtures at 100%; on real unambiguous cases, rule explanation cites all decisive VF/ETECFinder markers (no silent dropped decisive markers) | Frozen VF/ETECFinder outputs from ≥50 isolates spanning every decision-table row + synthetic fixtures for edge cases | Parser drift, DB-version drift, circular "agreement" claims | open | 2026-05-26 |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| Project framed as EP-4-candidate pathotype prediction (per the long-horizon roadmap) | 2026-05-26 | Continuation of "DNA Decode AI Project" session; ledger started for separated context |
| External CLI surface = FASTA-in (+ optional GFF3) → pathotype label + provenance + VirulenceFinder diff | 2026-05-26 | From originating goal; not yet decomposed to internal architecture |
| Substrate-first ordering: architecture-fork decision deferred until label-provenance audit | 2026-05-26 | Post-/brainstorm round 2; H1 is load-bearing — if labels are gene-rule-derived, classifier/rule-comparison is circular |
| v0 default architecture = Candidate 5 (deterministic multilabel cluster resolver + abstention); Candidate 2+3 demoted to v0.1+; Candidates 1+4 subsumed | 2026-05-26 | Per /brainstorm round 2; pending final lock at Short-term Action 4 conditional on substrate-survey findings |
| v0 ship target = 3 months (down from 12-month upper bound) | 2026-05-26 | No trained classifier in v0; deterministic resolver layer only |
| Provenance schema MUST capture VF software version + DB version + DB commit/checksum + parameters | 2026-05-26 | DB drift between revisions is real (455 alleles added 2024 ETECFinder revision) |
| **T1 RESOLVED: VirulenceFinder independence framing = option C** (VF as both internal caller AND comparison reference; non-independence flagged via `caller_is_independent_baseline: false`; diff value lives in derived-rule layer) | 2026-05-26 | Cheapest honest v0; avoids second-ontology mapping problem from option A; avoids gene-caller-reimplementation cost from option B |
| **T2 RESOLVED: v0 output surface = 11-class honest default + `--legacy-6class` flag** | 2026-05-26 | EHEC_COMPATIBLE / STEC_NON_LEE / tEPEC_COMPATIBLE / aEPEC_COMPATIBLE / ETEC_COMPATIBLE / EAEC_COMPATIBLE / UPEC_COMPATIBLE / HYBRID / AMBIGUOUS / UNCLASSIFIED / COMMENSAL_LOW_MARKER_BURDEN. Originating goal's 6-class promise preserved via flag |
| **Substrate-tier strategy LOCKED post-/research substrate survey** | 2026-05-27 | Tier-1 (H1-passing bedrock, ~100-200 strains): Whittam DECA + von Mentzer 2021 ETEC reference + N=5 prototype strains (clonal-lineage + literature-curated labels INDEPENDENT of v0 marker rules). Tier-2 (large pool, ~5,000+ strains pending per-record audit): Horesh 2021 curated 10,146-genome collection (hybrid label method — needs File F1 per-record audit to estimate independent-label fraction). Tier-3 (expansion for ExPEC/UPEC + commensal): NCBI Pathogen Detection `host_disease`-populated subset. **Tier-4 / out-of-scope-for-labels:** EnteroBase (pathotype labels STRUCTURALLY CIRCULAR via BlastFrost; usable for clonal-lineage + cohort assembly only, NOT as pathotype-label source). See `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` for full evidence. |
| **Architecture-fork FINAL LOCK 2026-05-27: Candidate 5 (deterministic multilabel cluster resolver + abstention) is the v0 default** | 2026-05-27 | Substrate survey validates the v0 choice — classifier track (Candidate 2) requires H1-passing labels, and the survey shows only Tier-1 (~100-200 strains) and partial-Tier-2 (Horesh ~5,000 pending audit) carry plausible-independent labels. Candidate 5 sidesteps the H1 risk by using rules+abstention on pinned VF/ETECFinder calls. Candidates 2+3 remain v0.1+ falsifier track conditional on H1 PASS on at least one Tier-2-or-better substrate. |
| **EcoCyc/Whittam naming correction** | 2026-05-27 | Originating goal + early ledger referenced "EcoCyc DEC reference panel." Substrate survey clarified: EcoCyc is the E. coli K-12 metabolic database, NOT a pathotype reference panel. The DEC reference panel lives at the Whittam STEC Center / MSU (~16,000 strains, 4 distributed sets including DEC + ECOR + O157 STEC + non-O157 STEC). Affects substrate-tier-1 sourcing language only; the originating-goal verbatim input is preserved unchanged in `## Project Context`. |
| **EP-1 SUSPEND-gate reuse feasibility verdict = ADAPTED_REUSE** | 2026-05-27 | Per `research_outputs/ep1-suspend-gate-pathotype-reuse-feasibility-2026-05-27.md`. Falsifier NOT triggered (7/7 noise-class symbols have direct pathotype analogs; threshold-shape reusable with retuning). Workhorse-side reuse: copy `scripts/drug_mechanism_phenotype_merge.py` as template → `scripts/pathotype_cluster_label_merge.py`; design `dna_decode/data/pathotype_tiers.py` analog to `mic_tiers.py` (workhorse-owned catalog enumeration); start with 0.70/0.40 verdict thresholds + 20/10/5 gate thresholds; retune on real cohort. Estimated 4-6 hr engineering post-Gate-A. Closes Action 5 in Short-term actions. |
| **Horesh build = BOUNDED VERTICAL SLICE, not full cohort** (two-machine disagreement adjudicated 2026-05-30 via /brainstorm) | 2026-05-30 | Per `research_outputs/pathotype_horesh_bounded_slice_decision_2026-05-30.md`. Laptop "off-critical-path" stance was too strong (ingestion plumbing IS mandatory); workhorse "build now" right on plumbing but overclaims if a 3-class result = "validated predictor" (resolver-vs-marker-derived-labels = rule fidelity, not skill). Decision: (1) re-pick smoke+first-cohort to direct-WGS-accession rows (Salipante-style; skip Kallonen lane-ID re-assembly); set exit criteria (accession yield / caller completion / independent-label count) before scaling. (2) Provenance-split eval: report resolver-conformance vs gene-rule labels AND external-validity vs independent labels SEPARATELY — never blend. (3) Reframe v0 as "auditable marker-based pathotype-COMPATIBILITY resolver with abstention," NOT "predictor." Supported classes = ExPEC/EPEC/ETEC; EAEC/commensal/clean-EHEC = documented scope-limit (abstain). Full materialize gated on slice-exit-criteria PASS + Gate B >=2 yes. |
| **v0 claim reframed: compatibility resolver, not predictor** | 2026-05-30 | Clean ExPEC labels are isolation-source-derived (blood/urine) → genuinely independent of markers → ExPEC subset IS a legitimate genotype→phenotype test; circularity bites only the dropped `(predicted)` 52.2%. v0 supported surface narrowed to the 3 H1-passing classes; honest-scope-limit discipline (abstention) carries the rest. |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| ~~T1 — VirulenceFinder independence framing~~ | claude (brainstorm 2026-05-26) | RESOLVED 2026-05-26 | **Resolved as option C** — see `## Decisions Made`. |
| ~~T2 — Output surface: 6-class vs 11-class~~ | claude (brainstorm 2026-05-26) | RESOLVED 2026-05-26 | **Resolved as 11-class honest default + `--legacy-6class` flag** — see `## Decisions Made`. |
| ~~Substrate source: EnteroBase / ECOR curated / NCBI Pathogen Detection / GenomeTrakr / EcoCyc DEC reference panel / hybrid~~ | claude | RESOLVED 2026-05-27 | **Resolved as 4-tier strategy** — see `## Decisions Made` ("Substrate-tier strategy LOCKED post-/research substrate survey"). |
| Hybrid + atypical pathotype handling: separate output class, multilabel flag, or SUSPEND-style audit verdict | claude | RESOLVED 2026-05-26 | HYBRID is a primary output class in the 11-class surface (when ≥2 DEC modules confidently present); SUSPEND-style audit verdict reserved for label-provenance-uncertain rows. |
| ~~DAEC and AIEC inclusion (out-of-scope v0?)~~ | user | RESOLVED 2026-05-26 (folded into T2) | Explicit out-of-scope for v0; emit UNCLASSIFIED (NOT commensal) when EIEC/DAEC/AIEC-like markers detected. |
| **NEW: Horesh 2021 File F1 supplementary audit** (per-record label-provenance distribution) | claude | requires download + parse of `https://pmc.ncbi.nlm.nih.gov/articles/PMC8208696/` File F1 supplementary | If isolation-source-only records ≥ 70%, Horesh becomes Tier-2 substrate. If gene-rule-only > 30%, Horesh demotes to Tier-3 (annotation source only). |
| **NEW: Whittam STEC Center direct contact** (DECA per-pathotype counts + MTA terms) | user | requires email to stec@cvm.msu.edu | Per-clone strain counts + per-pathotype breakdown of DECA collection not publicly documented. Confirms Tier-1 cohort feasibility at N≥50 per pathotype. |
| **NEW: NCBI Pathogen Detection host_disease facet query** (E. coli isolate-browser SOLR index) | claude (HONEST-GAP returned 2026-05-27) | programmatic facet access blocked (SOLR exposed only via JS-driven SPA shell; BigQuery requires GCP auth) | Per `research_outputs/ncbi-pathogen-detection-host-disease-facet-2026-05-27.md`. Pre-committed verdict bar: PASS ≥75 / PARTIAL 30-74 / FAIL <30 / HONEST-GAP if blocked. **Manual query paths documented in memo: (a) user-driven Isolates Browser UI facet retrieval (~10-30 min), (b) one-time GCP/BigQuery setup (~30 min + 1 SQL), (c) email pd-help@ncbi.nlm.nih.gov.** Decision NOT retired; awaiting one of the three paths. |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Ledger v5 (2026-05-30): Gate A PASSED (workhorse 2026-05-28). Horesh F1 audited (H1 FAIL on full 10,146; clean subset 2,077 = ExPEC/EPEC/ETEC only; EAEC/commensal/clean-EHEC HARD FAIL). Two-machine Horesh-build disagreement **adjudicated 2026-05-30**: **BOUNDED VERTICAL SLICE** (direct-WGS-accession rows, skip Kallonen re-assembly, exit-criteria before scaling) + **provenance-split eval** (resolver-conformance vs external-validity, never blended) + **v0 reframed = "compatibility resolver with abstention," NOT "predictor."** Supported v0 classes = ExPEC/EPEC/ETEC; rest = documented scope-limit. Decision memo: `research_outputs/pathotype_horesh_bounded_slice_decision_2026-05-30.md`. **SOLE binding next action = SEND Gate B outreach (still unsent, 0 replies);** kit at `research_outputs/pathotype_gate_b_send_kit_2026-05-29.md`. Full Horesh materialize gated on slice-exit-criteria PASS + Gate B ≥2 yes.

### Target state / terminal condition
Open-source CLI tool ships v0 release with: (1) FASTA-in → multilabel cluster profile + derived pathotype label + provenance + VirulenceFinder diff (11-class honest surface with optional `--legacy-6class` flag); (2) H2 + H3 + H4 ship-gate floors met on H1-passing substrate; (3) HYBRID + AMBIGUOUS + UNCLASSIFIED handled as explicit output classes (not buried under primary label); (4) abstention metrics reported alongside confident-call metrics; (5) pinned VF 3.2.0 + ETECFinder DB versions in every provenance record. v0.1+ extension target: classifier sidecar + disagreement-detector layer, gated on H1 PASS.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` count + `gates-passed` count (raw counts, unweighted). Current (2026-05-27): unknowns=9 (6 original + 3 new from substrate-survey) / 1 retired (EnteroBase resolution); hypotheses=5 / 0 tested but H1 now `under-investigation` with substrate-tier evidence; pending-decisions=8 / 5 resolved (T1+T2+substrate-source+hybrid-handling+DAEC-AIEC closed); Short-term actions = 9 total / 3 completed (substrate-survey + T1+T2 + architecture-fork-lock).
- **v0.2+:** weighted combination of unknowns-retired, gates-passed, evidence-confidence-improved, hypotheses-falsified (TBD via v0.2 design).

### Candidate next actions — REORDERED 2026-05-27 post-/research substrate-survey
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | **Horesh 2021 File F1 supplementary audit** (Short-term Action 1-new) | research | high (gates H1 status on Tier-2 substrate; determines whether Horesh is usable for the ~5,000-strain expansion or only Tier-3) | high (resolves U6-new + H1 status on the largest available substrate; per-pathotype × per-ST distribution feeds H3 fold construction) | low-medium (File F1 is publicly available at PMC8208696 supplementary; schema unknown until parsed) | hours |
| 2 | **Whittam STEC Center direct-contact follow-up** (Short-term Action 2-new) | ask-user | high (confirms Tier-1 cohort feasibility) | high (resolves DECA per-pathotype counts + MTA terms + non-export restrictions) | medium (depends on Whittam Center response time + access terms) | days-to-weeks (user-driven) |
| 3 | Confirm VirulenceFinder integration pathway + pin VF 3.2.0 + DB checksum (Short-term Action 3) | research | medium (Mid-term#3 prerequisite; H5 prerequisite) | medium (DB versioning + invocation contract clarified) | low (well-documented tool) | hours |
| 4 | **NCBI Pathogen Detection host_disease facet query** (Short-term Action 4-new) | research | medium (Tier-3 substrate density confirmation) | medium (resolves U8-new; supports ExPEC/UPEC + commensal cohort expansion) | medium (SOLR API / facet behavior on NCBI Pathogen Detection isolate browser) | hours |
| 5 | Cross-check audit-gate reuse + pathotype-specific opacity definition (Short-term Action 5) | research | medium (Mid-term#3 prerequisite; H5 prerequisite) | medium (resolves U4 + H5 partial) | low (EP-1 pattern documented in `scripts/cipro_mechanism_phenotype_merge.py`) | hours |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- **Default:** re-run `/project-state` after any action class fires (auto-append to Action Log triggers stale-state check).
- **Manual override:** user invokes `/project-state <slug>` at any time.
- **v0.2+:** automated trigger when N actions fire OR T days elapse OR a hypothesis falsifies (TBD).

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` — auto
- `research` — auto (delegate to /research / /athena-research / /research-verify)
- `write-plan` — auto (delegate to /technical-plan / /save-plan)
- `edit-local-code` — REQUIRES per-action human approval
- `run-tests` — auto if local + sandboxed
- `ask-user` — auto
- `stop` — auto

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-05-26 | propose | /idea-anchor invoked on pathotype CLI goal (continuation of "DNA Decode AI Project" session) | Formal Rephrase + 3 clarifications + 8 assumptions + blunt opinion produced; rule-vs-classifier fork surfaced as foundational; recommended /project-init next |
| 2 | 2026-05-26 | propose | /project-init invoked | ledger created; goal-normalization PASS/PASS/PASS; 4 architecture candidates surfaced; 5 short-term actions enumerated |
| 3 | 2026-05-26 | research | /brainstorm round 1+2 (Codex adversarial + generative-ideation v2.1 + sharpening) | 3 substantive issues (all accepted): action sequencing inverted; v0 architecture wrong; H1-H5 unsharp. Concrete v0 output contract delivered (≥23 marker clusters, 11-row decision table, provenance JSON schema). 5 rewritten falsifier hypotheses delivered. T1 + T2 trade-offs surfaced |
| 4 | 2026-05-26 | write-plan | Ledger updated to absorb /brainstorm round 2 | Refinement Candidates rewritten (Candidate 5 added; 1+4 subsumed; 2+3 demoted to v0.1+); Short-term actions reordered (substrate-first); H1-H5 rewritten; T1+T2 added to Pending Decisions; Decisions Made updated; Candidate next actions reordered; v0 Output Contract section appended |
| 5 | 2026-05-26 | ask-user | T1 + T2 resolved via AskUserQuestion | T1 = option C (VF as both + non-independence flag); T2 = 11-class honest default + `--legacy-6class` flag. Both recommendations accepted. Pending Decisions table updated; Open Questions Q1/T1/T2/Q5 closed; pending-decisions count 5 → 3 |
| 6 | 2026-05-27 | research | /research substrate-survey run executed (Mission Control L1; run-id 2026-05-27-0140-research-ecoli-pathotype-substrate) | 22 candidate rows → 20 supported (91% survival) + 2 unsupported; 5 new active candidates to followup queue; load-bearing finding: H1 (label-independence) at structural risk on EnteroBase + NCBI Pathogen Detection + Horesh 2021 hybrid method; only Whittam DECA + von Mentzer ETEC reference + N=5 prototype strains pass H1 cleanly. Wall-clock soft-breach (~50min vs 30min cap) driven by external retries. |
| 7 | 2026-05-27 | write-plan | Ledger updated to absorb /research substrate-survey findings | Evidence table +8 rows (4-11); Unknowns updated (closed EnteroBase resolution; added 3 new from survey); H1 status open → under-investigation; Decisions Made +3 (substrate-tier strategy LOCKED, architecture-fork FINAL LOCK Candidate 5, EcoCyc/Whittam naming correction); Pending Decisions table refreshed (substrate-source resolved; +3 new follow-up actions); Short-term actions reordered (3 completed, 5 open); Mid-term M1+M2 closed |
| 8 | 2026-05-27 | research | Discovery-machine closeout memos (Short-term Actions 4-new + 5) | Step 1 (NCBI host_disease facet) → HONEST-GAP; 3 manual query paths documented at `research_outputs/ncbi-pathogen-detection-host-disease-facet-2026-05-27.md`. Step 2 (EP-1 SUSPEND-gate reuse feasibility) → ADAPTED_REUSE; 7/7 noise-class symbols have direct pathotype analogs; workhorse-side reuse recommendation provided at `research_outputs/ep1-suspend-gate-pathotype-reuse-feasibility-2026-05-27.md`. Both memos committed at `272e90d`. |
| 9 | 2026-05-27 | write-plan | Ledger updated to absorb closeout-memo verdicts | Decisions Made +1 (EP-1 SUSPEND-gate ADAPTED_REUSE verdict); Pending Decisions: NCBI Pathogen Detection row NOT retired, HONEST-GAP note + 3 manual paths appended; Short-term Actions: #5 COMPLETED, #4-new PARTIAL/HONEST-GAP. NO Hypotheses table mid-row touch (H1/H5 status preserved). |
| 10 | 2026-05-29 | run-tests | Phase 4 pathotype Gate A executed on workhorse (Precision 7780); result bundle received on this laptop | PASS — 5/5 frozen panel correct; caller_runtime_ok / raw_outputs_usable / decision_table_computable all true. Gate A is the local-sanity pre-build gate; Gate B (adoption signal) is now the SOLE remaining pre-build gate. See Evidence rows 12-13. |
| 11 | 2026-05-29 | propose | Integrated workhorse Gate B bundle (10 files) + consolidated SEND KIT on this laptop | Full bundle local at reports/pathotype_gate_b/ (reports/ gitignored); 7 durable artifacts tracked in research_outputs/ (commit c9fba14). Canonical send kit = research_outputs/pathotype_gate_b_send_kit_2026-05-29.md (named 5-contact shortlist Seemann/Hazen/PHAC-NML/nf-core/PulseNet + post-Gate-A one-pager + outreach draft + rubric + CSV). Gate B PASS = 2+ credible yes-to-pilot within 60 days; binding action user-side. HOLD Horesh build until gate clears. |
| 12 | 2026-05-29 | research | Returned Horesh F1 metadata CSV + locator analysis to workhorse (2 bundles) | F1 CSV already on laptop (md5 1e76642e..., 10,146 rows) → returned + F1 audit. Locator finding: NO bulk Horesh assembly deposit; assemblies scattered per source-study. 1/5 smoke rows resolved (JSIS00000000 = E. coli upec-276, Salipante WGS); 4/5 are Kallonen Sanger lane IDs needing re-assembly. |
| 13 | 2026-05-30 | propose | Adjudicated Horesh-build two-machine disagreement via /brainstorm | Verdict: BOUNDED VERTICAL SLICE + provenance-split eval + v0 reframed as compatibility-resolver (not predictor). Decision memo at research_outputs/pathotype_horesh_bounded_slice_decision_2026-05-30.md. Binding next action remains user-side: SEND Gate B outreach (still unsent, 0 replies). Full materialize gated on slice-exit + Gate B >=2 yes. |
| 14 | 2026-05-30 | research | [Soraya advance run 2026-05-30-1200-ep4-pathotype] Resolved bounded-slice WGS-accession candidates (solo, laptop) | 260 Horesh clean rows carry direct WGS accessions (ExPEC 135 + EPEC 125 + ETEC 0). research_outputs/horesh_bounded_slice_wgs_accession_candidates_2026-05-30.csv. ETEC=0 → all 181 clean ETEC are Sanger lane IDs (need re-assembly). |
| 15 | 2026-05-30 | research | [Soraya advance] Resolved ETEC accession-bearing alternative = von Mentzer 2021 references | 8 strains / 7 lineages, BioProject PRJEB33365. research_outputs/etec_reference_vonmentzer_2026-05-30.{md,csv}. Slice can be 3-class by-accession after all. NUANCE: ETEC labels are toxin-typed = resolver markers → ETEC = near-conformance, NOT strong external-validity; ExPEC (isolation-site) remains cleanest arm. Per-strain GCA enumeration is the one remaining lookup. |
<!-- project-state:end:action-log -->

## Open Questions for User — UPDATED 2026-05-27 post-/research substrate-survey
- ~~**Q1 (architecture fork):**~~ **RESOLVED 2026-05-27** — Candidate 5 FINAL LOCK (deterministic multilabel cluster resolver + abstention).
- ~~**T1 (VirulenceFinder independence framing)**~~ **RESOLVED 2026-05-26:** option C (VF as both + `caller_is_independent_baseline: false`).
- ~~**T2 (6-class vs 11-class)**~~ **RESOLVED 2026-05-26:** 11-class honest default + `--legacy-6class` flag.
- ~~**Q2 (scope cap — folded into T2)**~~ **RESOLVED 2026-05-26:** DAEC + AIEC out-of-scope for v0; emit UNCLASSIFIED.
- **Q4 (audit gate):** Reuse EP-1 SUSPEND-gate pattern for pathotype-call opacity, or design a new audit verdict shape? Recommended: reuse EP-1 with pathotype-specific opacity (aggR-alone, partial LEE, contradictory toxin/regulator hits). Decided at Short-term Action 5.
- ~~**Q5 (substrate source)**~~ **RESOLVED 2026-05-27:** 4-tier strategy LOCKED per /research substrate-survey. Tier-1 DECA + von Mentzer + prototypes; Tier-2 Horesh 2021 pending File F1 audit; Tier-3 NCBI Pathogen Detection host_disease; Tier-4 EnteroBase (clonal-lineage only, NOT pathotype labels).
- **Q6 (PASS-path overlap):** Does this project consume Mash-cluster-N=147 phylogenetic infrastructure (built on Precision 7780)? Recommended: yes for clade-balanced cohort construction at Mid-term#4 (required by H3).
- **Q7 (NEW — direct-contact pacing):** Whittam STEC Center direct-contact response time + MTA terms + non-export restrictions are unknown; do you want to fire this contact NOW (parallel with Horesh File F1 audit) or HOLD until Horesh audit confirms substrate strategy? Recommended: fire NOW (response time may be days-to-weeks; no blocking dependency).
- **Q8 (NEW — Horesh File F1 audit ownership):** the supplementary file at PMC8208696 contains per-record metadata; want claude to do the parse + audit (Short-term Action 1-new) or user-driven? Recommended: claude (mechanical metadata parse + count aggregation).

## v0 Output Contract (DRAFT — added 2026-05-26 post-/brainstorm round 2)

> Draft target for the deterministic multilabel cluster resolver (Candidate 5). Subject to revision after substrate survey (Short-term Action 1) and T1/T2 resolution. NOT a ship spec — this is the design anchor that the substrate survey + audit-gate work feeds into.

### Marker clusters to track (minimum ≥23 markers)

| Cluster | Marker genes / aliases | Use |
|---|---|---|
| `STX1` | `stx1`, `stx1A`, `stx1B` | STEC/EHEC |
| `STX2` | `stx2`, `stx2A`, `stx2B` | STEC/EHEC |
| `LEE` | `eae` + optional `tir`, `escV`, `sepL` | EPEC/EHEC attaching-effacing |
| `BFP_EAF` | `bfpA`, `bfpB`, `perA`, EAF plasmid markers | typical EPEC |
| `LT` | `eltA`, `eltB` (+ LT-I/LT-II aliases) | ETEC |
| `ST` | `estA`, `sta`, `stb` (STh/STp where typed) | ETEC |
| `EAEC_REG` | `aggR` | EAEC regulator |
| `EAEC_TRANSPORT` | `aatA` (+ `aatPABCD`) | EAEC pAA transport |
| `EAEC_T6SS` | `aaiC` (+ `aaiA-Y`) | EAEC support marker |
| `AAF_I` | `aggA` | EAEC fimbrial variant |
| `AAF_II` | `aafA` | EAEC fimbrial variant |
| `AAF_III` | `agg3A` | EAEC fimbrial variant |
| `AAF_IV` | `agg4A` | EAEC fimbrial variant |
| `AAF_V` | `agg5A` | EAEC fimbrial variant |
| `P_FIMBRIAE` | `papC`, `papG` | ExPEC/UPEC-compatible |
| `S_FIMBRIAE` | `sfa`, `sfaS`, `sfa/focDE` | ExPEC/UPEC-compatible |
| `AFA_DRA` | `afa`, `afaD`, `draBC` | DAEC/ExPEC ambiguity |
| `HEMOLYSIN` | `hlyA` (+ `hlyCABD`) | ExPEC/UPEC-compatible |
| `CNF1` | `cnf1` | ExPEC/UPEC-compatible |
| `SIDEROPHORES` | `iutA`, `iroN`, `fyuA` | ExPEC support (not primary) |
| `CAPSULE_SERUM` | `kpsMII`, `traT` | ExPEC support (not primary) |
| `EIEC_FLAG` | `ipaH` | Out-of-scope EIEC/Shigella flag |
| `DAEC_FLAG` | `daaD`, `daaE`, `afa/dra-only pattern` | Out-of-scope DAEC flag |

### Derived-pathotype decision table (11 output classes)

| Output call | Required cluster combination | Blocking / downgrade |
|---|---|---|
| `EHEC_COMPATIBLE` | `STX1 ∨ STX2` ∧ `LEE` | hybrid if EAEC/ETEC module also present |
| `STEC_NON_LEE` | `STX1 ∨ STX2`, no `LEE` | NOT EHEC |
| `tEPEC_COMPATIBLE` | `LEE` ∧ `BFP_EAF`, no `STX` | hybrid if LT/ST or EAEC also present |
| `aEPEC_COMPATIBLE` | `LEE`, no `STX`, no `BFP_EAF` | lower confidence than tEPEC |
| `ETEC_COMPATIBLE` | `LT` ∨ `ST` | hybrid if STX/LEE/EAEC also present |
| `EAEC_COMPATIBLE` | `aggR` + (`AAF_*` ∨ `aatA` ∨ `aaiC`), OR `AAF_*` + (`aatA` ∨ `aaiC`) | `aggR`-alone or `aatA/aaiC`-alone → `AMBIGUOUS_EAEC` |
| `UPEC_COMPATIBLE` | ≥2 strong ExPEC markers (`papC/G`, `sfa/foc`, `afa/dra`, `hlyA`, `cnf1`) + optional support | prefer `ExPEC/UPEC-compatible` — definitive UPEC needs urine/infection metadata |
| `HYBRID` | ≥2 primary DEC modules confidently present | report all modules; no forced primary unless user flag requests priority |
| `AMBIGUOUS` | partial markers / weak ID-cov / contradictory / single weak EAEC support | covers `aggR`-alone, fragmented `LEE`, partial toxin hits |
| `UNCLASSIFIED` | EIEC/DAEC/AIEC-like markers OR unmapped virulence markers | NOT equivalent to commensal |
| `COMMENSAL_LOW_MARKER_BURDEN` | no DEC module + sub-threshold ExPEC + assembly QC PASS | absence of v0 markers, not biological commensal truth |

### Abstention rules

- Emit `AMBIGUOUS` when: any primary marker hit is below locked thresholds; only one weak EAEC support marker present; `LEE` is partial; rule depends on a plasmid marker fragmented across contigs.
- Emit `HYBRID` whenever ≥2 primary DEC modules are confidently present. Do not bury hybrids under the highest-priority label.
- Emit `UNCLASSIFIED` for EIEC/DAEC/AIEC-like evidence (out of v0 six-class surface).
- Emit `COMMENSAL_LOW_MARKER_BURDEN` only if assembly QC passes AND no v0 marker module fires. Otherwise use `AMBIGUOUS_LOW_QC`.

### Provenance JSON schema (per record)

```json
{
  "sample_id": "string",
  "analysis_date": "YYYY-MM-DD",
  "assembly": {
    "path": "input.fasta",
    "n_contigs": 0,
    "n50": 0,
    "total_bp": 0,
    "qc_verdict": "PASS|WARN|FAIL"
  },
  "caller": {
    "name": "VirulenceFinder|ETECFinder|supplemental-blast",
    "software_version": "pinned (e.g., VF 3.2.0)",
    "db_name": "string",
    "db_version": "pinned",
    "db_commit_or_checksum": "string",
    "parameters": {"min_identity": 0.90, "min_coverage": 0.60},
    "caller_is_independent_baseline": false
  },
  "marker_hits": [
    {
      "cluster": "STX2",
      "gene": "stx2A",
      "allele": "string|null",
      "accession": "string|null",
      "source_db": "string",
      "percent_identity": 99.4,
      "percent_coverage": 97.8,
      "contig": "contig_12",
      "start": 12345,
      "end": 13001,
      "strand": "+",
      "hit_status": "CONFIDENT|PARTIAL|LOW_ID|LOW_COV"
    }
  ],
  "cluster_profile": {
    "STX1": false, "STX2": true, "LEE": true, "BFP_EAF": false,
    "LT": false, "ST": false, "EAEC": false, "EXPEC_SCORE": 1
  },
  "derived_call": {
    "primary": "EHEC_COMPATIBLE",
    "secondary": ["STEC_COMPATIBLE"],
    "confidence_tier": "CONFIDENT|AMBIGUOUS|HYBRID|UNCLASSIFIED",
    "rule_id": "RULE_EHEC_001",
    "rule_version": "pathotype-rules-v0.1.0",
    "reason": "STX2 + LEE detected; no ETEC or EAEC module detected"
  }
}
```

### Stress-test slices (not ship gates, but design checks the v0 must survive)

- **EAEC stress test:** plasmid-associated + variable markers will likely cause undercalling or frequent abstention. Acceptable only if H4's abstention metrics are explicit.
- **STEC-non-LEE handling:** must NOT be silently collapsed into EHEC.
- **Hybrid prevalence:** require ≥5 known hybrid (STEC/ETEC, STEC/EAEC like O104:H4) in test set.
- **Atypical EPEC (LEE+ bfp-):** must be called as `aEPEC_COMPATIBLE`, not silently downgraded to `UNCLASSIFIED`.
- **ExPEC/UPEC genome-only call:** confidence must reflect the lack of site/context evidence; downgrade to `ExPEC/UPEC-compatible` is mandatory absent metadata.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-05-29
- **Progress signal:** Gate A PASS (5/5 frozen panel, real VF runtime, workhorse 2026-05-28); Gate B send kit consolidated + ready on this laptop (5 named contacts; commit c9fba14); Gate B (PASS=2+ yes/60d) is the sole remaining pre-build gate; Horesh heavy build held pending gate.
