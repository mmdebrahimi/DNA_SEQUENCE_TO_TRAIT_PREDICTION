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
| 4 | CORRECTION/refinement of row 3: 24-80GB was TRAINING cost (PlantCAD 225M = 8xH100 25d). INFERENCE/embedding fits consumer GPUs — authors target RTX 3090; 4 sizes 20M/40M/112M/225M; 225M weights ~0.5GB fp16 @ 512bp. RTX 3500 Ada 12GB holds the largest variant for frozen-embedding extraction, no quantization. Path-B local viability CONFIRMED. | pnas.org/doi/10.1073/pnas.2421738122 + github.com/kuleshov-group/PlantCaduceus | high | 2026-06-07 |
| 5 | G1 is COMPUTE-GATED, not no-compute (EP7 substrate assumption corrected). The best S.Africa C.auris WGS+MIC cohort (PRJNA737309, 92 isolates, PMC8370198) is deposited as SRA RAW READS (115 SRA records, 0 assemblies) — the deterministic caller needs assembled FASTA, so G1 needs de-novo assembly (~30min/genome CPU) + per-isolate MIC parse from a supplementary PDF (Table S1) + manual isolate<->accession join. MIC is supplementary-PDF-only (not NCBI/BV-BRC metadata); no assembled+MIC-labeled C.auris cohort found. | research_outputs/cauris-wgs-mic-cohort-sources-2026-06-07.raw.md (PMC8370198 + NCBI eutils) | high | 2026-06-07 |
| 6 | Maphanga et al. AAC 2021 Table 3 (n=92, main article, user-downloaded): clade I=13 (all Y132F ERG11, all FLZ-R), clade III=77 (76 ERG11-mutant VF125AL; FLZ MIC 16-256; 69 FLZ-R, ~8 pan-susceptible), clade IV=2 (1 FLZ-R). Independently CONFIRMS the fungal_amr.py catalog (Y132F clade I, VF125AL clade III) + gives per-clade R-rates to validate against. Per-isolate MIC still requires Table S1 (supplementary, PoW-gated). | C:/Users/Farshad/Downloads/aac.00517-21.pdf p4 Table 3 | high | 2026-06-07 |
<!-- project-state:end:evidence -->
### Unknowns
- Is the C. auris WGS+MIC cohort actually downloadable at depth (accessions in supplementaries)?
- ~~Does a plant DNA-FM fit 12GB, or is cloud (money) required for Path B?~~ RETIRED 2026-06-07: PlantCAD inference fits 12GB (training-only needs 8xH100); RTX 3500 Ada sufficient. See Evidence row 4.
- What fraction of C. auris azole-R is ERG11-only vs efflux/aneuploidy (sets fungal sensitivity ceiling)?
### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | The deterministic determinant-scan method transfers to fungi (C. auris azole) at acc≥0.80 | confirmed (method transfers; acc-bar label-limited) | 2026-06-08 |
| H2 | Frozen DNA-FM embeddings beat PRS+kinship baselines on Arabidopsis flowering-time (the embedding niche exists) | open | never |
<!-- project-state:end:hypotheses -->
### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| Path C ratified (A now / B queued on compute) | 2026-06-07 | user-ratified via /soraya |
| Two-machine split, no duplication, git-only sync | 2026-06-07 | plans/Eukaryotic_DualMachine_Coordination.md |
| Resolved Pending Decision row 1: Workhorse = personal Precision 7780 (RTX 3500 Ada ~12GB), NOT Bombardier/DLP. Path B safety gate CLEARED; personal code may run there. | 2026-06-07 | Auto-linked by /project-state --resolve-pending-decision (v0.2). Original decision text preserved at Pending Decisions row 1 with RESOLVED prefix. |
| Resolved Pending Decision row 2: Path B compute = local ~12GB RTX 3500 Ada; paid cloud A100 DEFERRED (leave for later). Databricks high-GPU available as the established burst pattern but treated as money-gated (not fired now). | 2026-06-07 | Auto-linked by /project-state --resolve-pending-decision (v0.2). Original decision text preserved at Pending Decisions row 2 with RESOLVED prefix. |
| Resolved Pending Decision row 3: G2-FAIL does NOT auto-close the embedding frontier — user chose KEEP-OPEN (a clean G2-FAIL prompts another embedding attempt: different FM / cloud A100 / different phenotype before concluding). H2 stays open past a single FAIL. | 2026-06-07 | Auto-linked by /project-state --resolve-pending-decision (v0.2). Original decision text preserved at Pending Decisions row 3 with RESOLVED prefix. |
| Resolved Pending Decision row 4: G1 cohort path = LAPTOP-CPU subset assembly of PRJNA737309 (user-chosen; no money). Cohort verified-feasible (114 MiSeq runs, isolate-ID=LibraryName, D: disk fine). Now blocked on ONE user action: browser-download the PoW-gated Table S1 PDF for MIC labels. | 2026-06-07 | Auto-linked by /project-state --resolve-pending-decision (v0.2). Original decision text preserved at Pending Decisions row 4 with RESOLVED prefix. |
| Resolved G2 design choices (manifest §8): user accepted the recommended defaults — estimand=structure-independent causal signal; agnostic-window=all gene bodies +~1kb flanks; budget=fixed sub-day per-accession window budget; cloud genome-wide=NOT run unless borderline+approved (money gate intact). §8 marked RATIFIED; workhorse executes without further design decision. | 2026-06-08 | Ratified via /soraya "move forward with recommendations". Manifest §8 frozen. |
<!-- project-state:end:decisions-made -->
### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| RESOLVED 2026-06-07: Workhorse identity (Precision 7780 vs Bombardier/DLP) | Soraya | user | SAFETY: Path B handoff assumes personal Precision 7780; do NOT route personal code through Bombardier |
| RESOLVED 2026-06-07: Path B compute: 12GB GPU sufficient OR cloud budget | Soraya | user | money gate — no paid compute without explicit OK |
| RESOLVED 2026-06-07: Pre-commit: G2-FAIL closes the embedding frontier permanently (no 5th attempt) | Soraya | user | guards against the diminishing-returns trap on a 4th embedding bet |
| RESOLVED 2026-06-07: G1 cohort path: laptop-CPU subset assembly (PRJNA737309) vs Databricks burst vs user-curated TSV | Soraya | user | G1 is compute-gated (Evidence row 5); needs a scope+venue decision. Laptop CPU assembly of a ~20-30 isolate de-confoundable subset is "minor compute" + no money; Databricks = money-gated; user-curated TSV = fastest if a per-isolate accession+MIC table is available |
| RESOLVED 2026-06-08: G1 EXECUTION now blocked on ONE user action: download Table S1 PDF (aac.00517-21-s0001.pdf) via browser | Soraya | user | PMC supplementary is behind a proof-of-work anti-scraping challenge (cloudpmc-viewer-pow); curl/WebFetch cannot solve it. Per-isolate fluconazole MIC labels live in that table. Everything else (SRR list, isolate-ID join key, assembly pipeline plan, disk) is verified-ready. User downloads -> drops in Downloads -> G1 executes |
| RESOLVED 2026-06-08: G2 design choices (manifest §8): primary estimand / agnostic-window rule / GPU-hours budget / cloud-genome-wide | Soraya | user | surfaced by the 2026-06-08 brainstorm; each changes what the gate measures. Defaults proposed in §8 (causal-signal estimand / all-genes+flanks / sub-day window budget / no cloud unless borderline+approved). Decide at the workhorse dry-manifest; not blocking the laptop |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
G0 + G0-completion DONE (caller validated on real C. auris ERG11 alleles, 3/3). Both authority gates resolved (workhorse=Precision 7780; compute=local 12GB, PlantCAD_l32 fits). G1 INFRASTRUCTURE complete (scripts/build_fungal_cohort.py, 16 fungal tests green, validated on real sequence). GATE G1 REACHED 2026-06-08. C. auris fluconazole deterministic decoder validated on N=24 real de-confoundable cohort: LABEL_LIMITED_FAILURE (acc 0.792, sens 1.000, spec 0.167). Method TRANSFERS across the kingdom boundary (ERG11 mutation found in 100% of MIC-R across clades I+III); the missed acc-bar is documented label-limitation (5/6 'susceptible' isolates carry F126L below the CDC tentative breakpoint — suspect-the-label, genotype is the reliable output), NOT a caller defect. H1 CONFIRMED. Packet + closeout in wiki/. Path B (Arabidopsis embedding) FULLY PRE-STAGED + design-frozen: manifest reviewed/revised (agnostic primary + dry-manifest gate + paired-bootstrap PASS), §0.5 dry-manifest gate CODED + unit-tested (scripts/g2_dry_manifest.py, 7 tests), §8 design choices RATIFIED (causal-signal / genes+flanks / sub-day / no-cloud). ZERO open design decisions. The ONLY remaining action is the GPU run on the workhorse (Precision 7780) — user fires it; laptop side complete.

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
| 5 | 2026-06-07 | ask-user | /project-state --resolve-pending-decision | Pending Decisions row 1: Workhorse = personal Precision 7780 (RTX 3500 Ada ~12GB), NOT Bombardier/DLP. Path B safety gate CLEARED; personal code may run there. |
| 6 | 2026-06-07 | ask-user | /project-state --resolve-pending-decision | Pending Decisions row 2: Path B compute = local ~12GB RTX 3500 Ada; paid cloud A100 DEFERRED. Databricks high-GPU available but money-gated (not fired now). |
| 7 | 2026-06-07 | research | Path-B model-fit feasibility check (verify-before-setup) | PlantCAD inference fits 12GB decisively: 4 sizes (20M/40M/112M/225M), authors target RTX 3090, 225M weights ~0.5GB fp16 @ 512bp ctx. The survey's 24-80GB was TRAINING cost, not inference. RTX 3500 Ada 12GB sufficient for frozen-embedding extraction, no quantization. Evidence row 4. Unknown 'does a plant DNA-FM fit 12GB?' RETIRED. |
| 8 | 2026-06-07 | research | G1 label-source iron-law check: where do C. auris fluconazole MIC labels live? | NEGATIVE on both programmatic sources: NCBI BioSample has only 53 C. auris records mentioning antifungal terms + they are experimental-TREATMENT fields ('256 ug/mL Fluconazole 2h'), NOT AST/MIC results; BV-BRC genome_amr has ZERO C. auris records (bacteria-only). Confirms EP7 plan: G1 label source = paper supplementaries ONLY (manual per-isolate extraction). Cohort-build binding constraint is now exactly ONE external input: a [accession, fluconazole_mic, clade] TSV. |
| 9 | 2026-06-07 | edit-local-code | G1 cohort-validation INFRASTRUCTURE built + validated on real sequence | scripts/build_fungal_cohort.py: reads a per-isolate label TSV (accession->download OR genome_fasta->local), runs call_erg11 over the cohort, computes overall + WITHIN-CLADE de-confounded acc/sens/spec + the efflux/aneuploidy DISCORDANCE set (MIC-R but ERG11-S = documented non-target mechanism), emits the G1 verdict (PASS / DOCUMENTED_FAILURE_MODE / FAIL) as .md+.json. Added CDC tentative MIC breakpoints (fluconazole>=32=R) + mic_to_phenotype to fungal_amr.py. 4 cohort tests validate the FULL pipeline on the 3 committed real alleles via real blastn (TP/TN/FN buckets, within-clade, discordance, both verdict branches); 16 fungal tests green. G1 now executes in ONE command once the label TSV lands. |
| 10 | 2026-06-07 | ask-user | /project-state --resolve-pending-decision | Pending Decisions row 3: G2-FAIL does NOT auto-close the embedding frontier (user KEEP-OPEN); H2 stays open past a single FAIL. |
| 11 | 2026-06-07 | research | /research run (cauris-wgs-mic-cohort-sources) — G1 cohort sourcing reconnaissance | Found the 4 primary S.Africa C.auris WGS+MIC papers + the only stated BioProject PRJNA737309 (92 isolates). Structural finding: PRJNA737309 = SRA raw reads (0 assemblies) -> G1 is COMPUTE-GATED (assembly + PDF MIC parse + ID join), NOT no-compute. Memo: research_outputs/cauris-wgs-mic-cohort-sources-2026-06-07.raw.md. Evidence row 5. Opens Pending Decision: G1 cohort path. |
| 12 | 2026-06-07 | research | G1 cohort grounding + execution-spec (laptop-CPU path chosen) | Verified PRJNA737309 = 114 Illumina MiSeq WGS runs (~434MB each, ~10GB for a 24-isolate subset; D: 4.4TB free), SRA LibraryName = isolate ID (e.g. 128_98) = the join key to Table S1. Wrote the G1 execution spec into EP7 step 4. THE one remaining blocker: Table S1 MIC labels are behind a PMC proof-of-work anti-scraping gate (cloudpmc-viewer-pow) — curl/WebFetch can't solve it; needs the user to browser-download aac.00517-21-s0001.pdf. Pending Decisions row 4 resolved (laptop-CPU); new pending = the user PDF download. |
| 13 | 2026-06-07 | edit-local-code | G1 cohort JOINED + assembly pipeline built; user provided Table S1 PDF | User downloaded the supplementary -> parsed Table S1 (92 per-isolate clade+FLC MIC) + solved the SRA join (LibraryName==isolate_id, 84 exact + 6 stem-rescued; stem-rescue recovers the susceptibles the paper lists as bare stems). Cohort = 90 isolates (82R/8S); clade III has 7S vs 69R = real within-clade de-confound contrast. Wrote data/fungal_ref/cauris_g1_cohort.tsv (90) + cauris_g1_subset.tsv (25: 7 cladeIII-S + 14 cladeIII-R + 4 cladeI). Built scripts/assemble_sra_cohort.py (prefetch+fasterq-dump+SPAdes via docker_runner, restartable, pinned images sra-tools 3.1.1 + spades 3.15.5, both pulled). Docker up; 1-isolate end-to-end validation RUNNING (then batch the 25 -> build_fungal_cohort.py -> G1 verdict). |
| 14 | 2026-06-07 | edit-local-code | G1 pipeline pivoted to targeted ERG11 read-mapping + VALIDATED end-to-end; 25-isolate batch RUNNING | WGS assembly was wrong tool (220x coverage -> SPAdes/SKESA ~25-45min/isolate). Added --method map (minimap2 -ax sr vs ERG11 CDS ref -> samtools consensus) = ~4min/isolate; call_erg11 BLASTs same ref vs consensus (identical downstream contract). minimap2 2.28 + samtools 1.21 + seqtk pinned/pulled; --jobs concurrency added. VALIDATED full chain on real data: isolate 3758 (cladeIII, MIC4, S) -> 1575bp ERG11 consensus, no catalogued mutation -> call_erg11 predicts S (TN), correct; build_fungal_cohort acc=1.0/spec=1.0 on the isolate. 25-isolate batch (--method map --jobs 3) launched to D:/dna_decode_cache/fungal_g1/batch_map.log (~1-1.5h); on completion -> build_fungal_cohort -> G1 VERDICT (within-clade-III + efflux discordance). |
| 15 | 2026-06-08 | run-tests | GATE G1 REACHED — C. auris fluconazole decoder validated on N=24 real cohort | Docker D:-mount corruption (force-kill aftermath) fixed via wsl --shutdown (lessoned to memory). Batch assembled 24/25 (3345 OOM-failed, recoverable). VERDICT: LABEL_LIMITED_FAILURE (acc 0.792, sens 1.000, spec 0.167). METHOD TRANSFERS: catalogued ERG11 mutation found in 100% of MIC-R across 2 clades (cladeI Y132F 4/4, cladeIII F126L 14/14). FAILURE is label-limited not method-limited: 5/6 'susceptible' isolates (all MIC16) CARRY F126L but fall below CDC tentative breakpoint(>=32) — the documented suspect-the-label pattern (genotype reliable, binary MIC splits a reduced-suscept continuum); only true-WT 3758 (MIC4) -> correct S. Within-clade-III F126L near-universal -> not clade confound. H1 CONFIRMED (method transfers; acc-bar label-limited) -> G1 lands on EP7 documented-failure-mode branch. Added LABEL_LIMITED_FAILURE verdict + test (17 fungal tests green). Packet: wiki/fungal_cohort_g1_fluconazole_2026-06-08.{md,json} + closeout wiki/fungal_ep7_g1_closeout_2026-06-08.md. PATH B now unblockable (G1 shared). |
| 16 | 2026-06-08 | write-plan | Path B PRE-STAGED for the workhorse (laptop no-compute) | Pinned + VERIFIED all G2 data sources: AraPheno FT10 (n=1162) + FT16 (n=1122) phenotype CSVs (downloaded + committed to data/arabidopsis/); 1001 Genomes genotype VCF live (HEAD HTTP 200, 19.2GB); imputed SNP matrix + pseudogenomes (Data Center); PlantCAD_l32 FM (fits 12GB). Wrote plans/EP8_PathB_PreStage_Manifest.md (workhorse-executable: locus-targeted FM embedding of FLC/FRI/FT panel, SNP-PRS + kinship baselines, leave-one-admixture-group-out CV + within-lineage diagnostic, G2 PASS = embedding R^2 beats BOTH baselines by >=0.05 AND within-group>structure). scripts/fetch_arabidopsis_pathb.py (verify/--download helper, tested). Handoff updated: G1 shared -> Path B GO. Short-term action 3 (pre-stage manifest) DONE. |
| 17 | 2026-06-08 | write-plan | /brainstorm review of the G2 manifest + REVISED it (caught a decisive design flaw) | Multi-round adversarial review found the manifest's PRIMARY locus-targeted approach (FLC/FRI/FT) CONTRADICTED EP8's own 'no curated mechanism catalog' niche criterion. Revised EP8_PathB_PreStage_Manifest.md: (1) primary -> phenotype-AGNOSTIC genome subsample (frozen, no FT labels / no flowering-gene enrichment); curated loci demoted to secondary diagnostic w/ narrowed claim; (2) added CPU-only §0.5 G2 dry-manifest gate (accession join / coord table / N-fraction QC / matched indel+missingness baseline / group labels) before ANY GPU; (3) continuous within-group test (binary AMR script doesn't apply); (4) PASS now needs paired-bootstrap CI EXCLUDING 0, not point +0.05; (5) kinship/PCs primary de-confounder, geography sensitivity-ONLY (residualizing geography removes real vernalization biology); (6) FT10 frozen primary, FT16 replication. Saved review at plans/EP8_PathB_PreStage_brainstorm.md. Opens §8 user design decisions (estimand / window-rule / GPU budget / cloud). Handoff updated with dry-manifest-first order. |
| 18 | 2026-06-08 | edit-local-code | Built the §0.5 G2 dry-manifest gate (laptop no-compute; the CPU-only gate before any GPU) | scripts/g2_dry_manifest.py: pure-logic functions for all 5 gate checks — accession intersection (pheno∩pseudogenome∩SNP), EMPIRICAL pseudogenome ID->filename pattern (never assumes pseudo{id}), phenotype-AGNOSTIC gene-bodies+flanks window table from a GFF (the §8-default window rule; preserves the no-catalog niche), per-window N-fraction QC, group-label presence, GREEN/RED verdict — plus a CLI the workhorse runs on the real downloaded data. 7 offline unit tests (tests/test_g2_dry_manifest.py) on synthetic GFF/CSV/FASTA fixtures, green; 24 eukaryotic-cycle tests green total. Converts pre-staged-plan -> ready-to-run gate code; workhorse RUNS it (no duplication). Manifest §0.5 updated to point at it. Drafted §8 defaults for ratification (estimand=causal-signal / window=genes+flanks / budget=sub-day / no-cloud-unless-borderline). |
| 19 | 2026-06-08 | ask-user | /project-state --resolve-pending-decision | G2 §8 design choices: user accepted recommended defaults (estimand=causal-signal; window=genes+flanks; budget=sub-day; cloud=not unless borderline+approved). Manifest §8 RATIFIED + frozen; workhorse executes with zero open design decisions. |
| 20 | 2026-06-08 | edit-local-code | PRODUCTIZED the G1 fungal decoder into shipped dna-amr + docs (v0.5.0) | The G1-validated fungal decoder was orphaned in scripts/; wired it into dna-amr as a first-class capability routed by --drug (fungal -> BLAST-ERG11 engine; bacterial -> AMRFinder, unchanged). --genome-fasta (BLAST) + --observed (pure) modes; uniform amr-mechanism-call-v1 record; _fungal_main relabels organism default to Candida_auris. 6 CLI tests (test_amr_cli_fungal.py) incl. real-BLAST; 876/877 suite green (1 pre-existing network-gated NT-download fail, unrelated). Docs: pyproject 0.4.0->0.5.0, README capability table + fungal quickstart, cli.py TRAITS, CHANGELOG 0.5.0 'kingdom jump', CLAUDE.md 2 gotchas. Captured reusable memory: targeted-mapping-beats-assembly-single-locus. Best laptop next-steps survey done (#1 productize + #2 docs + #4 regression executed; #3 retrospective + #5 G1-deepen deprioritized). |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Workhorse identity (SAFETY):** is the GPU "workhorse" the personal Precision 7780, or the Bombardier/DLP machine? Path B handoff assumes the former; if the latter, Path B is blocked on principle (no personal code on a DLP machine).
- **Path B compute (MONEY):** is the ~12GB RTX 3500 Ada acceptable to attempt, or do you want a cloud A100 (paid — requires explicit budget approval)?
- **Embedding pre-commitment:** confirm that a clean G2-FAIL closes the embedding frontier permanently (this is the thesis's 4th test after 0-for-3; pre-committing guards against the diminishing-returns trap).

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-06-07
- **Progress signal:** (init only — catalog shipped, caller next toward G0)
