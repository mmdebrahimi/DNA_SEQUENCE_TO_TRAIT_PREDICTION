# DNA Decode — AI/ML Project
<!-- project-schema: 0.1 -->

> Initialized 2026-05-11. Project ID: dna-decode-2026-05-11. Originating goal (verbatim user input): "use AI to decode and understand any DNA code which is a programming language written in the language of 3-5 proteins".

## Project Context
- **Project ID:** dna-decode-2026-05-11
- **Project root:** C:\Users\Farshad\PythonProjects\dna_decode
- **Captured:** 2026-05-11
- **Originating goal:** use AI to decode and understand any DNA code which is a programming language written in the language of 3-5 proteins
- **Refined goal (from 3c top-ranked candidate):** Build a tissue-specific gene-regulatory-element activity predictor for chromosome 21 with ≥0.7 AUROC, using a transformer-based DNA model trained on ENCODE / Roadmap Epigenomics multi-tissue data, by 2026-09-30.
- **Horizon (months):** 12
- **Schema:** project-schema 0.1

## Empirical Concerns
- **Verdict:** FAIL
- **Check status:** attempted
- **Provisional:** NO
- **Findings:**
  - **Nucleotide-count error:** DNA is composed of 4 nucleotides (adenine A, thymine T, cytosine C, guanine G — the "ACGT" four-letter alphabet), not 3-5. The "3-5" figure has no biological correspondence at the nucleotide level. Source: genome.gov ACGT glossary; Wikipedia Nucleotide-base.
  - **Protein-vs-DNA conflation:** Proteins are downstream translation products of DNA (DNA → mRNA → ribosomal translation → protein), built from 20 standard amino acids encoded by 3-nucleotide codons. Proteins are NOT the alphabet of DNA. There are ~20K-100K distinct human proteins, none of which fits the "3-5" figure either. Source: Wikipedia Genetic code; DNA-and-RNA codon tables.
  - The "programming language" analogy is a common intuition pump, but the "language" of DNA is the 4-nucleotide alphabet at the encoding level + 64-codon table mapping codons to 20 amino acids + start/stop signals + regulatory elements (promoters, enhancers, etc.). It is NOT "written in 3-5 proteins."
- **Corrected factual premise:** Use AI to model relationships between DNA sequences (composed of 4 nucleotides: A, T, C, G, encoded into 64 codons that translate via the genetic code to 20 standard amino acids) and downstream biological function — regulatory activity, protein expression, phenotype.

## Project vs Research-Program
- **Verdict:** FAIL
- **Provisional:** NO
- **Classification:** research-program
- **Rationale:** Originating goal phrases scope as "any DNA code" (unbounded — all DNA from all organisms). No terminal condition, no measurable success criterion, no scoping by organism / chromosome / functional element / phenotype. Matches the unbounded-goal pattern flagged by Step 3b ("understand any X", "decode all Y"). A bounded sub-goal is required for project-shape work; the unbounded framing is research-program territory and would not converge in any practical horizon. Successful DNA-AI projects (AlphaFold, Enformer, Basenji, scBERT, DeepBind) all operate on narrowly-scoped sub-problems.

## Refinement Candidates
- **Verdict:** FAIL
- **Provisional:** NO
- **Refined-from:** corrected-factual-premise
- **Candidates:**
  1. **Regulatory-element prediction on chr21** — Build a tissue-specific gene-regulatory-element activity predictor for human chromosome 21 with ≥0.7 AUROC, using a transformer-based DNA model fine-tuned on ENCODE/Roadmap Epigenomics multi-tissue ChIP-seq + ATAC-seq, by 2026-09-30. **Top-ranked — selected for HTN decomposition.**
  2. **Non-coding-variant pathogenicity classifier** — Build a baseline classifier on a ClinVar-style benchmark achieving ≥0.65 AUC vs. Enformer-style pre-trained baselines, focused on interpreting which DNA-sequence features drive classifier decisions, by Q3 2026.
  3. **AlphaFold reproduction on a small protein family** — Reproduce AlphaFold-style protein-structure prediction on a 50-protein homologous family; benchmark predicted vs experimental structure on TM-score, by Q3 2026.
  4. **LLM-based gene-function annotation tool** — Take gene name as input; produce natural-language summary of gene function; evaluate against UniProt-curated descriptions on 100 held-out genes with ≥0.7 BLEU or human-rating ≥3/5 by Q4 2026.
  5. **Sequence-to-phenotype model on constrained dataset** — Train on GWAS-significant SNPs for one disease (e.g., Type 2 diabetes); benchmark against existing polygenic-risk-score methods, by Q4 2026.

## Goal Hierarchy
### Long-term (12+ months tier)
Build a tissue-specific gene-regulatory-element activity predictor for chromosome 21 with ≥0.7 AUROC, then expand methodology to additional chromosomes and refine toward state-of-the-art performance benchmarks.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Literature review on DNA regulatory-element prediction | Annotated reading list of ≥10 primary papers (Enformer, Basenji, scBERT, DeepBind, Akita, etc.) with method-comparison table | 2026-06 (1 month) |
| 2 | Data pipeline for chr21 multi-tissue regulatory data | Working pipeline ingesting ENCODE/Roadmap chr21 data across ≥5 tissues with ChIP-seq + ATAC-seq tracks | 2026-07 (2 months) |
| 3 | Baseline model implementation | Enformer-style transformer fine-tuned on chr21 subset; AUROC measured on held-out tissue | 2026-09 (3 months) |
| 4 | Held-out evaluation + iteration to ≥0.7 AUROC | Test-set AUROC ≥0.7 with reproducible eval harness | 2026-12 (3 months) |
| 5 | Methodology write-up + decision: expand vs pivot | Internal write-up summarizing findings; go/no-go decision on chr1-22 expansion | 2027-01 (1 month) |

### Short-term (≤1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Confirm Refinement Candidate selection (default Candidate 1) | ask-user | user | 1 day |
| 2 | /research "DNA regulatory element prediction state of the art 2026" | research | user | 1 week |
| 3 | /research "ENCODE Roadmap Epigenomics data availability chromosome 21" | research | user | 1 week |
| 4 | Survey existing pre-trained models (Enformer weights, Basenji checkpoints) for fine-tunability | research | user | 1 week |
| 5 | Compute-budget decision: local GPU vs Colab vs cloud provider | ask-user | user | 1 week |

## State Snapshot
### Assumptions
- User has Python ML background sufficient to fine-tune transformer models (confidence: medium — confirmed by Athena/RCA-Engine work in PythonProjects/)
- ENCODE / Roadmap Epigenomics data is publicly available and downloadable for chr21 multi-tissue tracks (confidence: high — well-documented public datasets)
- Pre-trained DNA-sequence transformers (Enformer, Basenji) have permissive licenses allowing fine-tuning (confidence: high — Enformer is Apache-2.0)
- Local hardware OR cloud budget allows non-trivial GPU training (confidence: low — not yet confirmed; see short-term action 5)
- 12-month horizon is realistic for a single-person regulatory-element prediction project to reach 0.7 AUROC baseline (confidence: medium — depends heavily on prior experience + compute access)
- The "AI to decode DNA" intent is genuinely about ML modeling, not about a different framing (e.g., natural-language explanation of genomics; gene-therapy automation) (confidence: medium — user's verbatim goal is ambiguous)

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| 1 | DNA's alphabet is 4 nucleotides A/T/C/G | genome.gov + Wikipedia (WebSearch 2026-05-11) | high | 2026-05-11 |
| 2 | Proteins use 20 standard amino acids; encoded via 3-nucleotide codons | Wikipedia Genetic code (WebSearch 2026-05-11) | high | 2026-05-11 |
| 3 | Enformer is a published SOTA transformer for regulatory-element prediction | LLM training knowledge | high | 2026-05-11 |
| 4 | /idea-anchor pre-init pass produced a Formal Rephrase that deliberately stripped the "3-5 proteins" factual error from the goal text; the error was instead surfaced in /idea-anchor Assumptions and Blunt Opinion sections | (manual) | medium | 2026-05-11 |
| 5 | User refined goal 2026-05-11 post-brainstorm: G2P (genotype-to-phenotype) prediction platform starting with E. coli; specifically infer which genomic regions / patterns are responsible for particular traits, with biologically interpretable explanations. Long-term: multi-organism comparative genomics. Synthesizes Candidates 2+5+6 from initial Refinement Candidates list. | user post-brainstorm refinement | high | 2026-05-11 |
| 6 | E. coli is a prokaryote; Enformer / HyenaDNA / Caduceus (eukaryote-biased) are inappropriate. Prokaryote-aware foundation models: Evo (Together AI / Stanford, 7B params, microbial, 131K context — likely strongest), DNABERT-2 (BPE k-mer-free, multi-species), Nucleotide Transformer (174B nt training corpus, 850 species incl. bacterial), GENA-LM, DNABert-S | model landscape calibration | high | 2026-05-11 |
| 7 | Antibiotic resistance is the most tractable G2P phenotype for E. coli v1: large clean labeled corpora (CARD, ResFinder, NCBI AMRFinderPlus, BV-BRC AST), well-characterized resistance loci (gyrA, parC, blaCTX-M etc.) for ground-truth interpretability validation, binary or MIC continuous outputs, published baselines to beat | tractability calibration | high | 2026-05-11 |
| 8 | User refined spec via ChatGPT 2026-05-11; strong alignment with technical plan (~90% match). 5 substantive diffs identified: (A) 3-drug Phase 1 instead of 1 — accepted; (B) MIC regression Phase 1 vs Phase 2 — pushed back, Phase 2; (C) Comparative model benchmarking Phase 3 vs Phase 1 — push to Phase 1 (free given Step 7 wrapper); (D) SHAP+attention attribution Phase 1 vs Phase 2 — pushed back, Phase 2; (E) Pan-genome clustering Phase 1 vs Phase 2 — pushed back, Phase 2. | spec alignment review | high | 2026-05-11 |
| 9 | 4 failure modes neither spec version addresses: (1) phylogenetic confounding (model learns clade signature not resistance signal — mitigate via phylogeny-permutation negative control); (2) plasmid vs chromosomal encoding (β-lactamases on transferable plasmids; filter to complete-circle assemblies Phase 1); (3) annotation drift (Bakta vs Prokka gene calls differ — pin Bakta version); (4) BV-BRC AST measurement-method variance (filter to broth-microdilution Phase 1) | failure-mode audit | high | 2026-05-11 |
| 10 | v0.2 /project-state design requirements surfaced by post-/brainstorm-#3 readiness check: (1) --resolve-pending-decision flag for marking pending decisions as resolved (currently must manually edit rows); (2) --refresh-frame for updating Bellman target state + candidate-actions table when project scope shifts; (3) --update-last-evaluation for marking the Last Evaluation field on /project-step invocations; (4) --remove-row for cleanup of superseded rows. v0.1 limitations forced manual ledger edits for Candidate-1→Candidate-7 state drift cleanup 2026-05-11. | v0.2 design requirement | high | 2026-05-11 |
<!-- project-state:end:evidence -->
### Unknowns
- Specific tissue-set best suited for chr21 multi-tissue training (depends on data availability + biological relevance to user's interests)
- Whether 0.7 AUROC is achievable on chr21 alone with fine-tuning vs requires from-scratch training
- Whether user wants the project to bias toward novel ML methodology vs reproducing existing methods on a constrained dataset
- Compute-budget ceiling (local GPU available? cloud spend approved?)
- User's prior bioinformatics depth (Python ML known; genomics-specific tooling unknown)

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | Enformer pre-trained weights can be fine-tuned on chr21-only data to reach ≥0.7 AUROC on tissue-specific regulatory-element prediction | open | never |
| H2 | Public ENCODE/Roadmap chr21 multi-tissue data is sufficient (≥5 tissues with both ChIP-seq + ATAC-seq) for the proposed fine-tuning approach | open | never |
| H3 | User's existing compute access (local + any cloud credits) is sufficient for fine-tuning a single Enformer-scale model on chr21 subset within 3-month horizon | open | never |
| H4 | E. coli scale (~4.6 Mbp genome, ~5K genes, ~1000-5000 strain pan-genome) is tractable on a single mid-range GPU (RTX 3090/4090/A6000) when foundation-model embeddings are pre-computed once and cached | open | never |
| H5 | Frozen-embedding + XGBoost classifier baseline will reach within 5pp AUROC of LoRA fine-tuning on antibiotic-resistance prediction for E. coli, making it the right phase-1 default (frozen-first, fine-tune-later) | open | never |
| H6 | Pretrained Evo (7B params, microbial corpus, 131K context) will outperform DNABERT-2 and Nucleotide Transformer on E. coli antibiotic-resistance prediction by ≥3pp AUROC | open | never |
| H7 | REVISED 2026-05-11 post-tech-plan brainstorm: ISM-only attribution (gene-level + saturation; NO Captum in Phase 1 per C1) on Evo embeddings will overlap known resistance loci at top-K=20 with ≥0.6 precision for ciprofloxacin (gyrA, parC, qnr-family, aac(6')-Ib-cr) and ≥0.4 for ceftriaxone (CTX-M family, SHV family, AmpC overexpression loci). Captum IG with differentiable MLP head is Phase 2. | open | never |
| H8 | Phylogeny-permutation negative control (shuffle resistance labels within MLST clades) will reduce AUROC from ≥0.85 to ≤0.6, demonstrating the model learns mechanistic signal rather than clade signature | open | never |
| H9 | β-lactam attribution-precision will be lower than fluoroquinolone attribution-precision (target ≥0.4 vs ≥0.6) because β-lactam resistance involves distributed signal across plasmid β-lactamase genes + porin mutations + efflux pumps, not single point mutations | open | never |
| H10 | Drug-first cohort construction (per-drug ≥150 strains; broth-microdilution filter; assembly-quality threshold contig_count ≤500 + N50 ≥50K; AMRFinder plasmid/chromosome annotation as informational) retains ≥150 strains per drug for cipro/ceftriaxone/tet; 3-drug intersection ≥75 strains | open | never |
| H11 | Clade-only baseline classifier (Mash/ANI cluster identity features, no sequence) achieves ≤0.65 AUROC on held-out Mash clades for cipro resistance, vs Evo-embedding model ≥0.80 on the same clades; gap ≥0.15 validates mechanistic signal over clade signature on ≥75% of held-out clades | open | never |
<!-- project-state:end:hypotheses -->
### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| Refinement Candidate 1 (regulatory-element prediction on chr21) selected as actionable HTN target | 2026-05-11 | Default selection at init per /project-init Step 4 research-program-verdict branch; user may override via /project-state --append-decision |
| OVERRIDE: Selected Candidate 7 (synthesis of Candidates 2+5+6 reframed for E. coli) over default Candidate 1 | 2026-05-11 | User refined goal post-brainstorm: G2P platform starting with E. coli (prokaryote), antibiotic resistance as v1 phenotype, biologically interpretable predictions (NOT causal claims in v1), multi-organism comparative expansion long-term. Architecture sketched; full Technical Plan deferred to /technical-plan invocation. |
| Phase 1 Technical Plan locked: 16 steps, 7 waves, fluoroquinolone-resistance as the v1 drug target, MockFoundationModel as the CI-grade smoke-test substrate, Evo + DNABERT-2 as the real-run foundation models | 2026-05-11 | Written to plans/Ecoli_G2P_Platform_Technical_Plan.md after /technical-plan invocation. Next gates before /execute-plan: (1) git init + remote configuration; (2) compute strategy confirmation; (3) optional /brainstorm pass on the plan itself (post-tech-plan /brainstorm is load-bearing per user-memory feedback_post_tech_plan_brainstorm_loadbearing.md). |
| Phase 1 Technical Plan REVISED post-/brainstorm 2026-05-11: 17 steps in 8 waves; Captum REMOVED from Phase 1 (Step 12 deferred to Phase 2 with differentiable MLP head per C1); ISM extended with nucleotide-level saturation mutagenesis (Step 11); drug-first cohort with assembly-quality threshold + AMRFinder plasmid/chromosome metadata (Step 6 renamed); Mash/ANI phylogeny clustering + clade-only baseline + per-clade reporting (Step 10); Step 0.5 real-data pilot gate added; Step 17 comparative model benchmarking leaderboard added; Phase 1 drugs specified as ciprofloxacin + ceftriaxone + tetracycline (not "β-lactam" class) | 2026-05-11 | Captures Codex critique C1 (Captum-XGBoost incompatibility) + M1 (drug-first cohort) + M2 (stronger phylogeny controls) + T2 (pilot gate) + Adjustments A/C accepted/refined. β-lactam agent specified as ceftriaxone (CTX-M family attribution target; 3rd-gen cephalosporin clinically high-stakes). All changes preserved in plan diff at plans/Ecoli_G2P_Platform_Technical_Plan.md. Phase 1 ready for /execute-plan once git init + compute strategy confirmed. |
<!-- project-state:end:decisions-made -->
### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| RESOLVED 2026-05-11: Refinement Candidate selection → Candidate 7 (E. coli G2P platform) — see Decisions Made row 2 | /project-init | (resolved) | Kept for audit trail; v0.1 /project-state has no row-removal flag. v0.2 should add --resolve-pending-decision |
| RESOLVED 2026-05-11 (committed default): Compute strategy → single RTX 4090 with 4-bit Evo quantization via bitsandbytes; A100 fallback if 4-bit hurts attribution precision | /project-init | (resolved) | bitsandbytes added to Step 1 pyproject.toml deps; A100 rental on Lambda Labs as fallback documented in Risk Flags |
| RESOLVED 2026-05-11 (committed default): Target AUROC → SLO ≥0.80 (Phase 1 ships), target ≥0.85 (stretch); clade-only-baseline-gap ≥0.10 on ≥75% of held-out clades is the mechanistic-signal validation gate | /project-init | (resolved) | Defaults encoded in Verification section of plans/Ecoli_G2P_Platform_Technical_Plan.md |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

Per D5: surfaces the Bellman cognitive frame at the architecture level. v0.1 populates this at init; v0.2's `/project-step` reads from it to decide next actions. v0.1 does NOT iterate over the frame — it sets it up.

### Current state (one-line summary)
Project initialized + refined + planned 2026-05-11. Goal: E. coli G2P (Candidate 7). Phase 1 plan committed at plans/Ecoli_G2P_Platform_Technical_Plan.md (17 steps in 8 waves, post-/brainstorm revisions applied). Compute + AUROC defaults committed. Ready for git init + /execute-plan. No code written yet; no data ingested.

### Target state / terminal condition
Phase 1: Trained binary antibiotic-resistance classifier for ciprofloxacin + ceftriaxone + tetracycline on a drug-first E. coli pan-genome cohort (≥150 strains per drug; broth-microdilution AST; assembly-quality filter), evaluated via leave-one-Mash-clade-out CV with AUROC ≥0.80 SLO / ≥0.85 target + clade-only-baseline-gap ≥0.10 on ≥75% of held-out clades, with ISM attribution maps overlapping known resistance loci at top-K=20 precision ≥0.6 cipro / ≥0.4 ceftriaxone, all reproducible from a documented data + training pipeline. References plans/Ecoli_G2P_Platform_Technical_Plan.md Verification section.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` count + `gates-passed` count (raw counts, unweighted). At init: unknowns-retired = 0; gates-passed = 0 (3 of 3 goal-normalization sub-gates returned FAIL → no "passed" gates yet, the FAILs are intentional and successful catches).
- **v0.2+:** weighted combination of unknowns-retired, gates-passed, evidence-confidence-improved, hypotheses-falsified.

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | git init + remote configuration for dna_decode/ | propose | high (unblocks /execute-plan parallel mode) | low | low | minutes |
| 2 | /execute-plan plans/Ecoli_G2P_Platform_Technical_Plan.md → triggers Step 1 bootstrap (Wave 0) | write-plan | high (project scaffolding) | low | low | ~15 min |
| 3 | Step 0.5 real-data pilot gate (Wave 1; HARD gate) — verifies per-drug strain counts ≥150 before downstream ingestion fires | research | high (resolves H10) | high (real cohort numbers) | medium | ~1 hr (metadata-only HTTP) |
| 4 | Wave 1 parallel: Steps 2/3/4/5/7/10 ingest + parse + load deps | research+propose | high (unblocks all downstream) | medium | low | ~2-4 hrs |
| 5 | Wave 2: Step 6 cohort catalog + Step 8 embedding cache (depends on pilot pass) | research | high (drug-first cohort materialized) | medium | medium | hours |
| 6 | Wave 3: Step 9 baseline classifier + Step 11 ISM (gene-level + saturation) | run-tests | high (resolves H4, H5, H7 partial) | high | medium | hours |
| 7 | Wave 4-7: Steps 13 viz + 14 CLI + 15 smoke + 17 leaderboard + 16 docs | run-tests+write-plan | high (Phase 1 ships) | high (resolves H6, H7, H8, H11) | medium | days |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- **Default:** re-run `/project-state` after any action class fires (auto-append to Action Log triggers stale-state check)
- **Manual override:** user invokes `/project-state dna-decode-2026-05-11` at any time
- **v0.2+:** automated trigger when N actions fire OR T days elapse OR a hypothesis falsifies (TBD)

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
| 1 | 2026-05-11 | propose | /project-init invoked | ledger created |
| 2 | 2026-05-11 | research | WebSearch: "DNA nucleotide bases adenine thymine cytosine guanine four-letter alphabet" + "proteins amino acids count twenty standard DNA codon translation" | factual mismatches confirmed; Empirical Concerns FAIL |
| 3 | 2026-05-11 | propose | /project-state --append-observation | /idea-anchor pre-init pass produced a Formal Rephrase that deliberately stripped the "3-5 proteins" factual error from the goal text; the error was instead surfaced in /idea-anchor Assumptions and Blunt Opinion sections |
| 4 | 2026-05-11 | propose | /brainstorm (deep, generative-ideation v2.1) | 5 Refinement Candidates critiqued; Candidate 6 (DNA foundation-model benchmark + interpretability workbench) added as missing approach; ranking reassessed; cross-domain analogs from RCA Engine + Athena identified |
| 5 | 2026-05-11 | propose | /project-state --append-decision (×1) + --append-hypothesis (×4) + --append-observation (×3) | User refined goal to Candidate 7 (E. coli G2P platform); 3 new observations on prokaryote model landscape + tractability; 4 new hypotheses on compute scale, baseline performance, model ranking, and attribution overlap |
| 6 | 2026-05-11 | write-plan | /technical-plan | Phase 1 technical plan written at plans/Ecoli_G2P_Platform_Technical_Plan.md; 16 steps in 7 waves; max parallelism 6 (Wave 1); critical path Step 1 → 7 → 8 → 9 → 14 → 15 → 16; greenfield project; ready for /execute-plan once git init + remote configured |
| 7 | 2026-05-11 | propose | spec-alignment review against ChatGPT-refined user goal | 5 substantive diffs (A-E) + 4 failure modes captured as observations 8-9 + hypotheses 8-9; pending user decision on which adjustments to apply before /execute-plan |
| 8 | 2026-05-11 | propose | /brainstorm (deep, post-tech-plan, plan-review mode) against plans/Ecoli_G2P_Platform_Technical_Plan.md | Critical C1 (Captum-XGBoost incompatibility — load-bearing bug caught), Medium M1 (drug-first cohort) + M2 (phylogeny control insufficient), Tradeoff T2 (real-data pilot gate added). A/B/C/D/E adjustment verdicts: A partial-agree (refine to specific agents), B agree, C partial-agree (scope to frozen cohort), D agree-stronger (push Captum too), E agree (rename Step 6). |
| 9 | 2026-05-11 | write-plan | apply post-/brainstorm recommendations to plans/Ecoli_G2P_Platform_Technical_Plan.md | Step 12 removed; Step 11 extended with saturation ISM; Step 6 renamed + drug-first rewrite; Step 0.5 added; Step 10 strengthened with Mash/ANI + clade-only baseline + per-clade reporting; Step 17 leaderboard added; Step 13 dependency updated; Execution Preview rebuilt (8 waves; max parallelism 7); Risk Flags updated; Verification rewritten with ceftriaxone-specific targets; Problem Statement updated. |
| 10 | 2026-05-11 | propose | /project-state dna-decode-2026-05-11 (read-only) | Schema valid; 11 open hypotheses; 3 pending decisions; 9 action-log entries; 0 days since last evaluation. State observations: Pending Decision #1 stale (resolved); Bellman frame still references Candidate 1; v0.1 limitation flagged. |
| 11 | 2026-05-11 | propose | /brainstorm (deep, pre-/execute-plan readiness check, plan-review mode) | Caught C1 (Step 0.5 soft-vs-hard-gate ambiguity), M1 (Step 6 + Step 14 dependency drift), M2 (ledger state drift). Open-decision verdict: compute blocking-with-default (accept), AUROC non-blocking. Planning-paralysis verdict: proceed-to-execute after surgical fixes. |
| 12 | 2026-05-11 | write-plan | apply readiness-check fixes (C1 + M1 + M2 + M3) | Step 0.5 rewritten as HARD gate; Step 6 +Step 4 +Step 5 deps added; Step 14 +Step 10 dep added; Step 1 +bitsandbytes dep added; ledger state drift cleaned (Pending Decisions rows annotated RESOLVED with committed defaults; Bellman frame current+target state updated to E. coli G2P Phase 1; Candidate next actions table updated to current Phase 1 execution plan; Open Questions all marked resolved). |
| 13 | 2026-05-11 | propose | /probe — long-term-vision realignment check against ChatGPT's "DNA ↔ animal image" reframing | Verdict: Phase 1 aligns with "interpretable constrained G2P" goal but NOT with "DNA ↔ animal image" as a direct stepping stone. ChatGPT's 6-stage roadmap directionally OK but underestimates data/architecture discontinuity bacteria → eukaryotes → animals. Phase 1 transfers strongly to DNA→phenotype rigor; transfers weakly to multimodal DNA↔image. Multimodal pivot deferred — finish Phase 1 first; revisit multimodal track at Phase 2+ design time. |
| 14 | 2026-05-11 | write-plan | git init + initial commit + remote push | Initial commit (4 files: README.md, .gitignore, plans/Ecoli_G2P_Platform_Technical_Plan.md, project_state/dna-decode-2026-05-11.md). Remote: mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION (private). Default branch: main. Repo url: https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION. Ready for /execute-plan parallel mode (worktree branches will push to this remote). |
<!-- project-state:end:action-log -->

## Open Questions for User

All 5 init-time Open Questions resolved 2026-05-11:
- ~~Confirm Refinement Candidate selection~~ → RESOLVED: Candidate 7 (E. coli G2P; Decisions Made row 2).
- ~~12-month horizon realistic for Candidate 1?~~ → SUPERSEDED: Candidate 1 dropped. Phase 1 horizon for E. coli G2P is 3 months (per Phase-1 plan); 12-month horizon refers to Phase 1 + 2 + 3 completion.
- ~~Compute strategy?~~ → RESOLVED: single RTX 4090 + 4-bit Evo via bitsandbytes; A100 fallback (Pending Decisions row 2).
- ~~Target AUROC?~~ → RESOLVED: SLO ≥0.80 / target ≥0.85 / clade-baseline-gap ≥0.10 (Pending Decisions row 3).
- ~~ML modeling vs other intent?~~ → RESOLVED: ML modeling on DNA confirmed; biologically interpretable predictions, no causal claims in v1 (Decisions Made row 2).

No Open Questions blocking /execute-plan. Empirical questions remain (resolve via execution + pilot data):
- Will Step 0.5 pilot gate pass with ≥150 strains per drug? (Tests H10.)
- Will Evo 4-bit quantization preserve attribution precision? (Open assumption; revisit if H7 fails.)
- Will Mash CLI binary install cleanly on this Windows machine? (Risk Flag in plan.)

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-05-11
- **Progress signal:** (none yet — v0.1 init only)
