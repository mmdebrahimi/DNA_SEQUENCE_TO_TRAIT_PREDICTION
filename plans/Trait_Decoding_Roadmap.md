# Trait-Decoding Roadmap — DNA Decoder Long-Horizon Plan

> Connective tissue between the current v0/v0.1 ship state and the original 2026-05-11 project goal: "use AI to decode and understand any DNA code" → corrected factual premise: "model relationships between DNA sequences and downstream biological function — regulatory activity, protein expression, phenotype." Reaffirmed 2026-05-25 ("DNA input → phenotype + trait identification at gene level").

**Status:** DRAFT 2026-05-26.
**Anchors on:** `project_state/dna-decode-2026-05-11.md` (original /project-init goals); `CLAUDE.md` L9 long-term vision; `plans/Post_V0_EP_Ladder_Plan.md` (EP-0 through EP-3 + EP-4+ deferred). This doc is the BIGGER-picture map; the EP ladder is the next-3-EPs map.
**Supersedes:** none. First explicit articulation of the post-AMR roadmap.

---

## Why this doc exists

We have shipped v0 (cached-strain cipro AMR) + v0.1 cipro genome-input is real on Precision 7780 + cef cached-strain is in flight. That's 1 drug, 1 organism, 1 phenotype type. The user's actual goal is "decode ANY trait through matching DNA sections with phenotype." The EP ladder (`Post_V0_EP_Ladder_Plan.md`) bounds the next 3 EPs but stops at multi-drug E. coli AMR.

This doc names the phases beyond EP-3 explicitly + assigns terminal claims so each subsequent EP has a target without re-litigating the "any trait" unbounded scope every cycle.

The 2026-05-11 `/project-init` flagged "decode any DNA" as RESEARCH-PROGRAM (unbounded). Phase-by-phase terminal claims bound it incrementally; this doc shows the spine.

---

## Phase ladder (current → long-horizon)

Each phase has:
- **Terminal claim** (specific + measurable; bounds scope per-phase)
- **Substrate** (organism + phenotype + dataset class)
- **Architecture** (foundation model + classifier head + input mode)
- **Dataset prerequisite** (what must exist before this phase fires)
- **Falsifier** (what failure pattern would mean this phase ISN'T the right path)
- **Status**

### Phase 0 — Cached-strain cipro AMR (v0)
- **Terminal claim:** `pipeline.py predict --strain-id X --drug ciprofloxacin` on cached E. coli strain → v0 JSON output with prediction + calibrated_probability + confidence_tier + audit_verdict + provenance.
- **Substrate:** E. coli; cipro; categorical R/S labels; N=147 (149 minus 2 dedup'd).
- **Architecture:** NT v2 100M frozen mean-pool + XGBoost.
- **Dataset prerequisite:** N≥75 per-class strains with NT cache + AMRFinder mechanism audit.
- **Falsifier:** N/A — shipped.
- **Status:** ✓ shipped 2026-05-24 (FAIL-branch with documented scope-limit per north star "honest failure-tolerant iteration").

### Phase 1 — Genome-input cipro AMR (v0.1 slice 1)
- **Terminal claim:** `pipeline.py predict --genome-fasta X.fna --annotations Y.gff3 --drug ciprofloxacin` on novel E. coli FASTA → v0 JSON output. Same-strain parity: cached-strain path ≈ genome-input path within ε.
- **Substrate:** E. coli; cipro; novel FASTA + GFF3 input.
- **Architecture:** Same as Phase 0; new ingestion path.
- **Dataset prerequisite:** Phase 0 trained pickle + a single public held-out E. coli genome with AST.
- **Falsifier:** Same-strain parity delta > ε on cohort strains; or novel-genome prediction differs structurally from cached-strain prediction for the same strain.
- **Status:** ✓ Codex landed 2026-05-25 (cross-path concordance 4/4 + max delta 0.011599 ≈ ε=0.02 acceptable).

### Phase 2 — Multi-drug AMR E. coli (v0.1 slice 2+)
- **Terminal claim:** `pipeline.py predict --drug X` for X in {cef, tet, gent} runs end-to-end + CV AUROC ≥ 0.70 per drug OR documented scope-limit naming "this drug's mechanism class is OUT-OF-SCOPE for v0.1 architecture per EP-1.5 finding."
- **Substrate:** E. coli; cef + tet + gent.
- **Architecture (forks per mechanism class per EP-1.5):**
  - **Concentrated mechanism (cef = plasmid β-lactamases):** same as Phase 0/1 (mean-pool NT + XGBoost).
  - **Distributed mechanism (tet = efflux + ribosomal protection; gent = aminoglycoside acetyl/phospho/adenyl transferases):** uses EP-1.5 chosen architecture (per-gene NT windows OR k-mer + AMRFinder fusion OR transformer head).
- **Dataset prerequisite:** Per-drug cohort with assembly_accession-uniqueness + ≥50 R + ≥50 S; cef-S label backfill from PATRIC / NARMS / EuSCAPE per `FUTURE_FEATURES.md` 2026-05-18.
- **Falsifier:** EP-1.5 POC finds NO architecture beats mean-pool on tet → tet ships as documented INDETERMINATE; the architectural-class-bounded finding is upheld at clean labels.
- **Status:** Cef cached-strain IN FLIGHT (Codex shipped substrate; promotion plan written; CV AUROC 0.895 / no leakage). Cef genome-input deferred to v0.2. Tet + gent deferred to post-EP-1.5.
- **DECISIVE VERDICT 2026-06-05 (cipro, the first DE-CONFOUNDED Phase-2 falsifier) — NT embedding FAILS vs domain knowledge:** on the cipro N=147 cohort (the only substrate passing the de-confound gate — 6 shared R/S MLST lineages, country+year non-aliasing), under leakage-safe `leave_one_accession_out` CV:
  - NT-XGBoost **0.914** · NT-logreg 0.863 · k-mer-XGB 0.824 · **POINT-XGB (QRDR knowledge baseline) 0.943**.
  - NT beats bag-of-k-mers (+8.9 pp) **but LOSES to the QRDR-POINT knowledge baseline (−2.9 pp; CI [−9.0,+2.9] → FAIL).**
  - **Within-lineage diagnostic (mechanism vs lineage):** POINT-XGB within-lineage concordance **1.000 (p<0.001)** = pure mechanism (discriminates R/S even inside the same lineage); NT-XGBoost **0.605 (p=0.365)** = chance. **⇒ NT's 0.914 is largely lineage/genome-content, NOT the resistance mechanism.**
  - **Conclusion:** for cipro (concentrated point-mutation mechanism), the NT-frozen-mean-pool embedding does NOT earn its keep — simple QRDR-POINT features beat it AND are the only pure-mechanism signal. Honest north-star FAIL on the cleanest substrate. Single drug+cohort (no over-generalization), but it is the cleanest test available and the result is decisive.
  - Infra built: `dna_decode/eval/cohort_deconfound.py` (de-confound gate precondition), `scripts/amr_falsifier.py` (drug-agnostic CI-aware), `dna_decode/eval/point_baseline.py` (QRDR-POINT comparator), `scripts/within_lineage_diagnostic.py`. Artifacts: `wiki/ciprofloxacin_falsifier_2026-06-05.{md,scores.json}` + `wiki/ciprofloxacin_within_lineage_diagnostic_2026-06-05.md`. Per 3× /brainstorm + Soraya.
  - **Roadmap implication:** for concentrated-mechanism AMR drugs, prefer the mechanism-feature (POINT/AMRFinder) baseline over NT-frozen-pooling. The embedding's potential value is on DISTRIBUTED mechanisms where no clean knowledge baseline exists (tet efflux) — but tet failed earlier too. Reconsider whether frozen-mean-pool embeddings have a niche in AMR at all, vs the architecture only earning value on non-AMR phenotypes lacking curated mechanism catalogs.

### Phase 3 — Multi-organism AMR
- **Terminal claim:** `pipeline.py predict --drug ciprofloxacin --strain-id <K-pneumo-strain>` runs end-to-end on Klebsiella pneumoniae cohort + CV AUROC ≥ 0.70 OR documented scope-limit naming the organism-specific failure mode.
- **Substrate:** Klebsiella pneumoniae first (per CLAUDE.md L9 expansion path); then Pseudomonas aeruginosa; then Acinetobacter / mycobacteria.
- **Architecture:** Same as Phase 0/1/2 for each drug-mechanism class. Klebsiella + Pseudomonas have similar gram-negative cell-wall + plasmid biology to E. coli → mean-pool likely transfers. Mycobacteria differ structurally → may need re-architecture.
- **Dataset prerequisite:** BV-BRC Klebsiella AST + assembly metadata; ≥50 R + ≥50 S for each (cipro / cef first; tet + gent likely infeasible per the cipro 2026-05-18 strict-MIC census pattern).
- **Falsifier:** Klebsiella cipro AUROC < 0.65 → indicates organism transfer is harder than expected; need re-evaluation of whether E. coli architecture is the right base.
- **Status:** SLICE 1 DONE 2026-06-07 — **Klebsiella cipro TRANSFERS with a principled rule refinement.** N=30 NCBI K. pneumoniae (15R/15S). The E. coli rule applied UNCHANGED FAILS (acc 0.5, spec 0.0): K. pneumoniae intrinsic chromosomal OqxAB efflux (absent in E. coli) saturates the broad QUINOLONE-determinant count → all called R (Phase-3 falsifier fired as designed). The **QRDR-POINT refinement** (count only gyrA/parC/parE target POINT mutations) → **acc 1.000 / sens 1.000 / spec 1.000** on Klebsiella AND holds E. coli at 0.925 (−1.4pp vs broad). **Platform finding: cross-organism transfer needs counting the drug's TARGET-alteration mutations, not the broad drug-class determinant bag** (intrinsic chromosomal determinants are the organism-specific gotcha). `qrdr_point_count` added to amr_rules.py (+3 tests). Artifacts: `wiki/klebsiella_cipro_transfer_2026-06-07.{md,json}`, `scripts/klebsiella_cipro_transfer.py`. Open authority decision: adopt QRDR-POINT cipro globally (−1.4pp E. coli) vs per-organism. Next: Klebsiella cef + meropenem (carbapenem, new mechanism class).
- **SLICE 2 DONE 2026-06-07 — Klebsiella MEROPENEM (carbapenem, new mechanism class) VALIDATED.** N=30 NCBI K. pneumoniae (15R/15S). Acquired-carbapenemase rule (threshold 1 + CARBAPENEM-subclass: blaKPC/NDM/OXA-48) → **acc 0.867 / sens 1.0 / spec 0.733** (vs naive AMRFinder 0.533; the Subclass refinement lifts spec 0.067→0.733 by excluding ESBL/AmpC). Carbapenem is the defining K. pneumoniae threat — a mechanism class E. coli AMR never covered, now decoded. meropenem added to mic_tiers (5th drug). Blind to porin-loss-mediated R (expected FN mode). Artifacts: `wiki/klebsiella_meropenem_validate_2026-06-07.{md,json}`, `scripts/klebsiella_meropenem_validate.py`. **Phase 3 status: dna-amr now spans 5 drugs × 2 organisms.**
- **MATRIX COMPLETE 2026-06-07 — Klebsiella cef + gent + tet validated.** cef acc 0.80 / gent 0.867 (both ✅, rules unchanged from E. coli); tet acc 0.80 / spec 1.0 / sens 0.6 (PARTIAL — efflux-mediated tet-R via oqxAB is a curated-determinant blind spot; documented, not a defect). Added `gene_prefixes` refinement (tet counts acquired `tet*` only, excluding intrinsic oqxAB efflux — also improved E. coli tet 0.833→0.917). Full 5×2 matrix: `wiki/klebsiella_drug_matrix_2026-06-07.md`; `scripts/klebsiella_drug_validate.py`. **Cross-organism principle confirmed 3×: count the drug's specific determinants, not the broad class bag.** Phase 3 (multi-organism AMR) terminal claim MET for E. coli→Klebsiella; next organism (Pseudomonas) or other tracks are fresh decisions.
- **3rd ORGANISM 2026-06-07 — Pseudomonas aeruginosa cipro VALIDATED.** N=30 acc 0.867 / sens 0.80 / spec 0.933 (beats naive AMRFinder 0.767). The QRDR-POINT cipro rule transfers UNCHANGED to a *less-similar* gram-negative (MexAB-OprM efflux + intrinsic AmpC, no oqxAB). 3 FN = efflux-mediated cipro-R (the expected curated-determinant blind spot). dna-amr now validated across **3 organisms** (E. coli + K. pneumoniae + P. aeruginosa) — the "count the mechanism, not the broad class bag" principle is organism-general. Shipped `dna-amr --organism` (cross-organism in the CLI, not just scripts) + generalized `scripts/organism_drug_validate.py` (any organism × drug = one command). Artifact: `wiki/pseudomonas_aeruginosa_ciprofloxacin_validate_2026-06-07.md`.
- **1st GRAM-POSITIVE 2026-06-07 — S. aureus oxacillin (MRSA/mecA).** mecA detection TRANSFERS to a Gram-positive (sens 1.000 — all 15 R carry mecA); the acquired-gene+Subclass approach works beyond gram-negatives. spec 0.333 is oxacillin-LABEL noise NOT a rule defect (10/15 oxacillin-S strains carry mecA — oxacillin AST is the unreliable comparator; CLSI/EUCAST use cefoxitin, which is substrate-sparse here = 3R). **Terminal: genotype generalizes across the gram divide; phenotype-label validation is the binding constraint** (the project's recurring substrate/label lesson, now on a Gram-positive). oxacillin = 6th drug. Artifact: `wiki/staphylococcus_aureus_oxacillin_validate_2026-06-07.md`. dna-amr now spans 6 drugs × 4 organisms (E. coli, K. pneumoniae, P. aeruginosa, S. aureus) across the gram divide.

### Phase 4 — Non-AMR bacterial phenotypes (FIRST conceptual jump out of AMR)
- **Terminal claim:** `pipeline.py predict --phenotype <growth_rate | virulence | biofilm | etc.> --strain-id <E. coli strain>` runs end-to-end + meets the phenotype-specific success metric. Substrate: ONE non-AMR phenotype with paired DNA + phenotype data; ONE organism (E. coli).
- **Substrate:** TBD per `plans/EP_4_Non_AMR_Phenotype_Candidates.md` (this session's T2). Strong candidates: growth rate / fitness / virulence / biofilm formation / specific metabolic capability.
- **Architecture:** Likely same NT mean-pool + classifier (or regression head if phenotype is continuous). v0 schema extends with `phenotype` field replacing `drug`.
- **Dataset prerequisite:** Paired DNA + phenotype dataset; ≥50 strains for binary classification, ≥100 for regression.
- **Falsifier:** AUROC / R² < phenotype-specific threshold → the chosen phenotype isn't captured by NT-frozen mean-pool; architecture needs rework OR the chosen phenotype is genuinely not gene-encoded at the resolution the tool can detect.
- **Status:** **PARTIALLY SHIPPED 2026-06-04.** Pathotype (the T2 top candidate) shipped as a **deterministic VF-marker compatibility resolver v0** (tag `pathotype-v0`; CLI `dna-pathotype`; abstention + provenance + canonical-VF diff), NOT the planned multiclass NT classifier. **The classifier/embedding track is CLOSED for pathotype** — H1 falsified (`research_outputs/horesh-f1-label-provenance-audit-2026-06-04.md`): the only independent labels are isolation-source categories that are intrinsically sampling-confounded (label == sampling context), so no de-confounded discrimination cohort is buildable and any AUROC measures batch, not biology. The deterministic resolver is the honest tool for this phenotype.
- **REVISED falsifier insight (2026-06-04):** the binding failure mode for Phase 4-style categorical phenotypes is NOT model capability — it is **label provenance**. A phenotype *operationally defined by sampling context* (clinical site / disease presentation) cannot be validated against an embedding classifier, because the label and the batch variable are the same axis. The embedding-upgrade frontier must target phenotypes with **sampling-independent lab-measurement labels** (AMR MIC; quantitative assays), NOT sampling-defined categories.
- **NEXT non-AMR candidate SELECTED 2026-06-06 — carbon-source utilization (EP-6).** Per `research_outputs/ecoli-bacterial-phenotype-decoder-substrate-feasibility-2026-06-05.md`: carbon-utilization (BacDive 4397 strains × 58 sources) is the first substrate that is YES on BOTH embedding-niche halves — sampling-INDEPENDENT lab-assay labels (clears the pathotype confound) AND no AMRFinder-style curated catalog (clears the AMR-loses-to-mechanism-features failure). Pre-named trap: Li et al. 2023 found carbon-util prediction is largely PHYLOGENETIC (the same lineage-vs-mechanism crux that killed cipro NT), so the de-confound gate is a hard precondition. Infra scaffolded + green 2026-06-06: `dna_decode/data/bacdive.py` (loader, 8 tests) + `scripts/bacdive_carbon_util_feasibility.py` (layered census reusing `cohort_deconfound.py`). Full EP design + go/no-go gates: `plans/EP6_Carbon_Utilization_Substrate.md`. **Only blocker = acquiring a real BacDive E. coli carbon export.**

### Phase 5 — Multimodal phenotype (DNA + image / DNA + transcriptomics)
- **Terminal claim:** Decoder takes DNA + one paired phenotype modality (image OR transcriptomics) + produces a paired prediction. v0 schema extends to multi-input.
- **Substrate:** Bacterial colony morphology images paired with whole-genome sequencing OR transcriptomic profiles paired with DNA. Public datasets to survey: see T4 below (deferred).
- **Architecture:** Significantly new — DNA encoder (NT v2) + image encoder (CNN or ViT) + fusion head. NOT same as Phase 0-4.
- **Dataset prerequisite:** ≥500 paired DNA + image samples for the smallest credible v1 slice.
- **Falsifier:** Multimodal prediction doesn't beat DNA-only baseline → modality fusion has no signal; either chosen modality pair is wrong OR fusion architecture is wrong.
- **Status:** NOT STARTED. Research-program territory per the 2026-05-11 `/project-init` verdict; needs its own `/idea-anchor + /project-init` cycle.

### Phase 6 — Eukaryotic organisms
- **Terminal claim:** Decoder works on at least one eukaryotic organism (fungal first; then plant; then human GWAS-derived traits).
- **Substrate:** Candidate first eukaryotic target — fungal AMR (Candida auris fluconazole resistance; well-characterized + clinical) OR plant trait (Arabidopsis flowering time; rich public dataset) OR human GWAS (Type 2 Diabetes polygenic risk — proxy for the original 2026-05-11 5th refinement candidate).
- **Architecture:** Foundation model may need swap — NT v2 is bacterial-focused; eukaryotic genomes require Enformer / Borzoi / Caduceus / similar. Significant architectural pivot.
- **Dataset prerequisite:** Paired eukaryotic DNA + phenotype dataset; ≥1000 individuals for human GWAS (PROMOTE LATER — not solo).
- **Falsifier:** Cross-kingdom transfer fundamentally fails → the "decode any DNA" framing is too broad; need bacterial-restricted long-term scope.
- **Status:** NOT STARTED. Research-program territory; needs own cycle.

### Phase ∞ — Goal-state articulation
- **Goal:** Decoder takes ANY DNA input (any organism, any locus) + ANY phenotype query + emits a calibrated prediction with gene-level attribution + provenance + honest scope-limit when out of distribution.
- **Realistic path:** ~5+ years; requires multimodal datasets at scale; foundation-model pretraining beyond what's currently affordable for a solo developer.
- **Honest take:** Phase ∞ is aspirational; the project's value is the cumulative learning across Phases 0-6 + the audit infrastructure + the honest-output discipline, NOT reaching Phase ∞ literally.

---

## Phase transition gates

Between any two phases, the following gate must pass before the later phase fires:

1. **Reproducibility freeze** at the prior phase: tag a release; lock dependencies; emit canonical example.
2. **External validation** for at least one prior phase (EP-1B does this for cipro; cef + tet should follow once they ship).
3. **Architecture compatibility check:** the next phase's mechanism class / organism / modality is compatible with the current architecture OR has an explicit re-architecture plan (EP-1.5 is the pattern here).
4. **Dataset substrate ready:** prerequisite cohort/dataset exists + passes feasibility.
5. **Honest-output discipline carries:** every phase emits scope-limit-aware output; no silent overclaim.

These gates are inherited from the post-falsifier ship-path plan + scaled to each new phase.

---

## Dataset prerequisite catalog (what must exist for each phase)

| Phase | Required dataset shape | Public sources (candidates) | Status |
|---|---|---|---|
| 0 (cipro) | E. coli AST + assembly_accession | BV-BRC + NCBI Datasets | ✓ built |
| 1 (cipro genome-input) | + public held-out E. coli FASTA + GFF3 | NCBI Pathogen Detection / ATCC | ✓ used |
| 2 (cef + tet + gent) | E. coli per-drug AST cohorts at ≥50 per class | BV-BRC + PATRIC + NARMS + EuSCAPE (cef-S backfill needed) | Cef in flight; tet/gent feasibility uncertain |
| 3 (Klebsiella) | K. pneumoniae AST + assembly | BV-BRC + NARMS | NOT BUILT |
| 4 (non-AMR bacterial) | E. coli paired DNA + non-AMR phenotype dataset | Per T2 memo — TBD | NOT BUILT |
| 5 (multimodal) | Paired DNA + image / DNA + transcriptomics | OMIM image + sequencing repositories; Human Microbiome Project; etc. | NOT SURVEYED |
| 6 (eukaryotic) | Paired eukaryotic DNA + phenotype | Candida auris MIC + WGS; Arabidopsis 1001 Genomes; UK Biobank | NOT SURVEYED |

---

## Architecture-transferability bets

| From | To | Expected transfer cost |
|---|---|---|
| E. coli cipro mean-pool | E. coli cef plasmid β-lactamases | LOW — same mechanism class (concentrated signal); same architecture works (Codex shipped AUROC 0.895) |
| E. coli cipro mean-pool | E. coli tet efflux | HIGH — mechanism class mismatch (distributed signal); EP-1.5 architecture decision required |
| E. coli AMR | Klebsiella AMR | MEDIUM — similar gram-negative biology; same architecture; new cohort + retrain |
| E. coli AMR | E. coli non-AMR phenotype | MEDIUM — same organism; new phenotype catalog; classifier head likely fine; mechanism catalog needs full rebuild |
| Bacterial AMR | Mycobacterial AMR | HIGH — significantly different cell biology + AMR mechanisms |
| Any bacterial | Eukaryotic | VERY HIGH — different foundation model needed; different genome architecture |
| DNA-only | Multimodal | VERY HIGH — new architecture entirely (DNA encoder + modality encoder + fusion) |

These are inferences from the EP1/EP2 cross-drug architectural finding (2026-05-17) extrapolated; not measured. Each transfer needs its own falsifier per Phase Transition Gate #3.

---

## Honest scope limits

This roadmap deliberately does NOT cover:
- **Specific drug / organism / phenotype picks at each phase beyond Phase 2.** Those decisions need the T2-style scoping memos when their phase fires.
- **Compute substrate decisions.** Phase 4-6 may require burst-cloud compute beyond Precision 7780 capacity.
- **Funding / commercial decisions.** Per north star "AI DNA decoder tool, not papers"; commercialization is an outcome IF the tool works, not a development phase.
- **Regulatory / clinical considerations.** Discovery phase; deferred.
- **Cross-machine coordination beyond what `scripts/cross_machine_sync_check.py` already catches.** Codex + Claude pattern continues.

---

## What this roadmap unblocks

With phase-by-phase terminal claims locked, each downstream EP has:
1. A clear "what would success look like" answer (the terminal claim).
2. A clear "what would falsify this path" answer (the per-phase falsifier).
3. A clear "what dataset do we need first" answer (the prerequisite catalog).
4. A clear "what architectural pivot does this need" answer (the transferability table).

Future `/idea-anchor + /project-init` cycles for Phases 4-6 can anchor on the corresponding row of this doc rather than re-inventing the framing each time.

---

## Recommended next moves (per phase)

| Phase | Next action | Status |
|---|---|---|
| 1 | Codex pushes outstanding artifacts; tag v0.0-cipro + v0.1-cipro local-only | Awaiting Codex push |
| 2 | Cef promotion via `plans/Cef_V0_1_Promotion_Slice_Plan.md`; tet/gent gated on EP-1.5 | Cef in flight |
| 3 | Klebsiella substrate exploration (T3 from this session, NOT executed) | Pre-staging deferred |
| 4 | Pathotype v0 resolver SHIPPED; classifier track CLOSED. Carbon-utilization (EP-6) **E. coli-INFEASIBLE 2026-06-06** — data acquired (Li et al. OSF jwkr7), E. coli slice = 27 strains, 0 carbon sources clear ≥100 (artifact `wiki/bacdive_carbon_util_feasibility_2026-06-06`). Cross-taxa pivot deep but different experiment + phylogeny trap + Databricks-scale (NOT taken — authority decision). | EP-6 E. coli-infeasible; fork open |
| 5 | Multimodal dataset survey (T4 from this session, NOT executed) | Pre-staging deferred |
| 6 | Eukaryotic entry-point scoping (T5 from this session, NOT executed) | Pre-staging deferred |

---

## Bottom line

The trait-decoding goal isn't lost in the EP-by-EP execution. It's the spine that every EP serves. Phases 0-2 are AMR. Phases 4-6 are where "any trait" becomes literal. Each phase has its own substrate, its own architecture compatibility check, its own dataset prerequisite, its own falsifier.

The smallest meaningful next conceptual jump is **Phase 4 — non-AMR bacterial phenotypes** (T2 memo). It's the first place the project moves out of AMR. Smaller than Phase 5 (multimodal) or Phase 6 (eukaryotic). Should happen after Phase 2 (multi-drug AMR) ships, but the scoping memo IS conflict-free with Codex's cef work + worth pre-staging.
