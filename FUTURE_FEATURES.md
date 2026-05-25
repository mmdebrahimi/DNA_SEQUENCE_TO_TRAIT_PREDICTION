# Future Features

Capability ideas for Phase 2+. Captured from the technical plan's Phase 2 Backlog section + `/probe` + `/brainstorm` synthesis output. Phase 1 ships first; nothing below blocks Phase 1.

## High priority (Phase 2 — direct follow-ons)

### v0.1 cef decoder substrate via cef-S label backfill — 18/05/2026
BV-BRC strict-MIC 4-drug census 2026-05-18 found cef has 66 HIGH_R + 2 HIGH_S — heavily R-skewed but high assembly_accession retention (58% vs cipro's 4%). If cef-S labels can be backfilled from PATRIC / NARMS / EuSCAPE / other clinical surveillance sources, cef becomes the second viable decoder substrate (v0.1 in the locked v0 spec). Without S backfill, cef stays infeasible. Gating step: identify a source with sufficient HIGH-MIC cef-S strains + downloadable assemblies. Reuses the existing pipeline + `dna_decode/data/mic_tiers.py` per-drug catalog (cef breakpoints + AMRFinder class filter + loci-by-mechanism already configured) — only the cohort source changes.

### Audit-framework-as-product / standalone cohort-quality gate — 18/05/2026
The audit framework (mechanism × MIC × opacity merge with structurally-enforced SUSPEND gate, shipped Phase 1 closeout) was empirically validated as the project's most reusable artifact by the BV-BRC strict-MIC infeasibility finding 2026-05-18. Architectural insight: the audit framework RELAXES the cohort-construction bar — pre-filtering on strict-MIC duplicates the audit framework's job downstream AND discards usable training data. Productization idea: lift the audit framework to a standalone tool (`scripts/cohort_audit.py --drug X --cohort Y --ast-csv Z`) that ANY AMR researcher can run against their cohort + get a SUSPEND/MIXED/CLEAN verdict + per-strain noise classification. Reuses `dna_decode/data/mic_tiers.py` per-drug catalogs + the cipro_mic_audit / cipro_mechanism_audit / cipro_mechanism_phenotype_merge pipeline generalized via `mic_tiers.classify_gene_symbol`. Open-source release path = stronger portfolio artifact than another classifier.

### 4th-mechanism-class smoke test (partition-hypothesis falsifier) — 17/05/2026 (substrate updated 18/05/2026)
Symmetric with EP2 D5(b) for H17. The Phase 1 cross-drug architectural finding partitions NT-frozen-pooling behavior by mechanism-class data shape (concentrated PASS vs distributed-mobile-element FAIL). Validating the partition requires a 4th mechanism class. **Substrate update 2026-05-18**: gentamicin (originally proposed as the aminoglycoside 4th-class candidate) is INFEASIBLE at strict-MIC per the BV-BRC 4-drug census — only 2R / 132S strains pass all filters. New candidates: (a) **colistin via mcr-family plasmids** — concentrated plasmid acquired-gene signal; expected to PASS per the partition; (b) **trimethoprim via dfrA family** — concentrated plasmid + chromosomal; expected to PASS; (c) **fosfomycin via fosA family** — concentrated; expected to PASS. The partition is falsified if EITHER a concentrated-signal mechanism FAILS smoke OR a distributed mobile-element mechanism PASSES smoke. Reuses the existing 12-strain mini-cohort smoke runner pattern with `--drug X` + a freshly-built mini-cohort from BV-BRC AST. Gated on a BV-BRC strict-MIC feasibility check producing ≥6R/6S strains for the chosen 4th drug.

### Attribution Refinement Engine — 12/05/2026
Multi-stage local-refinement pipeline on top of Phase 1's gene-level ISM. Stages: sliding-window occlusion at sub-gene resolution, codon-aware mutagenesis (synonymous vs non-synonymous scoring), annotation overlap with RegulonDB / KEGG / EcoCyc, cohort variant association (Fisher exact / odds ratio), clade-controlled regression, bidirectional counterfactual editing (susceptible → introduce candidate → predict; resistant → revert → predict), multi-model attribution agreement scoring. Source: `/probe` + ChatGPT exchange 2026-05-12.

### Differentiable MLP head + Captum IG — 12/05/2026
Phase 1's XGBoost is non-differentiable so Captum IG was removed (post-tech-plan brainstorm C1). Phase 2 adds an alternate 2-layer MLP head trained on the same frozen embeddings + Captum IG attribution path. Cross-check against ISM on the same top-K genes; systematic disagreement = model exposing different signal under different probes.

### MIC regression head — 12/05/2026
Phase 1 is binary R/S only. Phase 2 adds continuous MIC regression head sharing the foundation-model embedding cache. Requires solving BV-BRC MIC sparsity + CLSI-vs-EUCAST breakpoint reconciliation. Source: original ChatGPT spec; deferred at Q2 resolution.

### Pan-genome graph layer (Panaroo + PyTorch Geometric) — 12/05/2026
Step 6's cohort catalog is drug-first selection. Phase 2 adds Panaroo/Roary pan-genome clustering → gene-family presence/absence matrix → PyG GNN over (gene-family nodes, co-occurrence/synteny edges). Captures epistasis + plasmid-context the per-strain embedding loses.

### Pygenometracks adapter — 12/05/2026
Phase 1 ships matplotlib + TSV viz (ship-path plan D4). Phase 2 adds publication-grade pygenometracks adapter for SME presentation work. Gate: project actually reaches a "show this to a biologist" phase.

## Medium priority (Phase 3 — architectural upgrades)

### Multi-task learning shared encoder — 12/05/2026
Currently each drug gets its own XGBoost classifier on frozen embeddings. Phase 3 shifts to a shared sequence encoder with multiple task heads (AMR per drug, MIC per drug, AMR-gene detection auxiliary, mutation-effect, plasmid-marker, calibration/OOD). Curriculum training: start with simple drugs + known markers; progressively add multi-drug + clade-aware + counterfactual edits.

### Haplotype Disambiguation Module — 12/05/2026
When SNP A + SNP B always travel together in resistant strains, both get high ISM/IG attribution but only one may be causal. Adds linkage-disequilibrium detection, conditional regression, haplotype-stratified tests, natural-recombination-event search.

### Annotation gene_symbol column — 12/05/2026
Root-cause cleanup for the Wave 3.5 C8 fix. Extend `parse_gff3` to extract the GFF3 `gene=` attribute as a separate `gene_symbol` column. Phase 1 workaround at the `build_attribution_report` layer is sufficient; Phase 2 cleans up upstream.

### Persisted MLST → clade_id mapping in cohort metadata — 14/05/2026
Replace runtime hash-based clade ID derivation (currently `zlib.crc32(scheme.st)` per `Sidework_Sequence_Ship_Path_Plan.md` D7) with an explicit `{mlst_string: int}` mapping serialized into cohort parquet metadata at cohort-build time. Eliminates ANY hashing concern (collisions, library upgrades, platform drift). Makes clade IDs auditable + greppable. "Cleanest scientific approach" per ChatGPT cross-engine review 2026-05-14. Trigger: when the deferred per_clade_baseline fix lands at Stage 1 and the coarser-binning decision (scheme-level vs Mash) is made — at that point, persist the mapping. Adds: schema field in cohort parquet, build-time numbering pass, serialization roundtrip test.

### Method-aware AST reweighting — 12/05/2026
Phase 1 filters BV-BRC AST to broth-microdilution only. Phase 2 adds method-aware label reweighting (handles disk-diffusion + E-test rows with appropriate uncertainty) + CLSI-vs-EUCAST breakpoint reconciliation. Source: failure-mode #4 from post-tech-plan brainstorm.

### Live BV-BRC + NCBI Datasets API integration — 12/05/2026
Replace the local-TSV + scaffolded `NotImplementedError` paths in `pilot.py` with real REST endpoint calls. Refresh-cadence-gated + retry-with-backoff per the pattern in `refseq.py`.

## Low priority (Phase 3-4 horizon)

### AlphaFold-inspired architecture — 12/05/2026
Long-term architectural shape: sequence encoder (existing) + pan-genome/evolutionary representation (Phase 2 graph) + pairwise interaction module (gene-gene + plasmid-gene + regulator-effector) + multi-task trait heads + confidence estimation + iterative refinement. **Caveat:** does NOT fit on consumer GPUs (project's actual GTX 860M 4 GiB is far too small; even a hypothetical RTX 4090 24 GiB would be tight); requires A100+ on rented compute. Source: `/probe` exchange — useful as architectural inspiration, not literal replication.

### Confidence estimation layer — 12/05/2026
Output prediction + confidence + reason, not just class label. Signals: probability calibration + ensemble agreement + OOD distance + clade-only-baseline gap + attribution stability + AMRFinder concordance.

## Nice to have (Phase 4+ horizon)

### Sequence-design search / RL — 12/05/2026
Once Phase 3 has a calibrated predictor, an "action-taking" agent can propose minimal sequence edits that flip predicted phenotype. Methods in increasing complexity: beam search → Bayesian optimization → genetic algorithms → RL. Per `/probe` exchange: RL is properly Phase 4+, NOT Phase 1 as initially intuited.

### Multimodal genotype-phenotype platform — 12/05/2026
Long-term vision: combine DNA + images + annotations + phenotype metadata for richer organisms. NOT a direct stepping stone from Phase 1 — would require a parallel multimodal track (contrastive learning, BioCLIP-style foundation model, image-paired genotype datasets). User original framing: "DNA → animal image / vice versa"; reframed via `/probe` to "constrained DNA → trait prediction → marker/lineage inference from observed phenotype."

### Cross-organism extension — 12/05/2026
Beyond E. coli: Mycobacterium tuberculosis (structurally similar; rich AST data) → yeast S. cerevisiae (first eukaryote; colony morphology phenotype) → plants (image + DNA paired datasets exist) → animals (constrained morphology prediction). Each stage is a substantial domain break, not a smooth progression.
