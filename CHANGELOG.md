# Changelog

All notable changes to `dna_decode`. Format loosely follows [Keep a Changelog](https://keepachangelog.com/);
this is a solo research-tool repo so the granularity is per-release-theme, not per-PR.

## [Unreleased]

- **`dna-decode fba --dead-ends` / `--gapfill-target` — find and repair the model's MISSING parts (design
  epoch, Track C).** The honest form of "predict what's absent": not nucleotide infilling, but the missing
  biochemistry, which is where a genome-scale model is actually wrong. Two capabilities kept apart by
  evidential weight — `--dead-ends` is a **structural fact** (metabolites produced but never consumed cannot
  carry steady-state flux; needs no labels, no donor, no network), while `--gapfill-target` is a
  **hypothesis** (which donor reactions would restore a wrongly-absent trait). **Worked end to end on a
  measured false negative:** iML1515 predicts no growth on sucrose, but BW25113 has a Sucrose carbon-source
  experiment in the Wetmore/Keio RB-TnSeq set (verified at source; that assay only runs sources the organism
  grows on). The diagnostic surfaces the gap **unprompted** — `suc6p_c` is one of 42 transport-fed dead ends,
  produced by `SUCptspp` and consumed by nothing, so the model carries a sucrose transporter that leads
  nowhere. Gap-filling against *Salmonella* iYS1720 proposes a single reaction, `FFSD`, and it takes sucrose
  growth **0.000 → 1.7798 /h** — about twice the glucose rate (0.877), as a disaccharide should — with
  glucose/xylose/cellobiose **unchanged**. Honesty rails in code, CLI output and tests: a gap is NOT
  automatically a defect (every proposal ships stamped `HYPOTHESIS`, since "repairing" a correct model
  fabricates biology), reversible reactions count as both producer and consumer (else the diagnostic invents
  dead ends), and `demand_reactions`/`exchange_reactions` are both off (inventing a sink is not a repair).
  Full write-up: `wiki/fba_gapfill_2026-08-07.md`.

- **`dna-decode fba --design-target … --milp` — MILP strain design; closes the enumeration's blind spot.**
  The bounded pair/triple enumeration is exhaustive only at depth 1 and cannot see a design whose members
  are individually unremarkable — which is most real ones. Two attempts to fix it by pre-ranking candidates
  both failed measurably (`competition_ranking`: first by ranking "is growth-associated" instead of
  "competes", flooding the top with ion transporters and BIOMASS; then, after matching growth, by leaving
  the product constraint too weak to expose `PFL`/`ALCD2x`), so the default stays the exhaustive strategy
  and is pinned by a test. The real fix is the bilevel MILP the literature uses: wired via **`straindesign`**
  (new optional `[design]` extra) as `find_coupled_designs_milp`. **It reaches `PFL + LDH_D + ALCD2x` and
  reports guaranteed flux 9.263835 — exactly matching the enumeration path**, because every MILP result is
  re-derived by `evaluate_knockouts` rather than trusted from the solver. Three diagnosable failures en
  route, now all pinned: a **PROTECT module is mandatory** (SUPPRESS alone returned 278 "designs" that
  simply kill the cell), **SCIP is required not optional** (GLPK lacks indicator constraints; its big-M
  fallback returned `unbounded`), and the **MILP floor is a fraction of WILD-TYPE growth, not the
  mutant-relative floor the enumeration uses** (at 0.9×WT the known design is `infeasible` by construction,
  since it grows at 52% of wild type). **Scope limit, measured:** the MILP finds what the pre-ranking cannot
  *within a scoped candidate set* (~23 s on the 11-reaction fermentation set) but does **not** make the
  unrestricted whole-model search tractable — the same formulation at the same working floor over all 2266
  gene-associated reactions hit `time_limit` with 0 solutions after 50 min on SCIP. Scoping to a pathway is a
  generic biological prior, not "knowing the answer", but `--milp` is not a point-it-at-the-whole-model
  button. Full write-up: `wiki/fba_strain_design_cell_2026-08-07.md`.

- **`dna-decode fba --design-target` — the DESIGN direction: product → edits (design epoch, Track A).**
  The FBA cell could answer *edit → trait*; this adds the inverse that strain engineering actually needs:
  given a product, which knockouts make producing it **necessary for growth**? Growth coupling is what
  makes an engineered strain stable — selection for growth becomes selection for production. Two-sided LP
  at a near-optimal growth floor: `min_flux > 0` = OBLIGATORY (a design), `min_flux ≈ 0 < max_flux` =
  POSSIBLE (the un-engineered case). `dna_decode/fba/design.py` + `scripts/fba_strain_design.py`.
  **Validated end-to-end: the search independently recovers the OptKnock-lineage anaerobic succinate
  design `PFL + LDH_D + ALCD2x`, guaranteed flux 9.26 vs a wild-type floor of 0.047** (a ~196× increase
  in the *guaranteed* floor), pinned by a slow test. Needs no labels and no new data — purely
  stoichiometric. Every design is stamped a **hypothesis for the bench, never a validated strain**.
  Four defects were caught by inspecting real output rather than by failing tests, and three of them
  failed in the same direction (return zero or garbage while appearing to work): counting `OBLIGATORY`
  instead of *improvement over wild type* reported **2096 of 2096** knockouts as designs; GPR isozymes
  blunt gene-level knockouts (default is now reaction-level); the growth floor must be **near-optimal and
  relative to each strain's own maximum** (the same design guarantees 0.0027 at 10% and 12.28 at 99%);
  and the top-ranked "design" was `ATPM`, the ATP-maintenance pseudo-reaction, now excluded by requiring
  a non-empty GPR. Full write-up: `wiki/fba_strain_design_cell_2026-08-07.md`.

## [0.12.1] — CORRECTNESS FIX: the FBA cross-organism registry shipped two wrong-organism models (2026-08-07)

**This supersedes 0.12.0 (tagged but never published). Anyone who ran `--organism saureus` or
`--organism paeruginosa` on 0.11.0 got a result for a DIFFERENT SPECIES.**

- **Fixed — `--organism` loaded the wrong organism's model for two of four aliases.** Verified against
  the BiGG Models API: `saureus` resolved to **iYS1720**, which is a ***Salmonella* pan-reactome**
  (1262/1707 gene ids carry the *S.* Typhimurium `STM` prefix), and `paeruginosa` resolved to
  **iJN1463**, which is ***Pseudomonas putida* KT2440**. Corrected registry:
  `saureus`→**iYS854** (*S. aureus* USA300_TCH1516), `salmonella`→iYS1720, `pputida`→iJN1463,
  `ecoli`→iML1515, `yeast`→iMM904 — every alias now loads a model of that actual organism.
- **Fixed — `--organism paeruginosa` now RAISES instead of silently substituting.** BiGG has no
  *P. aeruginosa* reconstruction at all; the alias reports that plainly and names the affected versions
  rather than handing back another species' model.
- **Fixed — emitted records hardcoded `"model": "iML1515"` / `"organism": "Escherichia coli K-12"`
  regardless of `--organism`.** Every `fba-metabolic-trait-v1` record (and the human-readable output +
  the `medium` field, which asserted E. coli's glucose-M9 default on any organism) now stamps provenance
  from the LOADED model via the new `model.MODEL_ORGANISM` / `organism_for()` source of truth.
- **Corrected framing — S. aureus / P. aeruginosa were never "LABEL_WALLED".** 0.12.0 reported both as
  blocked on a missing measured label. The real blocker was the wrong model. They are now separated:
  `LABEL_WALLED` (a model exists, no joinable label — S. aureus, Salmonella, P. putida) vs the new
  `MODEL_WALLED` (no model exists at all — P. aeruginosa; its gold standard *is* fetchable at PLOS
  `pcbi.1013945.s011`, sheets GOLD_84/GOLD_115, but is PA14-keyed with nothing to join to).
- **Regression guard** — `tests/test_fba.py` now asserts the invariant that would have caught this: the
  organism named by an alias must appear in the organism its resolved model actually reconstructs, every
  registered model must declare an organism, and `paeruginosa` must refuse. The previous tests *pinned
  the bug* (`assert resolve_model_id("saureus") == "iYS1720"`).

E. coli (iML1515, Keio-validated, MCC 0.652) and yeast (iMM904 vs SGD, MCC 0.252) results are
**UNAFFECTED** — both were correctly assigned. Frozen AMR/forward surfaces byte-unchanged.
Full evidence: `wiki/fba_wrong_organism_model_bug_2026-08-07.md`.

## [0.12.0] — the FBA cell matures: cross-organism generalization + honest per-organism validation, and the first NON-metabolic trait catalog (2026-08-03)

Purely **additive** — the frozen AMR/forward surfaces are **byte-unchanged from 0.11.0**. This release
hardens the v0.11.0 FBA metabolic cell with real per-organism + carbon-source validation (the cross-organism
claim is now *quantified*, not assumed), and opens a new axis with the first **non-metabolic** trait catalog
(motility). The recurring honest finding across all three: FBA/curated-catalog decoding is E. coli-strong and
label-walled elsewhere — the binding constraint is measured LABELS, not model capacity.

- **`dna-decode motility` — the first NON-metabolic trait catalog (flagellar swimming).** The AMR/metabolic
  determinant→phenotype paradigm applied to a physical BEHAVIOUR: gene presence → can the cell SWIM?
  MOTILE iff all 5 flagellar modules present (master `flhDC` → sigma-28 `fliA` → flagellin `fliC/fljB` →
  motor `motAB` → basal-body/export `fliF/fliG/flhA/fliI`); **chemotaxis (cheA/W/Y/Z) is reported SEPARATELY
  and does NOT gate motility** (a che-mutant still swims — gating on it would be a biology error).
  `dna-decode motility --genes flhD,flhC,... | --feature-table X.txt.gz`. KNOWLEDGE_BASELINE — curated
  flagellar catalog vs literature anchors (MG1655 + Salmonella motile; Shigella non-motile; flhDC/fliC/motAB
  KO non-motile); presence-based DIRECTION not speed; can't see a present-but-IS-disrupted flhD (v0.1
  sequence-mode). `dna_decode/motility/` + 12 tests; the four CLI-trait guards fire. Frozen AMR surface
  byte-unchanged. See `wiki/motility_catalog_v0_2026-08-03.md`.

- **FBA carbon-source growth validation — the QUANTITATIVE-growth axis (complement to essentiality).**
  `scripts/fba_carbon_growth_validate.py` + `dna_decode/fba/carbon_growth.py`: for a carbon source, swap it
  into iML1515's known-growing medium as the sole carbon (robust `model.medium` swap, not zero-then-reopen)
  → predicted growth RATE (/h). **RECALL 1.000 (21/21)** on measured-positive K-12 carbon sources
  (Keio/Wetmore), with a quantitative rate spread 0.21–0.94 /h (pentoses slower than glucose — biologically
  ordered). **Honest scope (R2 web probe):** a clean *measured growth-rate across carbon sources* dataset
  doesn't exist fetchably, and Biolog pos+neg for the K-12 strain is SI-locked (the 190-source Biolog set is
  E. coli Nissle — strain mismatch), so full specificity + a rate correlation are **EXTERNAL-walled**
  (named). Verify-in-batch surfaced a real finding: **BW25113 grows on sucrose but iML1515 has no sucrose
  transport → a false negative the validation exposes** (an honest model-gap, not a mapping gap). +3 tests +
  `wiki/fba_carbon_growth_validation_2026-08-03.md`. Frozen AMR/forward surfaces byte-unchanged.

- **Per-organism FBA essentiality validation — the cross-organism claim is now QUANTIFIED, not assumed.**
  `scripts/fba_essentiality_validate.py --organism <org>` generalizes the E. coli Keio validation: load the
  GEM → genome-wide single-gene-deletion essentiality → join a per-organism experimental gold standard →
  full metric panel (accuracy/MCC/precision/recall/ROC-AUC/PR-AUC). **Finding (honest): FBA essentiality does
  NOT transfer strongly from E. coli.** E. coli iML1515 = strong (MCC 0.652); yeast iMM904 vs SGD inviable-null
  (Giaever/SGD, keyed by systematic ORF) = **weak (accuracy 0.824 but MCC 0.252** — accuracy is flattered by
  the imbalanced majority, so MCC is the reported signal). S. aureus (iYS1720, Salmonella-style STM#### ids)
  + P. aeruginosa (iJN1463) are honest **LABEL_WALLED** (external — need a crosswalk / a fetchable Tn-seq set).
  > **⚠️ CORRECTED IN 0.12.1:** the S. aureus / P. aeruginosa half of this bullet is WRONG. iYS1720 is a
  > *Salmonella* model and iJN1463 is *P. putida* — the blocker was the wrong MODEL, not a missing label.
  > The E. coli + yeast numbers stand.
  `dna_decode/fba/essentiality_labels.py` (per-organism sources + pure parsers) + 4 tests + status roll-up
  `wiki/fba_per_organism_essentiality_2026-08-03.md`. Frozen AMR/forward surfaces byte-unchanged.

## [0.11.0] — the FBA metabolic-model cell: edit → quantitative cell-level trait (gene-KO, point-mutation, synthetic-lethality, cross-organism) (2026-08-03)

The **first GENERAL "gene edit → quantitative cell-level trait" cell** — mechanistic flux-balance analysis
over genome-scale metabolic models, in the OPEN mechanistic regime (computes from stoichiometry, sidesteps
the population-structure confounding that closed the learned organism-level track). The frozen AMR/forward
surfaces are **byte-unchanged**; this is purely additive. Ships the whole edit→cell-trait ladder: gene-KO,
point-mutation (composing `forward`), synthetic-lethality, and cross-organism generalization.

- **`dna-decode fba --organism` — the FBA cell generalizes across organisms.** The engine is organism-agnostic;
  `--organism ecoli|saureus|paeruginosa|yeast` (or `--model-id <BiGG id>`) loads the matching genome-scale
  model from BiGG (fetched + cached like iML1515). Verified on *S. aureus* iYS1720 (1707 genes, WT growth
  0.489/h, KO smoke).
  > **⚠️ CORRECTED IN 0.12.1:** `saureus` did NOT load an S. aureus model — iYS1720 is a *Salmonella*
  > pan-reactome, and `paeruginosa`'s iJN1463 is *P. putida*. The "verified on S. aureus iYS1720" claim
  > above is a **Salmonella** result. E. coli + yeast were correctly assigned. **Only E. coli is Keio-validated** — other organisms are v0 "the engine generalizes"
  with their own essentiality gold standard deferred (honest scope). `dna_decode/fba/model.py` `_BIGG_MODELS`
  registry + `resolve_model_id`.

- **`dna-decode fba --knockout A,B --synthetic-lethality` — two-gene edit → synthetic-lethality call.** Detects
  when a gene PAIR is lethal even though neither single is (each buffered by an isozyme/alternate route; the
  double breaks it) — how metabolic drug-target *pairs* are found. Verified on the real *dadX+alr* alanine-
  racemase isozyme pair (both singles viable, double lethal) vs a non-SL control (pgi+zwf). `model.synthetic_lethality`.

- **`dna-decode fba` point-mutation mode — composes `forward` + `fba` so a single MISSENSE edit → cell-level
  trait.** `dna-decode fba --gene gltA --mutation D362A --protein-seq S`: `forward` scores whether the
  missense breaks the enzyme (its own method-aware LOF call), and if damaging the gene is modelled as a
  knockout → FBA growth/essentiality; a *preserved* variant → wild-type (no metabolic change); an
  *uncertain* variant → the **conditional reported both-ways** (never a forced binary — the anti-theater
  rail). This is the rung above gene-KO: the user's original "minor edit to a genotype" (a point mutation),
  not just full knockouts. The chain **inherits** forward's DMS validation (missense→LOF) + fba's Keio
  validation (LOF→trait, accuracy 0.954); the only new piece is the ranker→LOF binarization, labelled a
  heuristic (forward's threshold), not a calibrated LOF probability. `--forward-method esm2/prosst/gemme/hybrid`
  upgrades the LOF call to the stronger DMS scorers. `dna_decode/fba/compose.py` + 5 tests; frozen AMR
  surface byte-unchanged (READ-only composition). See `wiki/fba_variant_compose_2026-08-03.md`.

- **`dna-decode fba` — the first GENERAL edit → quantitative cell-level-trait cell (mechanistic, not learned).**
  Knock out ANY of the 1515 genes in the **iML1515** E. coli genome-scale model → predicted **growth rate (/h)**
  + **essential/non-essential**, via flux-balance analysis (`cobrapy`). Unlike the learned organism-level regime
  (a closed negative — it learns population structure, not causation), FBA computes phenotype from stoichiometry
  + known biochemistry, so it sidesteps population-structure confounding by construction. `dna-decode fba
  --gene gltA` / `--knockout b0720,b0721` / `--wildtype`; `pip install dna-decode[fba]` (iML1515 auto-fetched
  from BiGG + cached; cobra bundles the GLPK solver). **Scope (honest): METABOLIC traits only** — growth /
  essentiality / secretion; NOT virulence / regulation. **Validated** genome-wide vs the free Keio-collection
  mutant-fitness gold standard (Bernstein 2023 method): **accuracy 0.954, MCC 0.652, ROC-AUC 0.863, PR-AUC 0.526**
  over 1339 genes (7.2% essential prevalence) — matching the published iML1515-vs-Keio literature (~0.93); 101
  FBA-essential genes corroborated by having no viable Keio mutant. In-distribution vs a knowledge baseline, not
  an independent-lab claim. `dna_decode/fba/` + `scripts/fba_keio_validate.py` + 10 tests; frozen AMR surface
  byte-unchanged. See `wiki/fba_keio_validation_2026-08-03.md`.

## [0.10.0] — the visible-trait animal fleet + more microbial/viral cells + the confound-free decoding-validation arm (2026-08-02)

Purely **additive** — the **frozen AMR decoder surface (`amr_rules.py` + `calibrated_amr_rules.json` +
`mic_tiers.py` + `shipped_decoder_surface.py`) is byte-unchanged from 0.9.0**, so every existing R/S call is
identical. This release expands the fleet toward the north-star "DNA → visible trait" direction and folds in
the confound-free genotype→phenotype validation arm.

**New visible-trait / physical-trait animal cells** (the north-star direction):

- `dna-coatcolor` — dog coat colour (E/K/A/B/D epistasis; Darwin's-Ark-validated) — first physical/visible
  animal-trait cell.
- `dna-morphology` — dog body size + ear type (pinned + Darwin's-Ark-validated catalog).
- `dna-horsecolor` — horse coat colour (E/A/CR/D/G epistasis).
- `dna-catcolor` — cat coat colour (W/O/A/B/D/C; X-linked orange → tortoiseshell).
- `dna-plumage` — chicken plumage (E/B/S/I/Bl/lav/c; Z-linked barring/silver).
- `dna-pigeoncolor` — pigeon plumage (B/TYRP1 Z-linked, SOX10, SLC45A2, NDP).
- A **shared mammalian coat-colour engine** spanning ~14 organisms (rabbit, mouse, goat, alpaca, guinea pig,
  fox, donkey, buffalo, camel, mink, roe-deer, …).

**More microbial + viral cells:**

- `dna-essentiality` — single-gene KO → essential/non-essential (conserved-core decoder).
- `dna-metabolic` — E. coli carbon-source utilization (uptake-gated catabolism decoder).
- `dna-kleb` — Klebsiella phage depolymerase → capsule KL-type (cross-organism, fetch-only).
- `dna-phage` — first-class bacteriophage genome/lineage → host-receptor class.
- HCMV v0.1 target-site (ganciclovir/cidofovir/foscarnet/letermovir) via `dna-amr --observed`.

**Confound-free decoding-validation arm** (research): genomic prediction (cv-ridge + gradient-boosted trees
+ permutation null) DECODES quantitative traits across **three kingdoms** on confound-free crosses — yeast
Bloom-2013 (12/12 traits), mouse BXD (brain weight r=0.57), Arabidopsis MAGIC (bolting r=0.57). Layer-2
finding: the nonlinear/epistasis advantage is **trait-architecture-dependent, not sample-size** (settled two
ways — a yeast power curve + a dense plant cross).

**FM-value regime map:** the definitive "when does a foundation model add value for genotype→phenotype"
synthesis — Regime A (reconstruction) = calibration not signal; Regime B (in-distribution) = the genotype is
a sufficient statistic → a simple model is near-optimal; Regime C (transfer / unseen variants) = the FM adds
real value (the DMS-validated `forward` cell, ProteinGym median Spearman ~0.49).

**Packaging:** sdist 8.4 MB → 588 KB (14×) via an explicit sdist include-list.

**Unchanged:** the frozen AMR/viral/fungal R/S surface, every evidence-contract, and all prior validation
numbers. No breaking changes.

## [0.9.0] — the usable-tool productization layer: input router + deployable strong-method variant-effect (2026-07-23)

Purely **additive** — the **frozen AMR decoder surface (`amr_rules.py` + `calibrated_amr_rules.json` +
`mic_tiers.py` + `shipped_decoder_surface.py`) is byte-unchanged from 0.8.1**, so every existing R/S call is
identical. This release turns the ~20-decoder fleet into a coherently *usable* tool and makes the
DMS-validated learned variant-effect frontier runnable, not just validated.

**Additions:**

- **`dna-decode decode <file>` — input-aware router.** Sniffs the input kind (nucleotide/protein FASTA or
  VCF) and lists every applicable decoder with its one-line claim, honest evidence tier (read live from the
  cell registry), and the exact command to run — the missing "I have a file, which decoders apply?" entry.
  It surfaces the 5 previously-orphaned cells (pgx/clinvar/hla/coloc/profile).
- **`dna-decode decode <file> --run`** — actually *runs* the auto-runnable decoders in one report (genome →
  the `profile` suite; protein → `inverse`), reporting the ones that need a specific parameter
  (`--mutation` / `--gene`) rather than guessing it.
- **Deployable strong-method variant-effect from the CLI.** `dna-decode forward` gains `--method
  esm2/prosst/gemme/hybrid/auto` + `--capabilities` (a preflight of which learned methods this host can
  run) + honest graceful degradation with provenance (`method_requested`/`method_used`/`degraded`). The
  strong methods (ESM2 / ProSST / GEMME / the validated ESM2+ProSST hybrid) were API-only.
- **Packaging extras for the learned frontier.** `pip install 'dna-decode[forward]'` (torch + transformers)
  for ESM2; `[forward,prosst]` (+ torch_geometric + biotite + pathos, plus a one-time AI4Protein/ProSST repo
  clone) for the structure quantizer. `[forward]` intentionally pins `transformers>=5` and is declared as a
  uv conflict with the closed-negative `[ml]` embedding extra (`transformers<5`). ProSST's local quantizer
  verified on this host: a raw AlphaFold PDB re-quantizes to the reference structure tokens 79/79.

**Unchanged:** the frozen AMR/viral/fungal R/S surface, every evidence-contract, and all prior validation
numbers. No breaking changes.

## [0.7.0] — the multi-kingdom decoder fleet: viral + human-PGx + typing cells + genome-map browser (2026-07-11)

First PyPI release since 0.6.4 (0.6.5 was an internal version bump, never published). The **frozen
bacterial/viral/fungal AMR decoder surface (`amr_rules.py` + `calibrated_amr_rules.json`) is byte-unchanged
from 0.6.4** — this release is purely additive, so every existing R/S call is identical. 310 commits since
0.6.4; themed per-release, not per-PR.

**Major additions since 0.6.4:**

- **Typing-cell fleet** — new `dna-decode` subcommands: `plasmid` (Inc-replicon), `serotype` (E. coli O:H),
  `resfinder` / `pointfinder` / `disinfinder` (independent cross-tool AMR checks), `mlst` (PubMLST ST),
  `ktype` (Klebsiella K-antigen), `salmserovar` (Salmonella Kauffmann-White), `pneumoserotype` (S. pneumoniae
  capsular — INDEPENDENT vs phenotypic Quellung, serogroup 0.939). Deterministic blastn callers, offline-safe degrade.
- **Human pharmacogenomics (`dna-pgx`)** — CYP2C19 / CYP2C9 / VKORC1 + CYP2D6 (real short-read WGS, 5 PGP-UK
  humans) + DPYD / NUDT15 / UGT1A1 / CYP4F2 / ABCG2; unified `dna-pgx --all` one-command 14-gene decode with drug
  annotations; PGx trust-surface report card; independent functional-evidence + trio co-segregation QC.
- **Viral expansion** — HIV-1 PI + INSTI + CAI classes (5 classes / 4 genes) with mutant-specific deconfounded
  v0.1 catalogs (NRTI / PI / INSTI); SARS-CoV-2 Mpro cell (Stanford CoV-RDB). Validated against free, independent,
  isolate-level wet-lab fold-change.
- **Genome-map graphical browser** (`dna_decode/genome_map/browser.py` + `scripts/genome_map_browser.py`) —
  self-contained HTML feature-track render with the evidence-tier honesty wall carried into the visual; a 5th
  `virulence-determinant` overlay tier.
- **Evidence-Contract Registry + certification capstone** — one test-enforced contract per shipped cell
  (67 cells / 5 tracks); no aggregate verdict (anti-"trust-theater" guardrail). Agent-era discoverability
  (`AGENTS.md`, quickstart / validation docs, CI green).
- **ABO blood-group decoder** (O/A/B/AB, free-label validated 0.902 on openSNP); photic-sneeze + asparagus
  single-locus falsification cells (both honestly fail vs the null baseline).
- **Validation infrastructure** — provenance-disjoint + external-cohort re-validation arms, clonality/lineage
  disclosure (Mash greedy-representative clustering), the prospective-lock accrual harness, a reproducibility
  freeze + negative-results map, and a fresh-cohort re-validation (the decoder reproduces on unseen genomes
  across all 10 SCORED cells).
- **Negative results (characterized, not hidden)** — a protein LM (ESM2) does NOT beat the curated catalog on
  antagonistically-selected AMR phenotypes (ESM2 peaks at 650M; below chance on HIV NNRTI), and the
  DNA-LLM-via-functional-alphabet path is a soft-negative. The binding constraint remains labels, not compute.

The 2026-06-26 infrastructure closeout that opened this release window (Evidence-Contract Registry +
certification capstone + agent-discoverability + the non-neural DNA-LLM probe) is detailed below.

- **Evidence-Contract Registry** (`dna_decode/data/cell_registry{,_vocab}.py`, `tests/test_cell_registry.py`):
  one checked-in, test-enforced contract per shipped cell so a new decoder cannot ship invisibly and
  abstention speaks ONE controlled vocabulary (not a confidence scale). **67 cells / 5 tracks** {amr 25
  (projected verbatim from the frozen `shipped_decoder_surface` via `canonical_cell_key`), viral 29 (HIV-1
  + SARS-CoV-2 — the CLI-routable-but-surface-absent gap), pgx 3, typing 6, finder 4}; `cli_routable_manifest()`
  derived LIVE from the CLI catalogs (drift-proof); `build_validation_report_card.py` now reads its AMR grid
  from the registry. NO numeric confidence / NO aggregate score (anti-"trust-layer-theater" guardrails).
- **Certification capstone** (`scripts/build_certification_capstone.py` → `wiki/certification_capstone.{md,json}`):
  a thin read-only presentation over the registry + per-domain report cards + the freeze/negative-results
  boundaries. NO aggregate boolean verdict (a 6-tier surface reduced to one bool would certify the weakest
  cell as strongly as the strongest) — `no_aggregate_verdict=True` pinned by test.
- **Agent-era discoverability** (`AGENTS.md`, `.claude/skills/dna-decode-demo/SKILL.md`, `Makefile`,
  `SECURITY.md`, `CITATION.cff`, `docs/quickstart.md`, `docs/validation.md`, `examples/README.md`,
  `.github/workflows/tests.yml`, richer `pyproject` keywords/urls): a `dna-decode list`-authoritative agent
  surface so coding agents can install/run/verify the tool. Research-use, not-clinical guardrails written
  down explicitly. CI verified green (1760 passed).
- **Non-neural functional-alphabet probe** (`dna_decode/eval/functional_tokens.py`,
  `scripts/functional_alphabet_probe.py`): the cheap CPU gate before any "DNA-LLM" GPU spend — a drug-general
  functional-determinant alphabet vs base-level k-mers, scored on the de-confounded within-lineage metric.
  Result across two mechanism regimes: cipro functional within-lineage **1.000** (concentrated QRDR), tet
  **0.963** (distributed efflux/ribosomal). The curated determinant alphabet separates R/S within-lineage on
  BOTH — **no demonstrated headroom for a learned model**; the DNA-LLM-via-alphabet path is a soft-negative.
  Re-confirms the project thesis: the binding constraint is labels, not compute. Closeout:
  `wiki/functional_alphabet_probe_closeout_2026-06-26.md`.
- FROZEN bacterial/viral/fungal AMR surface (`amr_rules.py` + `calibrated_amr_rules.json`) byte-unchanged
  throughout. 0 regressions (full suite 1760 passed; the 10 errors are the pre-existing xgboost-`[ml]`-extra
  predict-e2e, unrelated).

## [0.6.4] — PGx trust-surface report card + CYP2C9 sentinel v0.1 (2026-06-26)

Consolidation: a standing trust surface for the human-PGx phase + closing CYP2C9's non-core residual.

- **PGx trust-surface report card** (`scripts/build_pgx_report_card.py` → `wiki/pgx_report_card.{md,json}`):
  read-only roll-up (exit 0 always, a report not a gate) of every shipped PGx cell × validation axis
  (GeT-RM concordance / PharmCAT fixtures / functional-evidence verdicts / trio Mendelian QC), per the AMR
  report-card pattern. NO aggregate headline; each cell's honest tier stands alone; a missing sidecar renders
  NOT_RUN (never a fabricated number). 3 cells: CYP2C19 (72/72), CYP2C9 (73/73), VKORC1.
- **CYP2C9 sentinel v0.1** — honesty-parity with CYP2C19: the non-core SNP alleles **\*5/\*8/\*9/\*11**
  (rs28371686/rs7900194/rs2256871/rs28371685, GRCh38 coords grounded via Ensembl REST) are now **withheld**
  (any-non-ref-ALT-at-site sentinel) instead of silently mis-called \*1. GeT-RM re-validation: **core 73/73
  held; 10 non-core now correctly WITHHELD; genuine silent mis-call 14→4** (the remaining \*6-indel / \*61 /
  undefined tail — documented residual). The wildcard-ALT sentinel mode added to the shared caller.
- 3 new tests (`tests/test_pgx_report_card.py` + CYP2C9-sentinel assertions); 100 pgx-suite pass. FROZEN
  bacterial/viral/fungal AMR surface byte-unchanged.

## [0.6.3] — PGx independent functional-evidence + trio co-segregation layer (2026-06-25)

From `/hypothesise` → the two units that attack the PGx cells' one circular link (the per-allele FUNCTION
assignment is CPIC's own) + add a free internal QC. No new CLI surface (validation/evidence layer).

- **Unit A — independent functional-evidence cross-check** (`dna_decode/pgx/functional_evidence.py` +
  `scripts/pgx_functional_evidence.py` → `wiki/pgx_functional_evidence_2026-06-25.{md,json}`): per allele,
  a NON-CPIC signal + an AGREE/DISAGREE/FLAG/NO_SIGNAL verdict vs CPIC's function. Missense → Ensembl VEP
  predictors (live-fetched at curation: SIFT/PolyPhen); stop/splice → consequence class; regulatory →
  documented cis-regulatory expression effect (Sim 2006 / Rieder 2005). **6 alleles: AGREE 4 / DISAGREE 1 /
  FLAG 1.** The informative findings: **CYP2C9 *3 (I359L) DISAGREE** — predictors call it benign but CPIC =
  no-function (faithful-to-CPIC rests on clinical evidence the predictors miss); **CYP2C19 *2 FLAG** — VEP
  surface consequence is synonymous, the no-function mechanism is a documented cryptic splice. Honest tier:
  orthogonal cross-check, NOT ground truth; GTEx-eQTL confirmation deferred (didn't resolve via the v2 API).
- **Unit B — trio Mendelian co-segregation QC** (`scripts/pgx_trio_concordance.py` →
  `wiki/pgx_trio_mendelian_2026-06-25.{md,json}`): the 1000G 602 trios → **CYP2C19 602/602 + CYP2C9 602/602
  (100%) Mendelian-consistent, 0 violations** — an independent CALLING check (inheritance-physics axis,
  distinct from GeT-RM). One-pass VCF read (the 3202-column wide-line re-read footgun → fixed: 9-min stall
  → 0.9s).
- 13 new tests (`tests/test_pgx_functional_evidence.py` + `tests/test_pgx_trio.py`); 98 pgx-suite pass.
  FROZEN bacterial/viral/fungal AMR surface byte-unchanged.

## [0.6.2] — the warfarin pair: CYP2C9 (activity-score) + VKORC1 (2026-06-25)

Expands the human-PGx phase with the next two genes, reusing the CYP2C19 harness (the data-source research's
rec #1). The caller is now **gene-parameterized** (backward-compatible — CYP2C19 byte-behaviour unchanged,
72/72 GeT-RM regression held).

- **CYP2C9** (`dna_decode/pgx/cyp2c9_catalog.py`): core SNP-defined *2 (rs1799853, chr10:94942290) + *3
  (rs1057910, chr10:94981296) + *1; CPIC **ACTIVITY-SCORE** phenotype (*1=1.0/*2=0.5/*3=0.0 → AS 2=NM,
  1–1.5=IM, 0–0.5=PM). **Validated vs GeT-RM consensus on real 1000G: core diplotype 73/73 (100%)**
  (`scripts/pgx_getrm_concordance.py --gene cyp2c9` → `wiki/pgx_getrm_concordance_cyp2c9_2026-06-25.{md,json}`),
  caller independent of the consensus tools. 14/87 (16.1%) non-core residual (*5/*6/*8/*9/*11/*61 — no
  sentinel in v0; sentinel layer = v0.1, mirroring the CYP2C19 arc).
- **VKORC1** (`dna_decode/pgx/vkorc1.py`): single-SNP rs9923231 (c.-1639G>A) → warfarin sensitivity. Encodes
  the **minus-strand subtlety** (genomic chr16:31096368 C>T == cDNA G>A; genomic-T = the sensitive "A"
  allele) → G/G normal / G/A intermediate / A/A high-sensitivity (low dose). Absent record → assumed-ref,
  flagged.
- Wired into `dna-pgx --gene {cyp2c19,cyp2c9,vkorc1}` + the `dna-decode pgx` dispatch + `dna-decode list`.
  Coordinates grounded at dbSNP. 17 tests (`tests/test_pgx_cyp2c9.py`) + a CYP2C9 GeT-RM concordance test.
  FROZEN bacterial/viral/fungal AMR surface byte-unchanged.

## [0.6.1] — CYP2C19 v0.1 hardening: non-core sentinel withhold + honesty fields (2026-06-25)

Closes the safety gaps a review surfaced on the v0.6.0 cell (the *4b→*17 silent alias was a real
clinical-mis-call risk on real genomes — quantified at 0.16% *4-family / 1.37% *35 in the 1000G run):

- **Sentinel non-core layer** (`dna_decode/pgx/cyp2c19_catalog.py::SENTINELS`): rs28399504 (*4) +
  rs12769205 (*35). When a sentinel proves a non-core allele the single-SNP core proxy cannot resolve,
  the phenotype is **WITHHELD**, not mis-called. *35 fires only on an rs12769205 copy in EXCESS of the
  *2 (rs4244285) copies (so a plain *1/*2 with rs12769205 is NOT falsely withheld; the NA19122-style
  *2/*35 excess-copy case IS caught).
- **`phenotype_status` split from parse `status`** (`ok` / `phenotype_withheld` / `phase_ambiguous`);
  phenotype is `None` when withheld; `core_proxy_diplotype` always exposed. **CLI exits nonzero (3) when
  the phenotype is withheld** so automation can't consume a withheld call as valid.
- **Phase ambiguity surfaced:** ≥2 unphased het core sites whose two phase resolutions give different
  phenotypes → kept (standard trans call) but `phenotype_confidence=low` + `alternate_diplotype/phenotype`.
- **Provenance honesty fix:** the per-record `calling_is_independent_baseline=True` overclaim is replaced by
  `calling_independently_validatable` + `independent_validation_status` (`pending: faithful-to-PharmCAT
  6/6 done; GeT-RM not yet run`) + `is_core_marker_proxy=True` (NOT a full PharmVar star-allele caller).
- **`rs3758581` coordinate corrected** to chr10:**94842866** (GRCh38, verified at dbSNP; was 94852738).
- **Raises on a missing `--sample`** column (was a silent fall-back to the first sample).
- Re-validated: PharmCAT fixtures **core 6/6 hold** + the 2 non-core cases (`s1s35`, `s1s4b`) now **withheld**
  (2/2) instead of mis-called. Tests: `tests/test_pgx_cyp2c19.py` (29) + `tests/test_pgx_validate.py` (16).
  FROZEN bacterial/viral/fungal AMR surface byte-unchanged.
- **THE independent number landed — GeT-RM consensus concordance on real 1000G genomes**
  (`scripts/pgx_getrm_concordance.py` + `wiki/pgx_getrm_concordance_2026-06-25.{md,json}`): scored the caller
  vs the **GeT-RM NGS consensus** (Astrolabe+Stargazer+Aldy; Gaedigk 2022, via the ursaPGx benchmark table,
  vendored under `tests/data/pgx_getrm/`) on the **87** samples overlapping the 1000G 30× panel — genotypes
  from the public VCF (Docker bcftools), caller independent of the 3 consensus tools. **Core diplotype
  concordance 72/72 (100%)**; +7 `*38`==`*1` phenotype-equivalent (79/87 phenotype-correct); 2 non-core
  correctly **withheld** (`*4`/`*35`, incl. the real NA19122 `*2/*35`); **6/87 (6.9%) genuine non-core silent
  residual** (`*8`/`*13`/`*15`/`*39` — beyond the v0 SNP set + 2 sentinels, honestly disclosed). This is the
  strongest star-allele-CALLING validation tier available (vs the field's accepted consensus truth set). The
  per-record `independent_validation_status` is upgraded from "pending" to this achieved number. 4 tests
  (`tests/test_pgx_getrm.py`).

## [0.6.0] — the first HUMAN cell: CYP2C19 pharmacogenomics (2026-06-25)

- **NEW `dna-pgx`** (`dna_decode/pgx/`) — the first **human** decoder cell + the honest catalog-tractable form
  of the "higher organism" jump (the complex-trait/embedding path stays a closed 0-for-4 negative). A phased
  VCF (GRCh38) → CYP2C19 defining-SNP genotypes → **star-allele → diplotype → CPIC metabolizer phenotype**
  (PM/IM/NM/RM/UM). Pure-stdlib VCF parse (no pysam/Docker); handles phased/unphased, multiallelic, no-call,
  and absent-record (assumed-reference, **explicitly flagged**, never silent ref-by-absence).
- **Curated catalog** `dna_decode/pgx/cyp2c19_catalog.py` — core SNP-defined alleles **\*2 (rs4244285,
  chr10:94781859), \*3 (rs4986893, chr10:94780653), \*17 (rs12248560, chr10:94761900)** + \*1; GRCh38
  coordinates grounded vs PharmVar/dbSNP + Botton 2021 + Gaedigk 2022. CPIC standardized function +
  diplotype→phenotype table (Caudle 2020), incl. \*2/\*17 = provisional IM.
- **Honesty tier (load-bearing, two distinct claims):** star-allele **CALLING is independently validatable**
  vs the free GeT-RM consensus panel (`calling_is_independent_baseline=True`); the metabolizer **PHENOTYPE is
  FAITHFUL-TO-CPIC** (assigned, not measured; reference tool = PharmCAT;
  `phenotype_is_independent_baseline=False`). v0 covers the core SNP set; a non-core star allele is mis-called
  \*1 — a flagged blind spot. **NOT a clinical tool.**
- Wired into the `dna-decode` dispatcher (`dna-decode pgx`) + `dna-decode list`. 23 tests
  (`tests/test_pgx_cyp2c19.py`). FROZEN bacterial/viral/fungal AMR surface byte-unchanged.
- **VALIDATED** against the reference tool's own real VCF fixtures (`scripts/pgx_cyp2c19_validate.py` +
  `wiki/pgx_cyp2c19_report_card.{md,json}`): **core diplotype 6/6 + phenotype 6/6** on PharmCAT's CYP2C19
  test VCFs (`*1/*1, *1/*2, *1/*17, *2/*2, *2/*3` + a `*17`-site-missing no-call case). Honest tier =
  **FAITHFUL-TO-PHARMCAT** (in-distribution; the reference tool's expectations) — NOT yet the GeT-RM
  independent number. Two blind spots surfaced honestly: `*35` (rs12769205-defined, non-core)→ mis-called
  `*1/*1`; `*4b`→ aliases to `*1/*17` (shares the *17 SNP rs12248560 — clinically meaningful). 14 harness
  tests (`tests/test_pgx_validate.py`); fixtures vendored under `tests/data/pgx_cyp2c19/` (PharmCAT, MPL-2.0).
- **REAL-1000G run (`scripts/pgx_1000g_population.py` + `wiki/pgx_1000g_population_2026-06-25.{md,json}`):**
  ran the caller on the REAL 1000 Genomes 30× panel (**n=3202**), region fetched via Docker bcftools
  remote-querying the public VCF (the "no-htslib-on-Windows" wall dissolved — Docker provides it). Result:
  biologically-sane population distribution (NM 38.3% / IM 32.4% / RM 19.5% / PM 6.7% / UM 3.2%); blind-spot
  exposure quantified on real data (*35→*1 1.37%, *4-family 0.16%); grounded GeT-RM check NA19122 (consensus
  *2/*35 → v0 *1/*2, *35 haplotype confirmed invisible to v0). The full GeT-RM consensus concordance % is
  now a **data-access** step (labels in paper supplements), NOT a tooling wall; the harness consumes it via
  `--source getrm --expected-tsv`. Scoping: `plans/EP_PGx_CYP2C19_Cell_Scoping_2026-06-25.md`.

## [0.5.3] — salmserovar bugfix + pneumococcus AMR engines (2026-06-25)

- **FIX `dna-salmserovar` (was wrong on real genomes in 0.5.2):** `_best_per_axis` selected the H antigen by
  coverage only, but flagellin (fliC/fljB) alleles cross-hybridize at full coverage → it picked the WRONG
  antigen (S. Typhimurium LT2 gave 4:r:1,5,7 instead of 4:i:1,2). Now **identity-primary** → correct
  (LT2 → Typhimurium, all 100%). Regression test added.
- **NEW (library) `dna_decode/organism_rules/pneumo_betalactam.py`** — S. pneumoniae β-lactam PBP-type→MIC→R/S
  engine (CDC `Ref_PBPtype_MIC` lookup + context-aware breakpoints). Validated vs measured AST:
  penicillin@meningitis 0.974. Plus `organism_rules/pneumo_amr.py` (gene-presence macrolide/tet, 0.961/0.932,
  fully-independent via AMRFinder swap) + `data/pneumo_breakpoints.py`. Non-frozen overlays; FROZEN E. coli
  AMR surface byte-unchanged. (Engines are library modules; their DBs are gitignored, built via `scripts/`.)

## [0.5.2] — two new typing decoders (2026-06-25)

- **NEW `dna-salmserovar`** — Salmonella enterica serovar via the Kauffmann-White antigenic formula
  (O + H1=fliC + H2=fljB; SeqSero2-style antigen DB). Deterministic, offline-safe, faithful-to-tool.
- **NEW `dna-pneumo-serotype`** — S. pneumoniae capsular serotype via the cps-locus reference scheme
  (PneumoCaT/SeroBA-style). **INDEPENDENTLY validated vs phenotypic Quellung** (GPS Poland cohort, n=230):
  serogroup concordance 0.939 / exact 0.661 — the first independent measured-label validation for a non-AMR
  typing trait. See `wiki/pneumo_serotype_report_card.md`.
- Both wired into the `dna-decode` dispatcher + `dna-decode list`. The published 0.5.1 wheel predated them.
- Non-frozen pneumococcus AMR groundwork (library only, no new console script): `organism_rules/pneumo_amr.py`
  (gene-presence macrolide/tet rule; validated vs measured AST 0.961/0.932) + `data/pneumo_breakpoints.py`
  (context-keyed β-lactam breakpoints). FROZEN E. coli AMR surface byte-unchanged.

## [Unreleased] — Anchor-4: standing decoder-suite validation report card (2026-06-10)

- **NEW `dna_decode/eval/cohort_manifest.py`** — data-driven accession-manifest registry. `build_manifest()`
  scans EVERY `data/raw/*/selected.tsv` + EVERY `data/processed/*.parquet`; `prior_accessions(exclude_cohort
  =name)` excludes prior accessions by **EXACT-self cohort identity (not substring)**, and `Manifest.incomplete
  =True` on any load failure. REPLACES the hardcoded `_FLAGSHIP_PARQUET_COHORTS` list (3 cohorts) in
  `scripts/provenance_disjoint_validate.py` — leakage exclusion now covers all 8 parquet cohorts (744
  accessions vs ~175 before).
- **`scripts/provenance_disjoint_validate.py` now FAILS CLOSED on an incomplete manifest** (`INCOMPLETE_MANIFEST`,
  exit 2) unless `--allow-incomplete-manifest` is passed (which stamps degraded independence into the artifact).
  Artifact carries `manifest_complete` / `manifest_degraded`.
- **`scripts/ncbi_pd_provenance_census.py` now SELF-PERSISTS its powering verdict** to
  `wiki/provdisjoint_census_results.json` (group→organism normalizer; idempotent upsert per `(organism, drug)`;
  REFUSES to persist error/row-capped rows so a degraded run can't overwrite a good powering verdict).
  Previously stdout-only.
- **NEW `dna_decode/data/shipped_decoder_surface.py`** — the authoritative DEPLOYED-CLAIM surface registry
  (organism, drug, engine, organism_scope, phenotype_source_status, census_group).
- **`scripts/build_validation_report_card.py`** — rows now = shipped surface ∪ observed cells (a new decoder
  can't ship invisibly); added a 7th cell-state `LABEL_CONFOUNDED` (oxacillin×S.aureus, unreliable mecA
  surrogate); `NO_FREE_PHENOTYPE_SOURCE` is now surface-driven. Honest per-cell tier + no aggregate headline
  preserved. Writes `wiki/decoder_validation_report_card.{md,json}`. Current card: 25 cells (6 SCORED /
  4 NOT_CENSUSED / 1 UNDERPOWERED / 2 ABSTAINS_BY_DESIGN / 1 LABEL_CONFOUNDED / 11 NO_FREE_PHENOTYPE_SOURCE).

## [Unreleased] — intron-aware multi-HSP codon mapping (engine)

- **`observed_substitutions` (the shared target-site codon-mapper) is now INTRON-AWARE.** It stitches the
  query(CDS-reference)-position → subject-nucleotide map across ALL HSPs on the gene's best contig (not just
  the single best HSP). Because codon numbering is by *query* position (contiguous CDS), a codon that spans
  an exon-exon boundary — its 3 nts in two different HSPs (exons separated by an intron in the genome) — is
  still translated correctly. This generically improves **every** genome-mode caller (fungal ERG11, K13)
  for multi-exon / split-across-contigs genes, and **unblocks intron-containing targets** (pfcrt has 13
  exons; GenBank deposits are ~2471 bp genomic). Intronless genes are the single-HSP special case →
  identical prior behavior (guarded by the existing fungal/K13 tests; 967 passed, 0 regressions). Validated
  on the real (non-repetitive) 3D7 K13 CDS artificially split into exons: exon1 + deep-exon2 mutations both
  detected across the intron; a mid-codon boundary assembles to WT with no spurious call. 2 tests
  (`tests/test_intron_aware_mapping.py`). Known limit: a periodic CDS self-aligns at its period (use a real
  non-repetitive reference — real CDSs are).
- **pfcrt genome mode FLIPPED ON** (same run) — committed the 3D7 pfcrt CDS reference
  (`data/antimalarial_ref/Pf3D7_pfcrt_cds.fna`, NCBI RefSeq XM_001348968, 424aa, WT Lys@76) so
  `dna-amr --drug chloroquine --genome-fasta X.fna` works (`--pfcrt-ref` to override). **Validated on REAL
  13-exon genomic pfcrt alleles** (the 2471 bp GenBank field isolates): the intron-aware mapper recovered
  **K76T across the introns** + the full canonical CQ-R haplotype (A220S/Q271E/I356T/R371I = CVIET-type) on
  all 6 tested; WT 3D7 → S. Genome-mode ref is picked by target gene (K13 → `--k13-ref`, pfcrt →
  `--pfcrt-ref`); the intron guard is replaced by a ref-existence check. Real genomic fixture committed
  (`Pf_pfcrt_MN419894_K76T.fna`). This is the first REAL multi-exon end-to-end validation of the engine.

## [Unreleased] — antimalarial vertical: P. falciparum K13 (the 3rd kingdom, protozoan)

- **`dna_decode/data/antimalarial_amr.py`** + **`scripts/pf_kelch13_caller.py`** — extends the proven
  deterministic target-site method (bacterial AMRFinder → fungal ERG11) to *Plasmodium falciparum*
  artemisinin partial resistance via the WHO-validated **Pfkelch13 (K13)** propeller markers (C580Y, R561H,
  R539T, I543T, A675V, R622I, …). Hand-curated catalog (no AMRFinder-equivalent for Plasmodium); the caller
  REUSES the fungal caller's gene-generic `observed_substitutions` (BLAST K13-CDS-vs-genome → gap-aware
  codon-map → catalog), so the only new surface is the catalog + a thin wrapper. Artemisinin partial
  resistance is a clearance phenotype, not an MIC — the validated K13 marker IS the genotypic call; an S
  call surfaces non-K13 / partner-drug blind spots. Real 3D7 K13 reference committed
  (`data/antimalarial_ref/Pf3D7_K13_cds.fna`, NCBI XM_001350122.1, 726aa, WT Cys@580); G0-completion test
  validates C580Y numbering on the **real** reference. Offline-safe (absent BLAST+ → INDETERMINATE).
  6 tests. **3rd kingdom decoded** (bacteria → fungi → protozoa).
- **Wired into the unified `dna-amr` CLI** (`--drug artemisinin|artesunate|dihydroartemisinin`) — routes to
  the K13 engine (`--genome-fasta` real-BLAST, or `--observed K13:C580Y` wheel-only), mirroring the fungal
  productization; emits the same `amr-mechanism-call-v1` record; `--organism` default relabels to
  `Plasmodium_falciparum`. Shared `_target_site_record` + `_emit_target_site` now back both the fungal and
  antimalarial branches (extracted, not duplicated). 4 CLI tests. So `dna-amr` now spans **3 kingdoms**.
- **+ chloroquine (pfcrt K76T)** — extends the antimalarial vertical to the iconic chloroquine-resistance
  marker (`dna-amr --drug chloroquine --observed pfcrt:K76T` → R). `gene_for_drug` routes drug→target
  gene; K76T → CQ-R is unambiguous (pfmdr1 partner-drug markers deliberately omitted — their direction
  flips between amodiaquine and lumefantrine). Shipped first as `--observed`-only; **genome mode was then
  flipped on** once the intron-aware mapper + committed pfcrt CDS reference landed (see the intron-aware
  entry above) — `dna-amr --drug chloroquine --genome-fasta X.fna` now works on real 13-exon genomic pfcrt.
  5 tests.

## [Unreleased] — self-calibrating AMR rule (`calibrate_organism`)

- **`dna_decode/eval/calibrate_organism.py`** — auto-selects the per-organism AMR rule config from a
  ≥~15R/15S labeled cohort: chooses the determinant COUNTER (`qrdr_point` vs broad drug-class) and the
  count THRESHOLD by leave-one-out balanced accuracy, and auto-excludes INTRINSIC gene families (≥90% of
  both R and S, grouped at gene-family granularity so polymorphic intrinsics like blaOXA-51-family are
  caught). Returns a `CalibratedRule` (`.predict()` applies it); ABSTAINS with verdict `EXPRESSION_FLOOR`
  when no presence config clears the 0.70 LOO floor (expression-driven R that gene-presence cannot decode).
  Motivated by the wider-AMR boundary taxonomy (CONTENT/TUNING/EXPRESSION) — the counter, not just the
  threshold, is organism-specific (the Klebsiella-vs-Salmonella cipro contrast). Validated on cached cohorts
  (`wiki/calibrate_organism_validation_2026-06-08.md`): Campylobacter→1, Klebsiella→2 (+oqxAB excluded),
  Salmonella→broad@1 (deployed 0.567→1.0) all LOO 1.0; Acinetobacter + Pseudomonas meropenem → abstain.
  16 unit tests.
- **Wired into `call_resistance(..., organism=...)`** (opt-in) via a committed registry
  `dna_decode/data/calibrated_amr_rules.json` (built by `scripts/build_calibrated_registry.py`). When an
  organism is passed AND has a registry entry: a CALIBRATED entry applies its counter/threshold/intrinsic
  exclusions; an EXPRESSION_FLOOR entry returns **`prediction: "ABSTAIN"`** (refuses to predict an
  expression-driven organism×drug rather than over-call). `organism=None` (or unknown organism, or an
  explicit `resistance_threshold`) keeps the unchanged `DRUG_RULE` default — backward-compatible. Registry
  is IN-SAMPLE (N≈30) and opt-in by design; abstain verdicts are conservative. 9 wiring tests.
- **Design-review hardening (2026-06-09):** (a) `INSUFFICIENT_EVIDENCE` verdict + `MIN_CLASS_COUNT` guard —
  a one-class/under-powered cohort no longer yields a bogus `EXPRESSION_FLOOR` (fixed a degenerate
  Pseudomonas registry entry); (b) `loo_balanced_accuracy` now truly balanced (was plain accuracy), `None`
  when a class is absent; (c) deployed config is the deterministic full-cohort `_select_best_config` pick
  (removed the modal/tie-break ambiguity), LOO separately estimates the selection procedure; (d)
  `build_calibrated_registry` resolves runs via the validator `reuse_glob` (Pseudomonas under-load fix —
  all 5 entries now on valid 15R/15S); (e) promotion gate adds a **specificity floor** + min-10/class and
  treats config-match as a **flag, not a gate** (non-inferior OOS perf suffices) — all 3 cipro configs
  `promotion_eligible=True`. M2 (AMRFinder `Method`-column propagation) documented + deferred. 947 tests.

## [Unreleased] — cross-decoder analyses (concordance + profile + co-localization)

Three ANALYSES that compose the shipped decoders (no new DB) — variety roadmap Waves 1-2.

- **`dna-coloc`** (`dna-decode coloc`, Wave 2) — links each acquired resistance gene to plasmid replicon(s)
  on the SAME assembly contig → "is *this* AMR gene plasmid-borne?". Enabled by a new opt-in engine
  positions-mode (`call_alleles(..., with_positions=True)` returns each hit's subject contig+coords;
  default-off, every existing caller unchanged). Same-contig is suggestive, not proof (caveat shipped).
  Pure core + real-BLAST e2e (blaNDM-1 on plasmid-contig → plasmid-borne; sul1 on chrom-contig → not).

- **`dna-concordance`** (`dna-decode concordance`) — compares the two independent acquired-gene callers,
  AMRFinder (`dna-amr` main.tsv) vs ResFinder (`dna-resfinder` blastn), at the gene-family level (allele
  variant stripped; `sul1`≠`sul2`, `blaNDM-1`≈`blaNDM-19`) + Jaccard agreement. The cross-check `resfinder`
  was built to enable.
- **`dna-profile`** (`dna-decode profile`) — runs every assembly-FASTA decoder (pathotype + serotype +
  plasmid + resfinder) on one genome → a single unified report; each section degrades independently.
- Kept out of the `TRAITS` decoder registry (new `ANALYSES` dict; disjoint namespaces). 15 new tests.

## [Unreleased] — MLST sequence-typing decoder (the blocked one, unblocked via PubMLST)

- **`dna-mlst`** (`dna-decode mlst`) — 8th decoder. 7-gene MLST: exact per-locus allele (blastn 100/100 on
  the shared engine) → profile → Sequence Type via the PubMLST profiles table. v0: E. coli Achtman
  (adk/fumC/gyrB/icd/mdh/purA/recA, 16,242 STs). The earlier blocker was DB-discovery (CGE `mlst_db` raw
  paths 404); resolved by using **PubMLST's REST API** (`pubmlst_ecoli_achtman_seqdef` scheme 4).
  `dna-mlst --fetch-db` installs the scheme (DB on demand, gitignored). Novel/incomplete profiles report
  "ST not called", never guessed. **Validated end-to-end: K-12 MG1655 → ST10** (real DB + real genome) +
  synthetic real-BLAST e2e + pure-core tests. New shared `mlst/core.py` (profiles parse + ST lookup). 5 tests.

## [Unreleased] — DisinFinder biocide-resistance decoder (roadmap W4) + profile completion

- **`dna-disinfinder`** (`dna-decode disinfinder`) — 7th decoder. Acquired biocide/disinfectant resistance
  genes (quaternary-ammonium qac* + formaldehyde formA) via the DisinFinder DB on the shared engine (reuses
  resfinder's CGE gene parser). Hospital infection-control relevant; qac genes often share plasmids with AMR
  → pair with `dna-coloc`. Offline-safe; DB on demand. 3 tests.
- **`dna-profile` now also runs `pointfinder`** — the run-all covers all 5 assembly-FASTA decoders
  (pathotype/serotype/plasmid/resfinder/pointfinder).
- MLST (roadmap remainder) stays deferred: DB raw-paths 404 + it needs exact-allele/profile→ST semantics
  (a distinct batch, not the presence/codon pattern).

## [Unreleased] — PointFinder chromosomal point-mutation decoder (roadmap W3)

- **`dna-pointfinder`** (`dna-decode pointfinder`) — 6th decoder. Chromosomal AMR point mutations via the
  PointFinder DB: blastn each reference gene CDS vs the assembly, map the subject amino acid at each
  catalogued codon (new shared `typing/codon_map.py` — gap-aware codon→subject-AA, the proven fungal-ERG11
  pattern, now in-package), call a mutation when the subject AA matches a `Res_codon` in
  `resistens-overview.txt`. v0 scope: E. coli FQ QRDR (gyrA/parC/gyrB/parE). An INDEPENDENT point-mutation
  caller (`caller_is_independent_baseline: true`) complementing `amr` (AMRFinder POINT) + `resfinder`
  (acquired only). Epistasis (`Required_mut`) recorded, not enforced. Offline-safe; DB on demand.
  Validated on synthetic (S83L) + the real committed E. coli DB (gyrA codon83=S/87=D). 5 tests.

## [Unreleased] — typing-decoder family (plasmid + serotype + resfinder on one shared engine)

Three new deterministic curated-DB decoders — the tool grows from 2 traits to 5, all on one engine.

- **Shared engine** `dna_decode/typing/blast_caller.py` (`call_alleles`) — the generic best-HSP-per-allele
  blastn core (reuses pathotype/vf_runner's resolvers). The CGE curated-DB pattern (pathotype + plasmid)
  is now config-per-decoder, not a from-scratch build. Plasmid refactored onto it (DRY).
- **`dna-serotype`** (`dna-decode serotype`) — E. coli O:H serotyping via SerotypeFinder allele DB
  (best O-antigen + best H-antigen → `O25:H4`, partial → `O104:H?`). A genuinely new trait.
- **`dna-resfinder`** (`dna-decode resfinder`) — acquired-AMR-gene detection via ResFinder allele DB,
  per-class. Deliberately an **independent** caller vs the AMRFinder-based `dna-amr`
  (`caller_is_independent_baseline: true`) — the cross-tool concordance check the AMR decoder lacked.
- All offline-safe (status `unavailable`, exit 3); DBs downloaded on demand (not committed). ~16 new tests.

### plasmid replicon-typing decoder (earlier in this Unreleased cycle)

New deterministic trait decoder `dna-plasmid` (`dna-decode plasmid`) — the tool grows beyond AMR.

- **New capability:** plasmid incompatibility (Inc) replicon typing from a genome assembly via the curated
  PlasmidFinder allele DB + real blastn (identity 95 / coverage 60, PlasmidFinder defaults). Reports the Inc
  replicons present (IncF/IncH/IncI/IncX/IncN/…) — composing with `dna-amr`: AMR says *what* resistance,
  plasmid typing says whether it likely rides a known mobile element.
- Sibling architecture to `dna-pathotype` (curated-DB blastn caller; reuses `pathotype/vf_runner`'s blastn
  resolvers — DRY). Offline-safe: no blastn / no DB → `status: "unavailable"` (exit 3), never a crash.
  `caller_is_independent_baseline: false` (faithful to PlasmidFinder's own method, not an independent check).
- DB downloaded on demand (not committed), like the VirulenceFinder DB. 7 tests
  (`tests/test_plasmid_decoder.py`); cli-dispatch registry contract updated to the 3-decoder set.

## [0.5.0] — 2026-06-08 — Fungal AMR decoder (the kingdom jump)

`dna-amr` now decodes **fungal** azole/echinocandin resistance, not just bacterial — the determinant-scan
method validated across the bacteria→fungi kingdom boundary.

- **New capability:** `dna-amr --drug fluconazole|voriconazole|caspofungin|micafungin` routes to a
  BLAST-ERG11/FKS1 target-site engine (vs the AMRFinder engine for bacterial drugs — there is no
  AMRFinder-for-fungi). Two source modes: `--genome-fasta` (BLAST the committed C. auris ERG11 reference)
  and `--observed GENE:SUB[,...]` (pure, wheel-only, no BLAST). Emits the same `amr-mechanism-call-v1`
  record as the bacterial path (uniform tool surface); S calls surface the efflux/aneuploidy blind spots.
- **Validation (Gate G1, `wiki/fungal_ep7_g1_closeout_2026-06-08.md`):** on a de-confounded C. auris
  WGS+MIC cohort (S. Africa bloodstream, PRJNA737309 + AraPheno-style Table S1 MICs), the deterministic
  caller found the catalogued ERG11 mutation in **100% of fluconazole-MIC-R isolates across two clades**
  (clade I Y132F, clade III F126L/VF125AL) — sensitivity 1.0. Specificity is label-limited (reduced-
  susceptibility F126L carriers fall below the CDC tentative breakpoint), the documented "suspect the
  label" pattern; the genotype is the trustworthy output.
- **Infra shipped:** `dna_decode/data/fungal_amr.py` (hand-curated catalog + CDC tentative breakpoints),
  `scripts/fungal_erg11_caller.py` (BLAST→codon-map→catalog), `scripts/build_fungal_cohort.py` (cohort
  validation + within-clade de-confound + LABEL_LIMITED_FAILURE verdict), `scripts/assemble_sra_cohort.py`
  (targeted ERG11 read-mapping, ~4 min/isolate vs ~45 min full assembly). Committed real C. auris reference
  + 3 public allele fixtures (`data/fungal_ref/`).
- **Eukaryotic Path B (Arabidopsis flowering-time embedding test, Gate G2)** pre-staged + brainstorm-revised
  + CPU-only dry-manifest gate coded (`scripts/g2_dry_manifest.py`); GPU run deferred to the workhorse.
- ~43 new tests (fungal catalog/caller/cohort/CLI + dry-manifest). Bacterial path unchanged.

## [0.4.0] — 2026-06-07 — Multi-Organism AMR Decoder (capstone)

Milestone release consolidating the AMR arc. No new code — a capstone over v0.3.x.

`dna-amr` is a deterministic, interpretable AMR R/S decoder validated across **6 drugs × 4 organisms ×
4 mechanism classes, spanning the gram divide** (E. coli, K. pneumoniae, P. aeruginosa, S. aureus),
deployed as `dna-amr` / `dna-decode`. Every per-drug rule beats naive AMRFinder on independent data.

- **One engineering principle** held across every organism/mechanism: count the drug's SPECIFIC
  resistance determinants (target point-mutations / drug-specific Subclass / acquired gene-family), not
  the broad drug-class bag — intrinsic chromosomal genes (efflux) are the cross-organism gotcha.
- **Honest limits, named in every output:** `undetectable_mechanisms` (efflux/porin/regulatory expression
  phenotypes) + label-quality caveats (oxacillin → use cefoxitin). The recurring binding constraint is the
  de-confounded, reliably-labeled substrate — not the method.
- Capstone: `wiki/amr_multiorganism_capstone_2026-06-07.md`. 108 tests green.

**This is the milestone.** Further organism/drug breadth is diminishing-returns (re-confirms the same two
findings). The genuinely-different next leaps (cross-lab validation, a non-AMR sampling-independent
substrate, multimodal/eukaryotic) require resources beyond autonomous code work.

## [0.3.7] — 2026-06-07

1st Gram-positive: S. aureus oxacillin (MRSA/mecA) — genotype transfers; honest label finding.

### Added
- **oxacillin** (6th drug, 1st Gram-positive): mecA-based MRSA rule (threshold 1 + METHICILLIN-subclass,
  excludes blaZ penicillinase). Added to mic_tiers (breakpoints, classes, mec loci, primary mechanism) +
  DRUG_RULE. `supported_drugs()` now 6. +2 tests (104 → 106 → 108 green).
- **S. aureus oxacillin validation (4th organism, 1st Gram-positive):** N=30.
  `wiki/staphylococcus_aureus_oxacillin_validate_2026-06-07.md`.

### Finding (the honest result)
- **mecA genotype detection TRANSFERS to Gram-positive: sens 1.000** (all 15 R strains carry mecA). The
  acquired-gene + Subclass-refinement approach works on a Gram-positive, as on the gram-negatives.
- **spec 0.333 is oxacillin-LABEL noise, NOT a rule defect:** 10/15 "oxacillin-susceptible"-labeled strains
  carry full-length mecA — far above genuine OS-MRSA (<5%). Oxacillin direct AST is the documented unreliable
  comparator for mecA; CLSI/EUCAST recommend **cefoxitin** as the surrogate. The proper-label re-test is
  **substrate-blocked** (cefoxitin = only 3R on this NCBI dataset).
- **Terminal:** Gram-positive mecA detection generalizes; phenotype-label validation is the limit — the
  project's recurring "substrate/label is the binding constraint" lesson, now confirmed on a Gram-positive.

## [0.3.6] — 2026-06-07

3rd organism (Pseudomonas) + cross-organism shipped in the CLI.

### Added
- **`dna-amr --organism <O>`** (genome mode) — passes through to AMRFinder `-O` (organism-specific QRDR
  point-mutation detection); recorded in `provenance.amrfinder_organism`. Default Escherichia. Closes the
  gap where cross-organism support lived only in validation scripts, not the shipped CLI. +2 CLI tests.
- **Pseudomonas aeruginosa cipro VALIDATED** (3rd organism): N=30 acc 0.867 / sens 0.80 / spec 0.933
  (beats naive AMRFinder 0.767). The QRDR-POINT rule transfers UNCHANGED to a *less-similar* gram-negative
  (MexAB-OprM efflux, intrinsic AmpC — no oqxAB). 3 FN = efflux-mediated cipro-R (expected blind spot).
  `wiki/pseudomonas_aeruginosa_ciprofloxacin_validate_2026-06-07.md`.
- **`scripts/organism_drug_validate.py`** — generalized any-NCBI-organism × any-drug validator
  (auto-discovers latest PDG metadata, reuses cached runs). Every future organism is now a one-command run.

### Result
dna-amr validated across **3 organisms** (E. coli, K. pneumoniae, P. aeruginosa). The "count the drug's
specific determinants, not the broad class bag" principle holds across all three — strong evidence it is
organism-general, not E. coli-specific.

## [0.3.5] — 2026-06-07

Klebsiella full drug matrix complete — dna-amr validated 5 drugs × 2 organisms.

### Added
- **Klebsiella cef + gent + tet validated** (rules applied unchanged from E. coli):
  - ceftriaxone: acc 0.800 / sens 1.0 / spec 0.6 ✅
  - gentamicin: acc 0.867 / sens 0.867 / spec 0.867 ✅
  - tetracycline: acc 0.800 / spec 1.0 / **sens 0.600** ⚠️ PARTIAL (efflux blind spot — see below)
  - `scripts/klebsiella_drug_validate.py` (drug-agnostic; reuses cached Klebsiella runs across cohorts).
  - Consolidated: `wiki/klebsiella_drug_matrix_2026-06-07.md`.
- **`gene_prefixes` rule refinement** (`amr_rules.py`): tetracycline now counts only acquired `tet*` genes,
  excluding intrinsic K. pneumoniae OqxAB efflux (AMRFinder-tagged TETRACYCLINE but present in susceptible
  isolates). Same cross-organism principle as cipro QRDR-POINT. **Also improved E. coli tet 0.833 → 0.917.**

### Findings
- **tetracycline / Klebsiella is the honest PARTIAL:** the acquired-`tet*` rule is precise (spec 1.0) but
  sens 0.600 — 6/15 R strains are efflux-mediated (oqxAB overexpression), undetectable by ANY
  curated-determinant rule (an expression phenotype, surfaced in `undetectable_mechanisms`). A documented
  biological limit, not a rule defect.
- **Cross-organism principle confirmed 3× (cipro/tet/the gent+cef+mero Subclass refinements):** count the
  drug's specific resistance determinants, not the broad drug-class bag; intrinsic chromosomal determinants
  (efflux) are the organism-specific gotcha.

### Result
dna-amr spans **5 drugs × 2 organisms × 4 mechanism classes**; 4/5 Klebsiella drugs clear the bar
zero-tuning, all beat naive AMRFinder. +3 tet tests (104 green).

## [0.3.4] — 2026-06-07

Klebsiella meropenem — 2nd organism, NEW mechanism class (carbapenem). Phase 3 slice 2.

### Added
- **meropenem decoder** (5th drug): acquired-carbapenemase rule (threshold 1 + CARBAPENEM-subclass
  refinement — blaKPC/NDM/OXA-48). **Klebsiella N=30 acc 0.867 / sens 1.0 / spec 0.733** (vs naive
  AMRFinder 0.533; the CARBAPENEM-subclass refinement lifts spec 0.067→0.733 by excluding ESBL/AmpC that
  raise meropenem MIC without hydrolyzing it). `wiki/klebsiella_meropenem_validate_2026-06-07.md`.
  Carbapenem is the defining K. pneumoniae clinical threat — a mechanism class E. coli AMR never covered.
- meropenem added to `mic_tiers.py` (breakpoints CLSI R≥4/S≤1 + EUCAST; AMRFinder classes;
  carbapenemase loci catalog; primary mechanism). `supported_drugs()` now returns 5.
- `scripts/klebsiella_meropenem_validate.py` (reuses cached cipro-cohort AMRFinder runs on strain overlap).
- +2 tests (carbapenemase counted / ESBL excluded). 102 green.

### Honest scope
The meropenem rule is blind to porin-loss-mediated carbapenem resistance (no carbapenemase gene) — the
expected FN mode; surfaced in `undetectable_mechanisms`. 4 FP (susceptible strains carrying a carbapenemase
gene — likely low-expression / borderline MIC).

## [0.3.3] — 2026-06-07

Cross-organism: Klebsiella + the QRDR-POINT cipro rule (roadmap Phase 3, slice 1).

### Added
- **Cross-organism transfer (Klebsiella pneumoniae cipro):** N=30 NCBI cohort, **acc 1.000** with the
  deployed rule. The method generalizes across organisms. `scripts/klebsiella_cipro_transfer.py`,
  `wiki/klebsiella_cipro_transfer_2026-06-07.md`. `_run_amrfinder` gained an `organism` param
  (`-O Klebsiella_pneumoniae`).
- `qrdr_point_count` / `qrdr_point_determinants` in `amr_rules.py` — count fluoroquinolone QRDR
  target-alteration POINT mutations (gyrA/parC/parE) only.

### Changed (ratified — changes the deployed cipro number)
- **cipro `DRUG_RULE` switched to the QRDR-POINT rule globally** (`counter='qrdr_point'`): count QRDR
  target POINT mutations ≥2, not the broad QUINOLONE-class determinant bag. Rationale: the broad rule
  FAILED on Klebsiella (acc 0.5 — intrinsic chromosomal OqxAB efflux, absent in E. coli, saturates the
  count). The canonical target-mutation count is cross-organism-robust. Net effect:
  - Klebsiella cipro 0.5 → **1.000**
  - E. coli cipro in-cohort 0.939 → 0.925 (−1.4pp; dropped cases were qnr/efflux-mediated)
  - E. coli cipro **cross-source (NCBI) 0.955 → 1.000** (+4.5pp — QRDR-POINT generalizes better; the
    in-cohort −1.4pp was tuning-cohort overfit).
- **Platform finding:** cross-organism transfer requires counting the drug's TARGET-alteration mutations,
  not the broad drug-class bag; intrinsic chromosomal determinants are the organism-specific gotcha.
- Tests: cipro tests updated to QRDR-POINT (point-mutation rows); +4 new (qrdr helpers + intrinsic-exclusion);
  cohort regression re-pinned to 0.925/0.875/0.973. 24 total green.

## [0.3.2] — 2026-06-06

Trust-hardening + honesty taxonomy (adversarial-review follow-through). No new science — makes the
shipped decoders auditable + honest about their blind spots.

### Added
- **Blind-spot honesty:** every SUSCEPTIBLE `dna-amr` call now carries `undetectable_mechanisms`
  (`efflux` / `porin_loss` / `regulatory` — expression/regulatory resistance absent from AMRFinder's
  curated determinants) + a caveat that a negative means "no curated determinant found", not "definitely
  susceptible". `UNDETECTABLE_MECHANISMS` in `amr_rules.py`.
- **Discordance taxonomy:** `discordance_bucket(prediction, true_label)` + `evaluate_cohort` now emit a
  failure-mode breakdown — `FN_undetected_mechanism` (R missed → the blind spots) vs
  `FP_determinant_without_phenotype` (called R but susceptible → label noise / silent-or-low-expression /
  borderline MIC). The "failure-tolerant tool" deliverable: names where it fails.
- **Provenance pin:** output JSON `provenance.amrfinder_image` records the pinned AMRFinderPlus image
  (`ncbi/amr:4.2.7-2026-03-24.1`; tag encodes the DB date) so an R/S verdict is reproducible against a
  known determinant source. `AMRFINDER_IMAGE_PINNED` in `amr_rules.py`.
- **Value-add headline** in `wiki/dna_amr_multidrug_validation_2026-06-06.md`: explicit naive-AMRFinder
  vs dna-amr table (cef +0.28 acc/+0.50 spec; gent +0.43/+0.57) — proves the per-drug policy adds value
  over vanilla "any determinant → R", not just re-prints AMRFinder hits.
- 5 new tests (blind-spots on S/R calls, discordance taxonomy, cohort discordance breakdown). 20 total.

### Fixed
- Stale gentamicin "NOT yet cohort-validated" claims (`amr_rules.py` docstring + wiki caveat) reconciled
  with the N=128 acc 0.945 validation that DRUG_RULE/README/CHANGELOG already recorded (claim-hygiene).

### Validated (cross-source — closes the same-database gap)
- **Independent NCBI Pathogen Detection validation** (`scripts/xsource_validation.py`,
  `wiki/dna_amr_xsource_validation_2026-06-07.md`): 22 E. coli, balanced ~11R/11S per drug, **zero
  accession overlap** with the 176 BV-BRC cohort accessions (enforced at selection). Result: cipro 0.955,
  cef 0.864, gent 1.000, tet 0.909 — comparable to / better than in-cohort. **Answers the product
  question:** dna-amr beats naive "any-determinant→R" AMRFinder on UN-tuned data by +0.09 (cipro) /
  +0.23 (cef) / +0.41 (gent) / +0.0 (tet) accuracy — the per-drug threshold + Subclass refinement IS the
  value-add, not determinant discovery. Sensitivity 1.0 on all 4 drugs (FN=0); all errors are FP
  (determinant-present-but-susceptible). Honest scope: NCBI Pathogen Detection is a different
  source/curation, not a controlled different-lab study.

## [0.3.1] — 2026-06-06

The "one coherent tool" consolidation (after the v0.3.0 build settled the embedding question).

### Added
- **Unified `dna-decode` console entry** (`dna_decode/cli.py`) — single tool that dispatches to the
  trait decoders: `dna-decode amr …`, `dna-decode pathotype …`, `dna-decode list` (capability +
  validation status), `dna-decode --version`. Thin pass-through (argv delegated verbatim); the
  per-decoder entries (`dna-amr`, `dna-pathotype`) stay independently usable + unchanged. 6 dispatch tests.

### Changed
- Project ledger (`project_state/dna-decode-2026-05-11.md`) updated to record the strategic inflection:
  embedding frontier closed (0-for-3), deterministic mechanism-feature decoders are the product.

## [0.3.0] — 2026-06-06

The "deterministic decoders win" release. The frozen-genome-embedding (NT mean-pool) thesis was tested
to a decisive conclusion and the project committed to deterministic, interpretable mechanism-feature
decoders as the shipping path.

### Added
- **Multi-drug deterministic AMR caller** (`dna-amr`). Extended from cipro-only to **all 4 drugs**, with
  per-drug validated rules baked into `dna_decode/eval/amr_rules.py::DRUG_RULE`:
  - ciprofloxacin — threshold 2 (QRDR point-mutations) — N=147 acc 0.939.
  - ceftriaxone — threshold 1 + **extended-spectrum subclass refinement** (CEPHALOSPORIN/CARBAPENEM;
    excludes intrinsic blaTEM-1/blaEC that are ampicillin-R not ceftriaxone-R) — N=60 acc 0.933.
  - tetracycline — threshold 1 (acquired tet genes) — N=12 acc 0.833 (small N, provisional).
  - gentamicin — threshold 1 + **GENTAMICIN-subclass refinement** (excludes aph/aadA
    streptomycin-kanamycin genes that don't confer gentamicin-R) — N=128 acc 0.945.
  - General fix: a broad AMR class over-calls (cef spec 0.41, gent spec 0.39) because it counts genes for
    OTHER class members; AMRFinder's Subclass field is the drug-specific discriminator. One-line refinement.
  - `call_resistance(tsv, drug)` now auto-selects the per-drug rule; explicit threshold still overrides.
  - Validation: `wiki/dna_amr_multidrug_validation_2026-06-06.md`.
- **De-confound gate** (`dna_decode/eval/cohort_deconfound.py`) — within-lineage label-contrast
  precondition (3-state DE_CONFOUNDED/WARN/CONFOUNDED + promotability) for any embedding-vs-classical
  falsifier. The reusable guard against study==class confounding.
- **AMR embedding falsifier** (`scripts/amr_falsifier.py`) + **QRDR-POINT knowledge baseline**
  (`dna_decode/eval/point_baseline.py`) + within-lineage diagnostic.
- **dna-amr external validation** — held-out N=29 acc 0.862 / sens 0.882 / spec 0.833
  (`wiki/dna_amr_external_validation_2026-06-05.md`).
- **EP-6 carbon-utilization substrate infra** — `dna_decode/data/bacdive.py` loader +
  `scripts/bacdive_carbon_util_feasibility.py` census + `scripts/bacdive_li2023_to_long.py` adapter.

### Findings (recorded, not code)
- **AMR embedding thesis FALSIFIED on the cleanest substrate.** Cipro N=147 (de-confounded): NT-XGBoost
  0.914 beats k-mer (+8.9 pp) but **LOSES to QRDR-POINT 0.943**; NT within-lineage concordance = chance
  (0.605, p=0.365) → it learned lineage, not mechanism. Decision: `plans/AMR_embedding_niche_decision_2026-06-05.md`.
- **Carbon-utilization (EP-6) E. coli-INFEASIBLE.** Data acquired (Li et al. 2023 OSF jwkr7); E. coli
  slice = 27 strains, 0 carbon sources clear the ≥100 floor. The embedding-niche test needs a THIRD
  requirement beyond sampling-independent-label + no-catalog: **organism-specific depth at scale**.
  `wiki/bacdive_carbon_util_feasibility_2026-06-06.md`.

### Changed
- `dna-amr` CLI `--resistance-threshold` now defaults to the per-drug validated config (was hard-coded 2).
- README decoder table + `plans/Trait_Decoding_Roadmap.md` Phase 2/4 updated.

## [0.2.0] — 2026-06-05

- In-package `dna-amr` console entry (deterministic AMR mechanism caller, cipro-validated).
- `dna-pathotype` console entry (deterministic VirulenceFinder-marker pathotype resolver + abstention +
  canonical-VF diff; ExPEC recall 0.917).
- Packaging: both decoders ship in the wheel (`[project.scripts]`).

## [pathotype-v0] / [phase-1-shipped] — 2026-05/06

- Phase-1 closeout: NT-frozen-pooling characterized (passes concentrated-signal mechanisms, fails
  distributed). v0 cipro decoder + pathotype v0 resolver shipped. See git history + `wiki/`.

## [expression-context-v0] — 2026-06-10

- `dna_decode/eval/expression_context.py` — deterministic IS-element-upstream-of-target detector (blastn ISAba1+OXA-51 vs assembly; all-hits, no truncation; same-contig + strand-aware upstream proximity; offline-safe). Reads regulatory CONTEXT, not gene presence, to cross the EXPRESSION abstain floor.
- `amr_rules.call_resistance` gains an optional `genome_fasta` + a GATED EXPRESSION_FLOOR ABSTAIN->R override (registry `expression_context.enabled`; ships **off/experimental** — opt-in only, never default-on; default decoder behavior unchanged).
- Independent-cohort validation (15R/15S disjoint from the cached 30): **HOLD — UNDERPOWERED, not a falsification**. A mechanism audit found 14/15 R carry strong acquired carbapenemases (non-target); only 1 strain is intrinsic-only-R (the signal's target subset), detector-negative. The validator was stratified by acquired-carbapenemase and now gates on `target_R_rescues>=1 AND n_target_R>=10` — the cohort cannot test the signal (intrinsic-only carbapenem-R Acinetobacter is rare; acquired OXA-23 dominates). Override remains disabled (opt-in/off). 18 new tests, 0 regressions.
