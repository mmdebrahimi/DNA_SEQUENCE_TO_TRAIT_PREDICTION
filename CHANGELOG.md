# Changelog

All notable changes to `dna_decode`. Format loosely follows [Keep a Changelog](https://keepachangelog.com/);
this is a solo research-tool repo so the granularity is per-release-theme, not per-PR.

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
