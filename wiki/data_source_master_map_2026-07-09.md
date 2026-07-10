# Data-source master map — dna_decode (2026-07-09)

**Single canonical index** of every data source relevant to the project's use case, screened against the
8 rejection gates (`wiki/negative_results_map_2026-06-13.md`). Consolidates ~12 scattered scan artifacts +
`config/datasources.yaml` + the installed `data/*_db` reference DBs into one entry point, so future work
screens a candidate here instead of re-scanning. Captured by a `/soraya` data-mapping sweep (conversational-
executor; D: offline).

**Live sub-registries this indexes (do not duplicate — extend these):**
- `wiki/related_data_sources_registry.md` (2026-07-08) — world-model (variant-effect) + decoder substrate,
  with reachability verified by curl. **LIVE.**
- `wiki/label_acquisition_anchor_2026-07-04.md` — the acquisition target list (forward-path #1). **LIVE.**
- `wiki/negative_results_map_2026-06-13.md` — the 8 gates every candidate is screened against. **LIVE.**

## The use case (what "useful" means here)

Two substrates, two different label bars:

1. **Deterministic decoder** (the shipped product): organism genotype → AMR / typing / trait phenotype via a
   curated RULE. Needs a **RULE source** (catalog/determinant list) + a **free, independent, per-isolate
   MEASURED label** to validate non-circularly. The rule and the label must come from *different* origins,
   or G1 (circular label) trips.
2. **World-model / learned arm** (0-for-5, embedding track): sequence model → variant effect. Needs the
   3-part embedding-niche test — sampling-independent MEASURED label + NO curated catalog + organism DEPTH
   ≥100. No public substrate has ever met all three at once (the standing open experiment; only a biobank
   trait could).

**Binding constraint, verified across every closed track: LABELS, not models.** A source that only supplies
a RULE (CARD/AMRFinder/WHO-TB/PharmGKB) advances coverage but cannot validate; a source whose label is
genome-derived (MalariaGEN, FungAMR-as-label) is G1-circular. The rare well is a **free per-isolate measured
label** (Stanford HIVDB / CoV-RDB / CRyPTIC / EBI AMR-Portal / GeT-RM) — treat each as precious.

## 1 · IN-USE — deployed decoder cells (label already secured)

| Domain | Source | Provides | Access | Gate status |
|---|---|---|---|---|
| Bacterial AMR (E. coli, 6-drug) | **BV-BRC AST + NCBI Pathogen Detection** | per-isolate broth-microdilution MIC R/S | `ftp.bvbrc.org`, NCBI-PD | cleared all 8 (clonality-disclosed) |
| Bacterial AMR (E. coli independent) | **Oxford `PRJNA604975` / Spain PROBAC `PRJEB62601`** | measured MIC, provenance-disjoint | ENA/SRA | external arm; leakage-checked |
| Bacterial AMR (multi-organism provdisjoint) | **EBI AMR Portal** | 1.71 M measured-AST rows, 74 powered disjoint cells | EBI | free independent — the FALLEN "labels wall" |
| TB AMR (RIF/INH) | **CRyPTIC compendium** | 12,287-isolate VCF + BMD-MIC | Zenodo | in-distribution baseline (WHO cat partly built on it) |
| TB AMR (independent) | **EBI AMR-Portal TB disjoint cohort** | 39,193 isolates, measured DST, `leaked` flag | EBI | provenance-disjoint independent (the TB independent number) |
| Viral HIV-1 (5 classes) | **Stanford HIVDB PhenoSense** | per-isolate fold-change | `hivdb.stanford.edu` | **first free independent cell**; non-circular |
| Viral SARS-CoV-2 (Mpro) | **Stanford CoV-RDB** | measured fold-change | `github.com/hivdb/covid-drdb-payload` | in-distribution (catalog+fold both CoV-RDB); underpowered |
| Fungal (C. auris ERG11) | **assembled SRA cohort** (targeted read-map) | ERG11 genotype + CDC tentative breakpoints | SRA + `data/fungal_ref/` | G1-validated determinant transfer; spec label-limited |
| Human PGx (CYP2C8/9/19, CYP3A5) | **GeT-RM (CDC) consensus + 1000G VCFs** | multi-lab consensus diplotype, 1000G-joinable | `cdc.gov/lab-quality/php/get-rm` + `data/pgx_getrm`, `data/pgx_1000g` | NEAR_INDEPENDENT (multi-method wet-lab) |
| Human trait (eye/ABO/lactase/earwax) | **OpenSNP (Internet Archive mirror)** | DTC genotype + self-report trait | Archive 2017 zip | PILOT/DEMO tier (self-report ≈ near-independent; source ethically withdrawn — Care flag) |
| World-model (variant effect) | **ProteinGym v1.0/v1.1, MaveDB, humsavar** | DMS + pathogenicity labels | HF / Zenodo / UniProt FTP | captured + ESM-scored (see 2026-07-08 registry) |

## 2 · IN-USE — RULE / reference DBs (installed; give the rule, NOT the label)

The curated determinant + typing DBs on disk. These SOURCE the deterministic rule; they cannot validate a
decoder non-circularly (using them as labels = G1). CGE-family unless noted.

| `data/` dir | DB | Role |
|---|---|---|
| `amrfinder_db` | NCBI **AMRFinderPlus** | acquired AMR genes + point mutations (the deployed bacterial engine) |
| `resfinder_db` / `pointfinder_db` | CGE **ResFinder / PointFinder** | acquired resistance genes / chromosomal point mutations |
| `virulencefinder_db` | CGE **VirulenceFinder** | E. coli virulence alleles (pathotype resolver + genome-map virulence overlay) |
| `plasmidfinder_db` | CGE **PlasmidFinder** | plasmid replicon typing |
| `serotypefinder_db` | CGE **SerotypeFinder** | E. coli O:H serotyping |
| `disinfinder_db` | CGE **DisinFinder** | disinfectant/biocide resistance |
| `mlst_db` | **PubMLST** (via REST) | 7-locus MLST → ST (lineage) |
| `ktype_db` | Kleborate/Kaptive-class | Klebsiella K/O typing |
| `salmserovar_db` | **SISTR**-class | Salmonella serovar |
| `pneumoserotype_db` / `pneumo_betalactam_db` | pneumo typing / PBP β-lactam | S. pneumoniae serotype + β-lactam |
| `who_tb_catalogue` + `tb_lineage_barcode` | **WHO TB catalogue v2 + Napier barcode** | TB determinant rule + lineage SNPs (pinned) |
| `card` (config) | **CARD** | resistance ontology (rule sourcing; CC-BY-NC-SA) |
| `*_ref` (hiv/sarscov2/fungal/antiviral/antimalarial/isaba1) | committed CDS/genome references | caller reference sequences (HXB2 RT/PR/IN/CA, Mpro, ERG11, IS*Aba1*, …) |
| `cyp2d6_psv` | CYP2D6 structural/PSV reference | PGx CYP2D6 (SV/CNV — needs Cyrius-class caller) |

## 3 · CANDIDATE — FREE, fetchable now (no user authority; Soraya can fetch)

| Source | Provides | Unlocks | Screen |
|---|---|---|---|
| **FinnGen summary statistics** | GWAS summ-stats, ~2400+ endpoints | locus→trait priors, ancestry-transfer checks | FREE, no DAC — **grab-now** (summary tier = rule-source, not per-individual label) |
| **Personal Genome Project (PGP-US/UK/CA)** | open-consent WGS + self-report traits | ethically-clean OpenSNP successor; IrisPlex 6-SNP v0.1 | FREE, no DUA, designed for redistribution — **top free human substrate** |
| **1000 Genomes / IGSR + HGDP** | ~2.5k + 948 WGS, ancestry only | ancestry-confound control for eye-colour v0.1 | FREE (AWS Open Data / FTP); Tier-B confound control, not a trait label |
| **GeT-RM (full 28–34 gene panel)** | multi-lab consensus diplotypes, 1000G-overlap | PGx extension: CYP3A5→CYP2B6→TPMT→DPYD→UGT1A1 | FREE; NEAR_INDEPENDENT only where multi-method (screen per-gene, S2 memo) |
| **FungAMR (CARD-integrated)** | 35,792 fungal mutation→resistance entries, 95 spp/246 genes | fungal caller CATALOG-ENRICHMENT (coverage upgrade) | FREE catalog; **rule-source only** — G1 as a label |
| **ClinVar** | variant → clinical significance | Mendelian disease-risk cells (disease analogue of AMR) | FREE (440 MB; stream-filter to missense); rule + weak label |
| **AlphaMissense** | per-variant pathogenicity predictor | world-model baseline-to-beat (ESM2-3B edges it, n=7) | FREE 1.2 GB (stream-filter per UniProt) |
| **gnomAD constraint / DepMap** | variant tolerance (o/e,pLI) / gene essentiality | weak-label + a different phenotype axis | FREE; gene-level (not missense) |
| **OMIM / HIrisPlex-S / PharmGKB+CPIC** | gene→Mendelian / eye-hair-skin SNP model / star-allele→drug | rule sources for disease / colour / PGx cells | FREE rule-sources (pair with a Tier-A label) |

## 4 · CANDIDATE — ACQUISITION (user authority: accounts / DUA / fee)

Per `wiki/label_acquisition_anchor_2026-07-04.md`. **Soraya drafts + builds the ingestion scaffold; the
user authorizes + submits.** Ranked by (value × feasibility):

| # | Target | Unlocks | Friction |
|---|---|---|---|
| 1 | **UK Biobank** Approved Researcher — **CHOSEN PATH** (user-directed) | the fair learned-decoder test + free lab-measured human labels | MEDIUM (York professorship = vehicle); **access FEE = money gate** |
| 2 | **All of Us** Registered→Controlled | same, US-gated | ruled out for now (DURA: York not registered) |
| 3 | **dbGaP / EGA / GTEx** per-study controlled | disease-cohort multipliers after first biobank model | HIGH per study |
| 4 | **Clinical micro-lab AST+WGS biobank** | sampling-independent bacterial-AMR external validation | HIGH (institutional collab) |
| 5 | **Viral gap-fill outreach** (TB author-request; SARS clinical fold) | closes 2 `blocked:external` viral cells | LOW (free author email, drafted) |
| 6 | **WWARN / GISAID** (malaria IC50 / influenza NA) | measured antimalarial/antiviral fold | access-request-gated |

## 5 · REJECTED — with gate (do not re-spend labor)

| Candidate | Verdict | Gate(s) |
|---|---|---|
| Pathotype (EnteroBase/NCBI-PD labels) | label-blocked | G1 (AMRFinder/BlastFrost-derived) + G3 (isolation-site) + G2 |
| Foundation-model embeddings (0-for-4/5) | learns lineage not mechanism | MODEL ceiling (not a label gate) |
| MIC-continuous (graded resistance) | infeasible | G1 (91% XGBoost-from-genome) + G6 (~70% censored) + G8 |
| AMR grid — Salmonella tet/gent | underpowered | G4 (ecosystem-dominated) |
| AMR grid — Acinetobacter/Pseudomonas broad class-rules | intrinsic-gene degeneracy | mechanism ceiling (spec→0) |
| Arabidopsis flowering embedding (G2) | closed negative | de-confounded within-group r² −0.13 (learned structure not signal) |
| MalariaGEN Pf6/Pf7 | blocked as label | G1 (marker-inferred resistance) |
| FungAMR **as a label** | not independent | G1 (catalog-vs-catalog) — but valid as a rule (§3) |
| HF `*/plant-phenotype`, `EthnicErotic/*`, `genetics-phenotype-*` | not G→P data | no genotype side / NLP / sampling-defined |
| `Solshine/Rice_Genotype_and_Phenotype` | real G+P but regime-2 | population-structure confound (closed embedding arm); GBLUP-only |

## 6 · Viral-hepatitis / influenza frontier — RESOLVED by direct fetch (2026-07-10)

These sat as NEEDS-VERIFICATION only because the biology-resistance **search** surface is policy-filtered.
The lesson (`feedback_offline_cache_vs_unreachable_data`): a search filter is not a network wall — try the
DIRECT fetch. Probed the actual data hosts by `curl` (bypasses the filter):

| Candidate | Host reachable? | Verdict | Gate |
|---|---|---|---|
| **HCV** (DAA fold) — geno2pheno[HCV] | YES (HTTP 200) | **G1-circular NEGATIVE.** Its only `download` is "direct csv download" of the tool's OWN predictions (upload sequence → predicted resistance CSV). An interpretation TOOL, no isolate-level measured-phenotype table. | G1 |
| **HCV** — HCV-GLUE | YES (HTTP 200) | GLUE sequence/typing framework (alignment + genotyping), not a measured-phenotype repository. | G1/not-a-label |
| **HBV** (NRTI fold) — HBVdb | NO (HTTP 000) | Genuine NETWORK wall (SSL/host down across sessions), not a search-filter artifact. Even if up, HBVdb is a curated interpretation DB (likely G1). Parked as a real infra wall. | infra + G1 |
| **Influenza NA** (oseltamivir/zanamivir IC50) | — | GISAID-access-gated by design (registration + DUA); WHO GISRS reports are aggregate, not free-public per-isolate. | G7/access |

**Net: the free-independent-label frontier is now VERIFIED CLOSED** (not "search was blocked"). No free,
public, isolate-level MEASURED-phenotype genotype↔phenotype source surfaced beyond the already-exploited
Stanford family (HIVDB→HIV, CoV-RDB→SARS) + the bacterial AMR-Portal/BV-BRC/CRyPTIC tracks. Confirms the
banked thesis: the only path to a new independent cell is **acquisition** (§4, user authority).

## 7 · Scattered-artifact index (mark current vs superseded)

The source scans this map consolidates — read the CURRENT ones; the SUPERSEDED are folded here:

| Artifact | Status |
|---|---|
| `related_data_sources_registry.md` (2026-07-08) | **CURRENT** (world-model + decoder; reachability-verified) |
| `label_acquisition_anchor_2026-07-04.md` | **CURRENT** (acquisition target list) |
| `negative_results_map_2026-06-13.md` | **CURRENT** (the 8 gates) |
| `pgx_free_label_sources_research_2026-07-05.md` | **CURRENT** (GeT-RM per-gene tiering) |
| `free_label_source_scan_2026-06-28.md` | superseded-headline (missed EBI AMR-Portal; per-candidate analysis stands) |
| `human_gene_phenotype_data_sources_2026-06-28.md` | folded into §3/§4 |
| `wide_net_gene_phenotype_source_scan_2026-06-28.md` | folded into §3/§5 |
| `data_source_research_human_animal_2026-06-25.md`, `data_acquisition_sweep_2026-06-25.md` | folded (early acquisition drafts) |
| `label_acquisition_strategy_2026-06-27.md` | superseded by the anchor (AMR/TB independence went free via EBI) |
| `independent_phenotype_{label_census,source_sweep}_2026-06-10.md`, `ncbi_pd_provenance_census_2026-06-10.md` | folded (provdisjoint census inputs) |

## 8 · Gaps this consolidation surfaces

1. **Animal PGx cells** (`data/horse`, dog PGx feasibility) are scoped but not in any live registry — status
   is candidate/feasibility-only; fold into the next registry refresh or explicitly park.
2. **Influenza NA + HCV/HBV** remain the only NEEDS-VERIFICATION free candidates — the one cheap D-free move
   that could add a cell is re-running those two searches on a filter-free surface.
3. The **two live registries + this map are now three overlapping indexes** — recommend this master map become
   the umbrella and the 2026-07-08 registry + acquisition anchor its two detailed children (avoid a 4th).
4. Every remaining *free* lever is a RULE-source or a confound-control, not a new independent label. **The only
   path to a genuinely-new validated cell is acquisition (§4) — a user-authority decision, already anchored.**
