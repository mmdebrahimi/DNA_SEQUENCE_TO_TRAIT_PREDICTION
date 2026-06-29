# Human gene→phenotype data sources online — scan for the off-pathogen decoder (2026-06-28)

Extends the human-data branch of `wiki/wide_net_gene_phenotype_source_scan_2026-06-28.md`. Triggered by
"what other human-related data is available online" while the OpenSNP archive (eye-colour flagship) downloads.
**Filter = our paradigm:** a deterministic curated RULE + a FREE, INDEPENDENT, INDIVIDUAL-LEVEL phenotype label
(summary statistics and gene→trait association tables do NOT validate a per-genome decoder — they only source
the rule). Sorted by how directly each clears that bar.

## Tier A — free + individual-level + rich phenotype + NO application (the validation substrates)

These are the genuinely usable validation sources — the only two classes of human data we can score a decoder on
today without a Data Use Agreement.

| Source | What | Phenotype depth | Access | Verdict for us |
|---|---|---|---|---|
| **OpenSNP (Internet Archive mirror)** | DTC genotypes (23andMe/Ancestry) + self-reported traits | Eye/hair colour, height, handedness, lactose, many self-reports | Archive-only (live site **deleted 2025-04**; founders cited weaponization fear); single 21GB 2017 zip | **In flight** — the eye-colour flagship. Noisy self-report label = near-independent, non-circular. Ethically *withdrawn* by its creators (Care flag). |
| **Personal Genome Project (PGP)** | Open-consent WHOLE GENOMES + trait/health surveys; 5 centres (US/UK/Canada/Austria/China) | WGS + genome reports + self-reported age/sex/smoking/traits; PGP-UK adds methylome | **Fully open, no DUA, designed for public redistribution** (personalgenomes.org / .org.uk) | **THE ethically-clean OpenSNP successor.** Participants *consented* to public release — the clean path the eye_colour memo asked for. Caveat: voluntary phenotype reporting → sparse/uneven labels (P-index quantifies completeness). |

**Strategic read:** PGP is the recommended next human substrate after the OpenSNP flagship lands. Same caller
(`eye_colour.call_eye_colour` is format-agnostic on any DTC/VCF file), an ethically-uncomplicated source, and it
adds WGS → enables the IrisPlex 6-SNP v0.1 (OpenSNP DTC chips may not type all 6 SNPs; WGS does).

## Tier B — free + individual-level but THIN phenotype (confound control, not trait labels)

| Source | What | Why it matters to us |
|---|---|---|
| **1000 Genomes / IGSR** | ~2,500 individuals, 80M+ phased variants; open on AWS Open Data + FTP | Phenotype = population/ancestry/sex/family ONLY. **Not a trait substrate** — but it IS the **ancestry-confound control** the eye-colour v0.1 needs (rs12913832 is European-calibrated; a within-/across-ancestry split uses these allele frequencies). |
| **HGDP (Human Genome Diversity Project)** | 948 deep WGS across broader geography (Middle East, sub-Saharan Africa, Oceania) | Same role as 1000G; wider diversity → stronger confound stress-test. A harmonized 1kGP+HGDP public release exists (2023). |

## Tier C — curated RULE sources (give the rule, NOT the validation labels)

The human analogues of CARD/AMRFinder — they tell us *which variants → which trait*, the deterministic rule.
They do NOT contain per-individual phenotype labels, so they can't score the decoder; pair with Tier A/B.

| Source | Rule it provides | Decoder cell it enables |
|---|---|---|
| **HIrisPlex-S** | Deployed forensic eye/hair/skin-colour SNP model (incl. rs12913832) | eye/hair/skin colour (our flagship rule's upgrade path) |
| **ClinVar** | Variant → clinical significance (Mendelian disease), curated, free | Mendelian disease-risk cells (the disease analogue of AMR) |
| **OMIM** | Gene → Mendelian phenotype catalog | single-gene trait/disease rules |
| **PharmGKB + CPIC** | Star-allele → drug-response guidelines (PGx) | **pharmacogenomics** — note: the repo already has PGx wiki artifacts (`pgx_*`); a strong, deployed-rule cell |
| **GWAS Catalog (EBI)** | SNP → trait associations (SUMMARY STATS only) | rule-sourcing for polygenic traits — but polygenic = dead-embedding-arm territory (G-ARCH filter); prefer oligogenic |

## Tier D — rich phenotype but GATED (external wall: DUA / application / institution-gated)

Deep genotype+phenotype exists at scale, but behind a Data Use Agreement = the same "external wall" class as the
non-public AMR label sources. Named so we don't mistake them for free:

- **UK Biobank** — 500k individuals, deep phenotype + genotype/WGS. Application + fee + institutional approval.
- **All of Us (NIH)** — ~1M US participants, deep EHR phenotype. Registered-tier application; US-gated.
- **dbGaP** — thousands of US study cohorts, individual-level. Per-study controlled-access application.
- **GTEx** — genotype + tissue expression (eQTL). Genotypes are controlled-access (dbGaP).
- **EGA (European Genome-Phenome Archive)** — controlled-access cohorts; per-dataset DAC approval.

Conversion note: a single institutional DUA (e.g. UK Biobank) would convert this whole tier from external→code
wall — a USER acquisition decision (mirrors the reproducibility-freeze forward-path #1), not an executor action.

## Recommendation (best-judgement)

1. **Finish the OpenSNP eye-colour flagship** (downloading) — first real off-pathogen number.
2. **Then PGP** as the ethically-clean, free, individual-level human substrate — re-run the SAME eye-colour caller
   on PGP genome reports/VCFs (Tier A), and use 1000G/HGDP (Tier B) for the ancestry-confound split that upgrades
   eye-colour v0 → v0.1.
3. **Disease/PGx cells** (Tier C rule + Tier A label) are the natural human expansion beyond visible traits —
   ClinVar/OMIM Mendelian rules scored on PGP genomes, or the existing PGx rules. All single-/oligo-genic
   (deterministic-decoder-shaped); polygenic traits stay out (closed embedding arm).

## Sources
- PGP open-access multi-omics: [Nature Sci Data 2019](https://www.nature.com/articles/s41597-019-0205-4),
  [bioRxiv](https://www.biorxiv.org/content/10.1101/566711v1.full), portals personalgenomes.org / .org.uk
- P-index (phenotype completeness): [PMC5384003](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5384003/)
- 1000 Genomes open access: [IGSR data-access](https://www.internationalgenome.org/category/data-access/),
  [AWS Open Data](https://registry.opendata.aws/1000-genomes/)
- OpenSNP shutdown (Tier A caveat): see `wiki/eye_colour_cell_result_2026-06-28.md`
