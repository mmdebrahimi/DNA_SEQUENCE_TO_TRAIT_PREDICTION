# Data-source scouting — outside-the-box G2P label wells (2026-07-20)

**Tier: DISCOVERY (unverified).** Findings from 5 parallel web-research agents; label-provenance audited
but NOT independently confirmed by fetching/joining the data. Do NOT treat any row as a validated cohort
until its genome↔phenotype join is checked against the 8 rejection gates
(`wiki/negative_results_map_2026-06-13.md`) on real records.

## The binding pattern (re-confirmed)

The free + independent + **isolate/variant-level + WET-LAB-MEASURED** + provenance-separable quadrant is
narrow, and every source that clears it has the SAME shape: **national reference-lab AST programs**
(NARMS / GISP / CDC AR-Bank / CRyPTIC), **Stanford-style curated measured-fold-change DBs** (HIVDB /
CoV-RDB — the "labels not models" wins), and **DMS aggregators** (MaveDB / ProteinGym). Everything
tool-derived (G1), sampling-defined (G3), aggregate-only, or DUA-gated is out. Two genuinely
**outside-the-box** axes emerged: **phage susceptibility** (a new NON-AMR, non-circular phenotype) and
**MaveDB** (the deep well for the molecular forward/inverse cell).

## Tier 1 — free, isolate/variant-level, wet-lab measured, joinable (clears all 8 gates on paper)

| Source | Organism / trait | Label provenance | Scale | Why it's high-VOI now |
|---|---|---|---|---|
| **NARMS Now** (FDA/CDC/USDA) | Salmonella / Campylobacter / E. coli / Enterococcus AMR MIC | **WET-LAB** reference BMD (Sensititre); FDA-hosted MIC tables + NCBI genomes | thousands/yr | Biggest new BACTERIAL well; extends the deterministic decoder to 4 organisms. **G5 caveat: verify overlap vs BV-BRC/NCBI-PD** (some NARMS reaches you via PATRIC already) |
| **GISP + Euro-GASP** (PRJNA317462 + ENA) | *N. gonorrhoeae* AMR MIC (AZM/CRO/CIP/CFX/PEN/TET) | **WET-LAB** agar dilution / gradient | ~1,700 + ~1,479; Mortimer GWAS ~12k WGS | Highest-yield NEW organism; whole new species cell, provenance-clean |
| **MaveDB** | protein-variant effect (incl. TEM-1 β-lactam, DHFR trimethoprim) | **WET-LAB** multiplexed assays; open REST API | >7M variant effects, >200 score sets | Directly expands the `forward`/`inverse` cell (ProteinGym is a SUBSET). **De-dup vs ProteinGym first** |
| **CDC AR-Bank — non-Enterobacterales panels** | Pseudomonas ~55 / Acinetobacter ~51 / N. gonorrhoeae ~50 / S. aureus panels / **Aspergillus fumigatus ~24** / drug-R Candida ~45 | **WET-LAB** CLSI BMD; per-isolate MIC + BioSample link | tens/species | **SAME source + scraper already used** for the live Enterococcus + C. auris work → nearly-free expansion. Aspergillus/Candida extend the fungal arm |

## Tier 2 — free, wet-lab, but a per-cohort accession↔MIC join to build

- **Dedicated AMR-enzyme DMS** (feeds forward/inverse): TEM-1 β-lactamase — Stiffler/Firnberg (in ProteinGym) + a **55,296-variant / >8M-measurement** epistasis set (biorxiv 2025.07.08.663783, ampicillin+aztreonam); **DHFR/folA trimethoprim** — 1,536 homologs / ~759 species in-vivo TMP-resistance complementation (sciadv.adw9178). Drug-selection fitness, variant-level, wet-lab.
- **Pseudomonas aeruginosa** — Spanish nationwide n=1,445 BMD (13 agents, EUCAST); the ceftolozane-tazobactam study (PRJNA1220180) curates **6 BMD-only public Pseudomonas datasets** — a ready target list.
- **BSAC S. aureus** (agar-dilution MIC, 16 agents, WGS public) + **CRACKLE-2 Klebsiella** (reference BMD, ~593 assemblies) — filter Kp studies to reference-BMD (many use VITEK2).
- **Published C. auris WGS+MIC SRA cohorts** — SRX32028011–78 (n=22 measured MIC), PRJNA816104, PRJNA1264495 — extend the current FKS1/ERG11 fungal cell with MORE clades/isolates.

## Tier 3 — free, wet-lab, but manual transcription from paper tables (new viral cells)

- **WHO/ISIRV influenza NAI global-update** (annual, *Antiviral Research*, mostly OA): isolate-level **measured NA-inhibition IC50 fold-change** outliers + baloxavir PA/I38X. New influenza-NA cell (project already has an influenza-NA marker cell). Supplements, not an API.
- **Lurain & Chou 2010, *Clin Microbiol Rev*** (PMC2952978, OA): ~85 CMV UL97/UL54 mutations with **measured EC50 fold-change** (recombinant phenotyping). A new herpesvirus (CMV/ganciclovir) cell — the free measured herpesvirus source (CHARMD/MRA are interactive-only).

## Tier 4 — genuinely outside-the-box: NON-AMR, non-circular phenotype axes

- **★ Phage susceptibility (BASEL collection)** — *E. coli* host-range, **106 phages × enterobacterial strain panel**, phenotype = **wet-lab EOP plaque assay** (open: genomes PRJNA1207239, EOP matrix in S1 Data, Zenodo 10.5281/zenodo.14277981). Dodges G1 (not tool-derived) AND G3 (not sampling-defined) BY CONSTRUCTION — a physically-measured infection outcome. **Different modeling shape** (phage×host matrix, not per-genome R/S) — this is the strongest genuinely-new axis but needs a new cell design. Klebsiella (Beamud 2023, capsule-receptor) is a second panel.
- **Biolog/PMkbase + 96×10 carbon-utilization** (PMkbase; PLOS Comp Biol 96 isolates × 10 C-sources, all sequenced) — wet-lab metabolic phenotype. **Caveat: intraspecies E. coli depth is thin** (the known ~27-strain wall); the P. putida / cross-species sets are the usable-scale ones.
- **S. aureus biofilm GWAS** (ST-8 cohort, crystal-violet microtiter + WGS) — wet-lab virulence phenotype at GWAS scale.

## Disqualified (tripped a gate — recorded so we don't re-scout)

| Source | Gate tripped | Note |
|---|---|---|
| MalariaGEN Pf7/Pf8 (16k–33k genomes) | **G1** | "resistance" is rules-INFERRED from pfk13/pfcrt; no IC50/clearance. Genomes usable as corpus, label must be discarded |
| WWARN pooled IPD (~3,250 measured PC½) | **not free** | clean wet-lab clearance labels but DAC-gated (malariaDAC@iddo.org). ACQUISITION-decision territory |
| Pathogenwatch (CGPS) | **G1** | genome-DERIVED predicted AMR, not measured AST |
| HBVdb / geno2pheno[hbv,hcv] / FluSurver | **G1** | rules engines (Sensitive/Intermediate/Resistant from marker lists); no fold-change |
| FungAMR (35,792 entries) / MARDy | **G1 if used as labels** | curated CATEGORICAL confidence scores (1–8), not measured MIC — **but excellent as the ERG11/FKS1 determinant CATALOG** (FungAMR + its ChroQueTas scanner). Use as catalog, never as phenotype |
| EUCAST/ECOFF MIC distributions, CDC Reference AST | **aggregate-only** | measured but genome-UNLINKED counts-per-MIC-bin; breakpoint/calibration use only |
| NCTC 3000, FDA-ARGOS | **genome-only** | reference sequences, no paired AST |
| GISAID EpiFlu | **DUA-gated** | some measured NAI metadata but registration + Database Access Agreement |

## Antimalarial (a new kingdom, but thin)

**Cerqueira 2017** (Genome Biology, PRJNA262567) — ~150 P. falciparum isolates, open SRA + **measured
in-vivo artesunate clearance half-life**. The only free wet-lab antimalarial G2P cell found; modest N, one
region/drug, continuous phenotype. The HIVDB-analogue (free + large + isolate-level + measured) does NOT
exist openly for antimalarials — the large measured labels (WWARN IPD) are DAC-gated, the large free
resource (Pf7/Pf8) is label-inferred.

## Recommended next moves (ranked; all reversible, none committed here)

1. **NARMS overlap probe** — cheapest highest-VOI: does NARMS add measured-MIC E. coli/Salmonella/
   Campylobacter/Enterococcus isolates NOT already in the NCBI-PD / BV-BRC footprint? (a metadata-join
   count, no download). If yes → the biggest deterministic-decoder expansion available.
   **ACCESS CONFIRMED (2026-07-20 probe):** NARMS is reachable, NOT a paywall/DUA wall — WGS is public on
   ENA/NCBI (e.g. `PRJNA292668`, Campylobacter; more NARMS BioProjects exist), and an FDA CVM integrated
   NARMS MIC data file is directly downloadable (`https://www.fda.gov/media/79976/download`, HTTP 200). The
   remaining work is the actual **overlap join** (parse the FDA MIC file → isolate accessions → diff against
   the project's 744-accession `cohort_manifest.py` footprint → count net-new measured-MIC isolates per
   organism). **CAVEAT (the load-bearing risk):** NARMS isolates flow INTO NCBI-PD (the same source the 10
   frozen SCORED cells came from), so the overlap could be large — the net-new count is the whole question,
   and it needs the real join, not a guess. This is a clean scoped next-run task (fetch FDA file + manifest diff).
2. **CDC AR-Bank non-Enterobacterales panels** — reuse the existing AR-Bank scraper for the Pseudomonas /
   Acinetobacter / N. gonorrhoeae / Aspergillus / drug-R Candida panels (measured MIC + BioSample already
   scrapeable). Nearly-free given current infra.
3. **MaveDB de-dup + AMR-enzyme DMS pull** — expand the forward/inverse molecular cell beyond ProteinGym's
   217 assays; the TEM-1 8M-measurement + DHFR-TMP sets are the non-circular jewels.
4. **GISP gonococcus** — a whole new provenance-clean species cell (PRJNA317462 + open MIC supplements).
5. **Phage-susceptibility (BASEL)** — the genuinely-novel non-circular axis; needs a NEW cell design
   (phage×host matrix), so scope it as an idea-anchor candidate, not a drop-in.

**One authority fork to name:** #1–#4 are all free + reversible (Soraya's to execute). The WWARN IPD
(clean labels, DAC-gated) is the only ACQUISITION decision — a `malariaDAC@iddo.org` data request is a
user call, not an executor task.
