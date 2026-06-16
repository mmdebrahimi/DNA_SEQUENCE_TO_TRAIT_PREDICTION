# Acquirable label sources to reopen the decoder beyond free-public AMR — ranked shortlist (2026-06-16)

Research deliverable for the **acquisition** strategic fork (chosen 2026-06-16). Each candidate scored
pass/fail against the 4 non-negotiable gates: **G1 sampling-independent** (lab measurement, not sampling
context) · **G2 non-circular** (not tool-derived) · **G3 organism-depth** (≥~100-150 same-organism
isolates WITH downloadable assemblies) · **G4 provenance-disjoint-feasible** (≥20/class outside the
NARMS/CDC/FDA/GenomeTrakr surveillance ecosystem the project excludes, not clonally collapsed).

**Headline finding (reframes the fork):** the highest-value new label source is **FREE and immediately
fetchable (CRyPTIC), NOT a real-world MTA/contact action.** "Acquire a label source" → in practice a
download + a new-organism determinant catalog (modeling the executor can do), not a procurement step.

## Ranked shortlist

| # | Source | G1 | G2 | G3 | G4 | acquisition path | assemblies | N after filters |
|---|---|---|---|---|---|---|---|---|
| 1 | **CRyPTIC** M. tuberculosis | ✓ | ✓ | ✓✓ | ✓ (consortium, not US-surveillance; no TB tuning set exists → disjoint by construction) | **FREE FTP, immediate** | yes (ENA + VCF/compendium) | **12,289 × 13 drugs** |
| 2 | **Kayama 2023** Japan E. coli + K. pneumoniae | ✓ (reference BMD) | ✓ | ✓ | likely ✓ (Japanese national surveillance — geographically DISJOINT from the US NCBI-PD/NARMS tuning provenance) | free download IF genomes deposited (DDBJ) — **confirm data-availability** | likely (DDBJ) — verify | ~4,195 (E. coli + Kpn) |
| 3 | **Nguyen/Houston Methodist** K. pneumoniae | ✓ | ✓ | ✓ | **AT RISK** — genomes live in PATRIC/BV-BRC, the SAME source as the decoder's tuning set → must run the leakage check | free (PATRIC FTP) | yes (PATRIC IDs) | 1,667 × 20 drugs |
| 4 | **TransPred** E. coli growth-kinetics | ✓ (NEW phenotype: doubling-time/yield at sub-MIC) | ✓ | ✓ | ✓ (diverse non-surveillance: hospital/community/food/wild/wastewater) | free download (paper supplement) | yes | 1,407 |
| — | **Pfizer ATLAS** | ✓ | ✓ | **✗ NO paired genomes** | n/a | free (Wellcome/Micron) | **none** | 633k MICs but unusable (can't link to a genome) |
| — | **EUCAST MIC distributions** | ✓ | ✓ | **✗ aggregate ECOFF histograms, no per-isolate genomes** | n/a | free | none | unusable as training labels |
| — | **Whittam DECA / von Mentzer ETEC** | ✓ | ✓ | ✓ | ✓ | **MTA / direct-contact (real-world action)** | partial | reference-set scale |

## The recommendation: CRyPTIC is the move
- **Clears all 4 gates and needs no real-world acquisition** — `ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/` (reuse table CSV = genotype+phenotype for 12,288 isolates + VCF paths; reads on ENA). Reference-grade **measured** MIC (UKMYC broth-microdilution plates), 13 drugs.
- **Lands in the decoder's STRONGEST mechanism class.** TB resistance is overwhelmingly POINT-mutation-driven (rpoB/katG/gyrA/embB/pncA) — exactly the regime where the deterministic QRDR-style rule already wins (cipro 0.925, Klebsiella 1.0). The WHO TB mutation catalogue / AMRFinder TB give the determinant set.
- **Disjoint by construction:** the project has no TB tuning set, so leakage against the E. coli/Klebsiella accessions is impossible; only intra-CRyPTIC clonality needs the existing Mash-lineage check.
- **The "cost" is modeling, not procurement:** a TB determinant catalog + organism route — work the executor can do, the genome→phenotype kingdom-jump the project already proved on fungal C. auris.

## Honest caveats
- **CRyPTIC = new organism.** This is a coverage expansion to a 3rd kingdom-adjacent organism, not an E. coli deepening. If the goal is specifically a *deeper E. coli* substrate, **Kayama Japan (#2)** is the targeted pick — same organisms the decoder covers, reference BMD, plausibly provenance-disjoint by geography — but its genome data-availability/accessions must be confirmed first (couldn't fetch through the publisher auth wall).
- **Nguyen/Houston (#3)** is the cautionary one: large + paired, but PATRIC-sourced genomes risk overlapping the tuning set → it only counts after the accession-manifest leakage gate clears it; and BD-Phoenix is automated, not reference broth microdilution.
- **ATLAS + EUCAST are disqualified** as training substrates despite huge MIC volume — no per-isolate genomes to learn from (the project's recurring "96% assembly-availability drop", here 100%).
- **Whittam/von Mentzer** remain the only true MTA/contact options — pursue ONLY if the free genome-paired sources above prove insufficient.

## Sources
- Pfizer ATLAS (open raw MIC, 633k isolates, no WGS): https://www.nature.com/articles/s41467-022-30635-7
- Kayama 2023 Japan national genomic+BMD surveillance (4,195 E. coli/Kpn): https://www.nature.com/articles/s41467-023-43516-4
- Nguyen/Houston Methodist K. pneumoniae genome→MIC (1,667 × 20, PATRIC genomes): https://pmc.ncbi.nlm.nih.gov/articles/PMC11522393/
- CRyPTIC data compendium (12,289 M. tb × 13 drugs, UKMYC BMD, free FTP): https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001721
- EUCAST MIC/zone distributions (aggregate ECOFFs): https://www.eucast.org/bacteria/mic-and-zone-distributions-ecoffs/
- TransPred E. coli growth-kinetics under sub-MIC (1,407, non-surveillance): https://journals.asm.org/doi/10.1128/msystems.00346-21

## Recommended next step
A **CRyPTIC feasibility probe** (executor-doable, free, no MTA): fetch the reuse table, count isolates with
both measured MIC + downloadable genome per drug, score one POINT-driven drug (rifampicin/isoniazid) with
an AMRFinder/WHO-catalogue determinant rule, and check intra-cohort clonality. That converts the chosen
acquisition fork into a concrete first scored TB cell — without any real-world acquisition action.
