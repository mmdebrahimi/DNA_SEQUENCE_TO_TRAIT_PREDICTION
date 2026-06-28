# Wide-net gene→phenotype data-source scan — beyond AMR/pathogens (2026-06-28)

A deliberate, discipline-agnostic sweep for gene→phenotype data the decoder could consume — NOT confined to
antibiotics. The user's framing: "gene→eye colour / furry-or-not / skin colour … anything, think outside the
box." This catalogs candidate sources across 8 disciplines, gate-scored, with a ranked shortlist.

## The filter that makes this honest (read first)
The project has TWO products with OPPOSITE data needs:
- **The DETERMINISTIC decoder (the WORKING product):** a curated *small-gene-set rule* → trait. Wants
  **Mendelian / oligogenic** traits where a *curated determinant rule exists* + *free per-individual measured*
  labels exist. This is the AMR pattern (AMRFinder/WHO-catalogue rule + measured AST), generalized.
- **The LEARNED-embedding arm (CLOSED 0-for-4 negative):** wants polygenic/complex traits. **Do not feed it.**
  Arabidopsis *flowering time* (polygenic) already FALSIFIED the embedding bet (2026-06-12).

So the wide net is scored for the DECODER: a source is high-value iff its trait is **Mendelian/oligogenic
with a curated rule**, not merely "has genomes + phenotypes." A polygenic panel is data-rich but decoder-dead.

## Hard gates (per the negative-results map, generalized off-AMR)
G-IND per-individual (genome↔its OWN trait) · G-MEAS measured not genome-inferred (else circular) ·
G-FREE free/public/downloadable · G-POW enough trait-positive individuals · G-RULE a curated gene→trait
rule exists (decoder needs one) · G-ARCH trait architecture is Mendelian/oligogenic (decoder territory).

## The catalog (verified 2026-06-28)

| # | source (discipline) | what | G-IND | G-MEAS | G-FREE | G-RULE | G-ARCH | verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | **HIrisPlex-S + OpenSNP** (human visible traits) | HIrisPlex-S = deployed forensic RULE (eye/hair/skin colour, 41 SNPs, free webtool); OpenSNP = free per-person genotype + self-reported colour | ✅ | ◐ self-report | ✅ | ✅ (deployed rule) | ✅ oligogenic | **FLAGSHIP** |
| 2 | **Dog coat colour** (companion animal) | MC1R/ASIP/CBD103 curated rule (VGL UC Davis test); Darwin's Ark 3,287 dogs genotyped + owner-photo colour + Dog10K | ✅ | ◐ owner-photo | ✅ | ✅ | ✅ oligogenic | **STRONG** |
| 3 | **1011 Yeast Genomes** (microbial, non-AMR) | 1,011 isolates + 36 traits; the BIMODAL/Mendelian ones (CuSO₄, anisomycin resistance) are decoder-shaped; free portal | ✅ | ✅ lab-measured | ✅ | ◐ curatable | ◐ Mendelian subset | **STRONG** |
| 4 | **AraPheno / AraGWAS** (plant) | 462 phenotypes × 1,135 A. thaliana accessions, free CSV/PLINK; AraGWAS = curated determinant catalog. FILTER to Mendelian traits (seed/trichome) — flowering-time is polygenic = embedding-dead | ✅ | ✅ | ✅ | ✅ (AraGWAS) | ◐ Mendelian subset only | **MEDIUM** |
| 5 | **DGRP / DGRPool** (fly) | 205 inbred lines + 1,034 harmonized phenotypes, free, GWAS built-in | ✅ | ✅ | ✅ | ◐ | ✗ mostly polygenic | MEDIUM (Mendelian subset) |
| 6 | **CaeNDR** (worm) | C. elegans/briggsae/tropicalis wild strains + phenotype DB + GWAS, AWS-open, MIT | ✅ | ✅ | ✅ | ◐ | ✗ mostly quantitative | MEDIUM (Mendelian subset) |
| 7 | **3000 Rice Genomes / SNP-Seek** (crop) | 3K accessions + grain traits + API; some Mendelian (grain colour, aroma BADH2) | ✅ | ✅ | ✅ | ◐ | ◐ | MEDIUM |
| 8 | **Human Mendelian disease** (ClinVar/OMIM) | curated variant→disease (rule exists), but free PER-INDIVIDUAL measured-phenotype-linked genomes are scarce/gated (UK Biobank gated) | ◐ | ✅ | ✗ (labels gated) | ✅ | ✅ | LOW (label-gated) |

## The reframe (the answer to "are we stuck in AMR?")
**No — the deterministic decoder generalizes to ANY Mendelian/oligogenic trait with a curated rule + free
measured labels.** AMR was just the first instance of that pattern. The wide net found the SAME pattern in
non-pathogen biology, and one combo is a complete free cell available today:

**FLAGSHIP: the human visible-trait cell (HIrisPlex-S rule + OpenSNP labels).** It is the exact off-AMR
analog of the validated AMR cells — a curated, independently-developed gene→trait RULE (HIrisPlex-S, like
AMRFinder) scored against FREE, per-individual, independently-MEASURED labels (OpenSNP self-reported eye/hair
colour, like measured AST). It realizes the user's "gene→eye colour" example literally, and it is the
project's first jump OUT of the pathogen/AMR box. Honest caveats: OpenSNP labels are self-reported (noisier
than a lab assay — a `◐` on G-MEAS, but still independent of the genome → non-circular) and ancestry-confounded
(needs a within-ancestry check, the same de-confounding discipline as within-lineage AMR).

## Ranked shortlist (what to actually pursue, decoder-first)
1. **Human eye/hair colour — HIrisPlex-S vs OpenSNP** (flagship; free, today; non-pathogen; user's own example).
2. **Dog coat colour — MC1R/ASIP/CBD103 vs Darwin's Ark** (free citizen labels; the "furry/colour" example).
3. **Yeast Mendelian traits — CuSO₄/anisomycin resistance, 1011 panel** (free lab-measured; microbial non-AMR; cleanest measured labels of the lot).
4. (lower) Arabidopsis/rice/fly/worm — deep panels, but the decoder only wants their Mendelian SUBSETS; the polygenic bulk is embedding-arm territory (closed). Mine per-trait, don't adopt wholesale.

## Honest scope + next step
- This is a SOURCE SCAN (a plan input), not a built cell. Each shortlist item is a candidate `--until-mvp`
  cell: ingest the rule + the free measured labels, score within-ancestry/lineage (de-confound), report honestly.
- The G-ARCH filter is load-bearing: do NOT feed polygenic traits to the dead embedding arm; the wins are the
  Mendelian/oligogenic + curated-rule cells.
- For the visiting marine biologist (`wiki/collaborator_brief_marine_*` if saved): the same pattern in HIS
  domain = marine microbe / fish trait with a known major-gene + free measured per-individual data.

## Sources (verified 2026-06-28)
OpenSNP [opensnp.org](https://opensnp.org/) + [80-phenotype benchmark arXiv:2603.06768](https://arxiv.org/html/2603.06768v2);
HIrisPlex-S [hirisplex.erasmusmc.nl + FSI:Genetics 2018](https://www.fsigenetics.com/article/S1872-4973(18)30220-5/abstract);
Darwin's Ark + Dog10K [PNAS 2025](https://www.pnas.org/doi/10.1073/pnas.2421752122), [Dog10K Genome Biol 2023];
1011 Yeast [Nature 2018](https://www.nature.com/articles/s41586-018-0030-5) + [1002genomes.u-strasbg.fr];
AraPheno/AraGWAS [arapheno.1001genomes.org](https://arapheno.1001genomes.org);
DGRP/DGRPool [dgrpool.epfl.ch](https://elifesciences.org/articles/88981); CaeNDR [caendr.org](https://caendr.org/);
3000 Rice / SNP-Seek [snp-seek.irri.org](https://academic.oup.com/nar/article/43/D1/D1023/2435696).
