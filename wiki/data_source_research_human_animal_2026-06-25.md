# Wide data-source research — expanding the genotype→phenotype decoder (human + other animals), 2026-06-25

**Directive:** find free, machine-readable genotype→phenotype label + genotype sources to improve the
current decoder phase (the PGx / variant→catalog→call phase the CYP2C19 cell just validated at GeT-RM
core 72/72). **Tier: DISCOVERY** (web-sourced; surfaced for user acceptance, NOT written to the ledger as
fact). The lens is the project's standing constraint — *labels, not models* — + the 8 rejection gates.

## Headline
The exact pattern just shipped (PharmVar/ClinPGx defs → CPIC table → **GeT-RM consensus ⋈ 1000G genotypes**
via Docker bcftools) **generalizes to ~13 more human PGx genes for free**, reusing the same harness. AND a
genuinely new catalog-tractable kingdom opens up — **other animals** — with free labels (OMIA) + a free
open dog population VCF (Dog10K), including a direct veterinary-PGx analog (dog MDR1/ABCB1).

---

## A. Human PGx — the immediate, lowest-risk expansion (same pattern, same harness)
Each gene below has all three free deps the CYP2C19 cell used: **PharmVar/ClinPGx** allele defs +
**CPIC** diplotype→phenotype table + a **GeT-RM consensus** (or PharmCAT-on-1000G) truth set the existing
`scripts/pgx_getrm_concordance.py` already consumes (sample→diplotype TSV).

**Tractable NOW with the SNP-proxy caller (mostly simple SNP/indel star alleles):**
- **CYP2C9 + VKORC1** — the warfarin pair; both SNP-defined, both GeT-RM + CPIC. Natural next cell.
- **NUDT15, TPMT** — thiopurine dosing; small variant sets, high field concordance.
- **DPYD** — fluoropyrimidine toxicity; SNP-defined (truth via PharmCAT-on-1000G, defs updated post-GeT-RM).
- **SLCO1B1** — statin myopathy; SNP-defined.
- **CYP3A5, CYP4F2, UGT1A1, CYP2B6** — high reported concordance, low variant counts.
- *(some need an activity-score layer, not just function-pair — a small caller extension.)*

**Deferred — structurally hard (NOT the SNP pattern):**
- **CYP2D6** — CYP2D7 pseudogene + CNV + hybrids; the field's hardest (best tools ~87.5%). Needs
  structural/CNV calling. Out of the v0 SNP-caller scope (same reason the cell deferred it).
- **G6PD** (X-linked), **NAT2 / CYP1A2** (newer PharmVar additions).

**Honest tier (unchanged from CYP2C19):** GeT-RM consensus is itself caller-derived (Astrolabe/Stargazer/
Aldy), so it's "vs the field's accepted reference," not a wet-lab phenotype — calling independently
validatable, phenotype faithful-to-CPIC. The strongest-grounded expansion the project can make.

Sources: PharmVar (pharmvar.org/genes — free per-gene defs + API; ~15 genes incl. CYP1A2/2A6/2B6/2C8/2C9/
2C19/2D6/3A4/3A5/4F2/DPYD/NUDT15/NAT2/SLCO1B1); CPIC/ClinPGx (clinpgx.org — free diplotype→phenotype +
allele-function tables, downloadable); GeT-RM (Coriell + the J Mol Diagn characterization papers; the
ursaPGx benchmark table already vendored gives a machine-readable consensus for CYP2C8/2C9/2C19/2D6).

## B. Other animals — a new catalog-tractable kingdom (the directive's second half)
**Labels — OMIA** (omia.org): the "OMIM for animals." Free, curated, **downloadable causal-variant tables**
(~1000 likely-causal variants, OMIAvariant IDs since 2021) across 755 species. Mendelian-disease causal
variants: **dog 322 / cattle 189 / horse 73**, plus non-disease traits (coat colour, gait). Cross-linked to
NCBI/Ensembl for coordinates. This is the catalog layer — the animal analog of the AMR/PGx catalogs.

**Genotype source (the feasibility gate — can we TEST a call against real genomes?):**
| Species | Open population genotypes? | Resource |
|---|---|---|
| **Dog — STRONG** | ✅ fully open VCF | **Dog10K** (1987 dogs, 48M variants; kiddlabshare + Zenodo 8084059) + **DBVDC** (582 dogs + genotype freqs) + **Darwin's Ark** (3277 WGS, open, PLINK) |
| Horse — moderate | ✅ per-dataset (no single consortium VCF) | 88-horse WGS, 185 Thoroughbred, 6-breed tracks (ENA / animalgenome.org) |
| Cattle — gated | ⚠️ raw seq public, **merged VCF consortium-restricted** | 1000 Bull Genomes (SRA/ENA accessions) |

**Why dog is the cleanest "other animal" substrate:** free labels (OMIA) + free open population VCF
(Dog10K, like a canine 1000G) + a **demonstrated genotype↔phenotype validation precedent** (e.g. MC1R
c.916C>T e/e → solid-white coat, p=1.4e-26 across 26 dogs in DBVDC). The pattern is identical to our cells:
single causal variant (often one SNP/indel) → presence call → trait/disease call.

**The bridge cell — dog MDR1/ABCB1 (the standout):** the `nt230(del4)` / ABCB1-1Δ deletion → P-gp
loss → macrocyclic-lactone (ivermectin) neurotoxicity in collies. It is **veterinary pharmacogenomics** —
a direct continuation of the current PGx phase into another species, with a single well-characterized
causal variant (OMIA label) + Dog10K genotypes. Cleanest first animal cell.

**Animal gate caveat (must check per trait):** per-genome *phenotype* labels are thinner than human PGx —
Dog10K has breed metadata but not always the specific trait phenotype per dog. Validation denominator is
the gate: **dog coat-colour + MDR1 have the best phenotype availability**; rarer disease traits may be
underpowered. Same "suspect the label / check the denominator" discipline as the AMR cells.

## C. Cross-cutting truth set for the CALLING layer (defensive)
**GIAB CMRG benchmark** (NIST Genome-in-a-Bottle; HG002, 273 challenging medically-relevant genes incl.
pharmacogenes, GRCh37+38; wet-lab-grade haplotype-resolved assembly). Not a phenotype source — it validates
**variant DETECTION** on hard genes. Also flags reference false-duplication traps (CBS/CRYAA/KCNE1 — recall
8%→100% when masked) — a caution before adding any cell whose gene sits in a segdup.

## Recommended next moves (ranked; user picks — discovery-tier, nothing auto-started)
1. **Human PGx expansion — CYP2C9 + VKORC1 (warfarin pair).** Highest VOI, lowest risk: SNP-defined, both
   in GeT-RM + CPIC, reuses the entire harness + the 1000G fetch already on disk. Could add 3–5 genes fast.
2. **Veterinary-PGx bridge — dog MDR1/ABCB1 (ivermectin sensitivity).** The cleanest new-kingdom cell +
   thematic continuity (PGx → animal PGx); single causal deletion, OMIA label, open Dog10K genotypes.
3. **Animal trait cell — dog coat colour (MC1R/ASIP/TYRP1/MFSD12…).** Rich OMIA + Dog10K + a published
   genotype↔phenotype validation precedent; broader than #2 but well-supported.
4. **GIAB CMRG** wired as a calling-layer truth check for any future gene in a hard/segdup region.

## Gate summary
- Human PGx genes: clear the gates exactly as CYP2C19 did (free defs + free CPIC + free GeT-RM/1000G;
  not circular-to-a-tool; deep same-organism). **Green.**
- Dog OMIA cells: faithful-to-OMIA catalog tier; free labels + free open genotypes (Dog10K). Gate to watch
  = per-genome phenotype availability (coat-colour/MDR1 best). **Green-with-denominator-check.**
- Cattle: genotype VCF gated (consortium) → lower priority. Horse: workable but per-dataset assembly.

_Sources: pharmvar.org/genes; clinpgx.org (CPIC); Coriell GeT-RM + Gaedigk/Pratt J Mol Diagn papers;
omia.org (+ Nicholas/Tammen Animal Genetics); Dog10K (NAR 2025; Zenodo 8084059); DBVDC (Jagannathan 2019,
Animal Genetics, PMC6842318); Darwin's Ark (Science 2022; Dryad 10.5061/dryad.83bk3jb4r); 1000 Bull Genomes
(PRJEB56689); GIAB CMRG (Nature Biotech; NIST). All free except the 1000 Bull merged VCF._
