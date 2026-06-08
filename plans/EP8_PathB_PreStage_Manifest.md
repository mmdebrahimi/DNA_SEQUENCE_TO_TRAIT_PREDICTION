# EP-8 Path B — pre-stage manifest + baseline spec (workhorse-executable)

> Pre-staged 2026-06-08 on the laptop (no-compute). For the **Precision 7780 (RTX 3500 Ada ~12 GB)**.
> Purpose: turn Path B (Arabidopsis flowering-time embedding test, Gate G2) into pure execution — every
> download URL pinned + verified, the baseline + embedding + CV design fixed, the PASS/FAIL contract frozen.
> Anchors on `plans/EP8_Arabidopsis_Embedding_Test.md` + `wiki/HANDOFF_workhorse_eukaryotic_2026-06-07.md`.
> **G2 is the embedding thesis's decisive 4th test** (after 0-for-3 on bacterial AMR / pathotype / carbon-util).
> Per the ratified pre-commit: a clean G2-FAIL does NOT auto-close the frontier (user chose KEEP-OPEN).

## 0. Iron-law gate — all sources VERIFIED (2026-06-08)
| Asset | Endpoint | Verified | Notes |
|---|---|---|---|
| Phenotype FT10 (10 °C) | `https://arapheno.1001genomes.org/rest/phenotype/261/values.csv` | ✅ 1163 accessions | cols: phenotype_name, accession_id, accession_name, lon, lat, country, **phenotype_value** (days), obs_unit_id |
| Phenotype FT16 (16 °C) | `https://arapheno.1001genomes.org/rest/phenotype/262/values.csv` | ✅ 1123 accessions | same schema; AraPheno Study 12 (1001 Genomes Consortium 2016) |
| Genotype VCF (1135 acc) | `https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz` | search-surfaced (verify on download) | SNPs+short-indels, all 1135; PLINK-readable directly |
| Imputed SNP matrix (HDF5) | `1001_SNP_MATRIX.tar.gz` via `https://1001genomes.org/data-center.html` | search-surfaced | bi-allelic imputed matrix, 1135 acc (for baselines) |
| Pseudogenomes (per-acc FASTA) | 1001 Genomes Data Center → Pseudogenomes | search-surfaced | per-accession TAIR10+SNPs FASTA = the **FM sequence input** |
| Reference + annotation | TAIR10 genome + GFF (Araport11/TAIR10) | standard | to locate flowering-time loci by coordinate |
| Plant DNA-FM | HF `kuleshov-group/PlantCaduceus_l32` (225 M) | confirmed fits 12 GB | 512 bp context; fp16 weights ~0.5 GB; **frozen** (no fine-tune) |

**Join key:** AraPheno `accession_id` == 1001 Genomes accession id (the ecotype/CS id). FT10 ∩ genotype is the
analysis set (~1135 with phenotype; use FT10 primary, FT16 as a replication trait).

**First workhorse step (cheap):** re-pull the two AraPheno CSVs + `HEAD` the genotype URLs to confirm they
resolve before the big downloads. If a 1001genomes.org URL 404s (site reorg), the Data Center page lists the
current path — do NOT guess; surface it.

## 1. Sequence-input strategy for the FM (the crux)
Frozen-FM embedding needs each accession's **DNA sequence**. Two tiers, do (a) first:

**(a) Locus-targeted (primary, cheap, mirrors the bacterial "embed the genes" approach).**
Flowering time in A. thaliana is dominated by known loci — **FLC** (AT5G10140), **FRI** (AT4G00650), **FT**
(AT1G65480), plus **SOC1 / CO / VIN3 / FLM** as a panel. For each accession: pull the pseudogenome, extract
each locus ± ~2 kb (TAIR10 coords from the GFF), tile into 512 bp windows (stride 256), PlantCAD_l32 embed,
mean-pool per locus, concat loci → per-accession feature vector. This is ~8 loci × a few windows = cheap
(seconds/accession on the RTX 3500), and it's the fair test: "do FM embeddings of the causal loci beat a SNP
model of the same?"

**(b) Genome-wide mean-pool (secondary, only if (a) is inconclusive).** Tile the whole pseudogenome, embed,
mean-pool. Far more compute; defer unless (a) is borderline.

> Mirror the bacterial cache discipline: write per-accession-per-locus embeddings to an HDF5 cache once
> (`cache.populate`-style), `verify_complete` before the mean-pool consumer runs (half-flushed accession =
> valid-looking mean otherwise). Reuse `dna_decode/models/cache.py` patterns.

## 2. Baselines (the bars the embedding MUST beat) — no GPU
Run on the SAME accessions + SAME CV folds as the embedding.
1. **SNP-PRS (ridge / GBLUP).** Ridge regression (or genomic BLUP) on the SNP matrix → FT. Two variants:
   genome-wide SNPs, and FT-locus-only SNPs (the matched control for the locus-targeted embedding).
2. **Kinship / population-structure-only.** GBLUP with the genomic relationship matrix, OR ridge on the top
   ~10–20 genotype PCs → FT. This is the **de-confound bar** — the trap that sank cipro NT (the embedding
   must beat *pure relatedness/structure*, not just k-mer). See
   `[[feedback_embedding_vs_knowledge_baseline_and_within_lineage]]`.
3. (reference) **published PRS / heritability ceiling:** FT10/FT16 mapped loci explain ~55 %/48 % of genetic
   variance (Atwell/1001G) — a sanity ceiling, not a bar.

## 3. Embedding model
PlantCAD_l32 frozen embeddings (§1a) → ridge regression head (same α-search as the SNP-PRS, nested CV) → FT.
Report R² (and Spearman) out-of-fold.

## 4. Cross-validation = clade-stratified (the de-confound)
Cluster accessions into genetic groups (the 1001 Genomes **9 admixture groups**, or k-means on the top
genotype PCs if the admixture labels aren't readily joined) → **leave-one-group-out CV** (the Arabidopsis
analogue of LOSO/LOMO/leave-one-Mash-clade-out). Reuse `dna_decode/eval/cv.py` fold machinery.

**Within-lineage diagnostic (load-bearing).** Within each genetic group, recompute embedding-predicted vs
actual FT (R² / Spearman) and compare to the structure-only baseline within the same group. If the embedding
only predicts *between* groups (= population structure) but is at chance *within* a group, it learned lineage,
not flowering biology — the cipro-NT failure mode. Reuse the logic of `scripts/within_lineage_diagnostic.py`.

## 5. Gate G2 PASS/FAIL contract (frozen)
- **PASS** — embedding out-of-fold R² beats BOTH baselines (SNP-PRS AND structure-only) by a meaningful margin
  (pre-commit **≥ 0.05 R²** over the better baseline) under leave-one-group-out CV, **AND** within-group
  embedding R² > within-group structure R² (predicts mechanism, not just structure). First evidence the
  frozen-FM thesis earns its keep.
- **FAIL** — embedding ≤ baselines OR within-group = structure. Per the ratified pre-commit (ledger Decision
  row 3), a clean FAIL does **NOT** auto-close the embedding frontier — it prompts one more considered attempt
  (different FM / genome-wide §1b / different quantitative trait) before any closure. Document honestly either
  way (result packet `wiki/phase2_arabidopsis_result_<date>.md`).

## 6. Compute notes (RTX 3500 Ada 12 GB)
- PlantCAD_l32 inference fits 12 GB at fp16, batch the 512 bp windows (no quantization). Locus-targeted (§1a)
  is the cheap path — whole analysis set embeds in well under a day.
- Cache embeddings to local disk + `verify_complete` (per §1). Restartable; never re-embed a cached accession.
- Baselines (§2) + CV (§4) are CPU/sklearn — no GPU.
- Money gate: local 12 GB only. Do NOT provision cloud/Databricks GPU without explicit user approval
  (ledger Decision row 2).

## 7. Deliverable + handoff-back
`wiki/phase2_arabidopsis_result_<date>.md` = embedding R² vs both baselines + the within-group verdict = the
G2 PASS/FAIL. Push to origin + a one-line `project-state-row` (Path-B only). Then the laptop + workhorse decide
at G2 whether to iterate or conclude (frontier stays open on a single FAIL per the pre-commit).

## DO NOT (dual-machine hygiene)
- Do NOT route personal code through the Bombardier/DLP machine — Path B is the **personal Precision 7780** only.
- Do NOT touch Path-A fungal artifacts (laptop-owned). Append only Path-B rows + the Path-B result packet.
- Do NOT fire paid cloud compute without explicit approval (money gate).
