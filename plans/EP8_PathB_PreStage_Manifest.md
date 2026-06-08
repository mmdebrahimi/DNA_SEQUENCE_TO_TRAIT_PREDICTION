# EP-8 Path B — pre-stage manifest + baseline spec (workhorse-executable)

> Pre-staged 2026-06-08 on the laptop (no-compute). For the **Precision 7780 (RTX 3500 Ada ~12 GB)**.
> Purpose: turn Path B (Arabidopsis flowering-time embedding test, Gate G2) into pure execution — every
> download URL pinned + verified, the baseline + embedding + CV design fixed, the PASS/FAIL contract frozen.
> Anchors on `plans/EP8_Arabidopsis_Embedding_Test.md` + `wiki/HANDOFF_workhorse_eukaryotic_2026-06-07.md`.
> **G2 is the embedding thesis's decisive 4th test** (after 0-for-3 on bacterial AMR / pathotype / carbon-util).
> Per the ratified pre-commit: a clean G2-FAIL does NOT auto-close the frontier (user chose KEEP-OPEN).
> **REVISED 2026-06-08 (post-brainstorm):** primary sequence-input swapped from curated loci → phenotype-AGNOSTIC
> subsample (the original curated panel contradicted the "no curated catalog" niche criterion); added the
> CPU-only **§0.5 G2 dry-manifest** gate before any GPU; PASS now requires a paired-bootstrap CI excluding 0
> (not a point margin); kinship/PCs are the primary de-confounder (geography = sensitivity-only). Open
> design decisions for the user are listed at the end (§8).

## 0. Iron-law gate — all sources VERIFIED (2026-06-08)
| Asset | Endpoint | Verified | Notes |
|---|---|---|---|
| Phenotype FT10 (10 °C) | `https://arapheno.1001genomes.org/rest/phenotype/261/values.csv` | ✅ 1163 accessions | cols: phenotype_name, accession_id, accession_name, lon, lat, country, **phenotype_value** (days), obs_unit_id |
| Phenotype FT16 (16 °C) | `https://arapheno.1001genomes.org/rest/phenotype/262/values.csv` | ✅ 1123 accessions | same schema; AraPheno Study 12 (1001 Genomes Consortium 2016) |
| Genotype VCF (1135 acc) | `https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz` | ✅ HEAD HTTP 200, 19.2 GB | SNPs+short-indels, all 1135; PLINK-readable directly |
| Imputed SNP matrix (HDF5) | `1001_SNP_MATRIX.tar.gz` via `https://1001genomes.org/data-center.html` | search-surfaced | bi-allelic imputed matrix, 1135 acc (for baselines) |
| Pseudogenomes (per-acc FASTA) | 1001 Genomes Data Center → Pseudogenomes | search-surfaced | per-accession TAIR10+SNPs FASTA = the **FM sequence input** |
| Reference + annotation | TAIR10 genome + GFF (Araport11/TAIR10) | standard | to locate flowering-time loci by coordinate |
| Plant DNA-FM | HF `kuleshov-group/PlantCaduceus_l32` (225 M) | confirmed fits 12 GB | 512 bp context; fp16 weights ~0.5 GB; **frozen** (no fine-tune) |

**Join key:** AraPheno `accession_id` == 1001 Genomes accession id (the ecotype/CS id). FT10 ∩ genotype is the
analysis set. **FROZEN endpoint contract: FT10 is the sole PRIMARY endpoint; FT16 is replication ONLY** — both
frozen before any model comparison (no two-endpoint cherry-picking).

**First workhorse step (cheap):** re-pull the two AraPheno CSVs + `HEAD` the genotype URLs to confirm they
resolve before the big downloads. If a 1001genomes.org URL 404s (site reorg), the Data Center page lists the
current path — do NOT guess; surface it.

## 0.5 STEP 0 — G2 dry-manifest (CPU-only; MUST pass before ANY GPU embedding)
A pre-flight that proves the test is runnable before compute is spent. Produce + commit a dry-manifest report;
do not start §3 embedding until every check is green. Checks:
1. **Accession intersection** — AraPheno FT10 `accession_id` ∩ {pseudogenome files present} ∩ {SNP-matrix
   columns}. Report the final N and every dropped accession + reason. (Expect ~1122 with phenotype; the
   genotyped intersection is the real analysis N.)
2. **Pseudogenome ID join** — confirm the AraPheno `accession_id` → pseudogenome filename rule empirically
   (do NOT assume `pseudo{id}.fasta.gz`); validate by opening 3 files. Surface the real pattern.
3. **Window/coordinate table** — materialize the exact window set the FM will embed (see §1), as a BED/TSV
   with per-window TAIR10 coords. Frozen artifact; the SNP baseline (§2) reads the SAME table.
4. **Per-window/per-accession N-fraction QC** — fraction of `N`/uncalled bases per window per accession;
   flag windows/accessions above a threshold (uncalled regions silently corrupt both embedding + consensus).
5. **Matched variant feature matrix** — build the SNP+indel+missingness feature matrix over the SAME windows
   as the embedding, so §2 baselines are information-matched (not SNP-only vs sequence).
6. **Genetic-group labels present** — confirm the 9 admixture-group labels (or PC-derived groups) resolve for
   the analysis-N accessions (the CV needs them; see §4).
Output: `wiki/g2_dry_manifest_<date>.{md,json}`. Any red check → STOP + surface; do not burn GPU.

## 1. Sequence-input strategy for the FM (the crux — REVISED 2026-06-08 post-brainstorm)
> **Why revised:** the original primary (hand-picked FLC/FRI/FT panel) CONTRADICTS this substrate's own
> embedding-niche criterion (`EP8_Arabidopsis_Embedding_Test.md`: "No curated mechanism catalog"). Hand-picking
> the known causal loci injects the exact domain knowledge the niche claim says the test avoids — it silently
> becomes "can the FM represent these 8 curated loci", not "does a frozen DNA-FM find flowering biology
> unaided". Full genome-wide is the conceptually-clean alternative but is operationally hostile on 12 GB
> (~135 Mb → ~260k–520k windows/accession × ~1122 = a compute-endurance test) AND a naive global mean-pool
> dilutes sparse causal signal (a loss would falsify "global mean-pooling", not frozen embeddings).

**(a) PRIMARY — phenotype-AGNOSTIC genome subsample.** A predeclared window set chosen with **NO flowering-time
labels and NO known-flowering-gene enrichment** (this is what preserves the no-catalog niche). Pick ONE
selection rule + FREEZE it in the dry-manifest before any embedding:
- **all gene bodies + fixed flanks** (Araport11 genes ± ~1 kb), tiled to 512 bp; OR
- a **random stratified per-chromosome window sample** (fixed seed, fixed count/accession); OR
- **unsupervised top-variance windows** (rank windows by cross-accession sequence variance — uses genotype
  variation only, never FT). 
For each accession: pull the pseudogenome, extract the frozen window set, embed with frozen PlantCAD_l32,
mean-pool (and ALSO keep variance-pool — see note), → per-accession feature vector. Selection is frozen +
phenotype-blind → this is the niche-honoring test.

**(b) SECONDARY DIAGNOSTIC — known-flowering-locus representation (NOT the gate).** The original curated panel
(**FLC** AT5G10140, **FRI** AT4G00650, **FT** AT1G65480, + SOC1/CO/VIN3/FLM), same extract→embed→pool. Report
it with an **explicitly narrowed claim**: "does PlantCAD add value over matched local SNP/indel encodings AT
known causal loci" — an upper-bound diagnostic, NOT evidence for the no-catalog embedding thesis.

> **Pooling note:** mean-pool alone may be too blunt over a genome-wide window set; also compute a
> variance/max-pool (or per-chromosome block-pool) so a null result reflects the representation, not just the
> averaging. Decide the gate-pooling in the dry-manifest.
> **Cache discipline:** write per-accession-per-window embeddings to an HDF5 cache once (`cache.populate`-style),
> `verify_complete` before the pooling consumer runs (half-flushed accession = valid-looking pool otherwise).
> Reuse `dna_decode/models/cache.py` patterns.

## 2. Baselines (the bars the embedding MUST beat) — no GPU
Run on the SAME accessions + SAME CV folds as the embedding.
1. **SNP-PRS (ridge / GBLUP) — INFORMATION-MATCHED.** Ridge / genomic-BLUP on a variant feature matrix built
   over the **SAME frozen window set** the embedding sees (§0.5 step 5), and it MUST include **short indels +
   per-window missingness indicators**, not SNPs only — else a SNP-only baseline is unfairly weaker than a
   pseudogenome embedding that sees indels + `N` patterns. Two variants: the agnostic-window matrix (matched to
   §1a, the primary bar) and a genome-wide-all-SNP matrix (a broader reference bar).
2. **Kinship / population-structure-only (the PRIMARY de-confounder).** GBLUP with the genomic relationship
   matrix, OR ridge on the top ~10–20 genotype PCs → FT. The embedding must beat *pure relatedness/structure*,
   not just k-mer — the trap that sank cipro NT. See
   `[[feedback_embedding_vs_knowledge_baseline_and_within_lineage]]`.
3. (reference) **published PRS / heritability ceiling:** FT10/FT16 mapped loci explain ~55 %/48 % of genetic
   variance (Atwell/1001G) — a sanity ceiling, not a bar.

## 3. Embedding model
PlantCAD_l32 frozen embeddings (§1a) → ridge regression head (same α-search as the SNP-PRS, nested CV) → FT.
Report R² (and Spearman) out-of-fold.

## 4. Cross-validation + de-confounding (REVISED 2026-06-08 post-brainstorm)
**Folds:** cluster accessions into genetic groups (1001 Genomes **9 admixture groups**, or k-means on top
genotype PCs if labels aren't readily joined) → **leave-one-group-out CV**. Reuse `dna_decode/eval/cv.py`.

**Primary de-confounder = kinship / PCs, NOT geography.** Population structure (relatedness) is the confounder
the embedding must beat; the structure-only baseline (§2.2) + leave-one-group-out are the primary guards.
**Do NOT residualize the phenotype against geography in the gate** — latitude/climate are *partly causal* for
flowering (vernalization / local adaptation), so regressing them out subtracts real biology and then penalizes
the embedding for not recovering it. Geography residualization is a **sensitivity analysis only** ("does the
edge survive after removing geography-correlated signal?"), never a pass condition.

**Within-group diagnostic — CONTINUOUS, predeclared (NEW; the binary AMR script does NOT apply).**
`scripts/within_lineage_diagnostic.py` is binary within-MLST R-vs-S concordance — there is NO continuous
analogue, so a new one is required. Predeclared spec: take out-of-fold predictions; **center predictions and
phenotype within each genetic group**; compute **within-group R²** (pooled across groups, group-size-weighted)
for the embedding vs the structure-only baseline. A **minimum within-group N gate** applies — groups below it
return `indeterminate`, not pass/fail. If the embedding predicts only *between* groups but is at chance
*within* group, it learned lineage, not flowering biology (the cipro-NT failure mode).

## 5. Gate G2 PASS/FAIL contract (frozen — REVISED 2026-06-08 post-brainstorm)
Graded on the **§1a phenotype-agnostic embedding** (the niche-honoring primary), endpoint **FT10**, under
leave-one-group-out CV. The §1b curated-locus diagnostic is reported but does NOT decide the gate.
- **PASS** — the agnostic embedding out-of-fold R² beats BOTH baselines (information-matched SNP-PRS §2.1 AND
  structure-only §2.2), where the margin on `(embedding − best baseline) R²` is backed by a **paired bootstrap
  / grouped-permutation CI that EXCLUDES 0** (not a bare point estimate — with ~9 unequal groups, point margins
  are too brittle; the prior "≥0.05" is a guide, the CI-excludes-0 is the bar), **AND** within-group R²
  (§4 continuous) > within-group structure R². First real evidence the frozen-FM thesis earns its keep.
- **FAIL** — CI does not exclude 0 vs a baseline, OR within-group = structure. Per the ratified pre-commit
  (ledger Decision row 3), a clean FAIL does **NOT** auto-close the frontier — it prompts ONE more considered
  attempt (different FM / a stronger pooling or representation / paid-cloud true-genome-wide / a different
  quantitative trait) before any closure. Document honestly either way (`wiki/phase2_arabidopsis_result_<date>.md`).
- **INDETERMINATE** — dry-manifest red, or within-group N below the gate, or the agnostic-window run did not
  complete — report as such; never report a point margin as a verdict.

## 6. Compute notes (RTX 3500 Ada 12 GB)
- PlantCAD_l32 inference fits 12 GB at fp16, batch the 512 bp windows (no quantization). The §1a agnostic
  subsample is sized to a **fixed window budget per accession** (set in the dry-manifest) so the run is a
  bounded GPU-hours job, NOT a full-genome endurance test. Full genome-wide is explicitly out of local scope
  (money-gated cloud only).
- Cache embeddings to local disk (D:) + `verify_complete` (per §1). Restartable; never re-embed a cached
  accession. (Heed the Docker/USB-drive + `wsl --shutdown` lessons from Path A if Docker is in the loop.)
- Dry-manifest (§0.5), baselines (§2), CV (§4) are CPU/sklearn — no GPU.
- Money gate: local 12 GB only. Do NOT provision cloud/Databricks GPU without explicit user approval
  (ledger Decision row 2). The "true genome-wide" full-thesis run is the main thing that would need it.

## 7. Deliverable + handoff-back
`wiki/phase2_arabidopsis_result_<date>.md` = embedding R² vs both baselines + the within-group verdict = the
G2 PASS/FAIL. Push to origin + a one-line `project-state-row` (Path-B only). Then the laptop + workhorse decide
at G2 whether to iterate or conclude (frontier stays open on a single FAIL per the pre-commit).

## 8. Open design decisions for the user (surfaced by the 2026-06-08 brainstorm — decide before/at the dry-manifest)
These change what the gate measures; pick before the §0.5 dry-manifest freezes the run:
1. **Primary estimand** — cross-group *portability* vs structure-independent *causal signal* vs realistic-structure
   *prediction*. They grade differently; the chosen one sets which CV/baseline contrast is decisive. (Default if
   unspecified: structure-independent causal signal — it's the truest test of "the embedding learns biology.")
2. **Agnostic-window selection rule (§1a)** — "all gene bodies + flanks" vs "random stratified per-chromosome"
   vs "unsupervised top-variance". (Default: all gene bodies + flanks — interpretable + genome-representative.)
3. **GPU-hours / window budget (§6)** — fixed windows/accession, fixed GPU-hours, or full completion. Sets the
   subsample size. (Default: a fixed per-accession window budget sized to a sub-day local run.)
4. **True-genome-wide full-thesis run** — only via paid cloud (money gate). Default: NOT run unless the agnostic
   subsample is borderline AND the user approves spend.

## DO NOT (dual-machine hygiene)
- Do NOT route personal code through the Bombardier/DLP machine — Path B is the **personal Precision 7780** only.
- Do NOT touch Path-A fungal artifacts (laptop-owned). Append only Path-B rows + the Path-B result packet.
- Do NOT fire paid cloud compute without explicit approval (money gate).
