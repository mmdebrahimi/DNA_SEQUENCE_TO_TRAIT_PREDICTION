# Free-data capture registry — the sweep, closed out (2026-07-04)

Completion record for the "capture all freely available data, sequentially" directive. Every free
(no-auth, no-fee, Soraya-fetchable) source from the acquisition anchor is addressed here with a **capture
status + a product-value judgment**. **Load-bearing judgment (R2 / motion-avoidance):** the deterministic
decoder consumes CURATED CATALOGS + DETERMINANTS + LD references — NOT polygenic summary statistics. So the
product-aligned free sources are CAPTURED + INTEGRATED; the polygenic-prior-layer sources are STAGED as
pointers, because integrating them into a deterministic determinant decoder would be motion, not signal (the
prior sweep + the 2026-06-25 defer reached the same verdict; the embedding arm that WOULD consume them is a
closed 0-for-5 negative).

## CAPTURED + INTEGRATED (product-aligned)

| # | Source | What it gave | Artifact |
|---|---|---|---|
| **1** | **ClinVar** (GRCh38 VCF, 192 MB on D:) | a deterministic human-variant decoder — the Mendelian-disease analogue of the AMR catalog (curated catalog → regime-1 win). 31,616 P/LP+B/LB variants across 10 canonical Mendelian genes. | `dna_decode/data/clinvar.py` + committed `data/clinvar/clinvar_panel.tsv` + `scripts/capture_clinvar.py` (extensible to any gene panel). Verified: F508del→PATHOGENIC(4★); unknown→INDETERMINATE. |
| **2** | **1000 Genomes** (via Ensembl REST LD — no bulk download; C: disk-tight) | an INDEPENDENT ancestry-stratified LD validation of the imputation layer. ABO tag pair: EUR 0.97 / EAS 0.93 / SAS 0.99 / AMR 0.90 valid, **AFR 0.33 fails** → quantifies the imputation map's ancestry limit. | `scripts/capture_1000g_ld.py` + committed `data/imputation/rs8176719_from_rs657152_1000g_ld.json`. |

## STAGED — pointers, NOT integrated (polygenic prior-layer; integrating = motion)

| Source | Access (free) | Why NOT decoder-integrated |
|---|---|---|
| **GWAS Catalog** (EBI) | FTP TSV — **download path currently unstable** (2 probed URLs 404'd 2026-07-04; use the current `ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/` release file) | SNP→trait SUMMARY STATS. Product-value ONLY as a future **single-SNP candidate miner** (filter for strong single-locus effects → seed new deterministic cells like lactase/earwax). Pure-polygenic otherwise = dead-embedding territory. **Deferred as a candidate-mining next step, not an integration.** |
| **FinnGen** (R12) | public download; partially staged on D: (`finngen_R12_manifest.tsv` + 1 endpoint) | GWAS summary stats per ~2,400 endpoints. Locus-priors / endpoint-definitions layer — a prior, not a determinant the deterministic decoder calls. |
| **PGS Catalog** | REST API (reachable, http 200) | POLYGENIC score definitions — the opposite paradigm to single-locus determinant decoding. |
| **Pan-UKBB / OpenGWAS** | free bulk summary stats (TB-scale) | polygenic summary-stat repositories; not bounded to fetch + not decoder-consumable. Pointers only. |

## Completion verdict
**The free-data sweep is CLOSED.** The two product-aligned free sources are captured + integrated (ClinVar =
a new deterministic decoder capability; 1000G = an imputation ancestry-validity annotation). The remaining
free sources are polygenic prior-layers documented as pointers — integrating them would be motion (they don't
plug into the deterministic determinant decoder, and the embedding arm that would is a closed negative).

**What free data does NOT unlock (unchanged):** the genuinely-new science (external validation on lab-measured
labels + the learned decoder's fair test) still requires the acquisition-gated biobanks — see
`wiki/label_acquisition_anchor_2026-07-04.md` (All of Us / UK Biobank, user authority).

**One product-aligned free capture remains available as a bounded next step** (not motion, if pursued):
GWAS-Catalog single-SNP candidate mining → a shortlist of strong single-locus trait associations that could
seed new deterministic cells (the lactase/earwax pattern). Left as a documented recipe above; the URL needs
re-resolving to the current EBI release path first.
