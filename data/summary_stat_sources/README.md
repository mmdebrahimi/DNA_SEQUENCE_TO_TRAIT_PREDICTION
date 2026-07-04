# Summary-stat source indexes — captured 2026-07-04 (user-directed)

Committed compact INDEXES of the polygenic summary-stat sources. **Honesty caveat (unchanged):** these are
polygenic prior-layer sources the DETERMINISTIC determinant decoder does not directly consume (the embedding
arm that would is a closed 0-for-5 negative). They are captured as AVAILABLE, indexed resources per the
user's directive; the ONE product-aligned extract is the GWAS single-SNP candidate shortlist (seeds new
deterministic single-locus cells). The TB-scale per-phenotype bulk stays on D: / is fetched on demand via the
recipes below.

| File | Source | Rows | What it is | Bulk-fetch recipe |
|---|---|---|---|---|
| `pgs_catalog_index.tsv` | PGS Catalog | 5,385 | polygenic-score metadata (id, name, trait, n_variants, ftp) | per-score harmonized scoring file at each row's `ftp_scoring_file` |
| `pan_ukbb_phenotype_index.tsv` | Pan-UKBB | 7,223 | phenotype manifest (trait_type, phenocode, description, n) | per-pheno multi-ancestry sumstats via the manifest aws/gs paths (TB-scale) |
| `finngen_r12_endpoint_index.tsv` | FinnGen R12 | 2,469 | endpoint manifest (phenocode, phenotype, category, n_cases/controls) | `gs://finngen-public-data-r12/summary_stats/` per endpoint |
| `gwas_single_snp_candidates.tsv` | GWAS Catalog | 400 | **product-aligned**: strongest single-rsID, genome-wide-sig (P≤5e-8), OR 2–20 associations → candidate deterministic single-locus cells | full 716 MB associations TSV on `D:/dna_decode_cache/gwas_catalog/` |

## OpenGWAS — NOT captured (auth-gated)
OpenGWAS (MRC IEU) requires a **free JWT token** for most API requests since 2024-05. The study index cannot
be fetched without it. To capture: register at `https://api.opengwas.io/`, export the token, then
`GET https://api.opengwas.io/api/gwasinfo` with `Authorization: Bearer <token>`. This is a credential the
USER provides (not a Soraya-executable free fetch).

## Reproduce
```bash
uv run python scripts/capture_summary_stat_indexes.py          # PGS + Pan-UKBB + FinnGen indexes
uv run python scripts/mine_gwas_single_snp_candidates.py       # GWAS single-SNP candidate shortlist (needs the D: TSV)
```

## GWAS single-SNP candidate caveat
The GWAS Catalog `OR or BETA` field is a dirty free field (mixes OR / beta / malformed values up to ~1e14);
the miner caps at OR ≤ 20 to exclude data-entry artifacts (a verify-in-batch fix). The shortlist is a
candidate-DISCOVERY resource, not validated cells — each candidate still needs its own build (fetch a
joinable genotype×phenotype cohort, confirm the single-locus rule beats null), exactly like the
lactase/earwax/eye-colour cells.
