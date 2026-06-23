# The independent TB number — runbook (rec #2, Docker-ready host)

The saga's holy grail is now gated ONLY on compute, not on a label. The disjoint cohort + the scorer are
built; the one remaining step (assembly → CRyPTIC-style masked VCF) needs a Docker-ready host with the
clockwork/regeno variant-calling pipeline. **This host cannot run it** (no Docker; minimap2/bcftools/samtools/
datasets all absent), and a generic minimap2+bcftools VCF would NOT match the masked-VCF format
`tb_vcf.parse_masked_calls` expects (FILTER==PASS + GT alt; the MIN_DP/MIN_FRS/MIN_GCP floor) — so the proper
pipeline is required. The AMR Portal AMRFinder TB POINT calls are NOT a substitute (they would test AMRFinder's
TB caller, not the WHO catalogue — the documented confound).

## What's already done (free, code-side)
- **Cohort READY:** `data/raw/tb_goldset/amr_portal_tb_disjoint_cohort.tsv` (built by
  `scripts/amr_portal_tb_cohort.py`): **26,941 provenance-disjoint** M. tuberculosis isolates with measured
  RIF/INH (not in CRyPTIC, not in our cohorts); **2,845 carry a fetchable assembly accession** = the
  immediately-runnable independent set. Disjoint measured: RIF 26,645 / INH 25,773.
- **Scorer READY:** `scripts/score_tb_independent_goldset.py` + `scripts/build_tb_goldset_manifest.py` +
  `organism_rules/tb_amr` (the WHO-2023-catalogue rule). These already consume a goldset manifest of masked
  VCFs and emit the `INDEPENDENT_VALIDATION` arm (raw + lineage-collapsed sens/spec + Wilson CI, reusing the
  frozen `clonality.cluster_weighted_confusion`).

## The remaining step (Docker-ready host) — 4 commands
1. **Fetch the 2,845 assemblies** (the `assembly` column, GCA/ERZ accessions) via NCBI Datasets:
   ```bash
   # for each GCA in amr_portal_tb_disjoint_cohort.tsv (leaked==0 && assembly nonempty):
   datasets download genome accession <GCA> --include genome
   ```
2. **Produce CRyPTIC-style masked VCFs vs H37Rv NC_000962.3** (the COMPUTE-gated step) — the clockwork /
   gnomonicus / regeno pipeline (Docker), the same masked-VCF shape the project's `tb_vcf` already parses.
   One masked VCF per isolate (`<strain_id>.masked.vcf`).
3. **Build the leakage-checked goldset manifest** (re-confirms disjointness vs CRyPTIC by construction):
   ```bash
   # candidate TSV: strain_id, ena_accession/biosample, masked_vcf, regeno_vcf(optional), rif_label, inh_label
   uv run python -m scripts.build_tb_goldset_manifest --candidates <candidates_from_cohort>.tsv
   ```
4. **Score the WHO rule → the independent number:**
   ```bash
   uv run python -m scripts.score_tb_independent_goldset --drug rifampicin
   uv run python -m scripts.score_tb_independent_goldset --drug isoniazid
   ```
   → `INDEPENDENT_VALIDATION` (raw + lineage-collapsed sens/spec + Wilson CI), provenance-disjoint, measured
   DST — the first genuinely-independent TB number the gold-set saga (5 author/DUA/circular walls) could not
   get. FREE (no DUA, no author contact); only compute remained, and the LABEL is unblocked.

## Wall classification (honest)
**Code-closable but EXTERNAL-COMPUTE-gated.** Needs: a Docker-ready host + the clockwork/regeno masked-VCF
pipeline + ~2,845 genome fetches + variant-calling (hours-to-days; the project's known TB regeno cost). NOT
runnable on the current host (no Docker / no tools). Everything up to that line — the disjoint cohort, the
leakage gate, the scorer, the WHO rule — is built and tested. The independent TB number is one Docker-host
run away.

## Provenance
Cohort `scripts/amr_portal_tb_cohort.py` → `data/raw/tb_goldset/amr_portal_tb_disjoint_cohort.tsv` (+
`wiki/amr_portal_tb_cohort_stats.json`). Feasibility `wiki/amr_portal_feasibility_result_2026-06-23.md`;
independent-validation memo `wiki/amr_portal_independent_validation_2026-06-23.md`. Existing TB scorer
`scripts/score_tb_independent_goldset.py`.
