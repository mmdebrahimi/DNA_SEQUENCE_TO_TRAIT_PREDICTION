# DNA-LLM within-lineage cohorts — durability + provenance (2026-06-27)

How the multi-drug within-lineage probe data survives a disk loss, and how to reconstruct it.

## Where each artifact lives + its durability

| Artifact | Location | Git-tracked? | Durability |
|---|---|---|---|
| Verdict packets `wiki/functional_alphabet_probe_*_2026-06-2x.{md,json}` | `wiki/` | YES | **GitHub** (committed+pushed) — the irreplaceable conclusions |
| Cohort manifests `wiki/dna_llm_shared_lineage_manifest_*_2026-06-27.json` | `wiki/` | YES | **GitHub** — the accession lists (regeneration input) |
| Closeout + this note | `wiki/` | YES | **GitHub** |
| Cohort parquets `data/processed/shared_lineage_*_cohort.parquet` | `C:` | NO (gitignored `/data/`) | C: only + the D: backup tarball below |
| AMRFinder runs `data/amrfinder_runs/<acc>/{main,mutations}.tsv` (38M total) | `C:` | NO (gitignored) | C: only + the D: backup tarball below |
| Genomes `D:/dna_decode_cache/refseq/<acc>/genome.fna` | `D:` (external) | NO | D: only + re-downloadable from NCBI via the manifest |

## The durability guarantee

1. **Conclusions are on GitHub** (verdict packets + manifests + closeout) — survive total local-disk loss.
2. **Expensive intermediates have 2-location redundancy:** AMRFinder runs (38M) + cohort parquets are
   backed up to `D:/dna_decode_backups/dna_llm_cohorts_<date>.tar.gz` (so a C: failure doesn't force a
   ~3.5 hr Docker rebuild; a D: failure leaves them on C:).
3. **Full reconstruction is deterministic** even from total loss: manifest (GitHub) -> `download_genome`
   each accession (NCBI) -> `_run_amrfinder` each (Docker) -> `functional_alphabet_probe.py` -> verdict.
   Re-run: `scripts/build_drug_shared_lineage_cohort.py --manifest <wiki manifest> --drug <d>` then the probe.

## Regeneration command (per drug, fully offline once genomes are cached)

```bash
export UV_CACHE_DIR=C:/Users/Farshad/AppData/Local/uv_cache_c
export DNA_DECODE_AMRFINDER_DB=C:/Users/Farshad/dna_decode_stage2/amrfinder_db
uv run python scripts/build_drug_shared_lineage_cohort.py \
  --manifest wiki/dna_llm_shared_lineage_manifest_<drug>_2026-06-27.json \
  --drug <drug> --refseq-cache D:/dna_decode_cache/refseq \
  --out-parquet data/processed/shared_lineage_<drug>_cohort.parquet
uv run python scripts/functional_alphabet_probe.py \
  --cohort data/processed/shared_lineage_<drug>_cohort.parquet \
  --drug <drug> --refseq-cache D:/dna_decode_cache/refseq \
  --out wiki/functional_alphabet_probe_<drug>_n120_<date>
```

## Open option (user decision)
The cohort **parquets** are tiny (KB) — force-adding them to git (`git add -f`) would put the cohort
label-tables on GitHub too, leaving only the genomes (huge) + AMRFinder runs (38M, regenerable) off-git.
Not done unilaterally (the `/data/` gitignore is deliberate). Say the word to git-track the parquets.
