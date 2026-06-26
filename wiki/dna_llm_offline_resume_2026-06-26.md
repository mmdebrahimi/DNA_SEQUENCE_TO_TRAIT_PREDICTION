# DNA-LLM within-lineage levers — OFFLINE RESUME GUIDE (2026-06-26)

You're traveling with slow internet. This run **pre-downloaded the genomes** (the only internet-bound step)
to **D:** so everything downstream — AMRFinder (local Docker) + the within-lineage probe — runs **fully
offline**. Here is exactly what to do when you pick this back up.

## What was downloaded (internet, this session)
Shared-lineage cohort genomes for **5 AMR drugs**, each a 40-MLST-lineage / 400-strain manifest (200R/200S),
cached at `D:/dna_decode_cache/refseq/<accession>/genome.fna`:
- ciprofloxacin, tetracycline, gentamicin, ceftriaxone, meropenem.
- Manifests (the strain lists): `wiki/dna_llm_shared_lineage_manifest_<drug>_2026-06-26.json`.
- Cohort parquets (strain_id/assembly_accession/mlst/ast_<drug>): `data/processed/shared_lineage_<drug>_cohort.parquet`.

(Download was launched download-only; some drugs may not be 100% complete if internet dropped — the build
below is **restartable** and re-downloads only what's missing, so finishing it later costs only the gap.)

## Step 1 — finish AMRFinder per drug (OFFLINE, local Docker; ~35s/strain)
Run WITHOUT `--download-only`; it skips already-downloaded genomes + already-done AMRFinder runs:
```bash
DNA_DECODE_AMRFINDER_DB=C:/Users/Farshad/dna_decode_stage2/amrfinder_db \
uv run python scripts/build_drug_shared_lineage_cohort.py \
  --manifest wiki/dna_llm_shared_lineage_manifest_tetracycline_2026-06-26.json \
  --drug tetracycline --refseq-cache D:/dna_decode_cache/refseq \
  --out-parquet data/processed/shared_lineage_tetracycline_cohort.parquet
```
Repeat with `--drug ciprofloxacin|gentamicin|ceftriaxone|meropenem` + the matching manifest/parquet.
Needs only Docker (no internet). If Docker churn corrupts the WSL2 mount (125 mkdir errors), `wsl --shutdown`
+ relaunch Docker Desktop, then re-run (restartable).

## Step 2 — run the within-lineage probe per drug (OFFLINE, CPU-only; ~10min/drug for k-mer)
```bash
uv run python scripts/functional_alphabet_probe.py \
  --cohort data/processed/shared_lineage_tetracycline_cohort.parquet \
  --drug tetracycline --refseq-cache D:/dna_decode_cache/refseq
# -> wiki/functional_alphabet_probe_n147_<date>.{md,json}  (rename per drug)
```
The verdict bucket (BEATS_KMER / TIES / FAILS / UNDERPOWERED) + the within-lineage concordances + the paired
in-MLST permutation p are in the packet. With 40 shared lineages this is far more powered than the cipro
N=147 probe (6 lineages, p=0.0565).

## The scientific read (per the closeout `wiki/functional_alphabet_probe_closeout_2026-06-26.md`)
- **tet** is the headline test: if the curated-determinant alphabet FAILS/TIES within-lineage (tetA/B presence
  doesn't separate R/S *within* a lineage — resistance hinges on something AMRFinder misses), that is the first
  concrete **headroom** signal for a learned model (#3). If it cleanly separates (like cipro QRDR -> 1.0), it's
  another tautological determinant win, no headroom.
- The other 4 drugs are the comparison set: which mechanism regimes show the determinant alphabet failing
  within-lineage (= distributed/incomplete) vs winning (= concentrated/curated).

## State pointers
- Featurizer is drug-general (`dna_decode/eval/functional_tokens.py`, classify_gene_symbol + paren-strip).
- The cipro N=147 probe = TIES (functional 1.000 vs k-mer 0.721, p=0.0565) — the within-lineage win there is
  the curated QRDR determinants = the deterministic decoder, no headroom.
- Frozen AMR surface untouched throughout. No GPU needed for any of this.
