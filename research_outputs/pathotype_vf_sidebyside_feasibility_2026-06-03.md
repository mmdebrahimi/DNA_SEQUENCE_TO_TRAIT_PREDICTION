# VF side-by-side diff — feasibility (2026-06-03, Soraya overnight)

**Status: BLOCKED — needs external aligner tooling not installable unattended on this host. Deferred to user.**

The ledger requires a side-by-side diff of the v0 resolver vs CGE VirulenceFinder gene-calls. Real VirulenceFinder needs an external aligner — **KMA** (default) or **BLAST+** — to align the DB alleles to the assembly. Tonight's host check:

- `blastn` / `makeblastdb` / `kma` — **not on PATH**.
- `conda` / `mamba` — **absent** (no bioconda channel available).
- The uv-managed `.venv` has **no `pip`** (`uv pip` only); the `virulencefinder` PyPI package is a thin wrapper that still shells out to KMA/BLAST — pip-installing it does NOT provide the aligner.
- KMA is C source requiring compilation; BLAST+ on Windows is an interactive MSI (~200 MB). Neither is safe to attempt unattended, and a half-installed BLAST could break the environment.

**Decision (Soraya, money-only + judgment):** do NOT attempt the install unattended. Document + defer. This matches the recommendation's own "gate it" note and the dogfood finding that dep-installs should be gated.

## What the user can do (when ready)
1. Install BLAST+ for Windows (NCBI MSI) OR KMA (via WSL/conda), put it on PATH.
2. `uv pip install virulencefinder` (or use the Bitbucket repo) — DB already cached at `data/virulencefinder_db/`.
3. Run VF on a few genomes (e.g. `data/ena_wgs/AIEY01.fna`) → its results table.
4. Diff vs our resolver's `marker_hits` (already in `research_outputs/pathotype_cli_demo_*.json`).

## Honest interim position
Our v0 resolver uses the SAME VirulenceFinder allele DB, but with pure-Python k-mer-seed presence instead of KMA/BLAST alignment. So `caller_is_independent_baseline=false` is correctly set, and `percent_identity` is null (`method=kmer_seed_coverage_v0`). The two would agree on strong present/absent calls; they will differ on (a) borderline-identity alleles (BLAST %ID vs our coverage proxy) and (b) exact coordinates. A real diff quantifies that gap — valuable but not blocking for v0, which already validated (eae AUROC 1.0). Everything else in tonight's recommendation set (cohort eval, ExPEC calibration, caching) was executable and done.
