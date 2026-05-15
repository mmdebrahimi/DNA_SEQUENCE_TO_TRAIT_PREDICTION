# Stage 2 Toolchain Install Artifact — 2026-05-15

Reproducibility anchor for the Stage 2 Docker-tools install per `plans/Stage2_Docker_Tools_Install_Plan.md`. Records image SHAs, DB versions, smoke-test results, and known blockers.

## Host environment
- Windows 10 Home, build 10.0.19045
- Docker Desktop installed + running (`docker info` returns engine info)
- Git Bash shell on host (note: requires `MSYS_NO_PATHCONV=1` for direct `docker run` invocations; Python `subprocess.run` bypasses path munging — `tools/docker_runner.py` is therefore unaffected)
- C: free at install time: ~25 GB (chronically near-full per CLAUDE.md)
- D: free: 4.5 TB (Seagate Portable, USB-hiccup-prone)

## Pinned images + SHA digests (operational)
| Tool | Image | SHA digest |
|---|---|---|
| Mash 2.3 | `quay.io/biocontainers/mash:2.3--hb105d93_10` | `sha256:57e2af229b118c706b09d864f0a14c9fa99803c1d38e5f22ef3cbdc7c2ad3c1e` |
| AMRFinderPlus 4.2.7 | `ncbi/amr:4.2.7-2026-03-24.1` | `sha256:a42ead5b6fe439e399b2ceda212cd1dadb510116af8535c6be1cd1a5bc16b3f2` |

**Plan deviation:** Plan pinned `mash:2.3--he348c14_4` — that tag no longer exists on quay.io (manifest unknown). Switched to current `2.3--hb105d93_10` (Mash binary still 2.3).

## Failed-attempts (pre-blocker B1)
| Tool | Image | SHA digest | Result |
|---|---|---|---|
| Bakta v1.11.4 | `oschwengers/bakta:v1.11.4` | `sha256:cce1ab4c9f4eb1905f49f74ff15d4fa279af0c8db5c9956a49072d08dc635a61` | `bakta_db download` fails AMRFinderPlus dep-check (see B1) |
| Bakta v1.12.0 | `oschwengers/bakta:v1.12.0` | `sha256:4cb7c3f8327483d661013b5bcbb33f0dc71119d661c4c887ec1dce6d81542c60` | Same dep-check regression as v1.11.4 |

Both images remain in local Docker cache (~3.4 GB combined dead weight). Keep until B1 resolved — re-pulling burns network on each retry attempt.

## Databases
| DB | Path | Size | Version |
|---|---|---|---|
| AMRFinderPlus | `C:/Users/Farshad/dna_decode_stage2/amrfinder_db/latest/` | 240 MB | `2026-03-24.1` |
| Bakta light | (not installed — see Blocker B1) | — | — |

## Smoke-test results

### Mash — PASS
- Sketch + pairwise dist on 2 cipro-R N=40 strains (ST131 `GCA_000522345.1`, ST410 `GCA_003073955.1`).
- Self-distances: both 0.000 ✓
- Cross-strain (ST131 ↔ ST410): 0.0324 (sane for cross-lineage E. coli)
- Tool exited 0; output parseable.

### AMRFinderPlus — PASS (with cipro-R textbook validation)
- Cohort strain: `1328433.3` / `GCA_000522345.1` (ST131, cipro-R, 5 contigs, N50=5.08 Mbp) — pinned per `plans/Stage2_Docker_Tools_Install_Plan.md` D6 (replaces loose "any ST131 R" criterion).
- Command: `amrfinder -n /in/ST131_R.fna -O Escherichia --database /db/latest --mutation_all /out/mutations.tsv -o /out/main.tsv`
- Wall-time: 95 seconds (translated-nucleotide + mutation search; tblastn + blastn).
- `main.tsv`: 20 rows (including header).
- `mutations.tsv`: 668 rows (including header).
- **Textbook QRDR mutations confirmed:**
  - `gyrA_S83L` (QUINOLONE / QUINOLONE) ✓
  - `gyrA_D87N` (QUINOLONE / QUINOLONE) ✓
  - `parC_E84V` (QUINOLONE / QUINOLONE) ✓
  - (Note: ST131 commonly carries `parC_S80I`; this strain has `E84V` — still a QRDR mutation.)
- Closes Codex H14 SNP-baseline concern: AMRFinder POINT-row extraction works end-to-end on a known cipro-R substrate.

### Bakta — BLOCKED
- Tested both `v1.11.4` and `v1.12.0`. Both fail `bakta_db download` with:
  ```
  dependency not found! tool=AMRFinderPlus
  ERROR: AMRFinderPlus not found or not executable! Please make sure AMRFinderPlus is installed and executable
  ```
- Verified `amrfinder --version` returns `4.0.23` inside the v1.11.4 container (so the binary IS present); the dep check in `bakta_db` must be doing more than `which amrfinder`.
- Deferred. Bakta's consumer (gene-presence comparator) is itself blocked on Pending Decisions row 4 (RefSeq vs Bakta vs Roary annotation source decision); not a Stage 1 verdict blocker.

## Open blockers

### B1: Bakta `bakta_db download` AMRFinderPlus dep-check regression
Both v1.11.4 and v1.12.0 reject `bakta_db download` arguing AMRFinderPlus is unavailable, despite `/opt/conda/bin/amrfinder --version` returning a valid version string inside the same container. Root cause unknown — likely a brittle PATH check / version-format mismatch / `--threads` discovery. Workaround paths:
- (a) Direct Zenodo download of the v6 light DB tar.gz, skipping `bakta_db download`.
- (b) Try Bakta v1.10.x (before the AMRFinder dep was added).
- (c) Build a custom Docker image that satisfies the dep check.
- (d) Defer Bakta entirely until Pending Decisions row 4 picks an annotation source other than Bakta (Roary / Panaroo / Prokka alternatives exist for gene-presence comparator work).

Not a Stage 1 verdict blocker.

### B2: Direct `docker run` from Git Bash needs `MSYS_NO_PATHCONV=1`
Bind mounts via Git Bash silently fail to bind without `MSYS_NO_PATHCONV=1` (the `/db` container path gets munged to `C:/Program Files/Git/db`). Docker daemon doesn't error — it just creates an unbound ephemeral mount, downloads succeed inside the container, then disappear on `--rm`. Symptom: `du -sh <db_path>` returns 0 bytes after a "successful" download.

Mitigation: Python `subprocess.run` bypasses shell path-munging. `tools/docker_runner.py` is therefore unaffected; this only matters for ad-hoc shell invocations.

## What this artifact unblocks
- Phase B of `plans/Stage2_N150_Prep_Plan.md` — AMRFinderPlus POINT* extraction on the N=147 cohort.
- H14 SNP-baseline (Codex's "tighten classical comparator" concern from the 2026-05-15 brainstorm).
- Mash preflight at Stage 2 (leave-one-Mash-clade-out CV).

## Deferred from plan scope
- Step 6 (`phylogeny.py` Mash O(N²) refactor): separate commit, not part of this artifact.
- Bakta light DB download: blocked on B1.
- K-12 MG1655 (`GCF_000005845.2`) wild-type smoke: K-12 not in local refseq cache (cache populated from BV-BRC GCA accessions, not RefSeq GCF). Substituted with two existing N=40 cipro-R strains for Mash smoke; AMRFinder smoke run on the textbook-QRDR substrate (higher-value verification than wild-type).
