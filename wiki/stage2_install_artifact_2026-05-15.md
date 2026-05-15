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
| Bakta light | `C:/Users/Farshad/dna_decode_stage2/bakta_db/db-light/` | 4.0 GB (extracted from 1.34 GB tar.xz) | `6.0` (Zenodo DOI 10.5281/zenodo.14916843) |

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

### Bakta — PASS (after conda-init workaround)
- Initial attempts at `--entrypoint /opt/conda/bin/bakta_db download ...` failed with `dependency not found! tool=AMRFinderPlus` on both v1.11.4 and v1.12.0.
- Root cause: bypassing the image's default shell init left the conda env un-activated; `bakta_db`'s dep check (more than `which amrfinder`) found AMRFinder's PATH entry but failed something deeper.
- **Workaround:** invoke via `--entrypoint /bin/bash -c "bakta_db download --output /db --type light"`. The bash login-shell init activates conda properly. v1.11.4 then downloads + extracts the v6 light DB (Zenodo DOI 10.5281/zenodo.14916843) cleanly: 1.34 GB tar.xz → 4.0 GB extracted at `/db/db-light` + MD5 verified.
- Bakta annotation smoke (Step 4 plan §4) NOT yet run — annotation is CPU-heavy and would contend with the in-flight Stage 1 N=40 background job. Run post-Stage-1.

## Resolved blockers

### B1 RESOLVED: Bakta `bakta_db download` AMRFinderPlus dep-check (conda-init workaround)
Direct `--entrypoint /opt/conda/bin/bakta_db` invocation skips the image's default bash shell init → conda env un-activated → `bakta_db`'s AMRFinder dep check fails despite the binary being on PATH. Fix: wrap in `--entrypoint /bin/bash -c "..."`. Validated on v1.11.4. Same wrapper expected to work on v1.12.0; not re-tested since v1.11.4 + DB v6 is sufficient for Stage 2.

Wrapper pattern for `tools/docker_runner.py` callers: pass `image="oschwengers/bakta:v1.11.4"`, `args=["-c", "bakta_db download --output /db --type light"]`, and either expose an `entrypoint=` kwarg or use the existing API with a per-tool helper that prepends `/bin/bash -c`.

### B2: Direct `docker run` from Git Bash needs `MSYS_NO_PATHCONV=1`
Bind mounts via Git Bash silently fail to bind without `MSYS_NO_PATHCONV=1` (the `/db` container path gets munged to `C:/Program Files/Git/db`). Docker daemon doesn't error — it just creates an unbound ephemeral mount, downloads succeed inside the container, then disappear on `--rm`. Symptom: `du -sh <db_path>` returns 0 bytes after a "successful" download.

Mitigation: Python `subprocess.run` bypasses shell path-munging. `tools/docker_runner.py` is therefore unaffected; this only matters for ad-hoc shell invocations.

## What this artifact unblocks
- Phase B of `plans/Stage2_N150_Prep_Plan.md` — AMRFinderPlus POINT* extraction on the N=147 cohort.
- H14 SNP-baseline (Codex's "tighten classical comparator" concern from the 2026-05-15 brainstorm).
- Mash preflight at Stage 2 (leave-one-Mash-clade-out CV).
- Bakta gene-presence comparator path (DB ready; annotation smoke deferred until post-Stage-1 to avoid CPU contention with the in-flight N=40 background job).
- Mash O(N²) refactor of `dna_decode/eval/phylogeny.py` — 2 mash invocations across all strains instead of N*(N-1)/2 (10,731 calls at N=147 → 2). Routed through `tools/docker_runner.run` via the new `use_docker=True` kwarg.

## Deferred from plan scope
- Bakta annotation smoke on K-12 + cipro-R strain (Step 4 plan §4 partial): CPU-heavy, deferred to post-Stage-1 to avoid contention.
- K-12 MG1655 (`GCF_000005845.2`) wild-type substrate: not in local refseq cache (cache populated from BV-BRC GCA accessions, not RefSeq GCF). Mash smoke substituted with two N=40 cipro-R strains; AMRFinder smoke run on the textbook-QRDR substrate (higher-value verification than wild-type).
