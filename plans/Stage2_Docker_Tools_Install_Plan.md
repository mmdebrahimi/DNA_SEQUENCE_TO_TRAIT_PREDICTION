# Stage 2 Docker-Tools Install Plan

> After user starts Docker Desktop, install Mash + Bakta + AMRFinderPlus via pinned Docker images, write a single `tools/docker_runner.py` Python orchestration module (NOT `.sh` wrappers), and smoke-validate on K-12 + one cipro-R strain. Resolves the Phase A.1 / A.2 / A.5 install steps from `plans/Stage2_N150_Prep_Plan.md`.

---

## Problem Statement

Stage 2 prep needs three Linux-only bioinformatics CLI tools:
1. **Mash** — genome distance for leave-one-Mash-clade-out CV preflight at N=147 (used by `dna_decode/eval/phylogeny.py`).
2. **Bakta** — re-annotation for cross-strain stable `gene_symbol` (RefSeq's ~11% coverage is too low; `INDETERMINATE_IDENTIFIER_OOV` would fire on the gene-presence comparator otherwise).
3. **AMRFinderPlus** — POINT* SNP-table baseline for the gyrA/parC/parE textbook resistance signal that NT must beat.

Host is Windows 10. All three tools are Linux-only per their official channels (Bakta PyPI classifier `POSIX :: Linux`; AMRFinderPlus canonical bioconda Linux; Mash no Windows binary). Conda/mamba not installed. WSL2 has only Docker Desktop's internal distros (no Ubuntu). Docker Desktop IS installed.

The /brainstorm round (2026-05-14 PM) of an earlier proposal flagged four critical issues that this plan corrects:
- `.sh` wrapper scripts don't work with Python `subprocess.run([...])` on Windows + can't translate `C:\...` / `D:\...` host paths into Linux container paths.
- AMRFinderPlus command shape was wrong: `--database` not `--database_path`; `--mutation_all <file>` takes a path argument not a boolean.
- Image tags `:latest` are moving targets — Stage 2 results need reproducibility.
- D: drive concentration risk: putting Bakta DB + AMRFinder DB + outputs + active populate writes all on the same volatile USB drive compounds the failure mode that crashed an 8-hour populate today.

## Design Decisions

### D1: Docker (containers) — NOT WSL2 Ubuntu

**Decision:** Use Docker Desktop containers for all three tools.

**Rationale:** Docker Desktop is already installed on the host; user just needs to start it. WSL2 Ubuntu install is interactive (first-boot wizard wants username/password) and adds ~30 min of user-side setup before any tool can run. Docker images are pull-and-run; no host-environment changes beyond a `tools/docker_runner.py` Python module.

**Trade-off:** Container spin-up overhead ~3-5 sec per `docker run --rm` call. Matters at scale for Mash (11,175 pairwise calls under the existing nested loop in `phylogeny.py`). Mitigated by the companion `phylogeny.py` refactor (single `mash sketch` + single `mash dist sketch.msh sketch.msh` call) that turns 11,175 invocations into 2.

### D2: One Python `tools/docker_runner.py` module — NOT three `.sh` wrappers

**Decision:** Replace the originally-proposed `scripts/run_mash.sh` / `scripts/run_bakta.sh` / `scripts/run_amrfinder.sh` with a single `tools/docker_runner.py` Python module. The module exposes `run(tool: str, args: list[str], *, mounts: dict[str, str] = {}, **kwargs) -> CompletedProcess` and handles host→container path mapping (`D:/dna_decode_cache/...` → `/cache/...`, `C:/Users/Farshad/PythonProjects/dna_decode/...` → `/repo/...`) before invoking `subprocess.run([docker, run, ...])`.

**Rationale:** Per /brainstorm Round 1 finding:
- `dna_decode/eval/phylogeny.py` calls `subprocess.run([mash_binary, ...])` directly. A `.sh` wrapper is not directly executable via Windows `CreateProcess` — Python would fail to launch it without explicit `bash.exe` shim.
- The proposed `--volume <cwd>:<cwd>` mount uses Windows path syntax that doesn't translate inside Linux containers. Container target must be Linux-style.
- A Python module can normalize Windows paths, emit reproducible docker-run commands for logging, accept per-tool config, and slot directly into `subprocess.run([...])` from Python callers like `phylogeny.py`.

**Trade-off:** Slightly more code than three shell scripts (~200 LOC for one Python module vs ~20 LOC × 3 = 60 LOC of bash). Acceptable because the Python module is more testable (mockable subprocess), more maintainable (one place to add a new tool), and avoids the Windows-Python path-translation hazards.

### D3: Pin Docker image tags — NOT `:latest`

**Decision:** Pin concrete image tags:
- `ncbi/amr:4.2.7-2026-03-24.1` (current Docker Hub release; verified from `hub.docker.com/r/ncbi/amr`)
- `oschwengers/bakta:<concrete-version>` (TBD — query at install time)
- `quay.io/biocontainers/mash:2.3--he348c14_4` (already pinned in original plan, good)

Record image SHA digests in the Stage 2 install artifact at `wiki/stage2_install_artifact_<date>.md` alongside DB versions.

**Rationale:** Stage 2 results need reproducibility. `:latest` is a moving target; a re-pull weeks later could pull a different image version and silently shift POINT-row schemas or annotation outputs.

**Trade-off:** Slightly older versions in some cases; manual update step needed if/when a critical bugfix lands upstream. Acceptable because Stage 2 is a one-shot decision gate, not a continuously-updated pipeline.

### D4: Correct AMRFinderPlus invocation pattern

**Decision:** Use:
- `amrfinder_update --database /db` for DB setup (separate command, NOT `amrfinder --update`)
- Runtime: `amrfinder -n /cache/<acc>/genome.fna -O Escherichia --database /db/latest --mutation_all /out/mutations.tsv -o /out/amrfinder.tsv`

`--database` is the (singular) flag for custom DB path; `--mutation_all` takes a file path argument (NOT a boolean flag) — the point-mutation report writes to that path.

**Rationale:** Per /brainstorm Round 1 finding: the original proposal used `--database_path` (not a real flag) and treated `--mutation_all` as a boolean. The smoke test would have inspected only the main TSV (which wouldn't have POINT rows) and reported a false-clean install while the actual point-mutation artifact was never written.

**Trade-off:** None — this is a correctness fix, not a tradeoff.

### D5: Stage Bakta DB on C: drive if room; verify before flagging install complete

**Decision:**
- Bakta DB type = `light` (~5 GB) instead of `full` (~50 GB). Light gives gene-symbol assignment which is all Stage 2 needs.
- Download to C: drive first if ≥10 GB free (currently ~23 GB free per `wmic` check). Copy-then-verify to D: only if needed.
- After install, validate DB integrity via `bakta --version` + a smoke-parse on K-12.
- Record `bakta --version`, Bakta DB version + path, `amrfinder --version`, AMRFinder DB version, Mash version, image SHA digests in `wiki/stage2_install_artifact_<date>.md`.

**Rationale:** D: drive USB volatility already cost 8 hours today (HDF5 EOA-truncation from a USB hiccup). Putting Bakta DB + AMRFinder DB + outputs + active populate writes all on D: compounds the same failure mode. C: drive has enough room for the light DB; AMRFinderPlus DB is smaller (~few GB) and can also live on C: or D: depending on remaining space.

**Trade-off:** `light` Bakta DB misses some product annotations + signal peptides that `full` would include. Acceptable for Stage 2 (only gene_symbol matters); upgrade to full at Phase 3 if attribution work needs richer annotations.

### D6: Smoke-test on TWO strains, not one

**Decision:** Validate each tool on K-12 MG1655 (wild-type cipro-S; binary works check) AND one known cipro-R strain from the N=147 cohort (POINT-row parsing actually exercised). Specifically for AMRFinderPlus: assert the `mutations.tsv` output file is non-empty AND contains POINT-method rows AND parses to `(gene, codon, alt_residue)` tuples for the resistant strain.

**Rationale:** K-12 alone proves the binary launches; it does NOT prove POINT-row parsing handles real resistance mutations. A wild-type self-test would silently pass while the point-mutation pipeline downstream is broken.

**Trade-off:** Doubles the smoke-test runtime (still <10 min total). Acceptable.

## Implementation Plan

**Pre-requisite (USER):** Start Docker Desktop GUI app. Wait until system tray icon shows "Engine running." Tell me "go" or message back.

### Step 1: Write `tools/docker_runner.py` (NEW)
- New file. Single Python module exposing:
  ```python
  def run(
      image: str,                          # e.g., "quay.io/biocontainers/mash:2.3--he348c14_4"
      args: list[str],                     # CLI args inside the container
      *,
      mounts: dict[str, str] | None = None, # host_path -> container_path
      env: dict[str, str] | None = None,
      capture_output: bool = True,
      check: bool = True,
  ) -> subprocess.CompletedProcess: ...
  ```
- Normalizes Windows paths (`C:\...` → `/c/...` or `//c/...` depending on Docker volume syntax).
- Emits the exact `docker run` command via DEBUG log for reproducibility.
- 3 unit tests with mocked subprocess: argument-passing correctness, mount-translation correctness, env-passing correctness.

### Step 2: Pull and validate Docker images
- Pull each pinned image:
  - `docker pull quay.io/biocontainers/mash:2.3--he348c14_4`
  - `docker pull oschwengers/bakta:<concrete-version>` (query Docker Hub for current concrete release)
  - `docker pull ncbi/amr:4.2.7-2026-03-24.1`
- Record SHA digests: `docker images --digests | grep ...` → write to `wiki/stage2_install_artifact_<date>.md`.

### Step 3: Download tool databases (light/small variants)
- **Bakta DB (light, ~5 GB)** to `C:/bakta_db` if ≥10 GB free, else `D:/bakta_db`:
  - `docker run --rm -v <bakta_db_path>:/db oschwengers/bakta:<version> bakta_db download --output /db --type light`
- **AMRFinderPlus DB** to `C:/amrfinder_db` or `D:/amrfinder_db`:
  - `docker run --rm -v <amr_db_path>:/db ncbi/amr:4.2.7-2026-03-24.1 amrfinder_update --database /db`
- Mash needs no separate DB.
- Validate each DB: read a version/manifest file inside the DB dir; record version in install artifact.

### Step 4: Smoke-test on K-12 MG1655 (`GCF_000005845.2`)
- **Mash:** sketch + dist against itself; expect distance ~0.
- **Bakta:** run on `D:/dna_decode_cache/refseq/GCF_000005845.2/genome.fna`; output to a `/tmp` mount. Parse output GFF3 with `dna_decode.data.annotations.parse_gff3`. Assert CDS-row `gene_symbol` coverage ≥50% (per the parent→CDS propagation patch already landed at cb46b72).
- **AMRFinderPlus:** run with corrected invocation (`amrfinder -n ... -O Escherichia --database /db --mutation_all /out/mutations.tsv -o /out/main.tsv`). For K-12 (wild-type cipro-S), the mutations.tsv may be empty BUT the file must exist + parse to AMRFinderPlus's standard TSV schema.

### Step 5: Smoke-test on one cipro-R strain from N=147 cohort
- Pick a strain with known cipro resistance from `data/processed/stage2_n150_cipro_cohort.parquet` (e.g., a ST131 R strain).
- Re-run AMRFinderPlus. Assert:
  - `mutations.tsv` is non-empty.
  - Contains POINT-method rows (Method column = `POINT` / `POINTP` / `POINTX` / `POINTN`).
  - At least one row references `gyrA` or `parC` or `parE` (cipro-resistance textbook).
  - Parses to `(gene, codon_pos, alt_residue)` tuples without ValueError.

### Step 6: Mash O(N²) refactor in `phylogeny.py`
- Update `dna_decode/eval/phylogeny.py:compute_mash_distances` to use `mash sketch <fasta1> <fasta2> ... -o sketch.msh` (one call) + `mash dist sketch.msh sketch.msh` (one call) instead of pairwise nested-loop calls.
- Reduces 11,175 invocations at N=147 to 2. With Docker spin-up overhead, this saves ~hours of preflight runtime.
- Adjust to read through `tools/docker_runner.py` so calls slot cleanly.
- +1-2 regression tests: small-cohort distance matrix matches prior nested-loop output within tolerance.

### Step 7: Write the install artifact
- Create `wiki/stage2_install_artifact_<YYYY-MM-DD>.md` containing:
  - Mash version, Bakta version, AMRFinderPlus version
  - Bakta DB version + path
  - AMRFinder DB version + path
  - Image SHA digests for all 3 images
  - Smoke-test results (K-12 + cipro-R strain)
  - Date + machine identifier
- This becomes the audit-trail anchor for Stage 2 reproducibility.

## Verification

- Each tool returns `--version` successfully via `tools/docker_runner.py` invocation.
- Bakta smoke output on K-12: `parse_gff3` reports CDS-row `gene_symbol` coverage ≥50%.
- AMRFinderPlus smoke output on the cipro-R strain: `mutations.tsv` non-empty AND contains gyrA/parC/parE POINT-method rows.
- Mash smoke output: self-distance ~0; refactored `phylogeny.py` produces same matrix as nested-loop reference on a 4-strain subset (tolerance 1e-6).
- `wiki/stage2_install_artifact_<date>.md` exists with all 6 versions + 3 SHA digests + smoke-test PASS verdicts.
- All existing tests still pass: `uv run pytest tests/ -m "not slow" -q` returns ≥421 passed (current baseline) + new docker_runner tests + new phylogeny refactor tests.
- Unblocks Phase B of `plans/Stage2_N150_Prep_Plan.md`: Bakta scale-out + AMRFinderPlus extraction + Mash preflight on the N=147 cohort.

## Executed outcome (2026-05-15)

> **Status: PARTIAL.** Mash + AMRFinderPlus shipped end-to-end. Bakta blocked on a dep-check regression. Mash O(N²) refactor deferred. Full artifact at `wiki/stage2_install_artifact_2026-05-15.md`.

### What landed
- `tools/docker_runner.py` (Step 1) — 9 unit tests pass; supports `timeout`, `:ro` mount-mode pass-through, `FileNotFoundError` + `TimeoutExpired` wrapping as `DockerRunnerError`.
- Mash image pull (Step 2) — **plan deviation**: pinned tag `2.3--he348c14_4` no longer on quay.io; switched to current `2.3--hb105d93_10`. Same Mash binary version (2.3).
- AMRFinderPlus image pull (Step 2) — landed at pinned tag `4.2.7-2026-03-24.1`.
- AMRFinderPlus DB (Step 3) — 240 MB at `C:/Users/Farshad/dna_decode_stage2/amrfinder_db/latest/`. Version `2026-03-24.1`.
- Mash smoke (Step 4) — **plan deviation**: K-12 MG1655 (`GCF_000005845.2`) NOT in local refseq cache (cache is BV-BRC GCA accessions, not RefSeq GCF). Substituted with two existing N=40 cipro-R strains (`GCA_000522345.1` ST131 + `GCA_003073955.1` ST410). Self-distances 0; cross-distance 0.0324. PASS.
- AMRFinderPlus smoke (Step 5) — pinned cipro-R substrate = `1328433.3` / `GCA_000522345.1` (ST131, 5 contigs, N50=5.08 Mbp; D6-spec strain). `mutations.tsv` = 668 rows. **Textbook QRDR confirmed:** `gyrA_S83L` + `gyrA_D87N` + `parC_E84V` all QUINOLONE-class POINT rows. PASS. Closes the H14 SNP-baseline-extraction readiness question (separate from H14's analytical test, which fires at Stage 2 N=150).
- Install artifact (Step 7) — `wiki/stage2_install_artifact_2026-05-15.md`.

### What did NOT land
- **Bakta install (Step 3 partial)** — blocker B1: `bakta_db download` fails AMRFinderPlus dep-check in BOTH `oschwengers/bakta:v1.11.4` AND `v1.12.0`. `which amrfinder` returns `/opt/conda/bin/amrfinder` and `amrfinder --version` returns `4.0.23` inside the v1.11.4 container — the dep check in `bakta_db` does more than `which`. Root cause unknown. Workarounds enumerated in artifact §B1. Deferred — not a Stage 1 verdict blocker; Bakta's consumer is the gene-presence comparator (Pending Decisions row 4) which is itself deferred.
- **K-12 wild-type smoke** — substituted with the textbook-QRDR-positive ST131 cipro-R strain. Higher-value verification (a wild-type self-test would silently pass while a QRDR-detection failure would mean the install is unfit for purpose).
- **Step 6: Mash O(N²) refactor of `dna_decode/eval/phylogeny.py`** — algorithm-touching change; deserves its own commit + targeted regression tests. Deferred to follow-up.

### Operational deviations vs plan
- Plan §D5 specified Bakta DB at `C:/bakta_db` if ≥10 GB free. Actual install path: `C:/Users/Farshad/dna_decode_stage2/{bakta_db,amrfinder_db}` — `C:/<root>` requires admin; `C:/Users/Farshad/` is user-owned.
- Plan implicitly assumed `docker run -v <host>:<container>` works from any shell. Empirically: direct invocation from Git Bash requires `MSYS_NO_PATHCONV=1` env var (the `/db` container path gets munged to `C:/Program Files/Git/db`, silently producing an unbound ephemeral mount where downloads succeed inside the container then disappear on `--rm`). Python `subprocess.run` bypasses this; `tools/docker_runner.py` is unaffected.

### Verification scorecard
| Plan §10 success criterion | Status |
|---|---|
| Mash + Bakta + AMRFinder return `--version` via `tools/docker_runner.py` | PARTIAL: Mash + AMRFinder ✓; Bakta blocked |
| Bakta K-12 smoke: `gene_symbol` coverage ≥50% | BLOCKED on B1 |
| AMRFinderPlus cipro-R smoke: gyrA/parC/parE POINT rows present | ✓ PASS |
| Mash refactored phylogeny.py matches nested-loop reference | DEFERRED to follow-up commit |
| `wiki/stage2_install_artifact_<date>.md` exists with versions + SHAs + verdicts | ✓ PASS |
| All existing tests still pass (≥421 + new) | UNCHECKED in this session — Stage 1 background job holds the working tree; full suite run deferred to post-verdict |
