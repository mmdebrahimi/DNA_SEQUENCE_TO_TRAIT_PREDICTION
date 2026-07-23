# ProSST/GEMME full-hybrid resource wall — CLOSED (2026-07-23)

**Directive:** close the ProSST/GEMME full-hybrid packaging wall — the strong `--method prosst/hybrid` path
needed torch_geometric + a torch_scatter shim + biotite + pathos + the cloned AI4Protein/ProSST repo (heavy,
platform-specific, not cleanly pip-packageable), and C: is disk-tight. We have D: (4.3 TB free).

## Verdict: WALL_CLOSED

Real-surface proof on this exact Windows/CPU host (`scripts/prosst_hybrid_verify.py`, deps installed
ephemerally on top of `[forward]`, all caches on D: to spare the 97%-full C:):

| Check | Result |
|---|---|
| **Local quantizer reproduction** (torch_geometric GVP + k-means on a raw AlphaFold PDB → structure tokens, vs the committed manifest's reference tokens) | **79/79 = 1.0000** on CKS1B (P61024) |
| **ProSST transformer non-degenerate** (real forward on the locally-quantized tokens) | log-ratio **stdev 2.197** (degenerate failure is ~0; e.g. Q5A = −3.101 damaging) |

A perfect token reproduction proves the local torch_geometric quantizer reproduces the reference *exactly*
(the generalized "217/217"), and a non-degenerate ProSST table proves the transformer path runs end-to-end
from a raw structure — the two things the resource wall blocked.

## What was installed / where

- **`D:/prosst_repo`** — the cloned AI4Protein/ProSST repo (already present; bundles the GVP `AE.pt` +
  `{vocab}.joblib` k-means checkpoints). `PROSST_REPO` env points here.
- **torch_geometric + biotite + pathos** — the pip-able parts, now declared as the **`[prosst]` extra**
  (`uv sync --extra forward --extra prosst`). torch_scatter: recent torch_geometric uses torch-native
  scatter, so a pure-python shim suffices on Windows/CPU.
- Caches (ProSST-2048 LM, uv overlay) on **D:** via `HF_HOME` / `UV_CACHE_DIR` — C: untouched.

## Packaging / capability wiring

- **`[prosst]` extra** added to `pyproject.toml` (composes with `[forward]`; orthogonal to the `[ml]↔[forward]`
  transformers conflict).
- **Capability probe fixed**: it checked `find_spec("prosst")` (always False — ProSST is a *repo*, not a pip
  package). Now `probe_capabilities()` detects `torch_geometric` + the `$PROSST_REPO` entrypoint, so
  `dna-decode forward --capabilities` honestly reports prosst runnable once the extra + repo are present.
- **GEMME** needs no packaging — it is the Docker image `elodielaine/gemme:gemme` (validated in the 3-way
  last session, TEM-1 0.719); the probe already reports `gemme` runnable when Docker is present.

## Honest scope

- The `[forward]`+`[prosst]` deps are installed in an **ephemeral overlay** for the verification; a user who
  wants `--method prosst/hybrid` persistently runs `uv sync --extra forward --extra prosst` + clones the repo.
- The real ESM2+ProSST **hybrid combine** is a pure `rank_average` of two tables (already unit-tested) and
  ESM2 was separately proven (raw score −2.966); co-loading ESM2-650M *and* ProSST-2048 in one process OOMs
  this host's paging file, so a full one-process hybrid number needs a two-process (compute-each-table-then-
  combine) split — a runner detail, not a wall. The wall (the ProSST structure path) is closed.
