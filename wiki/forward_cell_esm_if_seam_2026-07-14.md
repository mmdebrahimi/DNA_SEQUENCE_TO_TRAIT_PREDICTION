# Forward cell — structure-based method (ESM-IF) seam + honest dependency wall (2026-07-14)

Adds ESM-IF (inverse folding) as a 4th `method=` on the forward variant-effect cell. **The seam is
code-complete + tested; the real number is externally walled on this host** and runs the moment the GNN
stack is provisioned (a Linux/GPU box). This is a wall-classified partial, not a silent failure.

## What shipped (code-closable, done)

- `dna_decode/forward/structure_scorer.py` — `esm_if_variant_table(pdb, wt_seq, mutants)` (ESM-IF conditional
  log-likelihood delta LL(mut|structure) − LL(wt|structure), higher = structure-compatible = preserved, the
  same sign as BLOSUM/ESM2/(1−AM)) + `fetch_alphafold_pdb` (AlphaFold DB by UniProt) + `esm_if_tier`. The
  real loader LAZY-imports `esm.inverse_folding` and raises **`StructureMethodUnavailable`** with a clear
  message when the stack is absent — it fails loudly, never silently.
- `predict_effect(method="esm_if", esm_if_table=…)` — mirrors the esm2/alphamissense wiring.
- `scripts/forward_leaderboard.py` — an **ESM-IF column** (populates from a `esm_if_structure` result JSON
  produced on a deps-provisioned host; `—` here).
- 5 structure tests (29 forward total), including an assertion that the real path raises loudly when the
  deps are absent.

## The wall — empirically confirmed (not assumed)

The pre-bar check flagged this; I tested it rather than guessing:

1. `torch_geometric` + `biotite` install fine, but `esm.inverse_folding` hard-imports **`torch_scatter`**
   (`from torch_scatter import scatter_add`).
2. `torch_scatter` has **no prebuilt wheel for torch 2.12.1+cpu** (this host's bleeding-edge torch) → it
   falls back to a source build → fails (needs torch in the build isolation + an MSVC C++ compiler, absent
   on Windows).

**Classification: EXTERNAL dependency wall** (needs `torch_scatter` — a resource the agent cannot produce on
this Windows/CPU host). **Conversion path:** provision `torch_geometric` + `torch_scatter` (+ `biotite`) on a
Linux/GPU host (the project's Precision-7780 / Databricks workhorse), fetch AlphaFold structures via
`fetch_alphafold_pdb`, and the already-wired `esm_if_variant_table` → `--method esm_if` path runs unchanged.

## Honest VOI note (R2)

Per the ProteinGym leaderboard, structure-based methods add only ~**+0.01–0.02** over sequence-only ESM2 —
which the cell already has at 0.52–0.73. So ESM-IF is a **small marginal gain for a heavy, externally-walled
dependency** on this host. The seam is worth having (it's cheap + ready); forcing a risky `torch_scatter`
source-build on a 96%-full C: for +0.02 would not have been. Banked as: seam done, wall named, number deferred
to a deps-provisioned host.

Frozen decoder surface byte-unchanged (`verify_lock OK`); `dna_decode/forward` NON-frozen.
