# ProSST-hybrid prospective validation on MaveDB — DATA WALL (2026-07-21)

**Verdict:** the deployed **ESM2+ProSST hybrid** CANNOT be prospectively validated on the MaveDB held-out
set. This is an EXTERNAL data wall (not code-closable), caught by a cheap pre-bar framing check BEFORE any
Kaggle spend. The sequence-tier prospective result (ESM2-650M 0.503, 90% paired BLOSUM beat) STANDS as the
deployable prospective number; the structure tier is data-walled on this substrate.

## The framing check (real MaveDB + UniProt API, read-only)
ProSST needs a per-protein structure whose token length EQUALS the target sequence length
(`prosst_variant_table` raises on mismatch). For a held-out MaveDB assay that means:
UniProt id (for AlphaFold) → AlphaFold model → structure length == target length. Measured over all 84
scored held-out assays:

| gate | count | note |
|---|---|---|
| has `uniprotIdFromMappedMetadata` | **42/84** | the other 42 have NO UniProt id → no AlphaFold fetch |
| target-len == UniProt-canonical-len | **12/84** | optimistic upper bound; even these need an AF model ≤2700 aa + the seq to BE canonical |

The 72 non-matching assays are **fragments / constructs** — the DMS target is a domain or an engineered
construct, not the full canonical protein (APP tgt-42 vs canon-770; ABL1 tgt-57 vs canon-1130; O94929 tgt-57
vs canon-683…). A full-length AlphaFold structure cannot align to a domain construct without residue-range
cropping + an offset MaveDB does not cleanly expose.

## Why the bar is refused (not looped)
- **Buildable-N ≈ 12 (14%)**, and optimistic. A 12-assay hybrid number does NOT pair with the 84-assay ESM2
  result (different subset) → it cannot be "the deployed scorer's true number ON the held-out set", only a
  tiny non-representative slice.
- **Build cost is large**: the ProSST quantizer needs `torch_geometric` + the cloned ProSST repo + the GVP
  `AE.pt` on Kaggle + AlphaFold fetch per protein + per-assay length reconciliation (the 5 documented
  GNN-stack API traps, memory `reference_kaggle_headless_gpu_kernels`). High cost × small non-comparable N =
  low VOI. This is the low-value build the R2 adversarial-pre-bar check exists to prevent.
- **The lift is already known in-distribution**: the hybrid's +0.067 vs ESM2 is established on ProteinGym
  (`wiki/forward_modality_hybrid_2026-07-17.md`, which ships pre-quantized structures). MaveDB cannot
  *re-measure it prospectively* because its targets lack matching structures.

## Wall classification
- **EXTERNAL (data), not code-closable.** More code cannot invent UniProt ids for the 42 that lack them, nor
  create full-length-matching structures for fragment/construct targets.
- **Conversion condition:** a held-out DMS substrate whose targets are FULL-LENGTH proteins with AlphaFold
  models (e.g. a curated set of full-length human clinical genes with DMS + structures). That is a new
  substrate-acquisition task, NOT this bar — surfaced for a user decision, not auto-started.
- **Small-N salvage (optional, user call):** a hybrid number on the ~12 length-matched assays is buildable
  but small-N + non-comparable to the 84. Not worth the Kaggle GNN-stack build autonomously; flagged, not done.

## What stands
The R2 molecular cell's prospective evidence remains: leakage-free 86-assay holdout + ESM2-650M **0.503**
(matches the field) + a **90% paired win over BLOSUM** (p=5e-15). That is the SEQUENCE-tier deployed path,
prospectively validated. The STRUCTURE tier's contribution is documented in-distribution (+0.067 on
ProteinGym) and is honestly out of prospective reach on MaveDB.
See `wiki/mavedb_holdout_esm2_2026-07-21.md` + `wiki/mavedb_esm_vs_blosum_paired_2026-07-21.md`.
