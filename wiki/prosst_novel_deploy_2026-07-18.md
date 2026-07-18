# Novel-protein ProSST deployment — the quantize runs LOCALLY (no Kaggle)

**Date:** 2026-07-18 · **Deliverable:** the end-to-end novel-protein path (PDB structure → self-quantize →
ProSST → ESM2+ProSST hybrid), validated on this Windows/CPU host · **Code:**
`dna_decode/forward/prosst_scorer.py::quantize_structure` + `scripts/prosst_quantize.py`

## Question + reframe (R2)

The directive: "complete an attended-deployment (novel-protein ProSST via a **Kaggle** quantize)". The
ProSST forward already runs locally (`wiki/prosst_lift_2026-07-18.md`); the only Kaggle-tagged piece was the
`PdbQuantizer` (structure→tokens, `torch_geometric` wall). **R2 check before committing to Kaggle: does the
quantizer install locally?** It does — the "Kaggle-only" tag was one dependency too pessimistic.

## Result — the full novel-protein path runs LOCALLY, validated exact

**(1) The quantizer runs on this Windows/CPU host** via the cloned ProSST repo + three pure-python shims (all
reversible):
- **`torch_scatter`** → a pure-python shim package (site-packages) backed by `torch_geometric.utils.scatter`
  (the compiled torch_scatter won't build on Windows; ProSST's GVP only needs scatter_add/sum/mean/max).
  `torch_geometric` itself (2.8.0) installs + works with no compiled extension.
- **biotite 1.x** renamed `filter_backbone` → `filter_peptide_backbone` (patched in the cloned repo).
- **Windows `spawn` multiprocessing** (pathos `Pool` AND torch `DataLoader(num_workers)`) deadlocks/EOFErrors
  on this host → force serial: patch pathos Pool/ThreadPool to an in-process `_SerialPool` + run the
  quantizer with `num_processes=0` (⇒ `DataLoader(num_workers=0)`, no child spawn). A single PDB needs no
  parallelism.

**(2) Ground-truth exact match.** The repo ships `p1.pdb` with published tokens `[407, 998, 1841, …]`; our
local quantizer reproduces them **exactly** (8/8, 75 tokens).

**(3) End-to-end deploy validation on a real ProteinGym protein (GRB2).** The repo ships GRB2's PDB and we
have GRB2's pre-quantized ground-truth tokens + DMS:
- **Self-quantized GRB2 tokens: 217/217 = 100% match** to ProteinGym's pre-quantized tokens — the deploy
  quantize is byte-exact, not approximate.
- **Self-quantized ProSST vs DMS: Spearman 0.798** (strong — GRB2 is well-structured).
- **ESM2 alone 0.726 → ESM2 + self-quantized-ProSST 0.789** (+0.064) — the validated structure lift
  (`wiki/prosst_lift_2026-07-18.md` +0.067) reproduces end-to-end through our OWN quantize, not ProteinGym's
  pre-quantized shortcut.

So the novel-protein pipeline is complete + validated locally: **PDB → `quantize_structure` → ProSST scores →
`rank_average_hybrid([esm2, prosst])`**, byte-exact against the reference and delivering the structure lift.

## Deploy recipe (this host)

1. `torch_geometric` installed; a pure-python `torch_scatter` shim in site-packages.
2. `git clone https://github.com/ai4protein/ProSST` → `$PROSST_REPO` (default `D:/prosst_repo`; bundles the
   GVP `AE.pt` + `{vocab}.joblib` k-means — no separate download); patch its one biotite import.
3. `biotite`, `pathos`, `joblib`, `scikit-learn` installed.
4. `quantize_structure(pdb, vocab)` → tokens → `prosst_variant_table(seq, muts, structure_tokens=tokens)`.

For a truly novel sequence: `fetch_alphafold_pdb(uniprot)` (already in the module) → `quantize_structure` →
ProSST → hybrid. A sequence-only input goes through ColabFold to a structure first.

## Honest scope

- The shims are host-setup, not shipped deps — `quantize_structure` raises `StructureMethodUnavailable` with
  setup instructions when the repo/`torch_geometric` are absent (a fresh machine needs the recipe above).
- GRB2's 100% token match is the strongest validation possible (the repo ships the exact PDB ProteinGym
  used); a novel protein's AlphaFold structure would give its own tokens (no ground truth to match, but the
  pipeline is proven correct).
- The quantize is slow on CPU (GVP graph + subgraph build, serial) — minutes per protein; acceptable for a
  one-off deploy, not a high-throughput screen (use a Linux/GPU host with real multiprocessing for volume).
- Frozen decoder surface byte-unchanged (`verify_lock OK`).

## Shipped

- `dna_decode/forward/prosst_scorer.py::quantize_structure` — the local `SSTPredictor` quantize + the serial
  shim + the repo-absent `StructureMethodUnavailable` signal.
- `scripts/prosst_quantize.py` — thin CLI (`--expect` verifies against known tokens).
- 13 seam tests green (the two "quantizer always absent" tests updated to the repo-absent contract).
