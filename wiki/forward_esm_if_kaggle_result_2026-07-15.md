# ESM-IF real number — the structure-method wall converted via Kaggle T4 (2026-07-15)

The ESM-IF seam (row 436) was externally walled on this Windows/CPU host (`torch_scatter` has no wheel for
torch 2.12.1+cpu). This run **converts that wall** on a Kaggle T4 GPU kernel and produces the real number the
seam was missing — closing the "run the real ESM-IF" loop.

## Result — PTEN, real ESM-IF, full 7,260 variants

| method | Spearman(pred, DMS) | source |
|---|---:|---|
| BLOSUM62 | 0.182 | deterministic |
| **ESM-IF (structure)** | **0.479** | `esm_if1_gvp4_t16_142M_UR50` on AlphaFold `AF-P60484-F1-model_v6`, Kaggle T4 |
| ESM2-650M (sequence) | 0.518 | zero-shot masked-marginal |
| AlphaMissense | 0.539 | learned, human-only |

n=7260, WT-mismatch=0, ~42 min GPU (2494 s). Running Spearman peaked ~0.62 mid-cohort and settled at **0.479**.

## The honest finding — structure does NOT beat sequence on PTEN

ESM-IF (0.479) lands **below** sequence-based ESM2 (0.518) and AlphaMissense (0.539). The expected
"+structure" margin (~+0.01–0.02 on the ProteinGym *average*) **does not materialize for PTEN** — a real,
per-protein result, consistent with ProteinGym reality where ESM-IF's median (~0.47) sits around ESM2's, and
the structure tier is led by ProSST/SaProt, not ESM-IF. So the 4th method is validated + real, but it is not
a lift here; it complements rather than tops the leaderboard on this protein.

## A latent sign bug the real run caught (verify-in-batch)

`score_sequence` returns the average per-token **log-likelihood** (higher = better). The local seam's
`ll()` **negated** it (`return -float(loss)`), which would have flipped the correlation sign (−0.479 instead
of +0.479). The Kaggle run's inline GPU-safe scorer had the sign right (positive running Spearman confirmed
it live); `dna_decode/forward/structure_scorer.py` is now fixed to match (no negation). The seam had never
run locally (deps absent), so this bug was latent until the real GPU run.

## Kaggle conversion — the reusable fixes (5 walls cleared in sequence)

The GNN stack that fails on Windows/CPU installs cleanly on Kaggle's Linux T4; the run then hit a chain of
API-version traps, each fixed:
1. **torch_scatter/torch_sparse** — install from the PyG wheel index matching Kaggle's torch (`-f
   https://data.pyg.org/whl/torch-{ver}+{cu}.html`); on Kaggle this is a prebuilt wheel (the Windows wall gone).
2. **biotite ≥1.0** renamed `filter_backbone` → `filter_peptide_backbone`; fair-esm imports the old name.
   Fix = keep latest biotite (numpy-2 compatible) + shim `biotite.structure.filter_backbone =
   filter_peptide_backbone` BEFORE importing esm (pinning `biotite<1.0` instead hits a numpy-1/2 ABI wall).
3. **AlphaFold PDB URL** — the `-model_v4` guess 404s; query `https://alphafold.ebi.ac.uk/api/prediction/
   {uniprot}` for the exact `pdbUrl` (P60484 is v6).
4. **ESM-IF loader name** — the installed fair-esm spells it `esm_if1_gvp4_t16_142M_UR50` (not
   `esm_if1_gvp_transformer_t16_142M_UR50`); resolve it by searching `dir(esm.pretrained)` for `if1`.
5. **GPU device mismatch** — the util `score_sequence` leaves tokens on CPU; call
   `CoordBatchConverter(alphabet)([(coords,None,seq)], device=model.device)` + `model.forward` manually.

Kernel: `emanueleebrahimi/esm-if-pten-forward` (private T4). Result JSON:
`wiki/tem1_forward_cell_pten_human_mighell_2018_esm_if_2026-07-15.json`. Leaderboard regenerated → PTEN is
now the full **4-method** row (BLOSUM/ESM2/AlphaMissense/ESM-IF). Frozen decoder surface byte-unchanged
(`verify_lock OK`); `dna_decode/forward` NON-frozen.
