"""Structure-based method for the forward variant-effect predictor — ESM-IF (inverse folding).

ESM-IF1 (Hsu et al. 2022) scores P(sequence | backbone structure): the conditional log-likelihood of each
residue given the 3D backbone. Its per-variant score = LL(mutant | structure) − LL(wild-type | structure)
(higher = more structure-compatible = more preserved), the same sign convention as BLOSUM / ESM2 / (1−AM).
On ProteinGym it is a TOP-TIER predictor (~structure adds a small margin over sequence-only ESM2).

DEPENDENCY REALITY (this host, 2026-07-14): ESM-IF needs `torch_geometric` (+ historically `torch_scatter` /
`torch_sparse`) for its GVP-GNN encoder + `biotite` for structure parsing — NONE installed, and these are
hard to build on Windows/CPU against a bleeding-edge torch. So `esm_if_variant_table` LAZY-imports and raises
`StructureMethodUnavailable` with a clear message when the stack is absent; the SEAM (predict_effect
method='esm_if' + the leaderboard column) is complete + mock-tested, and the real forward pass runs the
moment the deps are provisioned (a Linux/GPU host). Structures come from the AlphaFold DB by UniProt.
"""
from __future__ import annotations

from pathlib import Path


class StructureMethodUnavailable(RuntimeError):
    """Raised when the ESM-IF / structure stack (torch_geometric / biotite / esm.inverse_folding) is absent."""


# ESM-IF conditional-LL delta tiers (log-likelihood scale, like ESM2's masked-marginal).
_IF_PRESERVED = -1.0
_IF_DAMAGING = -3.0


def alphafold_pdb_url(uniprot: str) -> str:
    """AlphaFold DB predicted-structure PDB URL for a UniProt accession (v4)."""
    return f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"


def fetch_alphafold_pdb(uniprot: str, out_dir: Path) -> Path:
    """Download the AlphaFold predicted structure for `uniprot` -> a local .pdb path (cached)."""
    import urllib.request
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"AF-{uniprot}-F1-model_v4.pdb"
    if not dest.exists():
        urllib.request.urlretrieve(alphafold_pdb_url(uniprot), dest)
    return dest


def _load_esm_if():
    try:
        import esm  # noqa: F401
        import esm.inverse_folding  # noqa: F401
    except Exception as e:  # ModuleNotFoundError: torch_geometric / torch_scatter / biotite
        raise StructureMethodUnavailable(
            f"ESM-IF unavailable — the structure stack is not installed ({type(e).__name__}: {e}). "
            f"Install torch_geometric (+ torch_scatter/torch_sparse) + biotite on a compatible host "
            f"(Linux/GPU recommended), then re-run.") from e
    import esm
    model, alphabet = esm.pretrained.esm_if1_gvp_transformer_t16_142M_UR50()
    return model.eval(), alphabet


def esm_if_variant_table(pdb_path: Path, wt_seq: str, mutants, chain: str = "A") -> dict[str, float]:
    """{DMS mutation 'wt{pos}alt' -> ESM-IF conditional-LL delta (mutant − wild-type)} given the backbone.

    Raises StructureMethodUnavailable if the ESM-IF stack is absent (this host). Autoregressive scoring:
    LL(seq | structure) per full sequence; delta = LL(mut) − LL(wt) per single-point variant.
    """
    model, alphabet = _load_esm_if()                       # raises if deps missing
    import esm.inverse_folding as ifold
    structure = ifold.util.load_structure(str(pdb_path), chain)
    coords, native = ifold.util.extract_coords_from_structure(structure)

    def ll(seq: str) -> float:
        # average per-residue conditional log-likelihood of `seq` given the backbone
        loss, _ = ifold.util.score_sequence(model, alphabet, coords, seq)
        return -float(loss)

    base = ll(wt_seq)
    table: dict[str, float] = {}
    for m in mutants:
        m = m.strip()
        if ":" in m or len(m) < 3 or not m[1:-1].isdigit():
            continue
        wt, pos, alt = m[0], int(m[1:-1]), m[-1]
        if pos > len(wt_seq) or wt_seq[pos - 1] != wt:
            continue
        mut_seq = wt_seq[:pos - 1] + alt + wt_seq[pos:]
        table[m] = ll(mut_seq) - base
    return table


def esm_if_tier(delta: float) -> str:
    """ESM-IF conditional-LL delta -> forward-cell tier (higher = structure-compatible = preserved)."""
    if delta >= _IF_PRESERVED:
        return "preserved"
    if delta <= _IF_DAMAGING:
        return "damaging"
    return "uncertain"
