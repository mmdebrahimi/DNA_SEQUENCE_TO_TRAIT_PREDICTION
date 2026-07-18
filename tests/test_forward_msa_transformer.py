"""Offline tests for the MSA-Transformer module's pure helpers (dna_decode/forward/msa_transformer.py).

The model forward pass needs torch + fair-esm (heavy, GPU/CPU) and is NOT exercised in CI; the real run is
`scripts/msa_transformer_lift.py` (validated 2026-07-17). These pin the deterministic subsampling + the
focus-residue walk that the variant-position mapping depends on.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward.msa_transformer import _focus_residues, subsample_msa


def test_subsample_keeps_focus_and_is_deterministic():
    cols = [f"seq{i}" for i in range(1000)]
    a = subsample_msa(cols, 128)
    b = subsample_msa(cols, 128)
    assert a == b                       # deterministic (no RNG)
    assert len(a) == 128
    assert a[0] == cols[0]              # focus (row 0) always kept


def test_subsample_below_depth_returns_all():
    cols = ["a", "b", "c"]
    assert subsample_msa(cols, 128) == cols


def test_subsample_is_a_subset_in_order():
    cols = [str(i) for i in range(500)]
    s = subsample_msa(cols, 64)
    idxs = [cols.index(x) for x in s]
    assert idxs == sorted(idxs)         # strided -> increasing original order
    assert idxs[0] == 0


def test_focus_residues_counts_upper_and_lower_as_query_positions():
    # MKtAY: M,K match (pos1,2), t insert (pos3), A,Y match (pos4,5) -> all 5 are query residues
    r = _focus_residues("MKtAY")
    assert r == {1: "M", 2: "K", 3: "T", 4: "A", 5: "Y"}


def test_focus_residues_skips_alignment_gaps():
    # a '-' in the focus is not a query residue
    assert _focus_residues("M-KA") == {1: "M", 2: "K", 3: "A"}


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
