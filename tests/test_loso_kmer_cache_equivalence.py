"""Equivalence guard: the cached k-mer path in loso_kmer must be feature-level IDENTICAL to the
original build_kmer_vocabulary / kmers_to_feature_matrix (the perf optimization must not change
results). Pins the 2026-06-05 within-fold-rebuild speedup. Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.eval.loso_kmer import (
    _counts_by_strain, _matrix_from_cache, _vocab_from_cache,
)
from dna_decode.models.classical_baselines import (
    build_kmer_vocabulary, kmers_to_feature_matrix,
)

# synthetic genomes with realistic repetition + an N (ambiguous base, must be skipped)
_SEQS = {
    "s1": "ACGTACGTACGTTTGGCCAANNACGTACGT" * 10,
    "s2": "TTGGCCAATTGGCCAAGCGCGCGCACGTAC" * 10,
    "s3": "ACGTACGTGGGGCCCCTTTTAAAACGTACG" * 10,
    "s4": "GCGCGCGCGCGCATATATATACGTACGTAC" * 10,
}
STRAINS = ["s1", "s2", "s3", "s4"]
K = 8
TOP_N = 50


def test_vocab_cache_matches_original():
    cache = _counts_by_strain(_SEQS, STRAINS, K)
    for train in (STRAINS, STRAINS[:3], ["s2", "s4"], ["s3"]):
        orig = build_kmer_vocabulary([_SEQS[s] for s in train], k=K, top_n=TOP_N)
        cached = _vocab_from_cache(train, cache, TOP_N)
        assert cached == orig, f"vocab mismatch for {train}"


def test_matrix_cache_matches_original():
    cache = _counts_by_strain(_SEQS, STRAINS, K)
    vocab = build_kmer_vocabulary([_SEQS[s] for s in STRAINS], k=K, top_n=TOP_N)
    for subset in (STRAINS, ["s1"], ["s4", "s2"]):
        orig = kmers_to_feature_matrix([_SEQS[s] for s in subset], vocab, k=K)
        cached = _matrix_from_cache(subset, cache, vocab)
        assert np.array_equal(cached, orig), f"matrix mismatch for {subset}"


def test_full_loso_fold_equivalence():
    # reproduce one LOSO fold both ways → identical train+test matrices
    cache = _counts_by_strain(_SEQS, STRAINS, K)
    held = "s2"; train = [s for s in STRAINS if s != held]
    v_o = build_kmer_vocabulary([_SEQS[s] for s in train], k=K, top_n=TOP_N)
    v_c = _vocab_from_cache(train, cache, TOP_N)
    assert v_o == v_c
    assert np.array_equal(_matrix_from_cache(train, cache, v_c),
                          kmers_to_feature_matrix([_SEQS[s] for s in train], v_o, k=K))
    assert np.array_equal(_matrix_from_cache([held], cache, v_c),
                          kmers_to_feature_matrix([_SEQS[held]], v_o, k=K))


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
