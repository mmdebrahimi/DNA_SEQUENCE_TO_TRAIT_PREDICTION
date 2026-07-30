"""Offline tests for the dog body-SIZE catalog (dna_decode/pigment/dog_body_size.py).

Synthetic dosages only — NO genotype data / no network. Pins: the pinned+validated catalog integrity,
the additive polygenic score, the RELATIVE size-rank gradient, confidence tiers, and abstention on
v0-unmodelled structural loci. The measured functional_r values here are provenance constants asserted
against the committed catalog (the real validation ran on Darwin's Ark N=3276; see the wiki artifact).
Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.pigment.dog_body_size import (  # noqa: E402
    POLYGENIC_SCORE_R,
    SIZE_LOCI,
    SizeInputError,
    polygenic_size_score,
    reference_integrity_ok,
)


def test_reference_integrity():
    assert reference_integrity_ok() is True


def test_catalog_pinned_variants_wellformed():
    # exactly the 4 dominant loci, each with a canFam4 chrN:pos:ref:alt id whose alleles are big/small
    assert set(SIZE_LOCI) == {"IGF1", "HMGA2", "STC2", "GHR"}
    for L in SIZE_LOCI.values():
        chrom, pos, ref, alt = L.canfam4_variant.split(":")
        assert chrom.startswith("chr") and pos.isdigit()
        assert {L.big_allele, L.small_allele} == {ref, alt}
        assert 0.0 < L.functional_r < 1.0


def test_combined_beats_best_single_locus():
    # the whole point: the additive score correlates with height better than any single SNP
    assert POLYGENIC_SCORE_R > max(L.functional_r for L in SIZE_LOCI.values())


def test_all_small_is_toy():
    c = polygenic_size_score({k: 0 for k in SIZE_LOCI})
    assert c.polygenic_score == 0 and c.max_score == 8
    assert c.size_rank == "toy/small" and c.confidence == "high"


def test_all_big_is_large():
    c = polygenic_size_score({k: 2 for k in SIZE_LOCI})
    assert c.polygenic_score == 8 and c.size_rank == "large/giant"


def test_monotonic_gradient():
    scores = [polygenic_size_score({k: d for k in SIZE_LOCI}).polygenic_score for d in (0, 1, 2)]
    assert scores == [0, 4, 8]
    ranks = [polygenic_size_score({k: d for k in SIZE_LOCI}).size_rank for d in (0, 1, 2)]
    assert ranks[0] == "toy/small" and ranks[2] == "large/giant"


def test_partial_panel_shrinks_max_and_caps_confidence():
    # only the two secondary loci scored -> low confidence (no dominant IGF1/HMGA2)
    c = polygenic_size_score({"STC2": 2, "GHR": 2})
    assert c.n_loci_scored == 2 and c.max_score == 4 and c.polygenic_score == 4
    assert c.confidence == "low"


def test_one_dominant_locus_is_medium():
    c = polygenic_size_score({"IGF1": 1, "STC2": 1})
    assert c.confidence == "medium"


def test_unmodelled_structural_locus_abstains():
    c = polygenic_size_score({"IGF1": 2, "HMGA2": 2, "STC2": 2, "GHR": 2}, present_unmodelled=["FGF4"])
    assert any("FGF4" in a for a in c.abstains_on)
    assert c.confidence == "medium"  # capped down from high by the abstention


def test_unmodelled_locus_in_dosages_raises():
    with pytest.raises(SizeInputError):
        polygenic_size_score({"IGF1": 2, "FGF4": 1})


def test_bad_dosage_and_empty_raise():
    with pytest.raises(SizeInputError):
        polygenic_size_score({"IGF1": 3})
    with pytest.raises(SizeInputError):
        polygenic_size_score({})
    with pytest.raises(SizeInputError):
        polygenic_size_score({"NOSUCH": 1})


def test_as_dict_is_honest_about_relative():
    d = polygenic_size_score({"IGF1": 2, "HMGA2": 2, "STC2": 1, "GHR": 1}).as_dict()
    assert d["trait"] == "body_size" and d["regime"].startswith("A_curated_catalog_additive")
    assert "RELATIVE" in d["measure"]
    assert d["validation"]["polygenic_score_r"] == POLYGENIC_SCORE_R


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
            passed += 1
        except Exception as e:  # noqa: BLE001
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{passed}/{len(fns)} passed")
