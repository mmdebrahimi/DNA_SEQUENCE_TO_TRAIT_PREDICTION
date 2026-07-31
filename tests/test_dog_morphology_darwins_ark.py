"""Offline tests for the dog morphology validator + the pinned single-SNP MORPH_LOCI (ear).

No genotype data / no network — pins the pure `pearson` helper and the MORPH_LOCI['EAR'] catalog entry
(well-formed canFam4 id, high-allele is the ref/alt, functional_r in range, Q125 provenance). The real
r=+0.543 ear validation ran on Darwin's Ark N=2834; see the wiki artifact. Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.pigment.dog_body_size import MORPH_LOCI, reference_integrity_ok  # noqa: E402
from scripts.dog_morphology_darwins_ark_validate import pearson  # noqa: E402


def test_pearson_perfect_and_none():
    assert abs(pearson([0, 1, 2, 3], [0, 1, 2, 3]) - 1.0) < 1e-12
    assert abs(pearson([0, 1, 2, 3], [3, 2, 1, 0]) + 1.0) < 1e-12
    assert pearson([1, 1, 1], [1, 2, 3]) is None       # zero variance in x
    assert pearson([5], [5]) is None                   # too few points


def test_pearson_monotonic_dose_matches_measured_sign():
    # a synthetic monotonic dose->ordinal (like the ear locus) yields a strong positive r
    xs = [0] * 5 + [1] * 5 + [2] * 5
    ys = [-1.1] * 5 + [-0.3] * 5 + [0.4] * 5
    r = pearson(xs, ys)
    assert r is not None and r > 0.9


def test_ear_locus_catalog_wellformed():
    assert "EAR" in MORPH_LOCI
    ear = MORPH_LOCI["EAR"]
    chrom, pos, ref, alt = ear.canfam4_variant.split(":")
    assert chrom == "chr10" and pos == "8612500"
    assert ear.high_allele in (ref, alt)
    assert 0.0 < ear.functional_r < 1.0
    assert ear.darwins_ark_question == "Q125" and ear.gene == "MSRB3"


def test_reference_integrity_covers_morph_loci():
    # the size-catalog integrity guard now also validates MORPH_LOCI (ear) well-formedness
    assert reference_integrity_ok() is True


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
