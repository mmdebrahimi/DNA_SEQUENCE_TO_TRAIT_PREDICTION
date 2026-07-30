"""Tests for the general N-category multinomial pigment engine (dna_decode/pigment/multinomial.py).

The load-bearing test EQUIVALENCE: express the shipped IrisPlex EYE model as a PigmentModel and assert the
general engine reproduces `predict_eye_color`'s probabilities + call EXACTLY across a genotype grid. That
validates the engine against a known-correct, 1000G-population-validated model — so hair/skin only need
their coefficient tables filled (attended transcription), never fabricated. Runnable via pytest OR standalone.
"""
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.pigment.irisplex import _CONSTANT, IRISPLEX_SNPS, predict_eye_color  # noqa: E402
from dna_decode.pigment.multinomial import (  # noqa: E402
    MissingGenotypeError,
    PigmentModel,
    PigmentSNP,
    predict,
    reference_integrity_ok,
)


def _eye_model() -> PigmentModel:
    """The shipped IrisPlex eye model expressed as a general PigmentModel (blue = reference)."""
    snps = tuple(
        PigmentSNP(rsid=rsid, counted_allele=allele,
                   betas={"intermediate": b_int, "brown": b_brown},
                   required=(rsid == "rs12913832"))
        for rsid, allele, b_int, b_brown in IRISPLEX_SNPS
    )
    return PigmentModel(
        trait="eye_colour", categories=("blue", "intermediate", "brown"),
        intercepts={"intermediate": _CONSTANT[0], "brown": _CONSTANT[1]},
        snps=snps, source="IrisPlex Walsh 2011 (via irisplex.py)")


def _grid_genotypes():
    """Genotype dicts spanning x=0/1/2 at each of the 6 IrisPlex SNPs (a sample of the space)."""
    def gt(allele, x):
        other = next(b for b in "ACGT" if b != allele)
        return {0: other + other, 1: allele + other, 2: allele + allele}[x]
    # vary the two highest-impact SNPs fully, hold the rest at a few states
    hi = IRISPLEX_SNPS[0]  # HERC2
    for xs in itertools.product([0, 1, 2], repeat=3):  # HERC2 + next 2
        g = {}
        for i, (rsid, allele, _, _) in enumerate(IRISPLEX_SNPS):
            g[rsid] = gt(allele, xs[i] if i < 3 else 0)
        yield g


def test_engine_reproduces_shipped_eye_model():
    model = _eye_model()
    n = 0
    for g in _grid_genotypes():
        got = predict(model, g)
        ref = predict_eye_color(g)
        assert got.call == ref.call, (g, got.call, ref.call)
        assert abs(got.probabilities["blue"] - ref.p_blue) < 1e-9
        assert abs(got.probabilities["intermediate"] - ref.p_intermediate) < 1e-9
        assert abs(got.probabilities["brown"] - ref.p_brown) < 1e-9
        n += 1
    assert n == 27


def test_reference_integrity_eye_anchors():
    model = _eye_model()

    def g(herc2):
        return {rsid: (herc2 if rsid == "rs12913832" else next(b for b in "ACGT" if b != allele) * 2)
                for rsid, allele, _, _ in IRISPLEX_SNPS}
    assert reference_integrity_ok(model, [(g("GG"), "blue"), (g("AA"), "brown")]) is True
    # a corrupted anchor expectation fails loudly
    assert reference_integrity_ok(model, [(g("GG"), "brown")]) is False


def test_two_category_softmax():
    m = PigmentModel(trait="t", categories=("ref", "alt"), intercepts={"alt": 0.0},
                     snps=(PigmentSNP("rs1", "A", {"alt": 1.0}),), source="synthetic")
    # x=0 -> Z_alt=0 -> equal 0.5/0.5
    p0 = predict(m, {"rs1": "TT"}).probabilities
    assert abs(p0["ref"] - 0.5) < 1e-12 and abs(p0["alt"] - 0.5) < 1e-12
    # x=2 -> Z_alt=2 -> alt favored
    p2 = predict(m, {"rs1": "AA"}).probabilities
    assert p2["alt"] > p2["ref"] and abs(sum(p2.values()) - 1.0) < 1e-5


def test_probs_sum_to_one_four_categories():
    m = PigmentModel(trait="hair", categories=("black", "brown", "blond", "red"),
                     intercepts={"brown": 0.3, "blond": -0.2, "red": -1.0},
                     snps=(PigmentSNP("rs1", "A", {"brown": 0.5, "blond": 0.1, "red": 2.0}),
                           PigmentSNP("rs2", "T", {"brown": -0.4, "blond": 0.9, "red": -0.3})),
                     source="synthetic")
    p = predict(m, {"rs1": "AA", "rs2": "TG"}).probabilities
    assert abs(sum(p.values()) - 1.0) < 1e-5
    assert set(p) == {"black", "brown", "blond", "red"}


def test_pending_stub_refuses():
    stub = PigmentModel(trait="skin", categories=("very_pale", "pale", "intermediate", "dark", "dark_black"),
                        intercepts={}, snps=(), source="Walsh 2017 Table 2 (pending transcription)",
                        coefficients_pending=True)
    try:
        predict(stub, {"rs1426654": "AA"})
    except ValueError as e:
        assert "pending" in str(e) or "not yet populated" in str(e)
        return
    raise AssertionError("a coefficients_pending stub must refuse to predict")


def test_required_snp_and_allow_missing():
    m = _eye_model()
    # missing HERC2 (required) -> always raises
    g = {rsid: allele * 2 for rsid, allele, _, _ in IRISPLEX_SNPS if rsid != "rs12913832"}
    try:
        predict(m, g, allow_missing=True)
    except MissingGenotypeError:
        pass
    else:
        raise AssertionError("required HERC2 must raise even under allow_missing")
    # missing a NON-required SNP without allow_missing -> raises; with allow_missing -> imputes x=0, low conf
    g2 = {rsid: (allele * 2) for rsid, allele, _, _ in IRISPLEX_SNPS}
    del g2["rs1800407"]
    try:
        predict(m, g2)
    except MissingGenotypeError:
        pass
    else:
        raise AssertionError("missing non-required SNP without allow_missing must raise")
    assert predict(m, g2, allow_missing=True).confidence == "low"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
