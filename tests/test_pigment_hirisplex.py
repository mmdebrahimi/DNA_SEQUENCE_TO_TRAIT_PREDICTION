"""Tests for the recovered HIrisPlex-S hair/skin/eye models (dna_decode/pigment/hirisplex_models.py).

The load-bearing test PINS the loaded models against the deployed webtool's ACTUAL all-zero-genotype
outputs (the reference point the coefficients were recovered from) + the held-out validation the recovery
recorded. A corrupted coefficient JSON fails loudly. Offline / no network. Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.pigment.hirisplex_models import load_hirisplex_models, provenance  # noqa: E402
from dna_decode.pigment.multinomial import predict  # noqa: E402

# The webtool's ACTUAL predicted probabilities at the all-zero genotype (x=0 for all 41 SNPs),
# read from the extraction's out_full.csv ZERO row. These pin that the loaded softmax reproduces the
# deployed model at its baseline. (Category order = model category order.)
WEBTOOL_ZERO = {
    "eye_colour":  {"blue": 0.847810, "intermediate": 0.087663, "brown": 0.064527},
    "hair_colour": {"blond": 0.705643, "brown": 0.254482, "red": 0.002024, "black": 0.037850},
    "skin_colour": {"very_pale": 0.008511, "pale": 0.391839, "intermediate": 0.598458,
                    "dark": 0.001192, "dark_to_black": 0.0},
}


def _zero_genotype(model):
    """A genotype dict giving x=0 (counted allele absent) at every SNP of the model."""
    out = {}
    for snp in model.snps:
        other = next(b for b in "ACGT" if b != snp.counted_allele)
        out[snp.rsid] = other + other
    return out


def test_models_load_with_expected_shape():
    models = load_hirisplex_models()
    assert set(models) == {"eye_colour", "hair_colour", "skin_colour"}
    assert len(models["eye_colour"].categories) == 3
    assert len(models["hair_colour"].categories) == 4
    assert len(models["skin_colour"].categories) == 5
    # 41-SNP HIrisPlex-S panel on hair/skin; eye uses the same 41-col input (IrisPlex 6 are the drivers)
    assert len(models["skin_colour"].snps) == 41
    assert not models["skin_colour"].coefficients_pending


def test_reproduces_webtool_all_zero():
    """Loaded models must reproduce the deployed webtool's all-zero-genotype probabilities."""
    models = load_hirisplex_models()
    for trait, expected in WEBTOOL_ZERO.items():
        got = predict(models[trait], _zero_genotype(models[trait])).probabilities
        for cat, exp in expected.items():
            assert abs(got[cat] - exp) < 2e-3, (trait, cat, got[cat], exp)


def test_biology_direction_anchors():
    """Known pigment directions hold. rs12913832 (HERC2) dark allele shifts toward darker eye/hair."""
    models = load_hirisplex_models()
    eye = models["eye_colour"]
    # HERC2 counted allele homozygous (dark) -> brown-dominant; absent (blue allele) -> blue-dominant
    herc2 = next(s for s in eye.snps if s.rsid == "rs12913832")
    dark = _zero_genotype(eye); dark["rs12913832"] = herc2.counted_allele * 2
    assert predict(eye, dark).call == "brown"
    assert predict(eye, _zero_genotype(eye)).call == "blue"


def test_provenance_records_validation():
    p = provenance()
    assert "erasmusmc" in p["method"].lower() and "webtool" in p["source_tool"].lower()
    v = p["held_out_max_prob_error"]
    assert v["eye_colour"] < 1e-6 and v["hair_colour"] < 1e-6      # machine-precision
    assert v["skin_colour"] < 2e-2                                  # <2% (separation betas)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
