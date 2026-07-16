"""Offline tests for the IrisPlex eye-colour pigmentation cell. Pure math; no D:/no GPU/no network.

Pins: the reference-integrity biology contract (HERC2 GG->blue, AA->brown — the guard against fabricated
coefficients), probability normalization, missing-SNP fail-loud, genotype parsing, CLI dispatch.
"""
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pigment import (  # noqa: E402
    IRISPLEX_SNPS, MissingGenotypeError, predict_eye_color, reference_integrity_ok,
)
from dna_decode.pigment.cli import main as pigment_main  # noqa: E402
from dna_decode import cli as unified  # noqa: E402


def _all_snps(herc2_gt):
    """Genotype dict with a chosen HERC2 genotype + 0 counted alleles at the other 5."""
    out = {}
    for rsid, allele, _, _ in IRISPLEX_SNPS:
        if rsid == "rs12913832":
            out[rsid] = herc2_gt
        else:
            other = next(b for b in "ACGT" if b != allele)
            out[rsid] = other + other
    return out


def test_reference_integrity_biology_contract():
    # THE load-bearing guard: a fabricated/corrupted coefficient set fails this.
    assert reference_integrity_ok() is True


def test_herc2_gg_is_blue():
    r = predict_eye_color(_all_snps("GG"))
    assert r.call == "blue" and r.p_blue > 0.7 and r.confidence == "high"


def test_herc2_aa_is_brown():
    r = predict_eye_color(_all_snps("AA"))
    assert r.call == "brown" and r.p_brown > 0.7


def test_probs_sum_to_one():
    r = predict_eye_color(_all_snps("AG"))
    assert abs(r.p_blue + r.p_intermediate + r.p_brown - 1.0) < 1e-6


def test_missing_herc2_raises():
    g = _all_snps("GG")
    del g["rs12913832"]
    with pytest.raises(MissingGenotypeError):
        predict_eye_color(g, allow_missing=True)   # HERC2 required even under allow_missing


def test_missing_other_snp_requires_flag():
    g = _all_snps("GG")
    del g["rs1800407"]
    with pytest.raises(MissingGenotypeError):
        predict_eye_color(g)                       # default: fail loud
    r = predict_eye_color(g, allow_missing=True)   # imputed -> low confidence
    assert r.confidence == "low" and any("imputed" in n for n in r.notes)


def test_genotype_parsing_variants():
    for gt in ("A/G", "A|G", "AG"):
        g = _all_snps(gt)
        r = predict_eye_color(g)
        assert r.counted_alleles["rs12913832"] == 1


def test_cli_human_and_json(capsys):
    spec = "rs12913832=GG,rs1800407=TT,rs12896399=TT,rs16891982=GG,rs1393350=GG,rs12203592=GG"
    rc = pigment_main(["--genotypes", spec])
    out = capsys.readouterr().out
    assert rc == 0 and "eye colour" in out and "call:" in out.lower()
    rc = pigment_main(["--genotypes", spec, "--json"])
    d = json.loads(capsys.readouterr().out)
    assert rc == 0 and d["model"] == "IrisPlex" and d["call"] in ("blue", "intermediate", "brown")


def test_cli_bad_genotype_exits_2(capsys):
    rc = pigment_main(["--genotypes", "rs12913832"])   # no '=' -> parse error
    assert rc == 2 and "error" in capsys.readouterr().err


def test_dispatch_through_unified_cli(capsys):
    spec = "rs12913832=AA,rs1800407=TT,rs12896399=TT,rs16891982=GG,rs1393350=GG,rs12203592=GG"
    rc = unified.main(["pigment", "--genotypes", spec, "--json"])
    assert rc == 0
    assert json.loads(capsys.readouterr().out)["call"] == "brown"


def test_pigment_in_registry_and_list(capsys):
    assert "pigment" in unified.TRAITS
    rc = unified.main(["list"])
    out = capsys.readouterr().out
    assert rc == 0 and "pigment" in out and "eye colour" in out.lower()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
