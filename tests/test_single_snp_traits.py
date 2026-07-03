"""Offline pins for the single-SNP trait decoders + the generalized openSNP scorer.

Rule tests need no data. The end-to-end test builds a SYNTHETIC openSNP-shaped zip (phenotype CSV + per-user
genotype members) in tmp — no D: dump, no network."""
import io
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.single_snp_traits import (  # noqa: E402
    INDETERMINATE, TRAITS, bin_cilantro, bin_earwax, bin_lactase,
    call_cilantro, call_earwax, call_lactase,
)
from scripts.single_snp_opensnp_validate import run  # noqa: E402


def test_earwax_rule_strand_agnostic():
    assert call_earwax("TT") == "dry" and call_earwax("AA") == "dry"   # homozygous derived, either strand
    assert call_earwax("CT") == "wet" and call_earwax("CC") == "wet"   # ancestral C -> dominant wet
    assert call_earwax("AG") == "wet" and call_earwax("GG") == "wet"   # ancestral G -> dominant wet
    assert call_earwax("") == INDETERMINATE and call_earwax("--") == INDETERMINATE


def test_lactase_rule_strand_agnostic():
    assert call_lactase("AA") == "tolerant" and call_lactase("AG") == "tolerant"   # persistence allele
    assert call_lactase("CT") == "tolerant" and call_lactase("TT") == "tolerant"   # reverse strand T
    assert call_lactase("GG") == "intolerant" and call_lactase("CC") == "intolerant"
    assert call_lactase("") == INDETERMINATE


def test_cilantro_rule():
    assert call_cilantro("CC") == "soapy" and call_cilantro("CG") == "soapy"
    assert call_cilantro("AA") == "not-soapy" and call_cilantro("TT") == "not-soapy"
    assert call_cilantro("") == INDETERMINATE


def test_binners_polarity():
    assert bin_earwax("Dry / flaky") == "dry" and bin_earwax("Wet, sticky") == "wet"
    assert bin_earwax("rather not say") is None
    # lactase label POLARITY: the column is 'intolerance'
    assert bin_lactase("Lactose intolerant") == "intolerant"
    assert bin_lactase("not intolerant") == "tolerant"      # negation handled before 'intoler'
    assert bin_lactase("no") == "tolerant" and bin_lactase("tolerant") == "tolerant"
    assert bin_cilantro("tastes like soap") == "soapy" and bin_cilantro("no, tastes fine") == "not-soapy"


def test_registry_shape():
    assert set(TRAITS) == {"earwax", "lactase", "cilantro"}
    assert TRAITS["earwax"].rsid == "rs17822931" and TRAITS["earwax"].tier == "STRONG_MENDELIAN"
    assert TRAITS["cilantro"].tier == "WEAK_ASSOCIATION_CONTRAST"


def _synth_zip(tmp: Path) -> Path:
    """Build a minimal openSNP-shaped zip: phenotype CSV (';'-sep) + user genotype members."""
    # 6 users: 3 dry-earwax (TT) self-report dry, 3 wet (CC) self-report wet -> perfect separation
    pheno = "user_id;Earwax type\n"
    members = {}
    for uid, gt, rep in [("1", "TT", "dry"), ("2", "TT", "dry"), ("3", "AA", "dry"),
                         ("4", "CC", "wet"), ("5", "CT", "wet"), ("6", "GG", "wet")]:
        pheno += f"{uid};{rep}\n"
        members[f"user{uid}_file.23andme.txt"] = f"rsid\tchromosome\tposition\tgenotype\nrs17822931\t16\t48258198\t{gt}\n"
    p = tmp / "opensnp_test.zip"
    with zipfile.ZipFile(p, "w") as zf:
        zf.writestr("phenotypes_test.csv", pheno)
        for name, body in members.items():
            zf.writestr(name, body)
    return p


def test_end_to_end_synthetic_zip(tmp_path):
    z = _synth_zip(tmp_path)
    res = run(z, "earwax")
    assert res["status"] == "SCORED" and res["n_scored"] == 6
    # perfect synthetic separation: dry called dry, wet called wet
    conf = res["confusion_positive_dry"]
    assert conf["TP"] == 3 and conf["TN"] == 3 and conf["FP"] == 0 and conf["FN"] == 0
    assert res["accuracy"] == 1.0
