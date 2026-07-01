"""Pin the IrisPlex v0.1 OpenSNP validation HARNESS on a synthetic zip (no D:, no network).

The real-data run is gated on the external D: drive; this proves the harness logic — especially the
heterozygote RESCUE accounting (a user v0 abstains on that v0.1 resolves) — is correct independently.
"""
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.eye_colour_irisplex_validate import run  # noqa: E402

_HDR = ("rsid\tchromosome\tposition\tgenotype\n")


def _member(u: int, snps: dict[str, str]) -> tuple[str, str]:
    body = _HDR + "".join(f"{r}\t15\t1\t{g}\n" for r, g in snps.items())
    return (f"user{u}_file{u}_yearofbirth_1990_sex_XY.23andme.txt", body)


def _make_zip(tmp_path) -> Path:
    z = tmp_path / "osnp.zip"
    blue = {"rs12913832": "GG", "rs1800407": "CC", "rs12896399": "GG",
            "rs16891982": "GG", "rs1393350": "GG", "rs12203592": "CC"}
    brown = {"rs12913832": "AA", "rs1800407": "TT", "rs12896399": "TT",
             "rs16891982": "CC", "rs1393350": "AA", "rs12203592": "TT"}
    het = {**blue, "rs12913832": "AG"}     # v0 abstains; v0.1 -> brown (rescue)
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("phenotypes_x.csv",
                    "user_id;Eye color\n1;blue\n2;brown\n3;brown\n4;green\n")
        for u, snps in ((1, blue), (2, brown), (3, het)):
            name, body = _member(u, snps)
            zf.writestr(name, body)
        # user4 = green -> excluded; no genotype needed
    return z


def test_irisplex_harness_scores_and_rescues(tmp_path):
    res = run(_make_zip(tmp_path))
    assert res["status"] == "SCORED"
    assert res["n_complete_case_scored"] == 3          # users 1,2,3 (4 is green->other)
    assert res["n_other_excluded"] == 1
    # full-set v0.1 binary (argmax): blue->blue, brown->brown, het(brown)->brown = 3/3
    fb = res["v01_full_set_binary_argmax"]
    assert fb["n"] == 3 and fb["accuracy"] == 1.0
    # deployed 0.7-threshold mode present + never mis-calls (abstains below 0.7)
    thr = res["v01_deployed_threshold_0.7"]
    assert thr["FP"] == 0 and thr["FN"] == 0
    assert "n_undefined_or_intermediate_abstained" in thr
    # the rescue: user3 (AG) is abstained by v0, resolved correctly by v0.1
    r = res["heterozygote_rescue"]
    assert r["n_v0_abstained"] == 1
    assert r["n_v01_correct"] == 1
    assert r["v01_recall_on_rescue"] == 1.0
    # paired coverage: v0 binary-calls only homozygotes (2), v0.1 calls all 3
    assert res["paired_complete_case"]["v0_single_snp"]["n"] == 2
    assert res["paired_complete_case"]["v01_irisplex"]["n"] == 3


def test_irisplex_zip_absent_is_honest(tmp_path):
    assert run(tmp_path / "nope.zip")["status"] == "ZIP_NOT_PRESENT"
