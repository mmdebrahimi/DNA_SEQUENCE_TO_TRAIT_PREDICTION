"""Pin the OpenSNP-archive zip ingester on a SYNTHETIC OpenSNP-shaped zip (no 21 GB, no network).

Proves the zip-native ingest (phenotype-CSV parse + per-user genotype-member rs12913832 extraction +
blue/brown scoring) is correct BEFORE the real dump lands. Covers: ';'-delimited phenotype CSV, the
user<id>_ genotype-filename convention, 23andMe + AncestryDNA member shapes, and the strand-agnostic call.
"""
import io
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.eye_colour_opensnp_ingest import run  # noqa: E402


def _make_zip(tmp_path) -> Path:
    z = tmp_path / "osnp.zip"
    with zipfile.ZipFile(z, "w") as zf:
        # OpenSNP-style phenotype CSV: ';'-delimited, one row per user, an 'Eye color' column
        zf.writestr("phenotypes_201712080449.csv",
                    "user_id;date_of_birth;Eye color;Hair color\n"
                    "1;1980;blue;blonde\n"
                    "2;1975;brown;black\n"
                    "3;1990;Dark Brown;brown\n"
                    "4;1985;green;brown\n"            # 'other' -> excluded
                    "5;1970;blue;grey\n")             # user 5 has no genotype file -> n_no_genotype_file
        # 23andMe genotype member: user1 = GG (blue, correct)
        zf.writestr("user1_file11_yearofbirth_1980_sex_XY.23andme.txt",
                    "# 23andMe\nrsid\tchromosome\tposition\tgenotype\n"
                    "rs1\t1\t1\tAA\nrs12913832\t15\t28365618\tGG\n")
        # AncestryDNA member: user2 = T T (brown, correct) on the complement strand
        zf.writestr("user2_file22_yearofbirth_1975_sex_XX.ancestry.txt",
                    "#AncestryDNA\nrsid\tchromosome\tposition\tallele1\tallele2\n"
                    "rs12913832\t15\t28365618\tT\tT\n")
        # user3 dark-brown, genotype AG -> intermediate (a known single-locus limitation; not a binary hit)
        zf.writestr("user3_file33_yearofbirth_1990_sex_XY.23andme.txt",
                    "rsid\tchromosome\tposition\tgenotype\nrs12913832\t15\t28365618\tAG\n")
    return z


def test_ingest_scores_synthetic_zip(tmp_path):
    res = run(_make_zip(tmp_path))
    assert res["status"] == "SCORED"
    assert res["n_users_eye_labelled"] == 5          # users 1,2,3,4,5 have eye colour (4='green' counted then excluded)
    assert res["n_other_excluded"] == 1              # user4 green
    assert res["n_no_genotype_file"] == 1            # user5 no .txt
    # binary: user1 blue->blue (TN), user2 brown->brown (TP); user3 dark-brown->intermediate (not binary)
    c = res["confusion_brown_positive"]
    assert c["TP"] == 1 and c["TN"] == 1 and c["FP"] == 0 and c["FN"] == 0
    assert res["accuracy"] == 1.0
    # user3 lands in the brown/intermediate strata cell (single-locus limitation, surfaced not hidden)
    assert res["strata_pred_by_label"]["brown"]["intermediate"] == 1


def test_zip_absent_is_honest(tmp_path):
    res = run(tmp_path / "nope.zip")
    assert res["status"] == "ZIP_NOT_PRESENT"


def test_pick_eye_column_prefers_exact_over_decoy():
    """Real OpenSNP-2017 regression: 'Eye pigmentation' precedes 'Eye color'; a naive first-'eye'
    match grabs the decoy. Pin that we select the exact eye-colour column (caught by R3 inspect)."""
    from scripts.eye_colour_opensnp_ingest import _pick_eye_column
    hdr = ["user_id", "date_of_birth", "chrom_sex", "Retrognathia (Marfan Syndrome)",
           "eye pigmentation ", "vegetarianism", "form of foot", "eye color",
           "hair and eye color brown", "mother's eye color"]
    assert _pick_eye_column(hdr) == 7              # 'eye color', NOT index 4 'eye pigmentation'
    assert _pick_eye_column(["eye colour"]) == 0   # British spelling
    assert _pick_eye_column(["no relevant column"]) is None
