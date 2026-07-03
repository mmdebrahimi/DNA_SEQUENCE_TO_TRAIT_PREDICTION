"""Offline pin for the C. elegans ben-1 benzimidazole decoder on a synthetic isotype-variant table."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.celegans_ben1_decoder import ben1_determinant, run, score  # noqa: E402


def _fixture(path: Path):
    # 8 strains: 3 ben-1-variant (all resistant), 1 variant-but-sensitive (FP), 2 resistant-no-variant (FN),
    # 2 clean sensitive (TN). Continuous pheno higher = more resistant.
    rows = [
        ("A", "missense_F200Y", "NA", False, "Missense", "TRUE", 120.0),
        ("B", "NA", "Deletion", False, "Deletion", "TRUE", 90.0),
        ("C", "NA", "NA", True, "Low ben-1 expression", "TRUE", 80.0),
        ("D", "frameshift_31", "NA", False, "Frameshift", "FALSE", -60.0),   # determinant+ but S -> FP
        ("E", "NA", "NA", False, "No variant", "TRUE", 5.0),                  # R no determinant -> FN (marginal)
        ("F", "NA", "NA", False, "No variant", "TRUE", 8.0),                  # R no determinant -> FN
        ("G", "NA", "NA", False, "No variant", "FALSE", -100.0),              # TN
        ("H", "NA", "NA", False, "No variant", "FALSE", -90.0),              # TN
    ]
    df = pd.DataFrame(rows, columns=["strain", "ben-1_high_impact_SNV", "ben-1_SVs", "low_ben1_exp",
                                     "ben-1_clean_call", "abz_hta_norm_res", "abz_hta_norm_pheno"])
    df.to_csv(path, sep="\t", index=False)


def test_determinant_scan(tmp_path):
    p = tmp_path / "v.tsv"; _fixture(p)
    v = pd.read_csv(p, sep="\t")
    det = ben1_determinant(v)
    assert list(det) == [True, True, True, True, False, False, False, False]


def test_confusion_and_strong_vs_marginal(tmp_path):
    p = tmp_path / "v.tsv"; _fixture(p)
    v = pd.read_csv(p, sep="\t")
    det = ben1_determinant(v)
    r = score(v, det, "abz_hta_norm_res")
    assert (r["TP"], r["FP"], r["FN"], r["TN"]) == (3, 1, 2, 2)
    assert r["specificity"] == 0.667 and r["ppv"] == 0.75
    # determinant+ resistant strains are MUCH more strongly resistant than determinant- resistant (marginal)
    sm = r["strong_vs_marginal"]
    assert sm["mean_pheno_resistant_determinant_pos"] > sm["mean_pheno_resistant_determinant_neg"]


def test_run_end_to_end(tmp_path):
    p = tmp_path / "v.tsv"; _fixture(p)
    res = run(p)
    assert res["n_strains"] == 8 and res["n_ben1_determinant_pos"] == 4
    a = res["assays"][0]
    assert a["phenotype"] == "abz_hta_norm_res" and a["TP"] == 3
