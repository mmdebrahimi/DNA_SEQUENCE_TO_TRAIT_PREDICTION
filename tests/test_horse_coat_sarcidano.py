"""Regression: score the horse coat-colour rule on the committed Sarcidano-2022 real dataset (non-circular).

Data: data/horse/sarcidano_2022.tsv, reconstructed from Sarcidano 2022 (PMC9558981) Table 3 — functional
MC1R/ASIP genotype x VISUALLY-OBSERVED colour, N=70 (8 phenotypically-grey e/e horses excluded from the
base-3 test). Pins the real-data concordance so a rule change can't silently regress it. No network.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.horse_coat_validate import run  # noqa: E402

_TSV = Path(__file__).resolve().parent.parent / "data" / "horse" / "sarcidano_2022.tsv"


def test_sarcidano_real_data_concordance():
    if not _TSV.exists():
        import pytest
        pytest.skip("Sarcidano TSV absent")
    res = run(_TSV)
    assert res["status"] == "SCORED"
    assert res["n_scored"] == 62          # 4 bay + 19 black + 39 chestnut
    assert res["n_correct"] == 62
    assert res["concordance"] == 1.0      # deployed textbook rule on an isolated breed -> perfect base-3
    assert res["n_excluded_nonbase"] == 8 # phenotypically-grey e/e horses (separate masking locus)
    conf = res["confusion_observed_to_predicted"]
    assert conf["bay"] == {"bay": 4} and conf["black"] == {"black": 19} and conf["chestnut"] == {"chestnut": 39}
