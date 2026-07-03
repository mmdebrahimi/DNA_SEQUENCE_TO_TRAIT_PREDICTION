"""Offline pin for the GDSC fusion decoder — synthetic GDSC-shaped parquets + fusion csv, no network/D: data."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.gdsc_fusion_decoder import (  # noqa: E402
    _direction_ok,
    fusion_ids,
    load_gdsc,
    run,
    score_case,
)

LINS = ["leukemia", "lung", "skin"]
N = 30


def _fixture(d: Path, rng):
    d.mkdir(parents=True, exist_ok=True)
    cos = list(range(1000, 1000 + 3 * N))
    dm = [f"ACH-{i:06d}" for i in range(3 * N)]
    lin = [LINS[i % 3] for i in range(3 * N)]
    pd.DataFrame({"COSMIC_ID": cos, "DepMap_ID": dm, "CCLE_Name": dm, "lineage": lin}).to_csv(
        d / "sample_info.csv", index=False)
    # plant BCR-ABL1 in 6 leukemia lines; they get LOWER LN_IC50 for dasatinib (sensitize)
    bcr = [i for i in range(3 * N) if lin[i] == "leukemia"][:6]
    base = np.repeat(rng.normal(0, 1.0, 3), N)
    ic = base + rng.normal(0, 0.2, 3 * N)
    for i in bcr:
        ic[i] -= 5.0
    g2 = pd.DataFrame({"COSMIC_ID": cos, "DRUG_NAME": ["Dasatinib"] * (3 * N),
                       "LN_IC50": ic, "Z_SCORE": ic})
    g2.to_parquet(d / "gdsc2_fitted.parquet")
    # GDSC1 empty-ish (imatinib on a couple lines) so load_gdsc concatenates both
    pd.DataFrame({"COSMIC_ID": cos[:5], "DRUG_NAME": ["Imatinib"] * 5,
                  "LN_IC50": [0.0] * 5, "Z_SCORE": [0.0] * 5}).to_parquet(d / "gdsc1_fitted.parquet")
    # fusions
    pd.DataFrame({"DepMap_ID": [dm[i] for i in bcr], "#FusionName": ["BCR--ABL1"] * len(bcr)}).to_csv(
        d / "ccle_fusions.csv", index=False)
    return set(dm[i] for i in bcr)


def test_direction_helper():
    assert _direction_ok(-2.0, "sensitize") and not _direction_ok(2.0, "sensitize")


def test_fusion_ids(tmp_path):
    rng = np.random.default_rng(0)
    bcr = _fixture(tmp_path, rng)
    assert fusion_ids(tmp_path / "ccle_fusions.csv", "BCR--ABL1") == bcr


def test_load_gdsc_bridges_and_concats(tmp_path):
    rng = np.random.default_rng(1)
    _fixture(tmp_path, rng)
    g = load_gdsc(tmp_path)
    assert set(g["DATASET"]) == {"GDSC1", "GDSC2"}
    assert g["DepMap_ID"].notna().all() and "leukemia" in set(g["lineage"])


def test_score_case_recovers_powered_bcrabl_signal(tmp_path):
    rng = np.random.default_rng(2)
    bcr = _fixture(tmp_path, rng)
    g = load_gdsc(tmp_path)
    rec = score_case(g, bcr, "dasatinib", "GDSC2", "sensitize")
    assert rec["n_fusion_positive"] == 6
    assert rec["global_ln_ic50_delta"] < -2.0            # BCR-ABL+ much more sensitive
    assert rec["within_lineage_t"] < -1.0                # de-confounded within leukemia
    assert rec["direction_ok"]


def test_run_end_to_end(tmp_path):
    rng = np.random.default_rng(3)
    _fixture(tmp_path, rng)
    res = run(tmp_path, tmp_path / "ccle_fusions.csv")
    das = [c for c in res["cases"] if c["drug"] == "dasatinib"][0]
    assert das["biomarker_fusion"] == "BCR--ABL1" and das["direction_ok"]
