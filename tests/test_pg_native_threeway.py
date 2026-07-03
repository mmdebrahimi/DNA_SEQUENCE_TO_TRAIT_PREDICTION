"""Offline pin for the ProteinGym-native three-way join on synthetic assay/score/AM data."""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.pg_native_threeway import _blosum, assay_threeway, load_am, run  # noqa: E402


def _assay(dms_dir: Path, si_dir: Path, did: str):
    dms_dir.mkdir(parents=True, exist_ok=True); si_dir.mkdir(parents=True, exist_ok=True)
    # conservative subs -> high DMS_score (fitter); disruptive -> low; nonsense-like excluded
    rows = [("I10V", 0.9), ("L20M", 0.9), ("A30S", 0.8), ("V40I", 0.85), ("T50S", 0.82),
            ("G60W", 0.1), ("R70D", 0.15), ("P80W", 0.1), ("C90R", 0.2), ("W100G", 0.12),
            ("M110L", 0.88), ("F120Y", 0.86), ("K130R", 0.8), ("D140E", 0.83), ("S150T", 0.81),
            ("H160Q", 0.4), ("N170K", 0.45), ("E180D", 0.7), ("Q190H", 0.5), ("Y200F", 0.84),
            ("A210G", 0.3), ("L220P", 0.1), ("V230A", 0.6), ("I240T", 0.35), ("R250C", 0.18),
            ("G260A", 0.55), ("T270A", 0.5), ("S280A", 0.52), ("K290A", 0.48), ("D300A", 0.4)]
    pd.DataFrame(rows, columns=["mutant", "DMS_score"]).to_csv(dms_dir / f"{did}.csv", index=False)
    # Site-Independent: fitness-oriented (higher for conservative) -- correlated with DMS
    si = [(m, s * 2 + 0.1) for m, s in rows]
    pd.DataFrame(si, columns=["mutant", "Site_Independent_score"]).to_csv(si_dir / f"{did}.csv", index=False)
    return [m for m, _ in rows]


def test_load_am(tmp_path):
    f = tmp_path / "am.tsv"
    f.write_text("P1\tG60W\t0.95\tpath\nP1\tI10V\t0.05\tbenign\n", encoding="utf-8")
    lut = load_am(f)
    assert lut[("P1", "G60W")] == 0.95


def test_assay_threeway_alignment(tmp_path):
    muts = _assay(tmp_path / "dms", tmp_path / "si", "G1_HUMAN_Study_2020")
    # AM: pathogenicity -- HIGH for disruptive (low DMS), LOW for conservative
    am = {("P1", m): (0.9 if m in ("G60W", "R70D", "P80W", "L220P", "R250C") else 0.1) for m in muts}
    r = assay_threeway(tmp_path / "dms" / "G1_HUMAN_Study_2020.csv",
                       tmp_path / "si" / "G1_HUMAN_Study_2020.csv", "P1", _blosum(), am)
    # all three aligned positive = predicts fitness
    assert r["blosum_full"] > 0 and r["site_independent_full"] > 0.5 and r["alphamissense_full"] > 0.3
    assert r["n_intersection"] >= 20


def test_run_end_to_end(tmp_path):
    dms, si = tmp_path / "dms", tmp_path / "si"
    _assay(dms, si, "G1_HUMAN_Study_2020")
    ref = pd.DataFrame({"DMS_id": ["G1_HUMAN_Study_2020"], "UniProt_ID": ["G1_HUMAN"],
                        "coarse_selection_type": ["Activity"], "taxon": ["Human"]})
    rp = tmp_path / "ref.csv"; ref.to_csv(rp, index=False)
    am = tmp_path / "am.tsv"; am.write_text("P1\tI10V\t0.05\tb\n", encoding="utf-8")
    accp = tmp_path / "acc.json"; accp.write_text(json.dumps({"G1_HUMAN": "P1"}), encoding="utf-8")
    res = run(rp, dms, si, accp, am)
    assert res["identical_rows"] is True and res["n_assays_threeway"] >= 0
