"""Offline pin for the conservation (Site-Independent vs learned) benchmark on synthetic ProteinGym tables."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.dms_conservation_benchmark import load, run  # noqa: E402


def _tables(tmp: Path):
    # reference: DMS_id -> UniProt_ID, coarse_selection_type, taxon
    ref = pd.DataFrame({
        "DMS_id": [f"A{i}" for i in range(8)],
        "UniProt_ID": [f"G{i}_HUMAN" for i in range(8)],
        "coarse_selection_type": ["Activity", "Activity", "Activity",
                                   "Binding", "Binding", "Binding",
                                   "Stability", "Stability"],
        "taxon": ["Human"] * 8,
    })
    # spearman: Site-Independent competitive on Binding, trails on Activity/Stability
    sp = pd.DataFrame({
        "DMS ID": [f"A{i}" for i in range(8)],
        "Site-Independent": [0.42, 0.44, 0.40, 0.39, 0.41, 0.38, 0.30, 0.32],
        "ESM-1v (ensemble)": [0.52, 0.53, 0.51, 0.34, 0.35, 0.33, 0.52, 0.53],
        "EVE (ensemble)": [0.51, 0.52, 0.50, 0.40, 0.41, 0.39, 0.50, 0.51],
        "GEMME": [0.54, 0.55, 0.53, 0.35, 0.36, 0.34, 0.53, 0.54],
    })
    rp = tmp / "ref.csv"; ref.to_csv(rp, index=False)
    spp = tmp / "sp.csv"; sp.to_csv(spp, index=False)
    return rp, spp


def test_load_join(tmp_path):
    rp, spp = _tables(tmp_path)
    m = load(rp, spp)
    assert len(m) == 8 and "selection" in m.columns and m["Site-Independent"].notna().all()


def test_run_by_selection_and_verdict(tmp_path):
    rp, spp = _tables(tmp_path)
    res = run(rp, spp, blosum_json=None, am_json=None)
    bs = {r["selection"]: r for r in res["by_selection"]}
    # Binding: Site-Independent (median 0.39) beats ESM-1v (0.34) and GEMME (0.35) -> beats_a_learned_model
    assert bs["Binding"]["beats_a_learned_model"] is True
    # Activity: Site-Independent (0.42) trails all learned -> does not beat
    assert bs["Activity"]["beats_a_learned_model"] is False
    # Activity site-independent ~0.42 (>=0.35) + a modality beats a learned model -> LARGELY_COMPETES
    assert res["verdict"] == "DETERMINISTIC_CONSERVATION_LARGELY_COMPETES"
    assert res["function_site_independent"] == 0.42 and res["function_best_learned"] == 0.54
