"""Offline pin for the protein variant-effect cell: HGVS parsing + BLOSUM wiring + run() on a synthetic DMS."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.proteingym_variant_effect import parse_missense, run  # noqa: E402


def test_parse_missense():
    assert parse_missense("p.Tyr88Val") == ("Y", "V")
    assert parse_missense("p.Arg130Ter") is None      # nonsense
    assert parse_missense("p.Ala5Ala") is None         # synonymous
    assert parse_missense("p.Met1del") is None         # indel/unparseable
    assert parse_missense("garbage") is None


def test_run_on_synthetic(tmp_path):
    # plant a positive BLOSUM-vs-abundance relationship: conservative subs high score, disruptive low; nonsense low
    rows = [
        ("p.Ile10Val", 0.95),   # I->V conservative (BLOSUM +3) -> high abundance
        ("p.Leu20Met", 0.90),   # conservative
        ("p.Ala30Val", 0.80),
        ("p.Gly40Trp", 0.15),   # G->W very disruptive (BLOSUM -2) -> low
        ("p.Arg50Asp", 0.20),   # R->D disruptive
        ("p.Pro60Trp", 0.10),
        ("p.Ser70Ter", 0.05),   # nonsense -> low (counted separately)
        ("p.Cys80Ter", 0.02),
    ]
    csv = tmp_path / "dms.csv"
    pd.DataFrame(rows, columns=["hgvs_pro", "score"]).to_csv(csv, index=False)
    res = run(csv)
    assert res["n_missense_scored"] == 6
    assert res["spearman_blosum_vs_abundance"] > 0.5    # planted positive relationship recovered
    assert res["direction_sanity_ok"]                    # nonsense mean < missense mean
    assert res["n_nonsense"] == 2
