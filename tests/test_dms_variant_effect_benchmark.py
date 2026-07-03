"""Offline pin for the multi-gene DMS benchmark: polarity correction + aggregation on synthetic assays."""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.dms_variant_effect_benchmark import assay_result, run  # noqa: E402


def _assay(path: Path, damaging_high: bool):
    # conservative subs (high BLOSUM) preserve function; nonsense = maximally damaging.
    # functional_high: functional variants have HIGH score. damaging_high: damaging variants have HIGH score.
    rows = [("p.Ile10Val", 0.9), ("p.Leu20Met", 0.9), ("p.Ala30Ser", 0.8),
            ("p.Gly40Trp", 0.1), ("p.Arg50Asp", 0.15), ("p.Pro60Trp", 0.1)]
    ter = [("p.Ser70Ter", 0.05), ("p.Cys80Ter", 0.02), ("p.Trp90Ter", 0.03),
           ("p.Gln100Ter", 0.04), ("p.Arg110Ter", 0.01)]
    data = rows + ter
    if damaging_high:                      # flip: damaging variants get HIGH score
        data = [(h, 1.0 - s) for h, s in data]
    pd.DataFrame(data, columns=["hgvs_pro", "score"]).to_csv(path, index=False)


def test_polarity_correction(tmp_path):
    fh = tmp_path / "func_high.csv"; _assay(fh, damaging_high=False)
    dh = tmp_path / "dam_high.csv"; _assay(dh, damaging_high=True)
    rf = assay_result(fh, "G1", "abundance", min_n=6)
    rd = assay_result(dh, "G2", "function", min_n=6)
    assert rf["polarity"] == "functional_high" and rd["polarity"] == "damaging_high"
    # BOTH should give POSITIVE polarity-corrected Spearman (conservative -> preserved function)
    assert rf["spearman_polarity_corrected"] > 0.4
    assert rd["spearman_polarity_corrected"] > 0.4
    # raw signs are opposite
    assert rf["spearman_blosum62_raw"] > 0 and rd["spearman_blosum62_raw"] < 0


def test_run_aggregates(tmp_path):
    fh = tmp_path / "a.csv"; _assay(fh, damaging_high=False)
    dh = tmp_path / "b.csv"; _assay(dh, damaging_high=True)
    man = tmp_path / "man.json"
    man.write_text(json.dumps([
        {"urn": "u1", "gene": "G1", "modality": "abundance", "file": str(fh)},
        {"urn": "u2", "gene": "G2", "modality": "function", "file": str(dh)},
    ]), encoding="utf-8")
    res = run(man, min_n=6)
    assert res["n_assays"] == 2 and res["overall_median_spearman"] > 0.4
    assert "abundance" in res["by_modality"] and "function" in res["by_modality"]
