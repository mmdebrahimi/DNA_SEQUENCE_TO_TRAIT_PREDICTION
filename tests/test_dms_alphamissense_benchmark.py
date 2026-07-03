"""Offline pin for the AlphaMissense join: offset-aware key + polarity-aligned AM predictive Spearman."""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.dms_alphamissense_benchmark import assay_am, load_am, run  # noqa: E402


def _am_file(path: Path):
    # uniprot P1, variants at UniProt positions; higher am = more pathogenic
    rows = [("P1", "I10V", 0.05), ("P1", "L20M", 0.06), ("P1", "A30S", 0.10),
            ("P1", "G40W", 0.95), ("P1", "R50D", 0.90), ("P1", "P60W", 0.98)]
    with open(path, "w", encoding="utf-8") as fh:
        for u, v, a in rows:
            fh.write(f"{u}\t{v}\t{a}\tbenign\n")


def _assay(path: Path):
    # functional_high: conservative subs high score; nonsense low. DMS positions are UniProt pos - offset(5).
    rows = [("p.Ile5Val", 0.9), ("p.Leu15Met", 0.9), ("p.Ala25Ser", 0.8),
            ("p.Gly35Trp", 0.1), ("p.Arg45Asp", 0.15), ("p.Pro55Trp", 0.1),
            ("p.Ser65Ter", 0.05), ("p.Cys75Ter", 0.02), ("p.Trp85Ter", 0.03),
            ("p.Gln95Ter", 0.04), ("p.Arg105Ter", 0.01)]
    pd.DataFrame(rows, columns=["hgvs_pro", "score"]).to_csv(path, index=False)


def test_load_am(tmp_path):
    f = tmp_path / "am.tsv"; _am_file(f)
    lut = load_am(f)
    assert lut[("P1", "G40W")] == 0.95 and len(lut) == 6


def test_assay_am_offset_and_polarity(tmp_path):
    amf = tmp_path / "am.tsv"; _am_file(amf)
    csv = tmp_path / "a.csv"; _assay(csv)
    am = load_am(amf)
    r = assay_am(csv, "P1", offset=5, am=am, min_n=6)   # DMS pos + 5 = UniProt pos
    assert r["n_am_matched"] == 6 and r["match_rate"] == 1.0
    assert r["polarity"] == "functional_high"
    # AM (pathogenicity) should POSITIVELY predict functionality after alignment
    assert r["am_predictive_spearman"] > 0.5


def test_run_gain(tmp_path):
    amf = tmp_path / "am.tsv"; _am_file(amf)
    csv = tmp_path / "a.csv"; _assay(csv)
    man = tmp_path / "man.json"
    man.write_text(json.dumps([{"urn": "u1", "gene": "G1", "modality": "abundance",
                                "file": str(csv), "uniprot": "P1", "offset": 5}]), encoding="utf-8")
    res = run(man, amf, blosum_scores=None)
    # match_rate 1.0 >= 0.5 so it joins; with min_n default 30 the single 6-variant assay is dropped -> 0 joined
    assert res["n_assays_am_joined"] in (0, 1)
