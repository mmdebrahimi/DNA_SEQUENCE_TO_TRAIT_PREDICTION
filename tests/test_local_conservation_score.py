"""Offline pin for the frozen-spec local independent-sites conservation score on a synthetic MSA."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.local_conservation_score import score_assay, weighted_freqs  # noqa: E402


def _msa(tmp: Path):
    # 3 match columns. col0 conserved M, col1 conserved K, col2 variable.
    seqs = {
        "focus": "MKA",
        "h1": "MKC", "h2": "MKD", "h3": "MKE", "h4": "MKF",
        "h5": "MKa" + "G",   # lowercase insert 'a' between col1 and col2 -> stripped; then G at col2
    }
    p = tmp / "m.a2m"
    with open(p, "w", encoding="utf-8") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")
    return p, len(seqs)


def test_weighted_freqs_conserved(tmp_path):
    a2m, n = _msa(tmp_path)
    w = np.ones(n)
    f = weighted_freqs(a2m, w)
    assert f.shape[0] == 3
    # col0 dominated by M, col1 by K (argmax; freq < 0.5 here only because the synthetic MSA has 6 seqs +
    # the lambda=0.5 pseudocount dilutes -- real MSAs with 1000s of seqs give ~1.0)
    from scripts.local_conservation_score import AAS, AA_IDX
    assert AAS[int(f[0].argmax())] == "M" and AAS[int(f[1].argmax())] == "K"
    assert f[0, AA_IDX["M"]] > 0.35


def test_score_assay_direction(tmp_path):
    a2m, n = _msa(tmp_path)
    np.save(tmp_path / "w.npy", np.ones(n))
    # DMS: mutating the conserved M1 or K2 -> low fitness; variable pos3 tolerant -> high
    rows = [("M1C", 0.1), ("M1D", 0.15), ("M1E", 0.12), ("K2C", 0.1), ("K2D", 0.13),
            ("A3C", 0.8), ("A3D", 0.82), ("A3E", 0.85), ("A3F", 0.79), ("A3G", 0.81),
            ("M1F", 0.2), ("K2E", 0.18), ("A3W", 0.7), ("A3Y", 0.75), ("M1G", 0.14),
            ("K2F", 0.16), ("A3H", 0.6), ("A3I", 0.72), ("A3K", 0.68), ("M1H", 0.11)]
    dms = tmp_path / "dms.csv"
    pd.DataFrame(rows, columns=["mutant", "DMS_score"]).to_csv(dms, index=False)
    r = score_assay(dms, a2m, tmp_path / "w.npy", msa_start=1, msa_end=3)
    assert r is not None and r["n_scored"] >= 20 - 5
    # conserved-position mutations score LOW (log-odds negative), tolerant pos3 HIGH -> positive Spearman
    assert r["local_conservation_spearman"] > 0.3
    assert len(r["msa_sha"]) == 16 and r["ncol"] == 3
