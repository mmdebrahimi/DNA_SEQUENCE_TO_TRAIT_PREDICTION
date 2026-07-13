"""Kill-tests for the /innovate run toward the Genome world model (2026-07-13).

Anti-theater contract (falsification.py semantics): each test PASSES iff the candidate idea is
FALSIFIED. pytest-exit-0 (pass) => predicate TRUE => the idea is KILLED; pytest-exit-1 (fail) =>
predicate FALSE => the idea SURVIVED this executed falsification attempt. Every threshold is checked
against COMMITTED repo data — no idea is marked survived without an executed test that could have
killed it. Read-only over committed artifacts; frozen decoder surface untouched.
"""
from __future__ import annotations

import csv
import glob
import json
import statistics
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def _proteingym_median(col: str) -> float:
    with open(REPO / "wiki/refs/proteingym_v1.3_DMS_level_Spearman.csv", encoding="utf-8") as f:
        r = csv.reader(f); h = next(r); i = h.index(col)
        return statistics.median(float(row[i]) for row in r if row[i] not in ("", "NA"))


def test_g9_interventional_falsified():
    """g9-interventional DEAD iff learned models FAIL on clean INTERVENTIONAL (DMS) data.
    Precondition the idea rests on: on ProteinGym (edit->measure DMS), a learned scorer is well
    above chance. Falsified iff median ESM2-650M Spearman < 0.10 (chance). Committed: 0.484."""
    median = _proteingym_median("ESM2 (650M)")
    assert median < 0.10, (
        f"NOT falsified: ESM2-650M median Spearman={median:.3f} >> 0.10 chance on interventional DMS "
        "-> learned models DO work on edit->measure data (the idea's precondition holds)")


def test_g8_residual_detector_falsified():
    """g8-residual-detector DEAD iff there is NO mechanism signal beyond lineage.
    Precondition: the de-confounding (leave-one-clade-out) determinant AUC stays high after removing
    phylogeny. Falsified iff clade-grouped AUC <= 0.55 (signal was all lineage). Committed: 0.908."""
    d = json.loads((REPO / "wiki/crossaxis_lineage_deconfound_determinant_2026-07-12.json").read_text(encoding="utf-8"))
    clade_grouped = float(d["median_auc_clade_grouped"])
    assert clade_grouped <= 0.55, (
        f"NOT falsified: clade-grouped determinant AUC={clade_grouped:.3f} >> 0.55 -> residual signal "
        "beyond lineage EXISTS (the de-confounding machinery has something real to detect)")


def test_g5_catalog_extension_falsified():
    """g5-catalog-extension DEAD iff the curated catalog has NO blind spots to extend.
    Precondition: catalog-negative-but-resistant cases exist (something for a learned curator to find).
    Falsified iff the HIV NNRTI catalog-negative resistant count == 0. Committed: 53."""
    d = json.loads((REPO / "wiki/hiv_esm_vs_catalog_2026-07-09.json").read_text(encoding="utf-8"))
    blind_spot_resistant = int(d["n_resistant"])
    assert blind_spot_resistant == 0, (
        f"NOT falsified: {blind_spot_resistant} catalog-negative resistant isolates exist -> the catalog "
        "HAS blind spots a learned curator could extend (the idea's precondition holds)")


def test_g1_federation_falsified():
    """g1-federation DEAD iff there are too few independent SCORED cells to compose a federation.
    Precondition: >=3 provenance-disjoint independently-validated cells exist. Falsified iff < 3.
    Committed: 10 provenance-disjoint validation artifacts."""
    n_cells = len(glob.glob(str(REPO / "wiki/provenance_disjoint_validation_*.json")))
    assert n_cells < 3, (
        f"NOT falsified: {n_cells} provenance-disjoint SCORED cells exist -> a federation of "
        "independently-validated deterministic cells is constructible (the idea's precondition holds)")


def test_monolith_scaling_falsified():
    """CONTROL (a tempting-but-wrong idea): 'just scale the whole-genome embedding to 3B/15B to finally
    learn the mechanism'. DEAD iff bigger ESM does NOT beat ESM2-650M on the committed benchmark.
    Falsified iff median ESM2-15B <= median ESM2-650M. Committed: 0.438 <= 0.484 -> TRUE -> this test
    PASSES -> the scaling idea is correctly KILLED (demonstrates the pipeline discriminates)."""
    m650 = _proteingym_median("ESM2 (650M)")
    m15b = _proteingym_median("ESM2 (15B)")
    assert m15b <= m650, (
        f"NOT falsified: ESM2-15B median={m15b:.3f} > ESM2-650M median={m650:.3f} -> bigger WOULD help "
        "(this would REVIVE the scaling idea) — but committed data shows the opposite")
