"""Tests for multi-axis co-resistance (`scripts/coresistance_multiaxis.py`, Family C multi-axis).
Pure-function tests offline (numpy/sklearn core -> CI-safe); the real-data test SKIPS when the plasmid
cache (needs the local finder sweep) is absent.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import coresistance_multiaxis as ma  # noqa: E402


def test_build_joint_namespaces_replicons_and_filters_organism():
    acc_org = {f"g{i}": "klebsiella" for i in range(70)}
    acc_org["h0"] = "acinetobacter"                       # below MIN_GENOMES for its org -> dropped
    acc_dets = {a: {"sul1"} for a in acc_org}
    plasmid = {a: {"organism": acc_org[a], "replicons": ["IncFIB"]} for a in acc_org}
    per_org = ma.build_joint(acc_org, acc_dets, plasmid)
    assert "klebsiella" in per_org and "acinetobacter" not in per_org   # only >= MIN_GENOMES
    accs, X, names, nidx = per_org["klebsiella"]
    assert "sul1" in names and "plasmid:IncFIB" in names   # replicon namespaced with the prefix
    assert X.shape == (70, 2)


def test_rep_prefix_constant():
    assert ma.REP_PREFIX == "plasmid:" and ma.PASS_FRACTION == 0.5


def test_module_does_not_modify_frozen_surface():
    src = (REPO / "scripts" / "coresistance_multiaxis.py").read_text(encoding="utf-8")
    assert "calibrated_amr_rules" not in src and "shipped_decoder_surface" not in src


_CACHE = REPO / "data" / "processed" / "plasmid_axis_cache.json"


@pytest.mark.skipif(not _CACHE.exists(), reason="plasmid axis cache absent (needs the local finder sweep)")
def test_run_real_multiaxis():
    res = ma.run(cache=_CACHE, boot=100)
    assert res["verdict"] in ("PASS_MULTIAXIS_LINKAGE", "FAIL_AXES_INDEPENDENT", "NO_TESTABLE")
    assert res["n_plasmid_genomes"] >= 1
    # at least one organism should surface a determinant->plasmid lift entry
    assert any(r.get("determinant_to_plasmid_lift") for r in res["per_organism"].values())


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
