"""Offline tests for the forward capability preflight (dna_decode/forward/capabilities.py).

Pure: no torch import, no Docker call, no model download. Capabilities are probed via find_spec/which and
the runnable-method logic is tested against injected `caps` dicts so the result is host-independent.
"""
from __future__ import annotations

from dna_decode.forward.capabilities import (
    probe_capabilities,
    runnable_methods,
    strongest_runnable,
    render_capabilities,
    METHOD_STRENGTH,
)

ALL = {"torch": True, "transformers": True, "torch_geometric": True, "prosst_repo": True, "docker": True}
NONE = {"torch": False, "transformers": False, "torch_geometric": False, "prosst_repo": False, "docker": False}
ESM_ONLY = {"torch": True, "transformers": True, "torch_geometric": False, "prosst_repo": False, "docker": False}


def test_probe_returns_the_expected_keys():
    caps = probe_capabilities()
    assert set(caps) == {"torch", "transformers", "torch_geometric", "prosst_repo", "docker"}
    assert all(isinstance(v, bool) for v in caps.values())


def test_blosum_always_runnable():
    for caps in (ALL, NONE, ESM_ONLY):
        assert runnable_methods(caps)["blosum62"][0] is True


def test_esm2_needs_torch_and_transformers():
    assert runnable_methods(ESM_ONLY)["esm2"][0] is True
    assert runnable_methods({**NONE, "torch": True})["esm2"][0] is False   # transformers missing
    assert runnable_methods(NONE)["esm2"][0] is False


def test_prosst_needs_esm_deps_plus_geometric_and_repo():
    assert runnable_methods(ALL)["prosst"][0] is True
    assert runnable_methods(ESM_ONLY)["prosst"][0] is False                       # no torch_geometric/repo
    assert runnable_methods({**ALL, "prosst_repo": False})["prosst"][0] is False   # repo absent
    assert runnable_methods({**ALL, "torch_geometric": False})["prosst"][0] is False  # geometric absent


def test_gemme_needs_docker():
    assert runnable_methods({**NONE, "docker": True})["gemme"][0] is True
    assert runnable_methods(NONE)["gemme"][0] is False


def test_hybrid_needs_two_of_three():
    assert runnable_methods(ALL)["hybrid"][0] is True                       # esm2+prosst+gemme
    assert runnable_methods(ESM_ONLY)["hybrid"][0] is False                 # only esm2
    assert runnable_methods({**ESM_ONLY, "docker": True})["hybrid"][0] is True  # esm2+gemme


def test_strongest_runnable_orders_by_strength():
    assert strongest_runnable(ALL) == "hybrid"
    assert strongest_runnable(ESM_ONLY) == "esm2"
    assert strongest_runnable(NONE) == "blosum62"


def test_strength_ranking_is_ordered_and_complete():
    assert METHOD_STRENGTH[0] == "hybrid" and METHOD_STRENGTH[-1] == "blosum62"
    for m in ("hybrid", "prosst", "esm2", "gemme", "alphamissense", "blosum62"):
        assert m in METHOD_STRENGTH


def test_render_is_ascii_and_names_the_preflight():
    out = render_capabilities(ALL)
    assert "capability preflight" in out
    assert "runnable methods" in out
    assert all(ord(c) < 128 for c in out)   # no console mojibake on Windows cp1252
    # an unrunnable method shows an install hint
    out2 = render_capabilities(NONE)
    assert "pip install" in out2
