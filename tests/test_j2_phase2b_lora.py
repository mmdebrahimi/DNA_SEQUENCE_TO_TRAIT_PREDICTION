"""Offline tests for the J2 Phase 2b LoRA fine-tune scaffold (no torch, no peft, no GPU).

Covers the load-bearing pure logic: leakage-safe protein split, the ranking-pair sampler, protein-key
grouping, and variant extraction. The LoRA training/eval loop is GPU-run and not exercised here."""
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# phase2b imports the phase1 core; loading it exercises that wiring too
PB = _load(ROOT / "notebooks" / "j2_phase2b_lora_finetune.py", "j2b")


def test_protein_key_prefers_uniprot_then_prefix():
    assert PB.protein_key("A0A140D2T1_ZIKV_Sourisseau_2019") == "A0A140D2T1"
    assert PB.protein_key("A0A140D2T1_ZIKV_Sourisseau_2019", "P12345") == "P12345"
    assert PB.protein_key("X_Y_Z", "  ") == "X"          # blank uniprot -> fallback


def test_split_proteins_no_leakage_and_covers_all():
    # 3 proteins, 2 assays each (same protein must land in the same fold)
    assay_keys = {f"{prot}_assay{i}": prot for prot in ("P1", "P2", "P3", "P4", "P5") for i in (0, 1)}
    folds = PB.split_proteins(assay_keys, nfolds=3, seed=0)
    assert set(folds) == set(assay_keys)                 # every assay assigned
    # every protein's assays share ONE fold (no protein leaks across folds)
    by_prot = {}
    for a, f in folds.items():
        by_prot.setdefault(assay_keys[a], set()).add(f)
    assert all(len(fs) == 1 for fs in by_prot.values())


def test_split_proteins_deterministic_and_validates_nfolds():
    ak = {f"P{i}_a": f"P{i}" for i in range(10)}
    assert PB.split_proteins(ak, 5, 7) == PB.split_proteins(ak, 5, 7)     # same seed -> same split
    assert PB.split_proteins(ak, 5, 7) != PB.split_proteins(ak, 5, 8)     # different seed -> (very likely) differs
    with pytest.raises(ValueError):
        PB.split_proteins(ak, 1, 0)


def test_sample_pairs_bounds_and_determinism():
    assert PB.sample_pairs(1, 100, 0) == []
    assert PB.sample_pairs(4, 100, 0) == [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]  # all pairs (<=max)
    capped = PB.sample_pairs(50, 20, 0)
    assert len(capped) == 20 and len(set(capped)) == 20 and capped == sorted(capped)       # unique + sorted
    assert PB.sample_pairs(50, 20, 0) == PB.sample_pairs(50, 20, 0)                          # deterministic


def test_assay_variants_filters_mismatch_multi_and_oob():
    seq = "MKTAYIAK"                                     # positions 1..8
    dms = {"M1A": 1.0, "K2R": 2.0, "Q3D": 9.9,          # Q3 mismatches (seq[2]=='T') -> dropped
           "A99G": 8.8, "M1A:K2R": 7.7}                  # oob + multi-mutant -> dropped
    got = PB.assay_variants(seq, dms, maxlen=1022)
    assert sorted((w, p, m) for w, p, m, _y in got) == [("K", 2, "R"), ("M", 1, "A")]
    assert PB.assay_variants("M" * 2000, {"M1A": 1.0}, maxlen=1022) == []   # too long -> dropped
