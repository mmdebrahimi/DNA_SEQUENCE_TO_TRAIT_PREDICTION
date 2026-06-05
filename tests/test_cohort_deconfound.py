"""Tests for the cohort de-confound gate (dna_decode/eval/cohort_deconfound).

Pure synthetic cases pinning the verdict thresholds, incl. regression pins reproducing the two
real data points: cef gate_b (1 shared lineage -> CONFOUNDED) and cipro N=147 (6 shared, ~13%
minority -> DE_CONFOUNDED). Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.eval.cohort_deconfound import (
    CONFOUNDED, DE_CONFOUNDED, WARN, confound_report,
)


def test_confounded_disjoint_lineages():
    # R in lineages A..D, S in E..H -> 0 shared
    labels = [1, 1, 1, 1, 0, 0, 0, 0]
    lin = ["A", "B", "C", "D", "E", "F", "G", "H"]
    r = confound_report(labels, lin)
    assert r["verdict"] == CONFOUNDED and r["shared_lineages"] == 0


def test_confounded_one_shared_like_cef():
    # mirror cef: many distinct lineages, exactly 1 shared -> CONFOUNDED
    labels = [1] * 13 + [0] * 17
    lin = [f"R{i}" for i in range(12)] + ["SH"] + [f"S{i}" for i in range(16)] + ["SH"]
    r = confound_report(labels, lin)
    assert r["verdict"] == CONFOUNDED and r["shared_lineages"] == 1


def test_de_confounded_many_shared_like_cipro():
    # 6 shared lineages, each with both R and S; minority shared fraction >= 0.10
    shared = ["L1", "L2", "L3", "L4", "L5", "L6"]
    labels, lin = [], []
    for s in shared:                       # one R + one S per shared lineage
        labels += [1, 0]; lin += [s, s]
    labels += [1] * 6 + [0] * 40           # plus class-unique tails (keeps minority frac ~0.13)
    lin += [f"RU{i}" for i in range(6)] + [f"SU{i}" for i in range(40)]
    r = confound_report(labels, lin)
    assert r["verdict"] == DE_CONFOUNDED and r["shared_lineages"] == 6
    assert r["minority_shared_frac"] >= 0.10


def test_warn_borderline_shared():
    # exactly 3 shared lineages with good minority fraction -> WARN (not clean, not confounded)
    shared = ["L1", "L2", "L3"]
    labels, lin = [], []
    for s in shared:
        labels += [1, 0]; lin += [s, s]
    labels += [1, 0]; lin += ["RU", "SU"]
    r = confound_report(labels, lin)
    assert r["verdict"] == WARN and r["shared_lineages"] == 3


def test_geography_downgrades_clean_to_warn():
    # clean lineages (6 shared) but geography perfectly separable -> WARN
    shared = ["L1", "L2", "L3", "L4", "L5", "L6"]
    labels, lin, reg = [], [], []
    for s in shared:
        labels += [1, 0]; lin += [s, s]; reg += ["USA", "Kenya"]
    labels += [1] * 6 + [0] * 40
    lin += [f"RU{i}" for i in range(6)] + [f"SU{i}" for i in range(40)]
    reg += ["USA"] * 6 + ["Kenya"] * 40
    r = confound_report(labels, lin, reg)
    assert r["geo"]["separable"] is True and r["verdict"] == WARN


def test_degenerate_single_class():
    r = confound_report([1, 1, 1], ["A", "B", "C"])
    assert r["verdict"] == CONFOUNDED and "degenerate" in r["reason"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
