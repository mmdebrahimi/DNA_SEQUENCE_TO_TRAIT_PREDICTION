"""Tests for the hardened cohort de-confound gate (dna_decode/eval/cohort_deconfound).

v2 (2026-06-04): within-group label CONTRAST per provenance axis + matched-support floor + 3-state
promotability contract. Pins: BLOCK on no/thin within-lineage contrast, ADMIT only on real matched
support with no aliasing secondary axis, WARN (non-promotable) for borderline/secondary-aliased,
MLST-missing normalization. Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.eval.cohort_deconfound import (
    CONFOUNDED, DE_CONFOUNDED, WARN, confound_report,
)


def _interleave(n_lineages, r_per, s_per, prefix="L"):
    """n_lineages each carrying r_per R + s_per S → (labels, lineages)."""
    labels, lin = [], []
    for i in range(n_lineages):
        labels += [1] * r_per + [0] * s_per
        lin += [f"{prefix}{i}"] * (r_per + s_per)
    return labels, lin


def test_block_disjoint_lineages():
    labels = [1, 1, 1, 1, 0, 0, 0, 0]
    lin = ["A", "B", "C", "D", "E", "F", "G", "H"]
    r = confound_report(labels, lin)
    assert r["verdict"] == CONFOUNDED and r["promotable"] is False


def test_block_one_shared_like_cef():
    labels = [1] * 13 + [0] * 17
    lin = [f"R{i}" for i in range(12)] + ["SH"] + [f"S{i}" for i in range(16)] + ["SH"]
    r = confound_report(labels, lin)
    assert r["verdict"] == CONFOUNDED and r["primary"]["shared_groups"] == 1


def test_block_thin_matched_support():
    # 4 shared lineages (>= min 3) but each only 1R+1S → matched_minority 4 < 5 → BLOCK.
    labels, lin = _interleave(4, 1, 1)
    r = confound_report(labels, lin)
    assert r["verdict"] == CONFOUNDED and r["primary"]["matched_minority"] == 4


def test_admit_clean():
    # 6 shared lineages x 2R+2S → shared 6, matched 12; no secondary axis → ADMIT + promotable.
    labels, lin = _interleave(6, 2, 2)
    r = confound_report(labels, lin)
    assert r["verdict"] == DE_CONFOUNDED and r["promotable"] is True


def test_warn_secondary_axis_aliases():
    # primary clean, but country perfectly separable (all R=USA, all S=Kenya) → WARN, non-promotable.
    labels, lin = _interleave(6, 2, 2)
    regions = ["USA" if y == 1 else "Kenya" for y in labels]
    r = confound_report(labels, lin, regions)
    assert r["verdict"] == WARN and r["promotable"] is False


def test_warn_borderline_primary():
    # 3 shared lineages (>= min 3 but < clean 5) with good matched support → WARN.
    labels, lin = _interleave(3, 2, 2)
    r = confound_report(labels, lin)
    assert r["verdict"] == WARN and r["promotable"] is False


def test_missing_mlst_not_credited():
    # placeholder MLST values must NOT count as a shared lineage.
    labels = [1, 0, 1, 0]
    lin = ["unknown", "", None, "nan"]
    r = confound_report(labels, lin)
    assert r["verdict"] == CONFOUNDED and r["primary"]["shared_groups"] == 0


def test_admit_not_downgraded_by_mixed_secondary():
    # secondary country present AND mixed (each region has both classes) → stays ADMIT.
    labels, lin = _interleave(6, 2, 2)
    regions = ["USA", "Kenya", "USA", "Kenya"] * 6   # both regions carry both labels
    r = confound_report(labels, lin, regions)
    assert r["verdict"] == DE_CONFOUNDED and r["promotable"] is True


def test_degenerate_single_class():
    r = confound_report([1, 1, 1], ["A", "B", "C"])
    assert r["verdict"] == CONFOUNDED and "degenerate" in r["reason"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
