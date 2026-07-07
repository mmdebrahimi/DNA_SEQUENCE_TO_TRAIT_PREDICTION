"""Pins the CYP2D6-CYP2D7 PSV evidence table (Phase A — the read-level hybrid-identity falsifier).

The /brainstorm (2026-07-07) gated the hybrid-IDENTITY classifier behind a cheap decisive falsifier: does the
paired-coordinate PSV D6-fraction PROFILE reproduce the Cyrius hybrid signatures on a small labelled panel?
These tests pin the pure mpileup-base-count logic + the committed falsifier verdict (GO: *68/*36/*13 signalled,
non-hybrids flat). The real-CRAM pileups are gitignored (large); the committed artifact is the evidence JSON.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from cyp2d6_psv_evidence import _base_counts, load_psvs, sample_evidence  # noqa: E402

_EVIDENCE = REPO / "wiki" / "cyp2d6_psv_evidence.json"


# --- pure mpileup base-count parser ---
def test_base_counts_plain():
    assert _base_counts("GGGGaaCC") == {"A": 2, "C": 2, "G": 4, "T": 0}


def test_base_counts_strips_markers():
    # ^X = read-start (X=mapq char, skipped); $ = read-end; +2AC / -1G = indels (inner bases NOT counted)
    assert _base_counts("^GG$G") == {"A": 0, "C": 0, "G": 2, "T": 0}       # ^G skips the G(mapq), then G,G
    assert _base_counts("A+2CCA") == {"A": 2, "C": 0, "G": 0, "T": 0}      # the +2CC indel bases skipped
    assert _base_counts("T-1GT") == {"A": 0, "C": 0, "G": 0, "T": 2}       # the -1G indel base skipped


def test_load_real_psv_file():
    psvs = load_psvs()
    assert len(psvs) == 117                                      # the Cyrius 117 differentiating bases
    p = psvs[0]
    assert p["chrom"] == "chr22" and {"pos_d6", "base_d6", "pos_d7", "base_d7", "annotation"} <= set(p)


def test_sample_evidence_synthetic(tmp_path):
    """A synthetic 1-PSV pileup: all base_d6 reads -> D6-fraction 1.0."""
    psvs = load_psvs()[:1]
    p = psvs[0]
    (tmp_path / "S.d6.txt").write_text(f"{p['pos_d6']}\t20\t{p['base_d6'] * 20}\n", encoding="utf-8")
    (tmp_path / "S.d7.txt").write_text(f"{p['pos_d7']}\t20\t{p['base_d6'] * 20}\n", encoding="utf-8")
    ev = sample_evidence("S", tmp_path, psvs)
    hit = [x for x in ev["per_psv"] if x["d6_fraction"] is not None][0]
    assert hit["d6_fraction"] == 1.0 and hit["d6_like"] == 40 and hit["d7_like"] == 0


# --- committed falsifier verdict (real-CRAM panel) ---
def test_committed_falsifier_is_go():
    import json
    if not _EVIDENCE.exists():
        import pytest
        pytest.skip("evidence artifact not present")
    d = json.loads(_EVIDENCE.read_text(encoding="utf-8"))
    v = d["falsifier_verdict"]
    assert v["falsifier"] == "GO"
    assert v["n_hybrid_signalled"] == v["n_hybrid"] >= 3        # *68/*36/*13 all signalled
    assert v["n_nonhybrid_flat"] == v["n_nonhybrid"] >= 3       # normal/dup/deletion all flat
    # the three hybrid signatures are distinct
    sig = {s["label"]: s["identity_signal"] for s in d["samples"]}
    assert sig["hyb68"].startswith("directional_5p_high")
    assert sig["hyb13"].startswith("directional_5p_low")
    assert sig["hyb36"].startswith("exon9_tip_dip")


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
