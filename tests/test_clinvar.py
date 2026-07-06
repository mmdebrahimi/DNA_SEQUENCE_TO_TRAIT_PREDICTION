"""Offline pins for the deterministic ClinVar variant decoder (loads the committed panel; no VCF/network)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.clinvar import INDETERMINATE, ClinVarCall, ClinVarDecoder, _verdict  # noqa: E402


def test_verdict_mapping():
    assert _verdict("Pathogenic") == "PATHOGENIC" and _verdict("Likely_pathogenic") == "PATHOGENIC"
    assert _verdict("Benign") == "BENIGN" and _verdict("Likely_benign") == "BENIGN"
    assert _verdict("Uncertain_significance") == INDETERMINATE


def test_synthetic_decoder():
    dec = ClinVarDecoder({("7", "100", "A", "T"): {"significance": "Pathogenic",
                          "review_status": "reviewed_by_expert_panel", "gene": "X", "disease": "d", "clinvar_id": "1"}})
    c = dec.call("chr7", 100, "a", "t")            # chr-prefix + lowercase tolerated
    assert c.verdict == "PATHOGENIC" and c.stars == 3 and c.gene == "X"
    u = dec.call("7", 999, "A", "T")               # not in panel -> honest INDETERMINATE (not benign)
    assert u.verdict == INDETERMINATE and "not-in-panel" in u.provenance


def _panel():
    p = Path(__file__).resolve().parent.parent / "data" / "clinvar" / "clinvar_panel.tsv.gz"
    if not p.exists():
        import pytest
        pytest.skip("committed ClinVar panel not present")
    return ClinVarDecoder.from_tsv(p)


def test_committed_panel_f508del_pathogenic():
    dec = _panel()
    # canonical CFTR F508del (ClinVar 7105, GRCh38 chr7:117559590 ATCT>A, 4-star practice_guideline)
    c = dec.call("7", 117559590, "ATCT", "A")
    assert c.verdict == "PATHOGENIC" and c.stars == 4 and c.gene == "CFTR" and "7105" in c.provenance


def test_committed_panel_benign_and_unknown():
    dec = _panel()
    b = dec.call("11", 5225245, "C", "T")          # a benign HBB anchor
    assert b.verdict == "BENIGN" and b.gene == "HBB"
    assert dec.call("7", 1, "A", "C").verdict == INDETERMINATE    # unknown -> honest abstain
    assert len(dec.table) > 1000                   # the real panel is large (~31k P/LP+B/LB variants)
