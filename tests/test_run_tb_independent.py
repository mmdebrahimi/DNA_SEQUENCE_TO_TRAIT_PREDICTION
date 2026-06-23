"""Offline tests for the independent-TB runner pure logic (PASS-masking + aggregation)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.run_tb_independent_amr_portal import _pass_mask, aggregate  # noqa: E402


def test_pass_mask_only_filter_column():
    vcf = ("##fileformat=VCFv4.1\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
           "NC_000962.3\t761155\t.\tC\tT\t60\t.\tQNAME=x\tGT\t1/1\n")
    out = _pass_mask(vcf)
    body = [l for l in out.splitlines() if not l.startswith("#")][0].split("\t")
    assert body[6] == "PASS"               # FILTER '.' -> PASS
    assert body[1] == "761155" and body[3] == "C" and body[4] == "T"   # other columns untouched
    assert out.splitlines()[0] == "##fileformat=VCFv4.1"               # header lines untouched


def test_aggregate_confusion():
    rows = [
        {"rif_label": "R", "rif_pred": "R", "inh_label": "R", "inh_pred": "S"},  # rif TP, inh FN
        {"rif_label": "S", "rif_pred": "S", "inh_label": "S", "inh_pred": "S"},  # rif TN, inh TN
        {"rif_label": "R", "rif_pred": "S", "inh_label": "R", "inh_pred": "R"},  # rif FN, inh TP
        {"rif_label": "S", "rif_pred": "R", "inh_label": "", "inh_pred": "S"},   # rif FP, inh skip (no label)
    ]
    agg = aggregate(rows)
    assert (agg["rifampicin"]["tp"], agg["rifampicin"]["fp"], agg["rifampicin"]["tn"], agg["rifampicin"]["fn"]) == (1, 1, 1, 1)
    assert agg["rifampicin"]["sens"] == 0.5 and agg["rifampicin"]["spec"] == 0.5
    assert agg["isoniazid"]["n_R"] == 2 and agg["isoniazid"]["n_S"] == 1   # the blank-label row excluded


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
